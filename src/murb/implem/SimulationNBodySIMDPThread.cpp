#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <atomic>
#include <pthread.h>

#include "mipp.h"

#include "SimulationNBodySIMDPThread.hpp"

#define NB_THREADS 6
#define BATCH_SIZE 2

SimulationNBodySIMDPThread::SimulationNBodySIMDPThread(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    this->accelerations.ax.resize(this->getBodies().getN() + this->getBodies().getPadding());
    this->accelerations.ay.resize(this->getBodies().getN() + this->getBodies().getPadding());
    this->accelerations.az.resize(this->getBodies().getN() + this->getBodies().getPadding());

}

void SimulationNBodySIMDPThread::initIteration()
{
    std::fill(this->accelerations.ax.begin(),this->accelerations.ax.end(),0);
    std::fill(this->accelerations.ay.begin(),this->accelerations.ay.end(),0);
    std::fill(this->accelerations.az.begin(),this->accelerations.az.end(),0);
}


struct args {
    std::atomic<unsigned long>* current;
    SimulationNBodySIMDPThread *that;
    const float soft;
    const float G;
};
//Second version
void* worker(void *arg)
{
    struct args *arg_s = (struct args *) arg;
    const dataSoA_t<float> &d = arg_s->that->getBodies().getDataSoA();
    // compute e²
    const float softSquared = std::pow(arg_s->soft, 2); // 1 flops
    unsigned long n_bodies = arg_s->that->getBodies().getN();
    constexpr int N = mipp::N<float>();

    // flops = n² * 20

    unsigned long iBody_main = arg_s->current->fetch_add(N*BATCH_SIZE);

    // printf("thread %d got %ld\n", pthread_self(),iBody_main);


    while (iBody_main < n_bodies) {
        for (unsigned long iBody = iBody_main; iBody < iBody_main + N*BATCH_SIZE && iBody < n_bodies; iBody+=N) {
            // printf("thread %d : calculating %ld\n",pthread_self(),iBody);
            unsigned long jBody;
            mipp::Reg<float> i_qx = &d.qx[iBody];
            mipp::Reg<float> i_qy = &d.qy[iBody];
            mipp::Reg<float> i_qz = &d.qz[iBody];


            mipp::Reg<float> softSquared_v = softSquared;
            mipp::Reg<float> G_v = arg_s->G;

            mipp::Reg<float> ax = 0.0;
            mipp::Reg<float> ay = 0.0;
            mipp::Reg<float> az = 0.0;

            for (jBody = 0; jBody < n_bodies; jBody += 1) {
                mipp::Reg<float> j_qx = d.qx[jBody]; //We duplicate : the same j for multiple i.
                mipp::Reg<float> rijx = j_qx - i_qx;

                mipp::Reg<float> j_qy = d.qy[jBody];
                mipp::Reg<float> rijy = j_qy - i_qy;

                mipp::Reg<float> j_qz = d.qz[jBody];
                mipp::Reg<float> rijz = j_qz - i_qz;

                mipp::Reg<float> rijSquared = rijx * rijx + rijy * rijy + rijz * rijz;
                // mipp::Reg<float> rijSquared = rijx * rijx;  //fma doesn't seem to add to perf
                // rijSquared = mipp::fmadd(rijy,rijy,rijSquared);
                // rijSquared = mipp::fmadd(rijz,rijz,rijSquared);



                mipp::Reg<float> x = G_v / ((rijSquared + softSquared_v) * mipp::sqrt(rijSquared + softSquared_v));
                mipp::Reg<float> j_m = d.m[jBody];
                mipp::Reg<float> ai = x * j_m; // 1 flops


                mipp::Reg<float> i_ax = ai * rijx;
                ax += i_ax;
                mipp::Reg<float> i_ay = ai * rijy;
                ay += i_ay;
                mipp::Reg<float> i_az = ai * rijz;
                az += i_az;


            }
            

            ax.store(&arg_s->that->accelerations.ax[iBody]);
            ay.store(&arg_s->that->accelerations.ay[iBody]);
            az.store(&arg_s->that->accelerations.az[iBody]);
        }
        iBody_main = arg_s->current->fetch_add(N*BATCH_SIZE);
        // printf("thread %d got %ld\n", pthread_self(),iBody_main);
    }
    // printf("thread %d finished\n",pthread_self());
    return NULL;
}



void SimulationNBodySIMDPThread::computeBodiesAcceleration() {
    std::atomic<unsigned long> current;
    current = 0;
    pthread_t threads[NB_THREADS];
    struct args arg = {
        &current,
        this,
        this->soft,
        this->G
    };
    for (int i=0;i<NB_THREADS;i++) {
        pthread_create(&threads[i],NULL, &worker,(void*)&arg);
        // printf("thread %d created (%ld) \n",i,threads[i]);
    }
    for (int i=0;i<NB_THREADS;i++) {
        pthread_join(threads[i],NULL);
        // printf("thread %d joined\n",i);
    }
}


void SimulationNBodySIMDPThread::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
