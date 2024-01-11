#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "mipp.h"

#include "SimulationNBodySIMD_OMP.hpp"

SimulationNBodySIMD_OMP::SimulationNBodySIMD_OMP(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    this->accelerations.ax.resize(this->getBodies().getN() + this->getBodies().getPadding());
    this->accelerations.ay.resize(this->getBodies().getN() + this->getBodies().getPadding());
    this->accelerations.az.resize(this->getBodies().getN() + this->getBodies().getPadding());
}

void SimulationNBodySIMD_OMP::initIteration()
{
    std::fill(this->accelerations.ax.begin(),this->accelerations.ax.end(),0);
    std::fill(this->accelerations.ay.begin(),this->accelerations.ay.end(),0);
    std::fill(this->accelerations.az.begin(),this->accelerations.az.end(),0);
}

//Second version
void SimulationNBodySIMD_OMP::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();
    // compute e²
    const float softSquared = std::pow(this->soft, 2); // 1 flops
    unsigned long n_bodies = this->getBodies().getN();
    constexpr int N = mipp::N<float>();

    // flops = n² * 20
    unsigned long iBody;
    #pragma omp parallel \
            firstprivate(softSquared,n_bodies,N)
    {
    #pragma omp for schedule(guided,1)
    for (iBody = 0; iBody < n_bodies; iBody+=N) {
        unsigned long jBody;
        mipp::Reg<float> i_qx = &d.qx[iBody];
        mipp::Reg<float> i_qy = &d.qy[iBody];
        mipp::Reg<float> i_qz = &d.qz[iBody];
        //mipp::Reg<float> i_m = &d.m[iBody];


        mipp::Reg<float> softSquared_v = softSquared;
        mipp::Reg<float> G_v = this->G;

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


            // mipp::Reg<float> aj = x * i_m; // 1 flops


            // mipp::Reg<float> j_ax = aj * rijx;
            // mipp::Reg<float> j_ay = aj * rijy;
            // mipp::Reg<float> j_az = aj * rijz;

        }
        

        ax.store(&this->accelerations.ax[iBody]);
        ay.store(&this->accelerations.ay[iBody]);
        az.store(&this->accelerations.az[iBody]);

    }
    }
}



void SimulationNBodySIMD_OMP::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
