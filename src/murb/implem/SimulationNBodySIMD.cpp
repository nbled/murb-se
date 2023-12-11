#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "mipp.h"

#include "SimulationNBodySIMD.hpp"

SimulationNBodySIMD::SimulationNBodySIMD(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    this->accelerations.ax.resize(this->getBodies().getN());
    this->accelerations.ay.resize(this->getBodies().getN());
    this->accelerations.az.resize(this->getBodies().getN());

}

void SimulationNBodySIMD::initIteration()
{
    std::fill(this->accelerations.ax.begin(),this->accelerations.ax.end(),0);
    std::fill(this->accelerations.ay.begin(),this->accelerations.ay.end(),0);
    std::fill(this->accelerations.az.begin(),this->accelerations.az.end(),0);
}

void print_mipp_reg(mipp::Reg<float> reg) {
    constexpr int N = mipp::N<float>();
    float debug[N];
    reg.store(debug);
    printf("{");
    for (int i=0;i<N;i++) {
        printf("%f, ",debug[i]);
    }
    printf("\b}\n");
}


void SimulationNBodySIMD::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();

    // compute e²
    const float softSquared = std::pow(this->soft, 2); // 1 flops
    unsigned long n_bodies = this->getBodies().getN();
    constexpr int N = mipp::N<float>();
    unsigned long n_bodies_rounded = (n_bodies / N) * N;

    // flops = n² * 20
    for (unsigned long iBody = 0; iBody < n_bodies; iBody++) {
        unsigned long jBody;
        mipp::Reg<float> i_qx = d.qx[iBody]; //we duplicate for i, and take different values for j.
        mipp::Reg<float> i_qy = d.qy[iBody];
        mipp::Reg<float> i_qz = d.qz[iBody];
        mipp::Reg<float> i_m = d.m[iBody];


        mipp::Reg<float> softSquared_v = softSquared;
        mipp::Reg<float> G_v = this->G;

        for (jBody = 0; jBody <= n_bodies_rounded; jBody += N) {
            mipp::Reg<float> j_qx = &d.qx[jBody];
            mipp::Reg<float> rijx = j_qx - i_qx;

            mipp::Reg<float> j_qy = &d.qy[jBody];
            mipp::Reg<float> rijy = j_qy - i_qy;

            mipp::Reg<float> j_qz = &d.qz[jBody];
            mipp::Reg<float> rijz = j_qz - i_qz;

            mipp::Reg<float> rijSquared = rijx * rijx + rijy * rijy + rijz * rijz;
            // mipp::Reg<float> rijSquared = rijx * rijx;  //fma doesn't seem to add to perf
            // rijSquared = mipp::fmadd(rijy,rijy,rijSquared);
            // rijSquared = mipp::fmadd(rijz,rijz,rijSquared);



            mipp::Reg<float> x = G_v / ((rijSquared + softSquared_v) * mipp::sqrt(rijSquared + softSquared_v));
            mipp::Reg<float> j_m = &d.m[jBody];
            mipp::Reg<float> ai = x * j_m; // 1 flops


            mipp::Reg<float> i_ax = ai * rijx;
            this->accelerations.ax[iBody] += mipp::hadd(i_ax); //Sum all elements
            mipp::Reg<float> i_ay = ai * rijy;
            this->accelerations.ay[iBody] += mipp::hadd(i_ay);
            mipp::Reg<float> i_az = ai * rijz;
            this->accelerations.az[iBody] += mipp::hadd(i_az);


            // mipp::Reg<float> aj = x * i_m; // 1 flops


            // mipp::Reg<float> j_ax = aj * rijx;
            // mipp::Reg<float> j_ay = aj * rijy;
            // mipp::Reg<float> j_az = aj * rijz;

            // float tabx[N];
            // float taby[N];
            // float tabz[N];

            // j_ax.store(tabx);
            // j_ay.store(taby);
            // j_az.store(tabz);

            // for (int i=0;i<N;i++) {
            //     this->accelerations.ax[jBody+i] += tabx[i];
            //     this->accelerations.ay[jBody+i] += taby[i];
            //     this->accelerations.az[jBody+i] += tabz[i];
            // }




        }
        
        //We finish the remaining bodies using the code from optim
        for (; jBody < n_bodies; jBody++) {
            const float rijx = d.qx[jBody] - d.qx[iBody]; // 1 flop
            const float rijy = d.qy[jBody] - d.qy[iBody]; // 1 flop
            const float rijz = d.qz[jBody] - d.qz[iBody]; // 1 flop

            // compute the || rij ||² distance between body i and body j
            const float rijSquared = rijx * rijx + rijy * rijy + rijz * rijz; // 5 flops

            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            const float x = this->G / ((rijSquared + softSquared) * std::sqrt(rijSquared + softSquared));
            const float ai = x * d.m[jBody]; // 1 flops
//            const float aj = x * d.m[iBody]; // 1 flops

            // add the acceleration value into the acceleration vector: ai += || ai ||.rij
            this->accelerations.ax[iBody] += ai * rijx; // 2 flops
            this->accelerations.ay[iBody] += ai * rijy; // 2 flops
            this->accelerations.az[iBody] += ai * rijz; // 2 flops

            // this->accelerations.ax[jBody] += aj * -rijx; // 2 flops
            // this->accelerations.ay[jBody] += aj * -rijy; // 2 flops
            // this->accelerations.az[jBody] += aj * -rijz; // 2 flops

        }
    }
}

void SimulationNBodySIMD::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
