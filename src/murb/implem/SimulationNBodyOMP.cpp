#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyOMP.hpp"

SimulationNBodyOMP::SimulationNBodyOMP(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    this->accelerations.resize(this->getBodies().getN());
}

void SimulationNBodyOMP::initIteration()
{
    unsigned long nBodies = this->getBodies().getN();

    for (unsigned long iBody = 0; iBody < nBodies; iBody++) {
        this->accelerations[iBody].ax = 0.f;
        this->accelerations[iBody].ay = 0.f;
        this->accelerations[iBody].az = 0.f;
    }
}

void SimulationNBodyOMP::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();

    // compute e²
    const float softSquared = std::pow(this->soft, 2); // 1 flops
    unsigned long nBodies = this->getBodies().getN();

    // flops = n² * 20
#pragma omp parallel
{
#pragma omp for
    for (unsigned long iBody = 0; iBody < nBodies; iBody++) {
        // flops = n * 20
        for (unsigned long jBody = iBody + 1; jBody < nBodies; jBody++) {
            const float rijx = d.qx[jBody] - d.qx[iBody]; // 1 flop
            const float rijy = d.qy[jBody] - d.qy[iBody]; // 1 flop
            const float rijz = d.qz[jBody] - d.qz[iBody]; // 1 flop

            // compute the || rij ||² distance between body i and body j
            const float rijSquared = rijx * rijx + rijy * rijy + rijz * rijz; // 5 flops

            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            const float x = this->G / ((rijSquared + softSquared) * std::sqrt(rijSquared + softSquared));

            const float ai = x * d.m[jBody]; // 1 flops
            const float aj = x * d.m[iBody]; // 1 flops

            // add the acceleration value into the acceleration vector: ai += || ai ||.rij
            this->accelerations[iBody].ax += ai * rijx; // 2 flops
            this->accelerations[iBody].ay += ai * rijy; // 2 flops
            this->accelerations[iBody].az += ai * rijz; // 2 flops
            
            this->accelerations[jBody].ax += aj * -rijx; // 2 flops
            this->accelerations[jBody].ay += aj * -rijy; // 2 flops
            this->accelerations[jBody].az += aj * -rijz; // 2 flops
        }
    }
}

}

void SimulationNBodyOMP::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
