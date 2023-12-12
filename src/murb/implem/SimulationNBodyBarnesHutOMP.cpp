#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyBarnesHut.hpp"
#include "SimulationNBodyBarnesHutOMP.hpp"

SimulationNBodyBarnesHutOMP::SimulationNBodyBarnesHutOMP(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyBarnesHut(nBodies, scheme, soft, randInit)
{
}


void SimulationNBodyBarnesHutOMP::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();

    unsigned long n_bodies = this->getBodies().getN();
    #pragma omp parallel
    {
        #pragma omp for
        for (unsigned long iBody = 0; iBody < n_bodies; iBody++) {
            // printf("Computing for %e %e %e\n",d[iBody].qx,d[iBody].qy,d[iBody].qz);
            this->computeBodyAcceleration(this->tree,&d[iBody],&this->accelerations[iBody].ax,&this->accelerations[iBody].ay,&this->accelerations[iBody].az,0);
        }
    }
}
