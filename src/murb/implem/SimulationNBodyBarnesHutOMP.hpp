#ifndef SIMULATION_N_BODY_BARNES_HUT_OMP_HPP_
#define SIMULATION_N_BODY_BARNES_HUT_OMP_HPP_

#include <string>

#include "core/SimulationNBodyInterface.hpp"
#include "SimulationNBodyBarnesHut.hpp"
#include "core/Bodies.hpp"


class SimulationNBodyBarnesHutOMP : public SimulationNBodyBarnesHut {


  public:
    SimulationNBodyBarnesHutOMP(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodyBarnesHutOMP() = default;

  protected:
    void computeBodiesAcceleration() override;

};

#endif /* SIMULATION_N_BODY_OPTIM_HPP_ */
