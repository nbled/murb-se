#ifndef SIMULATION_N_BODY_SIMD_OMP_HPP_
#define SIMULATION_N_BODY_SIMD_OMP_HPP_

#include <string>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodySIMD_OMP : public SimulationNBodyInterface {
  protected:
    accSoA_t<float> accelerations;

  public:
    SimulationNBodySIMD_OMP(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodySIMD_OMP() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
};

#endif /* SIMULATION_N_BODY_OPTIM_HPP_ */
