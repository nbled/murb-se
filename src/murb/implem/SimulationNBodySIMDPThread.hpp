#ifndef SIMULATION_N_BODY_SIMD_PTHREAD_HPP_
#define SIMULATION_N_BODY_SIMD_PTHREAD_HPP_

#include <string>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodySIMDPThread : public SimulationNBodyInterface {
  public:
    accSoA_t<float> accelerations;

  public:
    SimulationNBodySIMDPThread(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodySIMDPThread() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
};

#endif /* SIMULATION_N_BODY_OPTIM_HPP_ */
