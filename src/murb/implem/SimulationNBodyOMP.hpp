#ifndef SIMULATION_N_BODY_OMP_HPP_
#define SIMULATION_N_BODY_OMP_HPP_

#include <string>
#include <omp.h>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodyOMP : public SimulationNBodyInterface {
  protected:
    std::vector<accAoS_t<float>> accelerations; /*!< Array of body acceleration structures. */

  public:
    SimulationNBodyOMP(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodyOMP() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
};

#endif /* SIMULATION_N_BODY_OMP_HPP_ */
