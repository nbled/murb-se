#ifndef SIMULATION_N_BODY_OPENCL_HPP_
#define SIMULATION_N_BODY_OPENCL_HPP_

#include <string>

/* TODO: check if the TX2 uses the same */
/* NOTE: there's a c++ version but easypap+course uses the c one */
#define CL_TARGET_OPENCL_VERSION 220
#include <CL/cl.h>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodyOpenCL : public SimulationNBodyInterface {
  protected:
    accSoA_t<cl_float> accelerations;
    cl_context context;
    cl_kernel kernel;
    cl_command_queue cmd_queue;

  public:
    SimulationNBodyOpenCL(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodyOpenCL() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
};

#endif /* SIMULATION_N_BODY_OPENCL_HPP_ */
