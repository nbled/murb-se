#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>

#include "SimulationNBodyOpenCL.hpp"

SimulationNBodyOpenCL::SimulationNBodyOpenCL(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    
    this->accelerations.ax.resize(this->getBodies().getN());
    this->accelerations.ay.resize(this->getBodies().getN());
    this->accelerations.az.resize(this->getBodies().getN());

    cl_uint num_platforms = 0;
    cl_int err = CL_SUCCESS;

    /* Fetch number of available platforms (platform = specific opencl implem) */
    if ((err = clGetPlatformIDs(0, NULL, &num_platforms)) != CL_SUCCESS) {
        std::cerr << "clGetPlatformIDs failed; error = " << err << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Fetch the platforms */
    cl_platform_id *platforms = new cl_platform_id [num_platforms] ();
    if ((err = clGetPlatformIDs(num_platforms, platforms, NULL)) != CL_SUCCESS) {
        std::cerr << "clGetPlatformIDs failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    /* Fetch number of devices assiociated with a platform (device = hardware linked to the opencl implem) */
    cl_uint num_devices;
    if ((err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices)) != CL_SUCCESS) {
        std::cerr << "clGetDeviceIDs failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    /* Fetch the devices */
    cl_device_id *devices = new cl_device_id [num_devices] ();
    if ((err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, num_devices, devices, NULL)) != CL_SUCCESS) {
        std::cerr << "clGetDeviceIDs failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    /* Create the context */
    this->context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "clCreateContext failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    /* Create the command queue */
    this->cmd_queue = clCreateCommandQueue(this->context, devices[0], 0, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "clCreateCommandQueue failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    /* Read kernel source */
    /* TODO: pass the source name as a constructor parameter */
    std::ifstream source_file("kernel/naive.cl");
    if (!source_file.is_open()) {
        std::cerr << "opening 'kernel/naive.cl' failed" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::stringstream source_stream;
    source_stream << source_file.rdbuf();
    std::string source = source_stream.str();

    /* Create & build the program */
    const char *c_source = source.c_str();
    printf("%s\n", c_source);
    const size_t lengths[] = { source.size() };

    cl_program program = clCreateProgramWithSource(this->context, 1, &c_source, lengths, &err);
    if (err != CL_SUCCESS) {
        std::cerr << "clCreateProgramWithSource failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }
    if ((err = clBuildProgram(program, 1, devices, NULL, NULL, NULL)) != CL_SUCCESS) {
        std::cerr << "clBuildProgram failed; error = " << err << std::endl;

        size_t log_size;
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        char *log = new char[log_size] ();
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
        std::cout << log << std::endl;

        exit(EXIT_FAILURE); 
    }

    /* At last, generate the kernel */
    this->kernel = clCreateKernel(program, "computeBodiesAcceleration", &err);
    if (err != CL_SUCCESS) {
        std::cerr << "clCreateKernel failed; error = " << err << std::endl;
        exit(EXIT_FAILURE); 
    }

    std::cout << "kernel built !" << std::endl;

    const dataSoA_t<cl_float> &d = this->getBodies().getDataSoA();


    this->in_buf_qx = clCreateBuffer(this->context, 
        CL_MEM_READ_ONLY , d.qx.size() * sizeof(d.qx[0]), NULL, &err);
    this->in_buf_qy = clCreateBuffer(this->context,
        CL_MEM_READ_ONLY , d.qy.size() * sizeof(d.qy[0]), NULL, &err);
    this->in_buf_qz = clCreateBuffer(this->context, 
        CL_MEM_READ_ONLY , d.qz.size() * sizeof(d.qz[0]), NULL, &err);
    this->in_buf_m = clCreateBuffer(this->context, 
        CL_MEM_READ_ONLY , d.m.size() * sizeof(d.m[0]), NULL, &err);

    this->out_buf_ax = clCreateBuffer(this->context,
        CL_MEM_READ_WRITE , this->accelerations.ax.size() * sizeof(this->accelerations.ax[0]), 
         NULL, &err);
    this->out_buf_ay = clCreateBuffer(this->context,
        CL_MEM_READ_WRITE , this->accelerations.ay.size() * sizeof(this->accelerations.ay[0]), 
        NULL, &err);
    this->out_buf_az = clCreateBuffer(this->context,
        CL_MEM_READ_WRITE , this->accelerations.az.size() * sizeof(this->accelerations.az[0]), 
         NULL, &err);

    std::cout << "kernel buffer created !" << std::endl;

}

void SimulationNBodyOpenCL::initIteration()
{
    unsigned long nBodies = this->getBodies().getN();
    for (unsigned long iBody = 0; iBody < nBodies; iBody++) {
        this->accelerations.ax[iBody] = 0.f;
        this->accelerations.ay[iBody] = 0.f;
        this->accelerations.az[iBody] = 0.f;
    }
}

void SimulationNBodyOpenCL::computeBodiesAcceleration()
{
    const dataSoA_t<cl_float> &d = this->getBodies().getDataSoA();
    const cl_float softSquared = this->soft * this->soft; // 1 flops
    cl_ulong n_bodies = this->getBodies().getN();

    /* CL_MEM_USE_HOST_PTR vs CL_MEM_COPY_HOST_PTR ? */
    /* if we use USE_HOST_PTR can we init only once ? */
    cl_int err = CL_SUCCESS;

    err = clEnqueueWriteBuffer(this->cmd_queue, in_buf_qx,
        CL_TRUE, 0,  d.qx.size() * sizeof(d.qx[0]), (void*) &d.qx[0],0,NULL,NULL);
    err = clEnqueueWriteBuffer(this->cmd_queue,in_buf_qy,
        CL_TRUE, 0,  d.qy.size() * sizeof(d.qy[0]), (void*) &d.qy[0],0,NULL,NULL);
    err = clEnqueueWriteBuffer(this->cmd_queue, in_buf_qz,
        CL_TRUE, 0,  d.qz.size() * sizeof(d.qz[0]), (void*) &d.qz[0],0,NULL,NULL);
    err = clEnqueueWriteBuffer(this->cmd_queue, in_buf_m,
        CL_TRUE, 0,  d.m.size() * sizeof(d.m[0]), (void*) &d.m[0],0,NULL,NULL);

    err = clEnqueueWriteBuffer(this->cmd_queue, out_buf_ax,
        CL_TRUE,0, this->accelerations.ax.size() * sizeof(this->accelerations.ax[0]), 
         (void*) &this->accelerations.ax[0], 0,NULL,NULL);
    err = clEnqueueWriteBuffer(this->cmd_queue, out_buf_ay,
        CL_TRUE,0, this->accelerations.ay.size() * sizeof(this->accelerations.ay[0]), 
         (void*) &this->accelerations.ay[0], 0,NULL,NULL);    
    err = clEnqueueWriteBuffer(this->cmd_queue, out_buf_az,
        CL_TRUE,0, this->accelerations.az.size() * sizeof(this->accelerations.az[0]), 
         (void*) &this->accelerations.az[0], 0,NULL,NULL);

    /* Setup the kernel */
    unsigned t = 0;
    clSetKernelArg(this->kernel, 0, sizeof(cl_mem), &in_buf_qx);
    clSetKernelArg(this->kernel, 1, sizeof(cl_mem), &in_buf_qy);
    clSetKernelArg(this->kernel, 2, sizeof(cl_mem), &in_buf_qz);
    clSetKernelArg(this->kernel, 3, sizeof(cl_mem), &in_buf_m);
    clSetKernelArg(this->kernel, 4, sizeof(cl_mem), &out_buf_ax);
    clSetKernelArg(this->kernel, 5, sizeof(cl_mem), &out_buf_ay);
    clSetKernelArg(this->kernel, 6, sizeof(cl_mem), &out_buf_az);
    clSetKernelArg(this->kernel, 7, sizeof(cl_ulong), &n_bodies);
    clSetKernelArg(this->kernel, 8, sizeof(cl_float), &softSquared);
    clSetKernelArg(this->kernel, 9, sizeof(cl_float), &this->G);

    /* Enqueue kernel */
    const size_t global = n_bodies;
    const size_t local = 8;

    clEnqueueNDRangeKernel(this->cmd_queue, this->kernel, 1, NULL, &global, &local, 0, NULL, NULL);
    
    clEnqueueReadBuffer(this->cmd_queue, out_buf_ax, CL_TRUE, 0, this->accelerations.ax.size() * sizeof(this->accelerations.ax[0]), (void*) &this->accelerations.ax[0],
        0, NULL, NULL);
    clEnqueueReadBuffer(this->cmd_queue, out_buf_ay, CL_TRUE, 0, this->accelerations.ay.size() * sizeof(this->accelerations.ay[0]), (void*) &this->accelerations.ay[0],
        0, NULL, NULL);
    clEnqueueReadBuffer(this->cmd_queue, out_buf_az, CL_TRUE, 0, this->accelerations.az.size() * sizeof(this->accelerations.az[0]), (void*) &this->accelerations.az[0],
        0, NULL, NULL);

    for (unsigned i = 0; i < n_bodies; i++) {
        printf("%e ", this->accelerations.ax[i]);
    }
    printf("\n");

    clFinish(cmd_queue);
}

void SimulationNBodyOpenCL::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
