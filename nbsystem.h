#ifndef NBSYSTEM_H
#define NBSYSTEM_H

#define __CL_ENABLE_EXCEPTIONS

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setprecision
#include <algorithm>
#include <CL/cl.hpp>
#include <libconfig.h++>


class NBSystem
{
private:
    // Number of particles
    unsigned int pnum;
    // Particles postions
    std::vector<cl_float4> pos;
    // Particle velocities
    std::vector<cl_float4> vel;
    // Vertices of boundaries triangles. Three consecutive elements make triangle.
    std::vector<cl_float4> tri;
    // Intitial total energy
    double enIn;
    // Kinetic energy
    double enKin;
    // Potential energy
    double enPot;
    // Current time
    double timeCur;
    // End time
    double timeEnd;
    // Time step
    double timeStep;
    double timeStepMax;
    // Time step etimation coefficient. Used in estimateTimeStep.
    double dtCoef;


    // Output precision
    int outPrec;
    // Output file name
    std::string outFileName;
    // Write time interval
    double timeWrite;

    // OpenCL device type (CL_DEVICE_TYPE_CPU / CL_DEVICE_TYPE_GPU)
    int platformId;
    cl_device_type deviceType;
    cl::Context context;
    cl::Device device;
    cl::Program program;
    cl::CommandQueue queue;

    // OpenCL kernels
    cl::Kernel stepInTime;
    cl::Kernel estimateDt;
    cl::Kernel calcAcc;
    cl::Kernel checkBoundaries;
    cl::Kernel calcEnKin;
    cl::Kernel calcEnPot;

    // Device buffers
    cl::Buffer d_pos;
    cl::Buffer d_vel;
    cl::Buffer d_newPos;
    cl::Buffer d_newVel;
    cl::Buffer d_acc;
    cl::Buffer d_est;
    cl::Buffer d_tri;
    cl::Buffer d_enKin;
    cl::Buffer d_enPot;

    cl_float *h_enKin;
    cl_float *h_enPot;
    cl_float *h_esTS;


    void initializeCLBuffers();
    void setupCLKernelsArgs();

    void addParticle(cl_float &m, cl_float4 &p, cl_float4 &v);
    void updateEn(cl::Event &event);

public:
    NBSystem();
    ~NBSystem();

    void setupCL(int platformId, cl_device_type deviceType);

    double getEnIn();
    double getEnTot();
    double getEnPot();
    double getEnKin();


    void setTimeEnd(double &time);
    void setTimeStepMax(double &time);

    void readSettings(std::string fileName);

    // Warnings with std::vector<cl_float> &mass,
    void addParticle(std::vector<float> &mass,
                     std::vector<cl_float4> &pos,
                     std::vector<cl_float4> &vel);
    void addParticle(std::string fileName);
    void removeParticle(unsigned int i);

    void addBoundary(std::vector<cl_float4> &vertices);
    void addBoundary(std::string fileName);

    double estimateTimeStep(cl::Event event);
    void evolve(double dt);

    std::string energyString();
    std::string stateString(int pId);
    std::string stateString();
    std::string boundariesString();
    std::string toString();
};

#endif // NBSYSTEM_H
