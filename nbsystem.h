#ifndef NBSYSTEM_H
#define NBSYSTEM_H

#define __CL_ENABLE_EXCEPTIONS

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> // std::setprecision()
#include <algorithm>
#include <cmath> // std::fabs()
#include <CL/cl.hpp>


class NBSystem
{
private:
    // Number of particles
    unsigned int pnum;
    // Particles postions. pos[i].s3 contains mass of the particle
    std::vector<cl_float4> pos;
    // Particle velocities
    std::vector<cl_float4> vel;

    unsigned int vnum;
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
    double dt;
    // Number of time steps
    unsigned int tsNum;
    // Time step etimation coefficient. Used in updateDt.
    double dtCoef;


    // Output precision
    unsigned int outPrec;
    // Output file name
    std::string outFileName;
    // Write time interval
    double timeWrite;

    void addParticle(cl_float4 &p, cl_float4 &v);

    // OpenCL device type (CL_DEVICE_TYPE_CPU / CL_DEVICE_TYPE_GPU)
    int platformId;
    cl_device_type deviceType;
    cl::Context context;
    cl::Device device;
    cl::Program program;
    cl::CommandQueue queue;
    cl::Event event;

    // OpenCL kernels
    cl::Kernel stepInTime;
    cl::Kernel estimateDt;
    cl::Kernel calcAcc;
    cl::Kernel checkBoundaries;
    cl::Kernel calcEnKin;
    cl::Kernel calcEnPot;

    // Device buffers
    size_t vectBuffSize; // Buffer size of vectors
    size_t scalBuffSize; // Buffer size of scalars
    unsigned int curBuffIndex;

    cl::Buffer d_pos[2];
    cl::Buffer d_vel;
    cl::Buffer d_acc;
    cl::Buffer d_tri;
    cl::Buffer d_enKin;
    cl::Buffer d_enPot;
    int d_est_len;     // To not to calculate over and over again
    size_t d_est_size; // same number then estDt() is called
    cl::Buffer d_est;

    cl_float *h_enKin;
    cl_float *h_enPot;
    cl_float *h_esTS;

    // This flag is set then pos and vel are the same as d_pos and d_vel
    bool syncFlag;



    void initializeCLBuffers();
    void setupCLKernelsArgs();

    void updateEn(cl::Event &event);
    double getEnKin(cl::Event &event);
    double getEnPot(cl::Event &event);
    void updateAcc(cl::Event &event);
    double estDt(cl::Event &event);
    void updatePos(double dt, cl::Event &event);
    void syncArrays();


public:
    NBSystem();
    ~NBSystem();

    void setupCL(int platformId, cl_device_type deviceType);

    double getEnIn();
    double getEnTot();
    double getEnPot();
    double getEnKin();
    double getEstDt();
    double getTimeCur();
    double getTSNum();
    unsigned int getPNum();
    cl_float4 getPos(unsigned int i);
    cl_float4 getVel(unsigned int i);

    //std::vector<cl_float4> data();

    void setDtCoef(double coeff);
    void setOutPrec(unsigned int precision);

    void addParticle(std::vector<cl_float4> &pos,
                     std::vector<cl_float4> &vel);
    void addParticle(std::string fileName);
    void removeParticle(unsigned int i);

    void addBoundary(std::vector<cl_float4> &vertices);
    void addBoundary(std::string fileName);

    void evolve(double timeStep);
    void evolveIn(double interval);

    std::string energyString();
    std::string stateString(int pId);
    std::string stateString();
    std::string boundariesString();
    std::string toString();
};

#endif // NBSYSTEM_H
