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
#include <random> // std::random_device
#include <CL/cl.hpp>
#include <CL/cl_platform.h>

#include "boundingbox.h"




class NBSystem
{
private:
    // Number of particles
    unsigned int pnum;
    // Particles postions. pos[i].s3 contains mass of the particle
    std::vector<cl_float4> pos;
    // Particle velocities
    std::vector<cl_float4> vel;

    // Vertices of boundaries triangles. Three consecutive elements make triangle.
    std::vector<cl_float4> tri;
    // Number of vertices
    unsigned int vnum;

    // Minimal bounding box
    BoundingBox<double, 3> box;

    // Time step etimation coefficient. Used in updateDt.
    double dtCoef;

//    // The depth of the potential well
//    double ljEps;
//    // The distance at which the potential reaches its minimum
//    double ljRmin;

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

    double getEnKin(cl::Event &event);
    double getEnPot(cl::Event &event);
    void updateAcc(cl::Event &event);
    double estDt(cl::Event &event);
    void updatePos(double dt, cl::Event &event);
    void syncArrays();

public:
    NBSystem();
    ~NBSystem();

    std::string setupCL(int platformId, cl_device_type deviceType);

    double getEnTot();
    double getEnPot();
    double getEnKin();
    double getEstDt();
    unsigned int getPNum() const;

    std::vector<std::array<double, 3> > getMinBox() const;

    const std::vector<cl_float4>* posData();
    const std::vector<cl_float4>* velData();
    const std::vector<cl_float4>* bndData();

    void setDtCoef(double coeff);



    void addParticle(std::vector<cl_float4> &pos,
                     std::vector<cl_float4> &vel);
    void addParticle(unsigned int num,
                     std::array<double, 2> mass,
                     BoundingBox<double, 3> posBox,
                     BoundingBox<double, 3> velBox);
    void removeParticle(unsigned int i);

    void addBoundary(std::vector<cl_float4> &vertices);

    void evolve(double timeStep);
    void evolveIn(double interval);

    std::string boundariesString(unsigned int prec);
};

#endif // NBSYSTEM_H
