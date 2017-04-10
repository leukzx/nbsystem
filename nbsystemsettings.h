#ifndef NBSYSTEMSETTINGS_H
#define NBSYSTEMSETTINGS_H

#include <iostream>
#include <libconfig.h++>
#include "nbsystem.h"
#include "boundingbox.h"
#include "CL/cl_platform.h"

class NBSystemSettings
{
public:
    NBSystemSettings(std::string fileName);

    libconfig::Config cfg;

    std::string outFileName;
    std::string particleFileName;
    std::string boundariesFileName;
    int platformId;
    int deviceId;
    cl_device_type deviceType;
    double dtCoef;
    int outPrec;
    double writeInterval;
    double timeStepMax;
    double timeStart;
    double timeEnd;
    bool fixedTimeStep;
    unsigned int randomParticlesNum;
    std::array<double, 2> randomMass;
    BoundingBox<double, 3> randomPosBox;
    BoundingBox<double, 3> randomVelBox;

    void addBoundary(const std::string& fileName, NBSystem& nbsystem);
    void addParticle(const std::string& fileName, NBSystem& nbsystem);
    void addParticleUniform(const unsigned int &num,
                     const std::array<double, 2>  &mass,
                     const BoundingBox<double, 3> &posBox,
                     const BoundingBox<double, 3> &velBox,
                            NBSystem& nbsystem);

};

#endif // NBSYSTEMSETTINGS_H
