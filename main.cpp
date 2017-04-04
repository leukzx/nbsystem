#include <iostream>
#include "nbsystem.h"
#include "nbsystemsettings.h"

std::string energyString(double enKin,
                         double enPot,
                         double enInit,
                         unsigned int prec);

std::string particleStateString(NBSystem &nbsystem,
                                const unsigned int& pId,
                                const unsigned int& prec);

std::string systemStateString(NBSystem &nbsystem,
                              const double& time,
                              const double& timeStep,
                              const double& enKin,
                              const double& enPot,
                              const double& enInit,
                              const unsigned int& prec);

void stateToFile(std::ofstream& outFileStream,
                 NBSystem& nbs,
                 const double& time,
                 const double& timeStep,
                 const double& enKin,
                 const double& enPot,
                 const double& enInit,
                 const unsigned int& prec);

std::string curentProgressString(const double& enTot,
                                 const double& timeStepAvg,
                                 const double& time,
                                 const double& enInit,
                                 const unsigned int& prec);


int main(int argc, char *argv[])
{
    if( argc < 1) {
      std::cerr << "Use: nbsystem <settings file>" << std::endl;
      exit(EXIT_FAILURE);
    }

    NBSystemSettings set(argv[1]);
    NBSystem nbs;

    std::cout << nbs.setupCL(set.platformId, set.deviceType) << std::endl;
    set.addParticle(set.particleFileName, nbs);
    set.addBoundary(set.boundariesFileName, nbs);
    nbs.setDtCoef(set.dtCoef);

    // Add random particles
    if (set.randomParticlesNum)
        nbs.addParticle(set.randomParticlesNum, set.randomMass,
                        set.randomPosBox, set.randomVelBox);

//    set.addParticleUniform(set.randomParticlesNum, set.randomMass,
//                    set.randomPosBox, set.randomVelBox, nbs);

    double time = set.timeStart;
    double timeStepMax = (set.timeStepMax < set.writeInterval) ?
                set.timeStepMax : set.writeInterval;
    double timeStepEst = timeStepMax;
    double timeStep = timeStepMax;

    double timeStepAvg = 0;
    double timeEnd = set.timeEnd;
    double writeInterval = set.writeInterval;
    double writeTime = 0.0; // Write time counter
    double eps = std::numeric_limits<double>::epsilon() * timeStep;
    unsigned int nStepsTotal = 0;
    unsigned int nStepsAvg = 0;
    double enKin = nbs.getEnKin();
    double enPot = nbs.getEnPot();
    double enTot = enKin + enPot;
    double enInit = enTot;
    unsigned int prec = set.outPrec;
    std::ofstream outFileStream;
    // Output, append
    outFileStream.open(set.outFileName, std::ofstream::out | std::ofstream::app);
    stateToFile(outFileStream, nbs, time, timeStep,
                enKin, enPot, enInit, prec);

    // Console output settings
    std::cout << std::scientific;
    std::cout << std::setprecision(prec);
    std::cout << std::right;

    // Two branches
    // First - time step estimates for every move of the system
    // Second - time step is fixed
    // It was choosen to not to check if time step is fixed in every iteration
    if (!set.fixedTimeStep)
        do {
            // Estimate time step
            timeStepEst = nbs.getEstDt();
            timeStep = (timeStepEst < timeStepMax) ? timeStepEst : timeStepMax;

            // Check timeStep to not overshoot write time
            if ((writeTime + timeStep) > writeInterval)
                timeStep = writeInterval - writeTime;

            // Check timeStep to not overshoot end time
            if ((time + timeStep) > timeEnd)
                timeStep = timeEnd - time;

            // Move particles
            nbs.evolve(timeStep);

            // Increase time counter
            time += timeStep;
            writeTime += timeStep;

            // Increase number of steps
            nStepsTotal++;
            nStepsAvg++;

            // Write data to file
            if (std::fabs(writeInterval - writeTime) < eps) {
                //|| std::fabs(time - timeEnd) < eps) {
                timeStepAvg = writeTime / nStepsAvg;
                nStepsAvg = 0;
                writeTime = 0;
                enKin = nbs.getEnKin();
                enPot = nbs.getEnPot();
                enTot = enKin + enPot;
                stateToFile(outFileStream, nbs, time, timeStepAvg,
                            enKin, enPot, enInit, prec);

                // Print current progress to console
                std::cout << curentProgressString(enTot, timeStepAvg, time,
                                                  enInit, prec)
                          << std::endl;
            }
        } while (std::fabs(timeEnd - time) > eps);
     else
        do {
            // Check timeStep to not overshoot write time
            if ((writeTime + timeStep) > writeInterval)
                timeStep = writeInterval - writeTime;

            // Check timeStep to not overshoot end time
            if ((time + timeStep) > timeEnd)
                timeStep = timeEnd - time;

            // Move particles
            nbs.evolve(timeStep);

            // Increase time counter
            time += timeStep;
            writeTime += timeStep;

            // Increase number of steps
            nStepsTotal++;
            nStepsAvg++;

            // Write data to file
            if (std::fabs(writeInterval - writeTime) < eps) {
                //|| std::fabs(time - timeEnd) < eps) {
                timeStepAvg = writeTime / nStepsAvg;

                enKin = nbs.getEnKin();
                enPot = nbs.getEnPot();
                enTot = enKin + enPot;

                stateToFile(outFileStream, nbs, time, timeStepAvg,
                            enKin, enPot, enInit, prec);

                nStepsAvg = 0;
                writeTime = 0;
                timeStep = set.timeStepMax;

                // Print current progress to console
                std::cout << curentProgressString(enTot, timeStepAvg, time,
                                                  enInit, prec)
                          << std::endl;
            }
        } while (std::fabs(timeEnd - time) > eps);


    // Write system's final state to file if it is not written already
    if (writeTime) {
        timeStepAvg = writeTime / nStepsAvg;
        nStepsAvg = 0;
        writeTime = 0;
        enKin = nbs.getEnKin();
        enPot = nbs.getEnPot();
        enTot = enKin + enPot;
        stateToFile(outFileStream, nbs, time, timeStepAvg,
                    enKin, enPot, enInit, prec);

        // Print current progress to console
        std::cout << curentProgressString(enTot, timeStepAvg, time,
                                          enInit, prec)
                  << std::endl;
    }

    outFileStream.close();

    std::cout << "Number of time steps: " << nStepsTotal << std::endl;
    std::cout << "Average time step: " << (time - set.timeStart) / nStepsTotal
              << std::endl;

    return 0;
}


std::string energyString(double enKin,
                         double enPot,
                         double enInit,
                         unsigned int prec)
{
    std::ostringstream out;
    out << std::scientific;
    out << std::setprecision(prec);

    std::string delim = " ";

    double enTot = enKin + enPot;
    out << "E_kin=";
    out << enKin;
    out << delim << "E_pot=";
    out << enPot;
    out << delim << "E_tot=";
    out << enTot;
    out << delim;
    if (enInit) {
        out << "(E_tot-E_init)/E_init=";
        out << (enTot - enInit)/enInit;
    } else {
        out << "E_init=0\t(E_tot-E_init)=";
        out << enTot - enInit;
    }

    return out.str();
}

std::string particleStateString(NBSystem &nbsystem,
                                const unsigned int& pId,
                                const unsigned int& prec)
{
    int cWidth = 8; // Minimum column width is 8
    std::string delim = " "; // data delimeter

    cWidth += prec;

    std::ostringstream out;

    out << std::scientific;
    out << std::setprecision(prec);
    std::cout << std::right;

    out << pId;
    for (unsigned int i = 0; i < 3; i++) {
        out << delim;
        out << std::setw(cWidth) << nbsystem.posData()->at(pId).s[i];
    }
    for (unsigned int i = 0; i < 3; i++) {
        out << delim;
        out << std::setw(cWidth) << nbsystem.velData()->at(pId).s[i];
    }
    out << delim;
    out << std::setw(cWidth) << nbsystem.posData()->at(pId).s[3];

    return out.str();
}

std::string systemStateString(NBSystem &nbsystem,
                        const double& time,
                        const double& timeStep,
                        const double& enKin,
                        const double& enPot,
                        const double& enInit,
                        const unsigned int& prec)
{
    std::ostringstream outStr;

    std::string cSymb = "#";

    outStr << cSymb << "Time = " << time; // Time
    outStr << "\t Avg time step = " << timeStep << std::endl;
    outStr << cSymb;
    outStr << energyString(enKin, enPot, enInit, prec) << std::endl;
    // All particles state
    unsigned int pnum = nbsystem.getPNum();
    for (int i = 0; i < pnum; ++i) {
        outStr << particleStateString(nbsystem, i, prec);
        outStr << std::endl;
    }

    return outStr.str().erase(outStr.str().size()-1);
}

void stateToFile(std::ofstream& outFileStream,
                                   NBSystem& nbs,
                                   const double& time,
                                   const double& timeStep,
                                   const double& enKin,
                                   const double& enPot,
                                   const double& enInit,
                                   const unsigned int& prec
                                   )
{
    if (outFileStream.is_open()) {
        outFileStream << systemStateString(nbs,
                                           time,
                                           timeStep,
                                           enKin,
                                           enPot,
                                           enInit,
                                           prec);
        // End of the data block
        outFileStream << std::endl << std::endl << std::endl;
    } else {
        throw std::runtime_error("Error at stateToFile()."
                                      " File is not opened.");
    }
}

std::string curentProgressString(const double& enTot,
                                 const double& timeStepAvg,
                                 const double& time,
                                 const double& enInit,
                                 const unsigned int& prec)
{
    // Current progress string
    std::ostringstream out;

    out << std::scientific;
    out << std::setprecision(prec);

    out << "Time: " << time << "  "
        << "Avg time step: " << timeStepAvg
        << "  ";
    if (enInit) {
        out << "(E_tot-E_init)/E_init=";
        out << (enTot - enInit)/enInit;
    } else {
        out << "E_init=0\t(E_tot-E_init)=";
        out << (enTot - enInit);
    }
    return out.str();
}
