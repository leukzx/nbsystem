#include <iostream>
#include "nbsystem.h"
#include <libconfig.h++>

class Settings
{
public:
    Settings(std::string fileName);

    libconfig::Config cfg;

    std::string outFileName;
    std::string particleFileName;
    std::string boundariesFileName;
    int platformId;
    cl_device_type deviceType;
    double dtCoef;
    int outPrec;
    double writeInterval;
    double timeStepMax;
    double timeEnd;
    unsigned int randomParticlesNum;
    bool fixedTimeStep;
};

Settings::Settings(std::string fileName)
{
    platformId = -1;
    deviceType = CL_DEVICE_TYPE_CPU;
    dtCoef = 0.01;
    outPrec = 6;
    writeInterval = std::numeric_limits<double>::max();
    timeStepMax = std::numeric_limits<double>::max();
    timeEnd = std::numeric_limits<double>::max();
    randomParticlesNum = 0;
    fixedTimeStep = false;

    // Read the file. If there is an error, report it and exit.
    try {
        cfg.readFile(fileName.data());
    }
    catch(const libconfig::FileIOException &fioex) {
        std::cerr << "I/O error while reading file." << std::endl;
        throw;
    }
    catch(const libconfig::ParseException &pex) {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                  << " - " << pex.getError() << std::endl;
        throw;
    }
    try {

        try {
            std::string outFileNameStr = cfg.lookup("Write_to_file");
            outFileName = outFileNameStr;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Write_to_file' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string particleFileNameStr = cfg.lookup("Particles_data_file");
            particleFileName = particleFileNameStr;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Particles_data_file' setting"
                         " in configuration file."
                      << std::endl;
            throw;
        }

        try {
            std::string boundariesFileNameStr
                    = (cfg.lookup("Boundaries_data_file"));
            boundariesFileName = boundariesFileNameStr;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Write_interval' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            platformId = cfg.lookup("OpenCL_platform_id");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'OpenCL_platform_id' setting."
                         " Using first available."
                      << std::endl;
        }

        try {
            std::string type = cfg.lookup("OpenCL_device_type");
            deviceType = type.compare("gpu")?
                                CL_DEVICE_TYPE_CPU:CL_DEVICE_TYPE_GPU;
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'OpenCL_device_type' setting. Using CPU."
                      << std::endl;
        }

        try {
            dtCoef = cfg.lookup("Time_step_estimation_coeff");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Time_step_estimation_coeff' setting"
                         " in configuration file."
                      << std::endl;
        }

        try {
            outPrec = cfg.lookup("Write_precision");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Write_precision' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            writeInterval =  cfg.lookup("Write_interval");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Write_interval' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            timeStepMax = cfg.lookup("Max_time_step");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Max_time_step' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            timeEnd = cfg.lookup("End_time");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'End_time' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            randomParticlesNum = cfg.lookup("Random_particles_number");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Random_particles_number'"
                         " setting in configuration file."
                      << std::endl;
        }

        try {
            fixedTimeStep = cfg.lookup("Fixed_time_step");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Fixed_time_step'"
                         " setting in configuration file."
                      << std::endl;
        }

    }
    catch (const libconfig::SettingTypeException &stex) {
        std::cout << "Something wrong in the settings file with '"
                  << stex.getPath() << "'." << std::endl;
        throw;
    }
}

void stateToFile(std::ofstream& outFileStream, NBSystem& nbs)
{
    if (outFileStream.is_open()){
        outFileStream << nbs.stateString();
//        std::string commentSymb = "#";
//        osFile << commentSymb << "Time = " << time; // Time
//        osFile << "\t Time step = " << timeStep << std::endl;
//        osFile << commentSymb << nbs.energyString() << std::endl; // Energy
//        // All particles state
//        unsigned int pnum = nbs.getPNum();
//        for (int i = 0; i < pnum; ++i) {
//            out << nbs.stateString(i);
//            out << std::endl;
//        }
        outFileStream << std::endl << std::endl; // End of the data block
    } else {
        throw std::runtime_error("Error at stateToFile()."
                                      " File is not opened.");
    }

}

int main(int argc, char *argv[])
{
    if( argc < 1) {
      std::cerr
        << "Use: nbsystem <settings file>"
        << std::endl;
      exit(EXIT_FAILURE);
    }

    Settings set(argv[1]);
    NBSystem nbs;

    nbs.setupCL(set.platformId, set.deviceType);
    nbs.addParticle(set.particleFileName);
    nbs.addBoundary(set.boundariesFileName);
    nbs.setDtCoef(set.dtCoef);
    nbs.setOutPrec(set.outPrec);

    double time = 0.0;
    double timeStep = set.timeStepMax;
    double timeEnd = set.timeEnd;
    double writeInterval = set.writeInterval;
    double writeTime = 0.0; // Write time counter
    double eps = std::numeric_limits<double>::epsilon() * timeStep;
    unsigned int nSteps;
    double eInit = nbs.getEnIn();
    double eCurrent = nbs.getEnTot();
    std::ofstream outFileStream;
    // Output, append
    outFileStream.open(set.outFileName, std::ofstream::out | std::ofstream::app);
    stateToFile(outFileStream, nbs);

    do {
        timeStep = nbs.getEstDt();

        // Check timeStep to not overshoot write time
        if ((writeTime + timeStep) > writeInterval)
            timeStep = writeInterval - writeTime;
        // Check dt to not overshoot end time
        if ((time + timeStep) > set.timeEnd)
            timeStep = timeEnd - time;


        // Move particles
        nbs.evolve(timeStep);

        // Increase time counter
        time += timeStep;

        // Increase number of steps
        nSteps++;

        // Write data to file
        writeTime += timeStep;
        if (std::fabs(writeInterval - writeTime) < eps
                || std::fabs(time - timeEnd) < eps) {
            stateToFile(outFileStream, nbs);
            writeTime = 0;

            // Print current progress
            std::cout << std::scientific;
            std::cout << std::setprecision(6);
            std::cout << std::right;
            std::cout << "Time: " << time << "  "
                      << "Time Step: " << timeStep
                      << "\t";

            if (eInit) {
                std::cout << "(E_tot-E_init)/E_init=";
                std::cout << (nbs.getEnTot() - nbs.getEnIn())/nbs.getEnIn();
            } else {
                std::cout << "E_init=0\t(E_tot-E_init)=";
                std::cout << (nbs.getEnTot() - nbs.getEnIn());
            }
            std::cout << std::endl;
        }

    } while (std::fabs(timeEnd - time) > eps);

    outFileStream.close();

    std::cout << "Time steps: " << nSteps << std::endl;

    return 0;
}
