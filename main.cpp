#include <iostream>
#include "nbsystem.h"
#include <libconfig.h++>

class NBSystemSettings
{
public:
    NBSystemSettings(std::string fileName);

    libconfig::Config cfg;

    std::string outFileName;
    std::string particleFileName;
    std::string boundariesFileName;
    int platformId;
    cl_device_type deviceType;
    double dtCoef;
    double ljEps;
    double ljRmin;
    int outPrec;
    double writeInterval;
    double timeStepMax;
    double timeEnd;
    bool fixedTimeStep;
    unsigned int randomParticlesNum;
    std::array<double, 2> randomMass;
    BoundingBox<double, 3> randomPosBox;
    BoundingBox<double, 3> randomVelBox;
};

NBSystemSettings::NBSystemSettings(std::string fileName)
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
            outFileName = cfg.lookup("Write_to_file").c_str();
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Write_to_file' setting in configuration file."
                      << std::endl;
            throw;
        }

        try {
            particleFileName = cfg.lookup("Particles_data_file").c_str();
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Particles_data_file' setting"
                         " in configuration file."
                      << std::endl;
            throw;
        }

        try {
            boundariesFileName = cfg.lookup("Boundaries_data_file").c_str();
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
            fixedTimeStep = cfg.lookup("Fixed_time_step");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Fixed_time_step'"
                         " setting in configuration file."
                      << std::endl;
        }

        try {
            randomParticlesNum = cfg.lookup("Random_particles.number");
            libconfig::Setting &min_position
                    = cfg.lookup("Random_particles.min_position");
            libconfig::Setting &max_position
                    = cfg.lookup("Random_particles.max_position");
            libconfig::Setting &min_velocity
                    = cfg.lookup("Random_particles.min_velocity");
            libconfig::Setting &max_velocity
                    = cfg.lookup("Random_particles.max_velocity");
            std::array<double, 3> randomPosMin;
            std::array<double, 3> randomPosMax;
            std::array<double, 3> randomVelMin;
            std::array<double, 3> randomVelMax;
            for (unsigned int i = 0; i < 3; i++) {
                randomPosMin.at(i) = min_position[i];
                randomPosMax.at(i) = max_position[i];
                randomVelMin.at(i) = min_velocity[i];
                randomVelMax.at(i) = max_velocity[i];
            }
            randomPosBox.setMinVertex(randomPosMin);
            randomPosBox.setMaxVertex(randomPosMax);
            randomVelBox.setMinVertex(randomVelMin);
            randomVelBox.setMaxVertex(randomVelMax);

            libconfig::Setting &mass
                    = cfg.lookup("Random_particles.mass");
            for (unsigned int i = 0; i < 2; i++)
                randomMass.at(i) = mass[i];
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Random_particles' setting in configuration file.";
        }

    }
    catch (const libconfig::SettingTypeException &stex) {
        std::cout << "Something wrong in the settings file with '"
                  << stex.getPath() << "'." << std::endl;
        throw;
    }
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
    outStr << "\t Time step = " << timeStep << std::endl;
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

void addParticle(const std::string& fileName, NBSystem& nbsystem)
{
    std::ifstream file;
    file.open(fileName);

    if (file.is_open()) {
        std::istream::pos_type blockPos = file.end;
        std::string blockMark = "#Time";
        std::stringstream sLine;
        std::string line;
        std::string word;

        while (std::getline(file, line)) {
            sLine = std::stringstream(line);
            sLine >> word;
            if (word == blockMark) blockPos = file.tellg();
        }

        if (blockPos == file.end)
            throw std::runtime_error("Can't find data "
                                     "block in particles file.");

        file.clear();
        // Setting position to the beginning of last time block
        file.seekg(blockPos);
        // Skipping lines with system's energy data
        std::getline(file, line);

        std::vector<cl_float4> pos;
        std::vector<cl_float4> vel;

        // Read positions and velocities until the end of file
        while (std::getline(file, line)) {
            if (line == "") break;
            sLine = std::stringstream(line);

            // Skip id of particle
            sLine >> word;

            cl_float4 p;
            for (int i = 0; i < 3; i++) {
                sLine >> word;
                p.s[i] = std::stof(word);
            }

            cl_float4 v;
            for (int i = 0; i < 3; i++) {
                sLine >> word;
                v.s[i] = std::stof(word);
            }
            vel.emplace_back(v);

            sLine >> word;
            p.s3 = std::stof(word);
            pos.emplace_back(p);
        }

        nbsystem.addParticle(pos, vel);
        file.close();
    } else throw std::runtime_error("Can't find particles data file.");
}

void addBoundary(const std::string& fileName, NBSystem& nbsystem)
{
    // This string should preceed the sequence of boundary points
    std::string blockMark = "#Boundary";

    std::ifstream file;
    file.open(fileName);
    if (file.is_open()){
        std::stringstream sLine;
        std::string line;
        std::string word;

        while (std::getline(file, line)) {
            sLine = std::stringstream(line);
            sLine >> word;
            if (word == blockMark) {
                std::vector<cl_float4> boundary;
                while (std::getline(file, line)){
                    if (line != "") {
                        cl_float4 vertice;
                        vertice.s3 = 0.0f;
                        sLine = std::stringstream(line);
                        sLine >> word;
                        vertice.s0 = std::stod(word);
                        sLine >> word;
                        vertice.s1 = std::stod(word);
                        sLine >> word;
                        vertice.s2 = std::stod(word);
                        boundary.emplace_back(vertice);
                    } else break;
                }

                nbsystem.addBoundary(boundary);
            }
        }
        file.close();
    } else {
        throw std::runtime_error("Can't find boundaries data file.");
    }
}



int main(int argc, char *argv[])
{
    if( argc < 1) {
      std::cerr << "Use: nbsystem <settings file>" << std::endl;
      exit(EXIT_FAILURE);
    }

    NBSystemSettings set(argv[1]);
    NBSystem nbs;

    std::cout << nbs.setupCL(set.platformId, set.deviceType) << std::endl;
    addParticle(set.particleFileName, nbs);
    addBoundary(set.boundariesFileName, nbs);
    nbs.setDtCoef(set.dtCoef);
    nbs.addParticle(set.randomParticlesNum, set.randomMass,
                    set.randomPosBox, set.randomVelBox);

    double time = 0.0;
    double timeStep = set.timeStepMax;
    double timeEnd = set.timeEnd;
    double writeInterval = set.writeInterval;
    double writeTime = 0.0; // Write time counter
    double eps = std::numeric_limits<double>::epsilon() * timeStep;
    unsigned int nSteps;
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

    do {
        timeStep = nbs.getEstDt();

        // Check timeStep to not overshoot write time
        if ((writeTime + timeStep) > writeInterval)
            timeStep = writeInterval - writeTime;
        // Check timeStep to not overshoot end time
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
            //stateToFile(outFileStream, nbs);
            enKin = nbs.getEnKin();
            enPot = nbs.getEnPot();
            enTot = enKin + enPot;
            stateToFile(outFileStream, nbs, time, timeStep,
                        enKin, enPot, enInit, prec);
            writeTime = 0;

            // Print current progress to console
            std::cout << "Time: " << time << "  "
                      << "Time Step: " << timeStep
                      << "  ";
            if (enInit) {
                std::cout << "(E_tot-E_init)/E_init=";
                std::cout << (enTot - enInit)/enInit;
            } else {
                std::cout << "E_init=0\t(E_tot-E_init)=";
                std::cout << (enTot - enInit);
            }
            std::cout << std::endl;
        }

    } while (std::fabs(timeEnd - time) > eps);

    outFileStream.close();

    std::cout << "Time steps: " << nSteps << std::endl;

    return 0;
}
