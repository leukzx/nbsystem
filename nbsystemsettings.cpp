#include "nbsystemsettings.h"

NBSystemSettings::NBSystemSettings(std::string fileName)
{
    platformId = -1;
    deviceId = 0;
    deviceType = CL_DEVICE_TYPE_CPU;
    dtCoef = 0.01;
    outPrec = 6;
    writeInterval = std::numeric_limits<double>::max();
    timeStepMax = std::numeric_limits<double>::max();
    timeStart = 0.0;
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
            std::cerr << "No 'Boundaries_data_file' setting"
                         " in configuration file."
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
            deviceId = cfg.lookup("OpenCL_device_id");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'OpenCL_device_id' setting."
                         " Using first available."
                      << std::endl;
        }

        try {
            std::string type = cfg.lookup("OpenCL_device_type").c_str();
            if (!type.compare("gpu")) deviceType = CL_DEVICE_TYPE_GPU;
            else if (!type.compare("cpu")) deviceType = CL_DEVICE_TYPE_CPU;
            else if (!type.compare("all")) deviceType = CL_DEVICE_TYPE_ALL;
            else {
                deviceType = CL_DEVICE_TYPE_CPU;
                throw;
            }
//            deviceType = type.compare("gpu")?
//                                CL_DEVICE_TYPE_CPU:CL_DEVICE_TYPE_GPU;
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
            timeStart = cfg.lookup("Start_time");
        }
        catch(const libconfig::SettingNotFoundException &nfex) {
            std::cerr << "No 'Start_time' setting in configuration file."
                      << std::endl;
        }

        try {
            timeEnd = cfg.lookup("End_time");
            if (timeEnd <= timeStart)
                throw std::runtime_error("Error while reading configuration "
                                         "file: 'End_time' <= 'Start_time.'");
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

void NBSystemSettings::addParticle(const std::string& fileName,
                                   NBSystem& nbsystem)
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
            v.s[3] = 0.0f;
            vel.emplace_back(v);

            sLine >> word;
            p.s[3] = std::stof(word);
            pos.emplace_back(p);
        }

        nbsystem.addParticle(pos, vel);
        file.close();
    } else throw std::runtime_error("Can't find particles data file.");
}

void NBSystemSettings::addBoundary(const std::string& fileName,
                                   NBSystem& nbsystem)
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
                        vertice.s[3] = 0.0f;
                        sLine = std::stringstream(line);
                        sLine >> word;
                        vertice.s[0] = std::stod(word);
                        sLine >> word;
                        vertice.s[1] = std::stod(word);
                        sLine >> word;
                        vertice.s[2] = std::stod(word);
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

void NBSystemSettings::addParticleUniform(const unsigned int &num,
                                          const std::array<double, 2>  &mass,
                                          const BoundingBox<double, 3> &posBox,
                                          const BoundingBox<double, 3> &velBox,
                                          NBSystem &nbsystem)
{
    // This function generates particles in the given area
    // with the given velocity and mass restrictions
    // and even spacing between them
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    std::vector<cl_float4> positions;
    std::vector<cl_float4> velocities;
    double sideN = 10;
    std::cout << sideN;
    std::array<double, 3> dist;

    for (unsigned int i = 0; i < 3; i++)
        dist.at(i) = (posBox.maxVertex()[i] - posBox.minVertex()[i]) / sideN;


    for (unsigned int i = 0; i < sideN; i++) {
        cl_float4 p; // Particle's position

        p.s[0] = posBox.minVertex()[0] + i * dist.at(0);
        for (unsigned int j = 0; j < sideN; j++) {

            p.s[1] = posBox.minVertex()[1] + j * dist.at(1);

            for (unsigned int k = 0; k < sideN; k++) {

                p.s[2] = posBox.minVertex()[2] + k * dist.at(2);

                cl_float4 v; // Particle's velocity
                for (unsigned int i = 0; i < 3; i++)
                    v.s[i] = 0.0;
//                    v.s[i] = velBox.minVertex()[i]
//                            + (velBox.maxVertex()[i] - velBox.minVertex()[i])
//                              * distribution(generator);
                p.s[3] = mass.at(0)
                        + (mass.at(1) - mass.at(0)) * distribution(generator);
                v.s[3] = 0;

                positions.emplace_back(p);
                velocities.emplace_back(v);
            }
        }
    }

    unsigned int remainderN = num - std::pow(sideN, 3);


    nbsystem.addParticle(positions, velocities);
}
