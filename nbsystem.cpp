#include "nbsystem.h"

NBSystem::NBSystem()
{
    pnum = 0;
    vnum = 0;
    enIn = 0;
    enKin = 0;
    enPot = 0;
    timeCur = 0;
    timeEnd = 0;
    dt = 0;
    tsNum = 0;
    outPrec = 6;
    dtCoef = 0.01;
    syncFlag = false;

    curBuffIndex = 0;

    // Default setting is to select first available cpu device
    setupCL(-1, CL_DEVICE_TYPE_CPU);

    h_enKin = new cl_float[pnum];
    h_enPot = new cl_float[pnum];
    h_esTS = new cl_float[pnum * pnum];
}

NBSystem::~NBSystem()
{
    delete[] h_enKin;
    delete[] h_enPot;
    delete[] h_esTS;
}

double NBSystem::getEnIn()
{
    return enIn;
}

double NBSystem::getEnPot()
{
    return getEnPot(event);
}

double NBSystem::getEnKin()
{
    return getEnKin(event);
}

double NBSystem::getEnTot()
{
    return getEnPot(event) + getEnKin(event);
}

void NBSystem::addParticle(cl_float4 &p, cl_float4 &v)
{
    pos.emplace_back(p);
    vel.emplace_back(v);
    pnum++;
}

void NBSystem::addParticle(std::vector<cl_float4> &pos,
                 std::vector<cl_float4> &vel)
{
    if (pos.size() != vel.size())
        throw std::runtime_error("Sizes of particle state vectors"
                                 " (mass, pos, vel) are not equal.");
    if (timeCur != 0) {
        queue.enqueueReadBuffer(d_pos[curBuffIndex],
                                CL_TRUE, 0, vectBuffSize, pos.data());
        queue.enqueueReadBuffer(d_vel,
                                CL_TRUE, 0, vectBuffSize, vel.data());
    }

    unsigned int num = pos.size();
    for (unsigned int i = 0; i < num; i++ )
        addParticle(pos.at(i), vel.at(i));

    // Update buffers
    initializeCLBuffers();
    setupCLKernelsArgs();

    // Update accelerations
    updateAcc(event);
    event.wait();

    // Update initial energy
    enIn = getEnTot();
}

void NBSystem::removeParticle(unsigned int i)
{
    if (i >= pnum) throw std::domain_error("No particle with given id.");

    if (timeCur != 0) {
        queue.enqueueReadBuffer(d_pos[curBuffIndex],
                                CL_TRUE, 0, vectBuffSize, pos.data());
        queue.enqueueReadBuffer(d_vel,
                                CL_TRUE, 0, vectBuffSize, vel.data());
    }

    pos.erase(pos.begin() + i);
    vel.erase(vel.begin() + i);
    pnum--;

    // Update buffers
    initializeCLBuffers();
    setupCLKernelsArgs();

    // Update accelerations
    updateAcc(event);
    event.wait();

    // Update initial energy
    enIn = getEnTot();
}

void NBSystem::addParticle(std::string fileName)
{
    std::ifstream file;

    file.open(fileName);
    if (file.is_open()){
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
        addParticle(pos, vel);
        file.close();
    } else throw std::runtime_error("Can't find particles data file.");

}

void NBSystem::addBoundary(std::vector<cl_float4> &vertices)
{
    if (vertices.size() < 3)
        throw std::runtime_error("Invalid boundary:"
                                 " number of vertices is less than 3.");
    // Should add check if vertices are coplanar
    if (vertices.size() > 3) {
        cl_float4 a;
        a.s0 = vertices.at(1).s0 - vertices.at(0).s0;
        a.s1 = vertices.at(1).s1 - vertices.at(0).s1;
        a.s2 = vertices.at(1).s2 - vertices.at(0).s2;
        a.s3 = 0.0f;
        cl_float4 b;
        b.s0 = vertices.at(2).s0 - vertices.at(0).s0;
        b.s1 = vertices.at(2).s1 - vertices.at(0).s1;
        b.s2 = vertices.at(2).s2 - vertices.at(0).s2;
        b.s3 = 0.0f;
        cl_float4 axb;
        axb.s0 = a.s1 * b.s2 - b.s1 * a.s2;
        axb.s1 = -(a.s0 * b.s2 - b.s0 * a.s2);
        axb.s2 = a.s0 * b.s1 - b.s0 * a.s1;
        axb.s3 = 0.0f;

        // Some number to compare result to zero
        float eps = 1e-5;

        for (unsigned int i = 3; i < vertices.size(); i++){
            cl_float4 c;
            c.s0 = vertices.at(i).s0 - vertices.at(0).s0;
            c.s1 = vertices.at(i).s1 - vertices.at(0).s1;
            c.s2 = vertices.at(i).s2 - vertices.at(0).s2;
            c.s3 = 0.0f;

            cl_float abc = axb.s0 * c.s0 + axb.s1 * c.s1 + axb.s2 * c.s2;
            if (abc > eps)
                throw std::runtime_error("Invalid boundary:"
                                         "vertices are not coplanar.");
        }
    }

    for (unsigned int i = 1 ; i < vertices.size() - 1; i++) {
        unsigned int l[] = {0, i, i + 1};
        for (unsigned int k : l) tri.emplace_back(vertices.at(k));
    }

    vnum = tri.size();
    initializeCLBuffers();
    setupCLKernelsArgs();
}

void NBSystem::addBoundary(std::string fileName)
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
                addBoundary(boundary);
            }
        }
        file.close();
    } else {
        throw std::runtime_error("Can't find boundaries data file.");
    }
}

std::string NBSystem::stateString(int pId)
{
    std::streamsize defaultPrecision = std::cout.precision();
    int cWidth = 8; // Minimum column width is 8
    std::string delim = " "; // data delimeter

    cWidth += outPrec;

    std::ostringstream out;

    out << std::scientific;
    out << std::setprecision(outPrec);
    std::cout << std::right;

    out << pId;
    out << delim;
    out << std::setw(cWidth) << pos[pId].s[0];
    out << delim;
    out << std::setw(cWidth) << pos[pId].s[1];
    out << delim;
    out << std::setw(cWidth) << pos[pId].s[2];
    out << delim;
    out << std::setw(cWidth) << vel[pId].s[0];
    out << delim;
    out << std::setw(cWidth) << vel[pId].s[1];
    out << delim;
    out << std::setw(cWidth) << vel[pId].s[2];
    out << delim;
    out << pos[pId].s[3];

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);
    return out.str();
}

std::string NBSystem::energyString()
{
    std::string delim = " ";
    std::ostringstream out;
    std::streamsize defaultPrecision = std::cout.precision();
    out << std::scientific;
    out << std::setprecision(outPrec);

    out << "E_kin=";
    out << getEnKin();
    out << delim << "E_pot=";
    out << getEnPot();
    out << delim << "E_tot=";
    out << getEnTot();
    out << delim;
    if (enIn) {
        out << "(E_tot-E_init)/E_init=";
        out << (getEnTot() - getEnIn())/getEnIn();
    } else {
        out << "E_init=0\t(E_tot-E_init)=";
        out << (getEnTot() - getEnIn());
    }

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);

    return out.str();
}

std::string NBSystem::stateString()
{
    syncArrays();

    std::string cSymb = "#";
    std::streamsize defaultPrecision = std::cout.precision();

    std::ostringstream out;

    out << cSymb << "Time = " << timeCur; // Time
    out << "\t Time step = " << dt << std::endl;
    out << cSymb << energyString() << std::endl; // Energy
    // All particles state
    for (int i = 0; i < pnum; ++i) {
        out << stateString(i);
        out << std::endl;
    }

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);

    return out.str().erase(out.str().size()-1);
}

std::string NBSystem::boundariesString()
{
    std::ostringstream out;
    std::streamsize defaultPrecision = std::cout.precision();
    out << std::scientific;
    out << std::setprecision(outPrec);

    std::string delim = " ";
    int cWidth = 8; // Minimum column width is 8
    cWidth += outPrec;

    for (unsigned int i = 0; i < tri.size(); i = i + 3) {
        out << i / 3;
        out << std::endl;
        for (unsigned int j = 0; j < 3; j++) {
            out << delim << std::setw(cWidth) << tri.at(i+j).s0;
            out << delim << std::setw(cWidth) << tri.at(i+j).s1;
            out << delim << std::setw(cWidth) << tri.at(i+j).s2;
            out << std::endl;
        }
    }

    std::cout.unsetf(std::ios::fixed | std::ios::scientific);
    std::cout << std::setprecision(defaultPrecision);

    return out.str().erase(out.str().size()-1);
}

void NBSystem::setupCL(int platformId, cl_device_type deviceType)
{
    try {
        // List available platforms and devices
        std::vector<cl::Platform> platforms;
        std::cout << "Available platforms:" << std::endl;
        cl::Platform::get(&platforms);
        for (auto& platform : platforms ) {
            std::string platformName;
            platform.getInfo((cl_platform_info)CL_PLATFORM_NAME,
                             &platformName);
            std::cout << " " << platformName << std::endl;
            try {
                std::vector<cl::Device> deviceList;
                platform.getDevices(CL_DEVICE_TYPE_ALL, &deviceList);
                for(auto& deviceItem : deviceList) {
                    std::string deviceName;
                    deviceItem.getInfo((cl_device_info)CL_DEVICE_NAME,
                                       &deviceName);
                    std::cout << "  " << deviceName << std::endl;
                }
            } catch (const cl::Error &error) {
                std::cerr << "  Error while getting device list: "
                          << error.what()
                          << '(' << error.err() << ')' << std::endl;
                //throw;
            }
        }

        //The constructor attempts to use the first platform that
        //has a device of the specified type.
        //cl::Context context(deviceType, contextProps);
        //cl::Context context(device);
        try {
            if (platformId < 0) context = cl::Context(deviceType);
            else {
                cl_context_properties cps[] =
                {CL_CONTEXT_PLATFORM,
                 (cl_context_properties)(platforms.at(platformId)()), 0};
                context = cl::Context(deviceType, cps);
            }

        } catch (const cl::Error &error) {
            std::cerr << "  Error while setting context: "
                      << error.what()
                      << '(' << error.err() << ')' << std::endl;
            throw;
        }

        // Selecting first available device of the context
        std::vector<cl::Device> devices;
        try {
            devices = context.getInfo<CL_CONTEXT_DEVICES>();
            device = devices.at(0);
            std::string deviceName;
            device.getInfo((cl_device_info)CL_DEVICE_NAME, &deviceName);
            std::cout << "Selected device: " << deviceName << std::endl;
        } catch (const cl::Error &error) {
            std::cerr << "  Error while selecting device: "
                      << error.what()
                      << '(' << error.err() << ')' << std::endl;
            throw;
        }

        // Kernel sources
        std::string kernelSourceFileName = "kernelsSource.cl";
        std::ifstream file(kernelSourceFileName);
        if (!file.is_open())
            throw std::runtime_error("Error reading kernels source file.");

        std::string kernelSourceStr(std::istreambuf_iterator<char>(file),
                                    (std::istreambuf_iterator<char>()));
        file.close();

        // Create vector of kernel sources. ?There are may be more then 1?
        cl::Program::Sources
                kernelSourceVec(1,
                                std::make_pair(kernelSourceStr.data(),
                                               kernelSourceStr.length()+1));
        // Create an OpenCL program object for the context
        // and load the source code specified by the text strings in each
        // element of the vector sources into the program object
        program = cl::Program(context, kernelSourceVec);
        try {
            program.build(devices);
        } catch (cl::Error& error) {
            std::cerr << "Program build from kernel sources failed! "
                      << error.what()
                      << '(' << error.err() << ')' << std::endl;
            std::cerr << "Retrieving  logs ...:" << std::endl;
            for (auto& dvce : devices) {
                std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dvce)
                          << std::endl;
            }
            throw;
        }

        // Create queue with profiling enabled
        // queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE);
        queue = cl::CommandQueue(context, device);

        // Create a kernel objects
        stepInTime = cl::Kernel(program, "tsForwardEulerCL");
        estimateDt = cl::Kernel(program, "estimateDtCL");
        calcAcc = cl::Kernel(program, "calcAccCL");
        checkBoundaries = cl::Kernel(program, "checkBoundariesCL");
        calcEnKin = cl::Kernel(program, "calcEkinCL");
        calcEnPot = cl::Kernel(program, "calcEpotLJCL");

        // Create event
        event = cl::Event();

    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
    initializeCLBuffers();
    setupCLKernelsArgs();
}

void NBSystem::initializeCLBuffers()
{
    // Create and initialize device buffers
    //cl::Event event;
    try {
        if (pnum) {
            vectBuffSize = pnum * sizeof(cl_float4);
            scalBuffSize = pnum * sizeof(cl_float);

            for (unsigned int i = 0; i < 2; i++) {
                d_pos[i] = cl::Buffer(context, CL_MEM_READ_WRITE, vectBuffSize);
            }
            d_vel = cl::Buffer(context, CL_MEM_READ_WRITE, vectBuffSize);
            d_acc = cl::Buffer(context, CL_MEM_READ_WRITE, vectBuffSize);
            d_enKin = cl::Buffer(context, CL_MEM_READ_WRITE,
                                 scalBuffSize);
            d_enPot = cl::Buffer(context, CL_MEM_READ_WRITE,
                                 scalBuffSize);

            d_est_len = pnum * pnum;
            d_est_size = d_est_len * sizeof(cl_float);
            d_est = cl::Buffer(context, CL_MEM_READ_WRITE,
                               d_est_size);

            queue.enqueueWriteBuffer(d_pos[curBuffIndex], CL_TRUE,
                                     0, vectBuffSize, pos.data(), NULL, &event);
            queue.enqueueWriteBuffer(d_vel, CL_TRUE,
                                     0, vectBuffSize, vel.data(), NULL, &event);
            event.wait();
            // Particles state at device and host are the same now
            syncFlag = true;

            // Create ancillary arrays
            delete[] h_enKin;
            delete[] h_enPot;
            delete[] h_esTS;

            h_enKin = new cl_float[pnum];
            h_enPot = new cl_float[pnum];
            h_esTS = new cl_float[pnum * pnum];
        }

        int vnum = tri.size(); // Number of vertices
        if (vnum) {
            size_t triBufSize = vnum * sizeof(cl_float4);
            d_tri = cl::Buffer(context, CL_MEM_READ_WRITE,
                               triBufSize);

            queue.enqueueWriteBuffer(d_tri, CL_TRUE,
                                     0, triBufSize,
                                     tri.data(), NULL, &event);
            event.wait();
        }
    }  catch (cl::Error error) {
        std::cerr << "Caught exception at initializeCLBuffers(): "
                  << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
}

void NBSystem::setupCLKernelsArgs()
{
    // Commented calls of setArgs are executed before kernels runs
    //int vnum = tri.size(); // Number of vertices
    try{
        // Set stepInTime kernel arguments
//      stepInTime.setArg(0, d_pos);
        stepInTime.setArg(1, d_vel);
        stepInTime.setArg(2, d_acc);
//      stepInTime.setArg(3, (cl_float) dt);
//      stepInTime.setArg(4, d_newPos);

        // Set estimateDt kernel arguments
//      estimateDt.setArg(0, d_pos);
        estimateDt.setArg(1, d_vel);
        estimateDt.setArg(2, d_acc);
        estimateDt.setArg(3, d_est);

        // Set calcAcc kernel arguments
//      calcAcc.setArg(0, d_pos);
        calcAcc.setArg(1, d_acc);

        // Set checkBoundaries kernel arguments
//      checkBoundaries.setArg(0, d_pos);
        checkBoundaries.setArg(1, d_vel);
//      checkBoundaries.setArg(2, d_newPos);
        checkBoundaries.setArg(3, d_tri);
        checkBoundaries.setArg(4, vnum);

        // Set calcEnKin kernel arguments
//      calcEnKin.setArg(0, d_pos);
        calcEnKin.setArg(1, d_vel);
        calcEnKin.setArg(2, d_enKin);

        // Set calcEnPot kernel arguments
//      calcEnPot.setArg(0, d_pos);
        calcEnPot.setArg(1, d_enPot);
    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
}

double NBSystem::getEnKin(cl::Event &event)
{
    // Calculate kinetic energy of the system
    try {
        calcEnKin.setArg(0, d_pos[curBuffIndex]);
        queue.enqueueNDRangeKernel(calcEnKin,
                                   cl::NullRange, cl::NDRange(pnum),
                                   cl::NullRange, NULL, &event);
        event.wait();

        queue.enqueueReadBuffer(d_enKin, CL_TRUE, 0,
                                scalBuffSize, h_enKin);
    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }

    double enKin = 0.0;
    for (size_t i = 0; i < pnum; ++i) {
        enKin += h_enKin[i];
    }

    return enKin;
}

double NBSystem::getEnPot(cl::Event &event)
{
    // Calculate potential energy of the system
    try {
        calcEnPot.setArg(0, d_pos[curBuffIndex]);
        queue.enqueueNDRangeKernel(calcEnPot,
                                   cl::NullRange, cl::NDRange(pnum),
                                   cl::NullRange, NULL, &event);
        event.wait();
        queue.enqueueReadBuffer(d_enPot, CL_TRUE, 0,
                                scalBuffSize, h_enPot);
    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
    throw;
    }

    double enPot = 0.0;
    for (size_t i = 0; i < pnum; ++i) {
        enPot += h_enPot[i];
    }
    enPot /= 2;

    return enPot;
}

double NBSystem::estDt(cl::Event &event)
{
    // Estimate time step for the next iteration
    try {
        estimateDt.setArg(0, d_pos[curBuffIndex]);

        queue.enqueueNDRangeKernel(estimateDt,
                                   cl::NullRange, cl::NDRange(pnum, pnum),
                                   cl::NullRange, NULL, &event);
        event.wait();

        queue.enqueueReadBuffer(d_est, CL_TRUE, 0,
                                d_est_size, h_esTS);

    } catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }

    return dtCoef * (*std::min_element(h_esTS, h_esTS + d_est_len));
}

void NBSystem::updateAcc(cl::Event &event)
{
    // Udpate accelerations for current state
    try {
        calcAcc.setArg(0, d_pos[curBuffIndex]);
        queue.enqueueNDRangeKernel(calcAcc,
                               cl::NullRange, cl::NDRange(pnum),
                               cl::NullRange, NULL, &event);
    }  catch (cl::Error error) {
        std::cerr << "Caught exception: " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
}

void NBSystem::updatePos(double dt, cl::Event &event)
{
    int currentIndex = curBuffIndex;
    int nextIndex = (curBuffIndex+1)%2;

    try {
        stepInTime.setArg(0, d_pos[currentIndex]);
        stepInTime.setArg(3, (cl_float) dt);
        stepInTime.setArg(4, d_pos[nextIndex]);

        queue.enqueueNDRangeKernel(stepInTime,
                                   cl::NullRange, cl::NDRange(pnum),
                                   cl::NullRange, NULL, &event);
        event.wait();

        // Check and resolve crossing of a boundaries
        checkBoundaries.setArg(0, d_pos[currentIndex]);
        checkBoundaries.setArg(2, d_pos[nextIndex]);

        queue.enqueueNDRangeKernel(checkBoundaries,
                                   cl::NullRange, cl::NDRange(pnum),
                                   cl::NullRange, NULL, &event);
        event.wait();

    }  catch (cl::Error error) {
        std::cerr << "Caught exception at updatePos(): " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
    curBuffIndex = nextIndex;
    // Particles state at device and host are not the same now
    syncFlag = false;
}

void NBSystem::setDtCoef(double coeff)
{
    if (coeff <=0)
        throw std::invalid_argument("Invalid argument to setDtCoef()");
    dtCoef = coeff;
}

double NBSystem::getEstDt()
{
    return estDt(event);
}


void NBSystem::evolve(double timeStep)
{
    if (timeStep <=0)
        throw std::invalid_argument("Invalid argument to evolve()");
    updatePos(timeStep, event);
    dt = timeStep;
    timeCur += timeStep;
    tsNum++;
    updateAcc(event);
}

void NBSystem::evolveIn(double interval)
{
    double time = 0;
    double eps = std::numeric_limits<double>::epsilon() * interval;

    do {
        dt = estDt(event);
        if ((time + dt) > interval) dt = interval - time;
        evolve(dt);
        time += dt;
    } while (std::fabs(interval - time) > eps);
    std::cout << tsNum << std::endl;
}

cl_float4 NBSystem::getVel(unsigned int i)
{
    syncArrays();
    return vel.at(i);
}

cl_float4 NBSystem::getPos(unsigned int i)
{
    syncArrays();
    return pos.at(i);
}

double NBSystem::getTimeCur()
{
    return timeCur;
}

double NBSystem::getTSNum()
{
    return tsNum;
}

void NBSystem::syncArrays()
{
    // Copy data from device to host
    if (!syncFlag) {
        try {
            queue.enqueueReadBuffer(d_pos[curBuffIndex],
                                    CL_TRUE, 0, vectBuffSize, pos.data());
            queue.enqueueReadBuffer(d_vel,
                                    CL_TRUE, 0, vectBuffSize, vel.data());

        } catch (cl::Error error) {
            std::cerr << "Caught exception at syncArrays(): " << error.what()
                      << '(' << error.err() << ')'
                      << std::endl;
            throw;
        }
        // Particles state at device and host are the same now
        syncFlag = true;
    }
}

void NBSystem::setOutPrec(unsigned int precision)
{
    outPrec = precision;
}
