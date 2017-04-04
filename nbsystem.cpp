#include "nbsystem.h"
#include "clkernels.h"

NBSystem::NBSystem()
{
    pnum = 0;
    vnum = 0;
    dtCoef = 0.01;
    bndNum = 0;
    syncFlag = false;

    curBuffIndex = 0;

    // Default setting is to select first available cpu device
    //setupCL(-1, CL_DEVICE_TYPE_CPU);

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
        throw std::runtime_error("Sizes of the particles' state vectors"
                                 " (mass, pos, vel) are not equal.");

    unsigned int num = pos.size();
    for (unsigned int i = 0; i < num; i++ )
        addParticle(pos.at(i), vel.at(i));

    // Update buffers
    initializeCLBuffers();
    setupCLKernelsArgs();

    // Update accelerations for next timestep
    updateAcc(event);
    event.wait();
}

void NBSystem::removeParticle(unsigned int i)
{
    if (i >= pnum) throw std::domain_error("No particle with given id.");

    pos.erase(pos.begin() + i);
    vel.erase(vel.begin() + i);
    pnum--;

    // Update buffers
    initializeCLBuffers();
    setupCLKernelsArgs();

    // Update accelerations for next timestep
    updateAcc(event);
    event.wait();
}

void NBSystem::addBoundary(std::vector<cl_float4> &vertices)
{
    if (vertices.size() < 3)
        throw std::runtime_error("Invalid boundary:"
                                 " number of vertices is less than 3.");
    // Should add check if vertices are coplanar
    if (vertices.size() > 3) {
        cl_float4 a;
        a.s[0] = vertices.at(1).s[0] - vertices.at(0).s[0];
        a.s[1] = vertices.at(1).s[1] - vertices.at(0).s[1];
        a.s[2] = vertices.at(1).s[2] - vertices.at(0).s[2];
        a.s[3] = 0.0f;
        cl_float4 b;
        b.s[0] = vertices.at(2).s[0] - vertices.at(0).s[0];
        b.s[1] = vertices.at(2).s[1] - vertices.at(0).s[1];
        b.s[2] = vertices.at(2).s[2] - vertices.at(0).s[2];
        b.s[3] = 0.0f;
        cl_float4 axb;
        axb.s[0] = a.s[1] * b.s[2] - b.s[1] * a.s[2];
        axb.s[1] = -(a.s[0] * b.s[2] - b.s[0] * a.s[2]);
        axb.s[2] = a.s[0] * b.s[1] - b.s[0] * a.s[1];
        axb.s[3] = 0.0f;

        // Some number to compare result to zero
        float eps = 1e-5;

        for (unsigned int i = 3; i < vertices.size(); i++){
            cl_float4 c;
            c.s[0] = vertices.at(i).s[0] - vertices.at(0).s[0];
            c.s[1] = vertices.at(i).s[1] - vertices.at(0).s[1];
            c.s[2] = vertices.at(i).s[2] - vertices.at(0).s[2];
            c.s[3] = 0.0f;

            cl_float abc =
                    axb.s[0] * c.s[0] + axb.s[1] * c.s[1] + axb.s[2] * c.s[2];
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

    box.setBoundingBox(tri);

    for (unsigned int i = 1; i < vertices.size() - 1; i++) {
        Boundary bnd;
        bnd.v0 = vertices.at(0);

        // Calculate edges of boundary
        for (unsigned int j = 0; j < 3; j++) {
            bnd.edge1.s[j] = vertices.at(i).s[j] - bnd.v0.s[j];
            bnd.edge2.s[j] = vertices.at(i+1).s[j] - bnd.v0.s[j];
        }
        bnd.edge1.s[3] = 0.0f;
        bnd.edge2.s[3] = 0.0f;

        // Calculate normalized normal vector to boundary
        bnd.n0.s[0] = bnd.edge1.s[1] * bnd.edge2.s[2]
                                        - bnd.edge2.s[1] * bnd.edge1.s[2];
        bnd.n0.s[1] = -(bnd.edge1.s[0] * bnd.edge2.s[2]
                                        - bnd.edge2.s[0] * bnd.edge1.s[2]);
        bnd.n0.s[2] = bnd.edge1.s[0] * bnd.edge2.s[1]
                                        - bnd.edge2.s[0] * bnd.edge1.s[1];
        bnd.n0.s[3] = 0.0f;

        float n0length = std::sqrt(bnd.n0.s[0] * bnd.n0.s[0]
                                 + bnd.n0.s[1] * bnd.n0.s[1]
                                 + bnd.n0.s[2] * bnd.n0.s[2]);
        bnd.n0.s[0] /= n0length;
        bnd.n0.s[1] /= n0length;
        bnd.n0.s[2] /= n0length;

        boundaries.emplace_back(bnd);
    }

    bndNum = boundaries.size();

    initializeCLBuffers();
    setupCLKernelsArgs();
}

std::string NBSystem::boundariesString(unsigned int prec)
{
    std::ostringstream out;

    out << std::scientific;
    out << std::setprecision(prec);

    std::string delim = " ";
    int cWidth = 8; // Minimum column width is 8
    cWidth += prec;

    for (unsigned int i = 0; i < tri.size(); i = i + 3) {
        out << i / 3;
        out << std::endl;
        for (unsigned int j = 0; j < 3; j++) {
            out << delim << std::setw(cWidth) << tri.at(i+j).s[0];
            out << delim << std::setw(cWidth) << tri.at(i+j).s[1];
            out << delim << std::setw(cWidth) << tri.at(i+j).s[2];
            out << std::endl;
        }
    }

    return out.str().erase(out.str().size()-1);
}

std::string NBSystem::setupCL(int platformId, cl_device_type deviceType)
{
    std::ostringstream out;
    try {
        // List available platforms and devices
        std::vector<cl::Platform> platforms;
        out << "Available platforms:" << std::endl;
        cl::Platform::get(&platforms);
        for (auto& platform : platforms ) {
            std::string platformName;
            platform.getInfo((cl_platform_info)CL_PLATFORM_NAME,
                             &platformName);
            out << " " << platformName << std::endl;
            try {
                std::vector<cl::Device> deviceList;
                platform.getDevices(CL_DEVICE_TYPE_ALL, &deviceList);
                for(auto& deviceItem : deviceList) {
                    std::string deviceName;
                    deviceItem.getInfo((cl_device_info)CL_DEVICE_NAME,
                                       &deviceName);
                    out << "  " << deviceName << std::endl;
                }
            } catch (const cl::Error &error) {
                std::cerr << "  Error while getting device list: "
                          << error.what()
                          << '(' << error.err() << ')' << std::endl;
                //throw;
            }
        }

        //The next constructor attempts to use the first platform that
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
            out << "Selected device: " << deviceName;
        } catch (const cl::Error &error) {
            std::cerr << "  Error while selecting device: "
                      << error.what()
                      << '(' << error.err() << ')' << std::endl;
            throw;
        }

        // Kernel sources
//        std::string kernelSourceFileName = "kernelsSource.cl";
//        std::ifstream file(kernelSourceFileName);
//        if (!file.is_open())
//            throw std::runtime_error("Error reading kernels source file.");

//        std::string kernelSourceStr(std::istreambuf_iterator<char>(file),
//                                    (std::istreambuf_iterator<char>()));
//        file.close();

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
        stepInTime = cl::Kernel(program, "timeStepCL");
        estimateDt = cl::Kernel(program, "estimateDtCL");
        calcAcc = cl::Kernel(program, "calcAccCL");
        checkBoundaries = cl::Kernel(program, "checkBoundariesCL");
        calcEnKin = cl::Kernel(program, "calcEkinCL");
        calcEnPot = cl::Kernel(program, "calcEpotCL");

        // Create event
        event = cl::Event();

    } catch (cl::Error error) {
        std::cerr << "Caught exception at setupCL(): " << error.what()
                  << '(' << error.err() << ')'
                  << std::endl;
        throw;
    }
    initializeCLBuffers();
    setupCLKernelsArgs();
    return out.str();
}

void NBSystem::initializeCLBuffers()
{
    // Create and initialize device buffers
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

        int bndNum = boundaries.size();
        if (bndNum) {
            size_t bndBufSize = bndNum * sizeof(Boundary);
            d_bnd = cl::Buffer(context, CL_MEM_READ_WRITE,
                               bndBufSize);

            queue.enqueueWriteBuffer(d_bnd, CL_TRUE,
                                     0, bndBufSize,
                                     boundaries.data(), NULL, &event);
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
        //checkBoundaries.setArg(3, d_tri);
        //checkBoundaries.setArg(4, vnum);
        checkBoundaries.setArg(3, d_bnd);
        checkBoundaries.setArg(4, bndNum);

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
    if (coeff <= 0)
        throw std::invalid_argument("Invalid argument to setDtCoef()");
    dtCoef = coeff;
    setupCLKernelsArgs();
}


//void NBSystem::setLJparameters(double eps, double rmin)
//{
//    if (eps < 0 || rmin < 0)
//        throw std::invalid_argument("Invalid argument to setLJparams()");
//    ljEps = eps;
//    ljRmin = rmin;
//    setupCLKernelsArgs();
//}

double NBSystem::getEstDt()
{
    return estDt(event);
}


void NBSystem::evolve(double timeStep)
{
    if (timeStep <=0)
        throw std::invalid_argument("Invalid argument to evolve()");
    updatePos(timeStep, event);
    updateAcc(event);
}

void NBSystem::evolveIn(double interval)
{
    double time = 0;
    double eps = std::numeric_limits<double>::epsilon() * interval;
    double dt = 0;
    do {
        dt = estDt(event);
        if ((time + dt) > interval) dt = interval - time;
        evolve(dt);
        time += dt;
    } while (std::fabs(interval - time) > eps);
}

const std::vector<cl_float4> *NBSystem::velData()
{
    syncArrays();
    return &vel;
}

const std::vector<cl_float4> *NBSystem::posData()
{
    syncArrays();
    return &pos;
}

const std::vector<cl_float4> *NBSystem::bndData()
{
    return &tri;
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

unsigned int NBSystem::getPNum() const
{
    return pnum;
}

std::vector<std::array<double, 3> > NBSystem::getMinBox() const
{
    return box.getBoundingBox();
}

void NBSystem::addParticle(const unsigned int &num,
                           const std::array<double, 2> &mass,
                           const BoundingBox<double, 3> &posBox,
                           const BoundingBox<double, 3> &velBox)
{
    // This function generates random particle in the given area
    // with the given velocity and mass restrictions
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    std::vector<cl_float4> randomPos;
    std::vector<cl_float4> randomVel;

    for (unsigned int n = 0; n < num; n++) {
        cl_float4 p; // Particle's position
        cl_float4 v; // Particle's velocity
        for (unsigned int i = 0; i < 3; i++) {
            p.s[i] = posBox.minVertex()[i]
                    + (posBox.maxVertex()[i] - posBox.minVertex()[i])
                      * distribution(generator);
            v.s[i] = velBox.minVertex()[i]
                    + (velBox.maxVertex()[i] - velBox.minVertex()[i])
                      * distribution(generator);
        }

        p.s[3] = mass.at(0)
                + (mass.at(1) - mass.at(0)) * distribution(generator);
        v.s[3] = 0;

        randomPos.emplace_back(p);
        randomVel.emplace_back(v);
    }

    this->addParticle(randomPos, randomVel);
}
