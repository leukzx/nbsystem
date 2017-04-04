project(nbsystem)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
set(SRC_LIST ${SRC_LIST} boundingbox.h boundary.h kernels.cl)
#find_package(OpenCL REQUIRED)
#find_package(libconfig REQUIRED)
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries (${PROJECT_NAME} OpenCL config++)
#set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fexceptions -fexpensive-optimizations -fopenmp -O3 -lconfig++ -lOpenCL")

# Convert kernels to header file before project build
set(FAKE_TARGET fakeTarget)
add_custom_target(fakeTarget
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/kernel2header.cmake
    ${CMAKE_CURRENT_SOURCE_DIR} boundary.h kernels.cl clkernels.h)
add_dependencies(${PROJECT_NAME} ${FAKE_TARGET})


#add_custom_command(TARGET ${PROJECT_NAME}
#                   PRE_BUILD
#                   COMMAND ${CMAKE_COMMAND}
#                   -P ${CMAKE_CURRENT_SOURCE_DIR}/kernel2header.cmake
#                   ${CMAKE_CURRENT_SOURCE_DIR} kernels.cl clkernels.h)
