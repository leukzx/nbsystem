#ifndef BOUNDARY_H
#define BOUNDARY_H

#ifndef __OPENCL_C_VERSION__

#ifndef __CL_PLATFORM_H
#include <CL/cl_platform.h>
#endif

typedef cl_float4 float4;
#endif

typedef struct
{
    // Triangle's zero vertex
    float4 v0;
    // Triangle's edges
    float4 edge1;
    float4 edge2;
    // Normal vector to triangle
    float4 n0;
} Boundary;

#endif // BOUNDARY_H
