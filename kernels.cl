#line 2 "kernel.cl" // This helps to find error in concatenated sources file

#define forwardeuler 1
#define leapfrog 2
#define rk4 3

#include "solverSettings.h"

#define lennardjones 0
#define gravity 1


inline float4 force_LJ_ptp(const float4 *pos1, const float4 *pos2)
{
    // Calculates force to pos1 from pos2

    // Depth of potential well
    float epsilon = EPSILON;

    // Distance at which the potential reaches its minimum (F=0)
    float rm = RMIN;
    float rm6 = pown(rm, 6);
    float rm12 = rm6 * rm6;

    //float rc = 2.22724679535085 * rm; // Truncation distance
    //float Fc = -4.37754333446651E-02 * epsilon; // Force shift value

    // Relative radius-vector of p2 to p1
    float4 relr;
    relr.xyz = (*pos2).xyz - (*pos1).xyz;

    // Distance from p2 to p1
    float d2 = relr.x * relr.x + relr.y * relr.y + relr.z * relr.z;
    float d = sqrt(d2);

    //if (d <= rc) {
        float d6 = pown(d2, 3);
        float d7 = d6 * sqrt(d2);
        float d13 = d7 * d6;
        float4 f; // Force
        f.xyz = 12
              * epsilon
              * (rm6 / d7 - rm12 / d13)
              * relr.xyz / d;
        f.w = 0.0f;
    //}

    return f;
}

inline float4 force_gravity_ptp(const float4 *pos1, const float4 *pos2)
{
    // Force to pos1 from pos2
    // Gravity constant
    float G = GRAVITY_CONSTANT;

    // Relative radius-vector of p2 to p1
    float4 relr;
    relr.xyz = (*pos2).xyz - (*pos1).xyz;

    // Squared distance from p2 to p1
    float d2 = relr.x * relr.x + relr.y * relr.y + relr.z * relr.z;
    float d3 = d2 * sqrt(d2);

    float4 f; // Force
    f.xyz = G * (*pos1).w * (*pos2).w * relr.xyz / d3;
    f.w = 0.0f;

    return f;
}

inline float4 acc_total(const size_t id,
                        const float4 *r,
                        const __global float4 *pos,
                        const size_t pnum)
{
    //float4 r = pos[id]; // Position of the particle
    float4 a = (float4) 0.0f; // Acceleration

    // Accleration of the i-th particle
    float4 pos2;
    for (int j = 0; j < id; ++j) {
        pos2 = pos[j];
#if POTENTIAL == lennardjones
        a += force_LJ_ptp(r, &pos2);
#elif POTENTIAL == gravity
        a += force_gravity_ptp(r, &pos2);
#else
#error Error! Interaction potential is not set.
#endif
    }

    for (int j = id + 1; j < pnum; ++j) {
        pos2 = pos[j];
#if POTENTIAL == lennardjones
        a += force_LJ_ptp(r, &pos2);
#elif POTENTIAL == gravity
        a += force_gravity_ptp(r, &pos2);
#else
#error Error! Interaction potential is not set.
#endif
    }

    a /= (*r).w;
    a.w = 0.0f;

    return a FORCE_FIELD_ACC;
}

inline float energy_LJ_ptp(const float4 *pos1, const float4 *pos2)
{
    // Function calculates potential energy of a particle
    // in L-J potential field of other particle.

    // Depth of potential well
    float epsilon = EPSILON;
    // Distance at which the potential reaches its minimum (F=0)
    float rm = RMIN;

    //float rc = 2.22724679535085 * rm; // Truncation distance
    //float Vc = -1.63168911360000E-02 * epsilon; // Potential shift value

    float4 relr; // Relative radius-vector of p1 to p2
    relr.xyz = (*pos1).xyz - (*pos2).xyz;

    // Distance from p1 to p2
    float d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
    //float d = length(relr);
    //if (d <= rc) {
        float rm_over_d6 = pown(rm / d, 6);
        float rm_over_d12 = rm_over_d6 * rm_over_d6;
        float en = epsilon * (rm_over_d12 - 2 * rm_over_d6); //- Vc;
    //}

    return en;
}

inline float energy_gravity_ptp(const float4 *pos1, const float4 *pos2)
{
    // Function calculates potential energy of a particle
    // in gravity field of other particle.

    // Gravity constant
    float G = GRAVITY_CONSTANT;

    // Relative radius-vector of p2 to p1
    float4 relr;
    relr.xyz = (*pos2).xyz - (*pos1).xyz;

    // Distance from p2 to p1
    float d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);

    float en = - G * (*pos1).w * (*pos2).w / d;

    return en;
}

__kernel void calcEpotCL(const __global float4 *pos, __global float *epot)
{
    // Kernel calculates total potential energy of each particle

    size_t id = get_global_id(0); // Id of particle
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[id]; // Position of the particle

    float en = 0.0f; // Potential energy of the particle

    float4 pos2; // Position of other particle

    // Potential energy of the particle
    for (int j = 0; j < id; ++j) {
        pos2 = pos[j];
#if POTENTIAL == lennardjones
        en += energy_LJ_ptp(&r, &pos2);
#elif POTENTIAL == gravity
        en += energy_gravity_ptp(&r, &pos2);
#else
#error Error! Interaction potential is not set.
#endif
    }

    for (int j = id + 1; j < pnum; ++j) {
        pos2 = pos[j];
#if POTENTIAL == lennardjones
        en += energy_LJ_ptp(&r, &pos2);
#elif POTENTIAL == gravity
        en += energy_gravity_ptp(&r, &pos2);
#else
#error Error! Interaction potential is not set.
#endif
    }

    epot[id] = en FORCE_FIELD_ENR;
}

__kernel void calcEkinCL(const __global float4 *pos,
                         const __global float4 *vel,
                               __global float *ekin)
{
    // Kernel calculates kinetic energy of each particle
    size_t id = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)

    float m = pos[id].w; // Mass of the partice
    float4 v = vel[id]; // Velocity of the particle

    // Kinetic energy the particle
    float en = 0.5f * m * (v.x * v.x + v.y * v.y + v.z * v.z);

    ekin[id] = en;
}

__kernel void timeStepCL(const __global float4 *pos,
                               __global float4 *vel,
                               __global float4 *acc,
                         const float dt,
                               __global float4 *new_pos)
{
    // This function moves particles to new positions
#if INTEGRATION_METHOD == leapfrog
    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position
    float4 v = vel[gid]; // Velocity
    float4 a = acc[gid]; // Acceleration

    // New position and velocity
    // Kick
    v.xyz += a.xyz * dt / 2;
    // Drift
    r.xyz += v.xyz * dt;
    // Kick
    a = acc_total(gid, &r, pos, pnum);
    v.xyz += a.xyz * dt / 2;

    new_pos[gid] = r;
    vel[gid] = v;

#elif INTEGRATION_METHOD == forwardeuler
    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position
    float4 v = vel[gid]; // Velocity
    float4 a = acc[gid]; // Acceleration

    // New position and velocity
    r.xyz += v.xyz * dt;
    v.xyz += a.xyz * dt;

    new_pos[gid] = r;
    vel[gid] = v;

#elif INTEGRATION_METHOD == rk4
    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r0 = pos[gid]; // Initial position
    float4 v = vel[gid]; // Velocity
    float4 a0 = acc[gid]; // Initial Acceleration

    float dtdt = dt * dt;

    float4 r = r0 + v * dt / 2 + a0 * dtdt / 8;
    //new_pos[gid] = r;
    //barrier(CLK_GLOBAL_MEM_FENCE);

    float4 a1 = acc_total(gid, &r, pos, pnum);
    r = r0 + v * dt + a1 * dtdt / 2;
    //new_pos[gid] = r;
    //barrier(CLK_GLOBAL_MEM_FENCE);

    float4 a2 = acc_total(gid, &r, pos, pnum);
    r = r0 + v * dt + (a0 + 2 * a1) * dtdt /6;
    v += (a0 + 4 * a1 + a2) * dt / 6;
    new_pos[gid] = r;
    vel[gid] = v;
#else
#error Error! No integration method is defined.
#endif
}

__kernel void calcAccCL(const __global float4 *pos, __global float4 *acc)
{
    // Kernel This function calculates acceleration of particles in current position.

    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position of the particle

    float4 a = acc_total(gid, &r, pos, pnum);

    acc[gid] = a;
}


__kernel void estimateDtCL(const __global float4 *pos,
                           const __global float4 *vel,
                           const __global float4 *acc,
                                 __global float *est)
{
    // Function estimates collision time of particle id0 with
    // particle id1.
    size_t id0 = get_global_id(0);
    size_t id1 = get_global_id(1);
    float t = FLT_MAX;

    if (id0 != id1) {

        float est0 = FLT_MAX;
        float est1 = FLT_MAX;

        float4 relPos = pos[id0] - pos[id1];
        relPos.w = 0;
        float d =  sqrt(relPos.x * relPos.x
                + relPos.y * relPos.y
                + relPos.z * relPos.z);

        float4 relVel = vel[id0] - vel[id1];
        //if (dot(relPos, relVel) < 0) {
            float v =  sqrt(relVel.x * relVel.x
                    + relVel.y * relVel.y
                    + relVel.z * relVel.z);
            est0 = d / v;
            if (isnan(est0)) {est0 = FLT_MAX;}
        //}

        float4 relAcc = acc[id0] - acc[id1];
        //if (dot(relPos, relAcc) < 0) {
            float a =  sqrt(relAcc.x * relAcc.x
                    + relAcc.y * relAcc.y
                    + relAcc.z * relAcc.z);
            est1 = sqrt(2 * d / a);
            if (isnan(est1)) {est1 = FLT_MAX;}
        //}
        //est0 = (v != 0)?(d / v):FLT_MAX;
        //est1 = (a != 0)?(sqrt(d / a)):FLT_MAX;

        t = (est0 < est1)?est0:est1;
    }

    size_t width = get_global_size(0);

    est[id0 * width + id1] = t;
}

__kernel void checkBoundariesCL(const __global float4 *pos,
                                      __global float4 *vel,
                                      __global float4 *new_pos,
                                const __global Boundary *bnd,
                                const  int bnd_num)
{
//
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
//

    size_t gid = get_global_id(0);
    size_t bnum = bnd_num;

    float4 r1 = new_pos[gid]; // New position
    float4 r0 = pos[gid]; // Previous position
    float4 v = vel[gid];


    // Some variables to get intersection point of displacement and boundary
    float4 tuv = (float4) 0.0f;
    float  det = 0.0f;
    float4 dir, edge1, edge2 , pvec, tvec, qvec;
    dir = edge1 = edge2 = pvec = tvec = qvec = (float4) 0.0f;
    float eps = 0.00001f;//10.0f * FLT_MIN; // Some very small number to compare derterminant with
    float4 n0 = (float4) 0.0f; // Normalized plane normal
    float4 xpoint = (float4) 0.0f; // Intersection point
    float4 dr_out = (float4) 0.0f; // Out of space displacement of particle
    // dir * vpar is the relative position to r0 of intersection point
    float vpar = FLT_MAX;

    // Boundary crossing flag
    // Each time there is a crossing, particle moves to new position
    // Boundaries are checked several times until there are no new crossings
    bool flag = 0;

    do {
        // Particle displacement or direction vector, or segment
        dir = r1 - r0;
        //dir.w = 0;
        if (length(dir) == 0) break;
        // Set intersection flag to no crossing
        flag = 0;
        // Process all boundaries
        //__attribute__((opencl_unroll_hint(BOUNDARIES_NUM))).
        for (size_t i = 0; i < bnum; i++) {
            edge1 = bnd[i].edge1;
            edge2 = bnd[i].edge2;

            pvec = cross(dir, edge2);
            det = dot(edge1, pvec);
            // If determinant is 0, displacement lies in plane of triangle
            // If determinat is negative, triangle is crossed from outside.
            // (Vertices of boundaries should be listed in counter clockwise order
            // from inside of area.)
            // In both cases there is no boundary crossing. Goes to next triangle
            if (det < eps) continue;

            // Calculate vector from triangle's zero vertice to displacement origin
            tvec = r0 - bnd[i].v0;
            //tvec.w = 0;

            // Calculate U parameter and test bounds
            tuv.y = dot(tvec, pvec);
            if (tuv.y < 0.0f || tuv.y > det) continue;

            // Prepare to test V parameter
            qvec = cross(tvec, edge1);

            // Calculate V parameter and test bounds
            tuv.z = dot(dir, qvec);
            if (tuv.z < 0.0f || (tuv.y + tuv.z) > det) continue;

            // Calculate T parameter and test bounds
            // tuv.x = 0.0f if origin of displacement lies in the triangle's plane
            // This may happen only when there was boundary crossing and
            // it was aleady processed.
            // This crossing should be excluded from further consderation.
            tuv.x  = dot (edge2, qvec);
            if (tuv.x <= 0.0f || tuv.x > det) continue;

            // Segment intersects triangle
            // Set intersection flag
            flag = 1;

            // Scale parameters
            tuv /= det;

            // Find closest intersection point
            if (vpar > tuv.x) {
                vpar = tuv.x;
                // Normal vector to crossed boundary
                n0 = bnd[i].n0;
            }
        }

        if (flag) {
            // Calculate intersection point
            xpoint = r0 + dir * vpar;

            // New velocity
            v = v - 2 * dot(v, n0) * n0; // Mirrors v against boundary

            // Out of plane displacement
            dr_out = r1 - xpoint;
            // Mirror dr against boundary
            dr_out = dr_out - 2 * dot(dr_out, n0) * n0;
            // Save old position.
            //r0 = r1;
            // Crossing point becomes origin of the displacement
            // to process all futher possible boundary crossings
            r0 = xpoint;

            //Calculate new position
            r1 = xpoint + dr_out;

            // Save new values to array
            new_pos[gid] = r1;
            vel[gid] = v;
        }
    } while (flag); // Do loop until there are no new intersections
}
