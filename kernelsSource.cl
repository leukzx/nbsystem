// EPSILON is the depth of the potential well
#define EPSILON 1.0
// RMIN is the distance at which the potential reaches its minimum
#define RMIN 1.0

// FORCE_FIELD_ACC is the mathematical expression of external field acceleration
//#define FORCE_FIELD_ACC +(float4)(0.0f, 0.0f, -10.0f, 0.0f)
#define FORCE_FIELD_ACC

// FORCE_FIELD_ENR is the mathematical expression of external field energy
//#define FORCE_FIELD_ENR +r.w*dot(-2*(float4)(0.0f, 0.0f, -10.0f, 0.0f), r)
#define FORCE_FIELD_ENR

#include "boundary.h"

inline float4 get_acc_LJ(__global float4 *pos, size_t id, size_t pnum)
{
    float4 r = pos[id]; // Position of the particle
    float4 a = (float4) 0.0f; // Acceleration

    float epsilon = EPSILON;

    float rm = RMIN; // Distance at which the potential reaches its minimum (F=0)
    float rm6 = pown(rm, 6);
    float rm12 = pown(rm6, 2);

    float rc = 2.22724679535085 * rm; // Truncation distance
    float Fc = -4.37754333446651E-02 * epsilon; // Force shift value

    float4 relr = (float4) 0.0f; // Relative radius-vector of p2 to p1
    float d = 0.0f; // Distance from p2 to p1
    float d7 = 0.0f;
    float d13 = 0.0f;

    // Accleration of the i-th particle

    for (int j = 0; j < id; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
        if (d <= rc) {
            d7 = pown(d, 7);
            d13 = pown(d7, 2) / d;
            a.xyz += (12
                   * epsilon
                   * (rm6 / d7 - rm12 / d13) - Fc)
                   * relr.xyz / d;
        }
    }

    for (int j = id + 1; j < pnum; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
        if (d <= rc) {
            d7 = pown(d, 7);
            d13 = pown(d7, 2) / d;
            a.xyz += (12
                   * epsilon
                   * (rm6 / d7 - rm12 / d13) + Fc)
                   * relr.xyz / d;
        }
    }

    a.w = 0.0f;
    a /= r.w;

    return a;
}

inline float get_energy_LJ_ptp(float4 *pos1, __global float4 *pos2)
{
    // Function calculates potential energy of a particle
    // in L-J potential field of other particle.

    float en = 0.0f; // Potential energy

    float epsilon = EPSILON;

    float rm = RMIN; // Distance at which the potential reaches its minimum (F=0)
    float rm_over_d6 = 0.0f;
    float rm_over_d12 = 0.0f;

    float rc = 2.22724679535085 * rm; // Truncation distance
    float Vc = -1.63168911360000E-02 * epsilon; // Potential shift value

    float4 relr = (float4) 0.0f; // Relative radius-vector of p1 to p2
    float d = 0.0f; // Distance from p1 to p2

    relr.xyz = (*pos1).xyz - (*pos2).xyz;
    d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
    if (d <= rc) {
        rm_over_d6 = pown(rm / d, 6);
        rm_over_d12 = pown(rm_over_d6, 2);
        en = epsilon * (rm_over_d12 - 2 * rm_over_d6) - Vc;
    }

    return en;
}


__kernel void calcEpotLJCL(__global float4 *pos, __global float *epot)
{
    // Kernel calculates total potential energy of each particle

    size_t id = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[id]; // Position of the particle

    float en = 0.0f; // Potential energy

    // Potential energy of the particle
    for (int j = 0; j < id; ++j) {
        en += get_energy_LJ_ptp(&r,  &pos[j]);
   }

    for (int j = id + 1; j < pnum; ++j) {
        en += get_energy_LJ_ptp(&r, &pos[j]);
    }

    epot[id] = en FORCE_FIELD_ENR;
}

__kernel void calcEkinCL(__global float4 *pos,
                       __global float4 *vel,
                       __global float *ekin)
{
    // Kernel calculates kinetic energy of each particle
    size_t id = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)

    float m = pos[id].w; // Mass of the partice
    float4 v = vel[id]; // Velocity of the particle

    float en = 0.0f; // Kinetic energy the particle

    en = 0.5f * m * (v.x * v.x + v.y * v.y + v.z * v.z);

    ekin[id] = en;
}

__kernel void tsForwardEulerCL(__global float4 *pos,
                         __global float4 *vel,
                         __global float4 *acc,
                         float dt,
                         __global float4 *newPos)
{
    // This function moves particles in new position.
    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position
    float4 v = vel[gid]; // Velocity
    //float4 a = (float4) 0.0f ; // Acceleration
    float4 a = acc[gid];

    // Particle's acceleration
    //a = get_acc_LJ(pos, gid, pnum);

    // New position and velocity
    r.xyz += v.xyz * dt;
    v.xyz += a.xyz * dt;

    newPos[gid] = r;
    vel[gid] = v;
}

__kernel void calcAccCL(__global float4 *pos, __global float4 *acc)
{
    // Kernel This function calculates acceleration of particles in current position.

    size_t gid = get_global_id(0);
    size_t pnum = get_global_size(0); // Number of particles (work-items)
    float4 r = pos[gid]; // Position of the particle

    float4 a = (float4) 0.0f; // Acceleration

    a = get_acc_LJ(pos, gid, pnum);

    acc[gid] = a FORCE_FIELD_ACC;
}


__kernel void estimateDtCL(__global float4 *pos,
                           __global float4 *vel,
                           __global float4 *acc,
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

__kernel void checkBoundariesCL(__global float4 *pos,
                                __global float4 *vel,
                                __global float4 *newPos,
                                __global Boundary *bnd,
                                            int bnd_num)
{
//
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
//
    size_t gid = get_global_id(0);
    size_t bnum = bnd_num;

    float4 r1 = newPos[gid]; // New position
    float4 r0 = pos[gid]; // Previous position
    float4 v = vel[gid];

    // Border vertices
    float4 vertex[3] = {(float4) 0.0f, (float4) 0.0f, (float4) 0.0f};

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
            newPos[gid] = r1;
            vel[gid] = v;
        }
    } while (flag); // Do loop until there are no new intersections
}
