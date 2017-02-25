float4 get_acc_LJ(__global float4 *pos, size_t id, size_t pnum)
{
    float4 r = pos[id]; // Position of the particle
    float4 a = (float4) 0.0f; // Acceleration

    float epsilon = 1.0f;

    float rm = 1.0f; // Distance at which the potential reaches its minimum (F=0)
    float rm6 = pown(rm, 6);
    float rm12 = pown(rm6, 2);

    float4 relr = (float4) 0.0f; // Relative radius-vector of p2 to p1
    float d = 0.0f; // Distance from p2 to p1
    float d7 = 0.0f;
    float d13 = 0.0f;

    // Accleration of the particle
    for (int j = 0; j < id; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
        d7 = pown(d, 7);
        d13 = pow (d7, 2) / d;
        a.xyz += 12
               * epsilon
               * (rm6 / d7 - rm12 / d13)
               * relr.xyz / d;
    }

    for (int j = id + 1; j < pnum; ++j) {
        relr.xyz = pos[j].xyz - r.xyz;
        d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
        d7 = pown(d, 7);
        d13 = pow (d7, 2) / d;
        a.xyz += 12
               * epsilon
               * (rm6 / d7 - rm12 / d13)
               * relr.xyz / d;
    }

    a.w = 0.0f;
    a /= r.w;

    return a;
}

float get_energy_LJ_ptp(float4 *pos1, __global float4 *pos2)
{
    // Function calculates potential energy of a particle
    // in L-J potential field of other particle.

    float en = 0.0f; // Potential energy

    float epsilon = 1.0f;

    float rm = 1.0f; // Distance at which the potential reaches its minimum (F=0)
    float rm_over_d6 = 0.0f;
    float rm_over_d12 = 0.0f;

    float4 relr = (float4) 0.0f; // Relative radius-vector of p1 to p2
    float d = 0.0f; // Distance from p1 to p2

    relr.xyz = (*pos1).xyz - (*pos2).xyz;
    d = sqrt(relr.x * relr.x + relr.y * relr.y + relr.z * relr.z);
    rm_over_d6 = pown(rm / d, 6);
    rm_over_d12 = pow (rm_over_d6, 2);
    en = epsilon * (rm_over_d12 - 2 * rm_over_d6);

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

    epot[id] = en;
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

    acc[gid] = a;
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

    float4 relPos;
    float4 relVel;
    float4 relAcc;

    float d, v, a, t, est0, est1;


    if (id0 != id1) {
        est0 = FLT_MAX;
        est1 = FLT_MAX;

        relPos = pos[id0] - pos[id1];
        relPos.w = 0;
        d =  sqrt(relPos.x * relPos.x
                + relPos.y * relPos.y
                + relPos.z * relPos.z);

        relVel = vel[id0] - vel[id1];
        //if (dot(relPos, relVel) < 0) {
            v =  sqrt(relVel.x * relVel.x
                    + relVel.y * relVel.y
                    + relVel.z * relVel.z);
            est0 = d / v;
            if (isnan(est0)) {est0 = FLT_MAX;}
        //}

        relAcc = acc[id0] - acc[id1];
        //if (dot(relPos, relAcc) < 0) {
            a =  sqrt(relAcc.x * relAcc.x
                    + relAcc.y * relAcc.y
                    + relAcc.z * relAcc.z);
            est1 = sqrt(2 * d / a);
            if (isnan(est1)) {est1 = FLT_MAX;}
        //}
        //est0 = (v != 0)?(d / v):FLT_MAX;
        //est1 = (a != 0)?(sqrt(d / a)):FLT_MAX;

        t = (est0 < est1)?est0:est1;
    } else {
        t = FLT_MAX;
    }

    size_t width = get_global_size(0);

    est[id0 * width + id1] = t;
}

__kernel void checkBoundariesCL(__global float4 *pos,
                                __global float4 *vel,
                                __global float4 *newPos,
                                __global float4 *verts,
                                            int vNum)
{
//
//https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
//
    size_t gid = get_global_id(0);
    size_t vnum = vNum;

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
    float epsilon = 0.00001f;//10.0f * FLT_MIN; // Some very small number to compare derterminant with
    float4 n0 = (float4) 0.0f; // Normalized plane normal
    float4 xpoint = (float4) 0.0f; // Intersection point
    float4 dr_out = (float4) 0.0f; // Out of space displacement of particle

    // Boundary crossing flag
    // Each time there is a crossing, particle moves to new position
    // Boundaries are checked several times until there are no new crossings
    bool flag = 0;

    do{
    flag = 0;

    // Particle displacement or direction vector, or segment
    dir = r1 - r0;
    //dir.w = 0;

    if (length(dir) == 0) {break;}

    for (size_t i = 0; i < vnum; i += 3) {
        // Get another 3 vertices from vertices array
        for (size_t k = 0; k < 3; ++k) {vertex[k] = verts[i+k];}

        // ?Maybe it is better to pre process edges
        edge1 = vertex[1] - vertex[0];
        edge1.w = 0;
        edge2 = vertex[2] - vertex[0];
        edge2.w = 0;

        pvec = cross(dir, edge2);
        det = dot(edge1, pvec);
        // If determinant is 0, displacement lies in plane of triangle
        // If determinat is negative, triangle is crossed from outside.
        // (Vertices of boundaries should be listed in counter clockwise order
        // from inside of area.)
        // In both cases there is no boundary crossing. Goes to next triangle
        if (det < epsilon) {continue;}

        // Calculate vector from triangle's zero vertice to displacement origin
        tvec = r0 - vertex[0];
        tvec.w = 0;

        // Calculate U parameter and test bounds
        tuv.y = dot(tvec, pvec);
        if (tuv.y < 0.0f || tuv.y > det) {continue;}

        // Prepare to test V parameter
        qvec = cross(tvec, edge1);

        // Calculate V parameter and test bounds
        tuv.z = dot(dir, qvec);
        if (tuv.z < 0.0f || (tuv.y + tuv.z) > det) {continue;}

        // Calculate T parameter and test bounds
        tuv.x  = dot (edge2, qvec);
        if (tuv.x < 0.0f || tuv.x > det) {continue;}

        // Segment intersects triangle
        // Scale parameters
        tuv /= det;

        // Set intersection flag
        flag = 1;
        // Calculate intersection point
        xpoint = r0 + dir * tuv.x;
        // Calculate plane normal of triangle
        n0 = normalize(cross(edge1, edge2));

        // New velocity
        v = v - 2 * dot(v, n0) * n0; // Mirrors v against boundary

        // Out of plane displacement
        dr_out = r1 - xpoint;
        dr_out = dr_out - 2 * dot(dr_out, n0) * n0; // Mirrors dr against boundary
        // Save old position. Calculate new position
        r0 = r1;
        r1 = xpoint + dr_out;

        // Save new values to array
        newPos[gid] = r1;
        vel[gid] = v;
        // Stop check other triangles. Do all checks again. (?)
        break;
    }

    } while (flag); // Do loop until there are no new intersections
}
