// Possible integration methods: forwardeuler, leapfrog, rk4
#define INTEGRATION_METHOD rk4

// Particle to particle interaction potential
#define POTENTIAL gravity

// FORCE_FIELD_ACC is the mathematical expression of external field acceleration
//#define FORCE_FIELD_ACC +(float4)(0.0f, 0.0f, -10.0f, 0.0f)
#define FORCE_FIELD_ACC

// FORCE_FIELD_ENR is the mathematical expression of external field energy
// (r.w is the mass of a particle)
//#define FORCE_FIELD_ENR +r.w*dot(-2*(float4)(0.0f, 0.0f, -10.0f, 0.0f), r)
#define FORCE_FIELD_ENR

// GRAVITY_CONSTANT is the constant in Newton's gravity law
// (for gravity potential)
#define GRAVITY_CONSTANT 1

// EPSILON is the depth of the potential well (for L-J potential)
#define EPSILON 1.0f
// RMIN is the distance at which the potential reaches its minimum
// (for L-J potential)
#define RMIN 1.0f
