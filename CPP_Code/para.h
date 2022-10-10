#include <string>
#include <math.h>

using namespace std;


const uint32_t N = 9; // Number of particles
const float T = 0.306; // Temperature of particles i.e used to set velocities
const float dt = 0.005; // time-step of simulation
const float tf = 100; // total time of simulation
const uint8_t d = 2; // dimension of simulation
const float rc = 2.5; // cut-off radius of Lennard-Jones Potential
const float rho = 0.00304; // density of particles

// directory used to store and input data
const string dir_name = "sim_data";
const string positions_input_file = "postions_N1000_d2_rho0.00304_sqaregrid.dat";
const string velocities_input_file = "velocities_N1000_d2_T0.306.dat";

//-----------Not to edit------------

float L = pow(N / rho, 1.0 / d);                       // size of box
float PE_cut = 4 * (1 / pow(rc, 12) - 1 / pow(rc, 6)); // lj potential at rc
int N_times_d = d * N;
