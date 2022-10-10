#include <string>
#include <math.h>

using namespace std;


const uint32_t N = 9; // Number of particles
const float T = 0.306; // Temperature of particles i.e used to set velocities
const float dt = 0.005; // time-step of simulation
const float tf = 100; // total time of simulation
const uint32_t d = 2; // dimension of simulation
const float rc = 2.5; // cut-off radius of Lennard-Jones Potential
const float rho = 0.00304; // density of particles


const string dir_name = "sim_data"; // directory used to store data
const string positions_input_type = "squaregrid"; // squaregrid or random
const string velocities_intput_type = "random"; // random or mbdist
// positions input file name which needs to be present in input_files folder
const string positions_input_file = "postions_N"+to_string(N) + "_d"+ to_string(d) +"_rho"+to_string(rho)+"_"+positions_input_type+".dat";
// velocities input file name which needs to be present in input_files folder
const string velocities_input_file = "velocities_N"+to_string(N) + "_d"+ to_string(d) +"_T"+ to_string(T)+"_"+velocities_intput_type+".dat";



//-----------Not to edit------------
float L = pow(N / rho, 1.0 / d);                       // size of box
float PE_cut = 4 * (1 / pow(rc, 12) - 1 / pow(rc, 6)); // lj potential at rc
uint32_t N_times_d = d * N;

float rho_max = 1 / pow(rc+1, d);

unsigned long long variables_datasize = N_times_d * sizeof(float);

