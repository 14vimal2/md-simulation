// para.h

#include <string>
#include <math.h>

#ifndef _para_
#define _para_


#define NUMBER_OF_PARTICLES 900 // Number of particles
#define DIMENSION 2 // dimension
#define TEMPERATURE  10 // Temperature of particles i.e used to set velocities
#define TIMESTEP  0.001 // time-step of simulation
#define FINAL_TIME 20 // total time of simulation
#define CUTOFF_DISTANCE  2.5 // cut-off radius of Lennard-Jones Potential
#define DENSITY 0.1 // density of particles


//-----------Not to edit------------
extern double LENGTH;// = std::pow(NUMBER_OF_PARTICLES / DENSITY, 1.0 / DIMENSION);                       // size of box
extern double CUTOFF_POTENTIAL;// = 4 * (1 / std::pow(CUTOFF_DISTANCE, 12) - 1 / std::pow(CUTOFF_DISTANCE, 6)); // lj potential at CUTOFF_DISTANCE
extern uint32_t ND;// =DIMENSION* NUMBER_OF_PARTICLES; // number of variables

extern double DENSITY_MAX;// = 1 / std::pow(CUTOFF_DISTANCE+1, DIMENSION);

#endif