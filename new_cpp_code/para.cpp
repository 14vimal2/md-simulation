// para.cpp

#include "para.h"



//-----------Not to edit------------
double LENGTH = std::pow(NUMBER_OF_PARTICLES / DENSITY, 1.0 / DIMENSION);                       // size of box
double CUTOFF_POTENTIAL = 4 * (1 / std::pow(CUTOFF_DISTANCE, 12) - 1 / std::pow(CUTOFF_DISTANCE, 6)); // lj potential at CUTOFF_DISTANCE
uint32_t ND =DIMENSION* NUMBER_OF_PARTICLES; // number of variables

double DENSITY_MAX = 1 / std::pow(CUTOFF_DISTANCE+1, DIMENSION);

unsigned long long variables_datasize = ND * sizeof(double);

