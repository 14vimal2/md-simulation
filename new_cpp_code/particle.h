// particle.h

#include "para.h"
#include <ctime>
#include <random>
#include <iostream>
#include "kd_tree.h"

#ifndef _particle_
#define _particle_

// Path: particle.h


class Particle
{
private:
    // all vector types data all stored in 1d array like the following sequence
    // r_0_0, r_0_1, r_0_2, r_1_0, r_1_1, ...
    // where r_i_j, i represent ith particle while j represent jth index of the ith particle

    double *r; // position
    double *rp; // previous position
    double *v; // velocity
    double *f; // force
    double *t; // temporary variable
    double T; // temperature
    // double **f_val; // force matrix for validation
    kd_tree *tree;
public:
    Particle();
    ~Particle();
    void initPosition();
    void initVelocity();
    void initPrevPosition();
    void displayParticle();
    void calculateForce();
    void timeAdvanceSingleStep();
};

#endif
