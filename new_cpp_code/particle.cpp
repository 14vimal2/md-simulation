// particle.cpp

#include "particle.h"
#include "kd_tree.h"
#include <vector>
#include <map>
#include <set>
#include <algorithm>

Particle::Particle()
{
    r = new double[ND];
    rp = new double[ND];
    v = new double[ND];
    f = new double[ND];
    t = new double[ND];
    // make kd tree
    tree = new kd_tree();

    initPosition();
    initVelocity();
    initPrevPosition();
}

// initiate position of all particles into lattice points
void Particle::initPosition() {
    int sn = ceil(std::pow(NUMBER_OF_PARTICLES, 1.0/DIMENSION));
    for (size_t i = 0; i < ND; i++)
    {
        r[i] = (( (i / DIMENSION) / (int)(std::pow(sn, i % DIMENSION))) % sn + (1.0 / 2.0)) * LENGTH / sn;
    }    

}

// initiate velocity of all particles such that it matches with velocities and vcm = 0
void Particle::initVelocity() {
    double *vcm = new double[DIMENSION]; // velocity of center of mass
    std::fill_n(vcm, DIMENSION, 0);

    for (size_t i = 0; i < ND; i++)
    {
        v[i] = (rand() / (double)(RAND_MAX)) - 1 / 2.0;
        vcm[i % DIMENSION] += v[i];
    }

    for (size_t j = 0; j < DIMENSION; j++)
    {
        vcm[j] /= NUMBER_OF_PARTICLES;
    }

    double ke = 0;

    // subtract velocity of center of mass from each particle
    for (size_t i = 0; i < ND; i++)
    {
        v[i] -= vcm[i % DIMENSION];
        ke += v[i]*v[i];
    }    

    double scale_factor = sqrt(TEMPERATURE * ND / ke);
    ke = 0;

    for (size_t i = 0; i < ND; i++)
    {
        v[i] *= scale_factor;
        ke += v[i]*v[i];
    }

    delete[] vcm;
}

// initiate previous position using present position and velocity
void Particle::initPrevPosition() {
    for (size_t i = 0; i < ND; i++)
    {
        rp[i] = r[i] - v[i] * TIMESTEP;
    }
}

// calculate force on each particle
void Particle::calculateForce() {
    // initiate force to 0
    std::fill_n(f, ND, 0);

    // build tree
    tree->build(r);

    std::vector<int> indices;
    std::map<int, std::unordered_set<int>> processed;

    // calculate and update force for each particle
    for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        tree->searchNeighbour(i, r, indices);
    }

    int s = indices.size()/2;

    int duplicates = 0;


    for (size_t i = 0; i < s; i++)
    {
        double r2 = 0.;
        int pi = indices[i * 2], pj = indices[i * 2 + 1];

        if (processed[pi].count(pj) > 0 ) // if already processed then ignore otherwise add to processed pair
        {
            duplicates++;
            continue;
        } else {
            processed[pi].insert(pj);
        }
        


        for (size_t j = 0; j < DIMENSION; j++)
        {
            r2 += (r[ pi * DIMENSION  + j ] - r[ pj * DIMENSION  + j ])
             * (r[ pi * DIMENSION  + j ] - r[ pj * DIMENSION  + j ]);
        }

        if (r2 < CUTOFF_DISTANCE * CUTOFF_DISTANCE)
        {
            double r2i = 1.0 / r2;
            double r6i = r2i*r2i*r2i;
            double ff = 48.0*r2i*r6i*(r6i-0.5);
            double force_mag = 0;
            for (size_t j = 0; j < DIMENSION; j++)
            {
                f[pi * DIMENSION + j] += ff * (r[ pi * DIMENSION  + j ] - r[ pj * DIMENSION  + j ]);
                f[pj * DIMENSION + j] -= ff * (r[ pi * DIMENSION  + j ] - r[ pj * DIMENSION  + j ]);
            }
        }    

    }

    if (s/2 != duplicates )
        std::cout << "s = " << s << " and d = " << duplicates << std::endl;
        
}


// move particles by one time step and update velocity
void Particle::timeAdvanceSingleStep() {
    for (size_t i = 0; i < ND; i++)
    {
        t[i] = r[i];
        r[i] = 2 * r[i] - rp[i] + f[i] * TIMESTEP * TIMESTEP;
        rp[i] = t[i];
        v[i] = (r[i] - t[i])/(2 * TIMESTEP);

        if (r[i] < 0)
        {
            r[i] = -r[i];
            rp[i] = -rp[i];
            v[i] = -v[i];
        }
        else if (r[i] > LENGTH) {
            r[i] = 2*LENGTH - r[i];
            rp[i] = 2*LENGTH -rp[i];
            v[i] = -v[i];
        }
    }
}

// prints all particles position and velocities
void Particle::displayParticle() {

    for (size_t i = 0; i < NUMBER_OF_PARTICLES; i++)
    {
        for (size_t j = 0; j < DIMENSION; j++)
        {
            std::cout << r[i*DIMENSION + j] << " ";
        }
        for (size_t j = 0; j < DIMENSION-1; j++)
        {
            std::cout << v[i*DIMENSION + j] << " ";
        }

        std::cout << v[i*DIMENSION + DIMENSION-1] << "\n";
    }
    
}

Particle::~Particle()
{
    delete[] r;
    delete[] rp;
    delete[] v;
    delete[] f;
    delete[] t;
    delete tree;
}