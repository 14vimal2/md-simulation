# md-simulation
This repository contains codes for md-simulation and its application

Right now I am developing this, but it has bugs. You are welcome to contribute.

# Download

```bash
git clone https://github.com/14vimal2/md-simulation
```

# Usage

```bash
cd md-simulation
```

Make a file lets [testrun.cpp] and copy this code

```CPP
#include <iostream>
#include "CPP_Code/particles.h"
using namespace std;

int main()
{
    Particles P = Particles<uint32_t, float>();
    int tf_ = tf / dt;
    for (int i = 0; i < tf_; i++)
    {
        P.computeForce();
        P.time_advance_single_step(true);
        cout << endl << P.get_KE() << endl;
    }

    P.DisplayPositions();
    P.DisplayVelocities();
    cout << P.get_KE();
}
```