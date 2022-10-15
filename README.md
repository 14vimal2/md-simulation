# md-simulation
This repository contains codes for md-simulation and its application

Right now I am developing this, but it has bugs. You are welcome to contribute.

## Download

```bash
git clone https://github.com/14vimal2/md-simulation
```


## Usage

```bash
cd md-simulation
```
### Quickstart
Download this repo first and go in to `md-simulation\` folder first. Make a file say `testrun.cpp` and copy this code and run.

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

## Project Structure

```
│   README.md
│   test.cpp
│   test.exe
│
├───CPP_Code
│   │   bindata.h
│   │   create_input_files.cpp
│   │   create_input_files.exe
│   │   kdtree.exe
│   │   kdtree.h
│   │   para.exe
│   │   para.h
│   │   particles.h
│   │   test.exe
│   │
│   └───input_files
│           positions_N100_d2_rho0.003040_squaregrid.dat
│           positions_N9_d2_rho0.003040_squaregrid.dat
│           velocities_N100_d2_T0.306000_random.dat
│           velocities_N9_d2_T0.306000_random.dat
│
└───Python_Code
        fns.py
        main.py
        para.py
        particles.py
        plot.py
        plottercopy.py
```
## Contributing
Pull requests are welcomed. For major changes, please open an issue first to discuss what you would like to change.