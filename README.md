# molecular-dynamics
The following code may be used to perform a full simulation of molecular dynamics (MD) for a set of 500 particles in a cubic box with a side length of 10 units, using periodic boundary conditions. The set of programs is capable of constructing the initial configuration of the system and predicting its evolution for any arbitrary time, influenced by a Lennard-Jones potential. The code keeps track of the coordinates, velocities, and accelerations of every particle. Using this information, thermodynamic variables are calculated. A detailed explanation of the physics and analysis techniques is provided in the [Explicación](Explicación.pdf).

## Installation and using of the code

The different files are mainly written in Fortran, and for compilation, CMake is used, so both are needed. Moreover, the installation guide is designed for a UNIX system, but in the end, it will be shown how to use this code on a Windows system

### Download the repository

The simplest way to download the code is using git, 

```git clone https://github.com/pgrusell/molecular-dynamics.git```

but, it can also be donwloaded directly as a ```.dot``` file.

### Installation, compilation and execution

As usually with CMake compilations, the steps are.

```
from molecular-dynamics directory
$ mkdir build
$ cd build
$ cmake ..
$ make
```

this process will generate in the build directory an executable file named ```crea_red```. To simplify and automate this, the file ```inicializacion.sh```, on UNIX systems, has been created. After grating permision with, 

```$ chmod +x inicializacion.sh```

run

```$ ./inicializacion.sh```

This script compiles and runs the programs, generating a new folder named  ```resultados```, which will contain the coordinates, energies, and more.

## Modification of the simulation

The simulation data source is the file ```configuracion.txt``` wich contains the size of the box, the constant energy of the system and the number of iterations of the simulation. 





