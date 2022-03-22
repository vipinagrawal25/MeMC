# MeMC
A Monte-Carlo code to simulate of fluctuating enclosed membranes. Such system can be
of relevance for the study of viruses, exosomes, e.t.c.

## Prerequisites

MeMC requires following libraries:

1) A C++ compiler. We have tested the code against gnu [g++](https://gcc.gnu.org/) versions 5.4.0 and above on x86_64 CPU.

2) [Hdf5](https://www.hdfgroup.org/solutions/hdf5) libraries for reading and writing data.  

3) [Python](https://www.python.org/) version with [scipy](https://www.scipy.org), [numpy](https://www.numpy.org), [h5py](https://www.h5py.org) and [numpy-quartenion](https://https://pypi.org/project/numpy-quaternion/) installed. We have tested the code against python version 3.8 and above.

For details of the installation for different packages, check the instructions on the official web page of the packages. If the operating system is Ubuntu then, g++ and hdf5 can be installed using the package manager apt,

```bash
apt install g++ hdf5
```

In order to install the required python libraries we suggest using the standard
python package manager pip, 
```bash
pip install scipy numpy-quaternion
```

## Installation

Once the prerequisites are satisfied, the installation of the MeMC amounts executing
the following bash commands: 

```bash
git clone https://github.com/vipinagrawal25/MeMC
cd MeMC
make
```
If successful, one should find an extra "bin" directory in the folder. The bin
directory will contain the binaries "exe_start" and "exe_memc".

---
**NOTE**
If the installation fails and shows the error

```c
In file included from hdf5.h:22:0,
                 from src/hdf5_io.cpp:1:
H5public.h:68:17: fatal error: mpi.h: No such file or directory
compilation terminated.
```
then we suggest the user to compile the code as:

```bash 
make CC=/path/to/h5c++.
```



# Using the Package:

The package is built to understand the physics of Exosomes and viruses where the
thermal noise are relevant. Therefore, we expect the user to have a prior
understanding of Metropolis algorithm of the Monte Carlo simulation and basic
concept of Elasticity. Both the concept is described in detail in the document doc/paper.pdf. 
That being said, we shall now dive deeper and explain how different part of the
code functions in the sections that follow. 

## Constructing the membrane

To begin the simulation we generate a equilibrated randomized position. For this
purpose use the executable "bin/exe_start". The executable takes three arguments:
Number of random points to be equilibrated, geometry of the surface "sph" for a
surface of sphere and "cart" for flat plane and folder name to write the output
files. The following command in the root directory of the repository will generate a
configuration of 1024 points in surface of sphere. All the outputs will be dumped
inside data_sph.

NOTE :: The simulation will run for 6e4 monte-carlo steps, which amounts to some
time in intel core. To go ahead with the other part it is sufficient to kill after
1e3 monte-carlo steps  

```bash
./bin/exe_start 1024 sph data_sph

```

 Before equilibration      |  After equilibration
:-------------------------:|:-------------------------:
![](./doc/figs/surf_mc_random.png)   |  ![](./doc/figs/surf_mc_lattice.png)

## The turn (Generating the connections)

For MeMC, one needs to triangulate the points generated as described in para above.
We rely on the python libraries to do so. The code to do the trick is  utils/gen_memc_conf.py, which take the data to be triangulated as an input.  Executing the following

```bash
python utils/gen_memc_conf.py data_sph/part_pos0003.bin 
```
will write the input file with all the connections in conf/dmemc_conf.h5.  


## The prestige  
Once the connections are set one can do the final simulation with bin/exe_memc. The
binary takes two arguments parameter file and the folder to write the simulation
data. One example of the parameter file is given in `Examples/para_file.in` where we give input to the simulation in plain text ascii format. For more details about the various parameters used, see the advanced document of the code.  Along with arguments, it is also expected to have a directory conf with file
"dmemc_conf.h5" inside it. Once all is set, run the code code with

```bash
./bin/exe_memc Examples/para_file.in out_memc

```

# Output's
Apart from the snapshot of the position, the code outputs % of accepted moves in
second column and total energy of the configuration in the third column in "mc_log".
The log file is written inside the specified folder.

# Examples

An example shell script to conduct all the steps stated above is stored in
"Examples/execute.sh".  Change the directory to Examples and run

```bash
sh execute.sh
```

The code will finish in about 30 minutes on Intel(R) Core(TM) i5-8265U CPU. You can
relax till then. Once done we have provided a gnuplot plot script  and sample
histogram "hist_start.dat" and "hist_memc.dat" of energies from the same simulation
we had conducted. Execution of script can be done as:
```bash 
gnuplot plot.gnu
```
Two plot window should appear. The histogram of energy for the
data generated locally is plotted using lines and points, whereas the case we have
done is in points. To make histogram we have used gsl-histogram library. One can
also use numpy or other standard libraries to do the same. 

# Visualization
Visualization can be done using [visit](https://visit.org) or [paraview](https://www.paraview.org). We provide a python utility "utils/viz_memc.py". The program takes three arguments: output from the exe_memc stored in the specified folder, the connections set by utils/gen_memc_conf.py, and the output file with extension .vtk. A sample execution will look like
```bash
python utils/viz_memc.py out_memc/snap_0004.h5 conf/dmemc_conf.h5 check_viz.vtk
```

If the execution is successful, the file "check_viz.vtk" will be written in the root
directory. In visit load the .vtk file and select `subset->domains` or `mesh->mesh`
to see the result.

