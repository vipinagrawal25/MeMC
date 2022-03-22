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

The package is built to understand the physics of Exosomes or viruses where the
thermal noise are relevant. These nano-vesicles are studied experimentally using
Atomic Force Microscopy. This package tends to simulate the experiment in the best
possible way. Therefore, we expect the user to have a prior understanding of
the Metropolis algorithm for the Monte Carlo simulation and a basic concept of Elasticity
pertaining to membranes. Both the concept is described in detail in the document
doc/paper.pdf.  That being said, we shall now dive deeper and explain how different
part of the code functions in the following sections:

## Constructing the membrane

To begin the simulation we generate a equilibrated randomized position on a surface
of a sphere.  For details we refer the reader to section xx.xx of doc/paper.pdf. The
main code for this purpose is given in `main/start.cpp` and the relevant  executable is "bin/exe_start". 
The executable takes three arguments:
1) Number of random points to be equilibrated
2) Geometry of the surface "sph" for a
surface of sphere and "cart" for flat plane
3) Folder name to write the output files. 

Rest of the parameter are hard coded. For instance, the radius of the sphere on
which we generate randomize positions is unity. 

In order to generate a randomized configuration of 1024 points in surface of sphere with radius one, the user can execute

```bash
./bin/exe_start 1024 sph data_sph

```
All the outputs will be dumped inside data_sph.

---
**NOTE**

The simulation will run for 60000 Monte Carlo steps, which will take more than 2 hours
to complete. To proceed with the other sections it is sufficient to kill the
simulation after 1000 Monte Carlo steps. 

## Triangulating the surface

The surface of the membrane will respond to all the internal and external forces
present. The internal ones are due to bending or stretching. If the user has read
the section xx.xx of the document, then evaluation of these forces requires the
knowledge of all the neighbours of a point. For this purpose we rely on the python
wrapper for [qhull](qhull.org), known as ConvexHull.  We provide a utility script,
`utils/gen_memc_conf.py` which takes the output from the `exe_start`, set all the
connections and write the output as `conf/dmemc_conf.h5`. The user can copy paste
the following command to set the connection for dump `data_sph/part_pos0003.bin`. 

```bash
python utils/gen_memc_conf.py data_sph/part_pos0003.bin.
```

## MeMC  

Once the connections are set, the executable `bin/exe_memc` can be used to
numerically study the nano-vesicles. The
binary takes two arguments:
1) Name of the parameter file from which all the physical parameters is read.
2) The folder name to store all the data.

The para file is not a format-free input. The variables has to be in the specified
format to be read correctly. We urge the user to copy the following block and paste
them identically to a plain text file for input. 

```text
## Membrane parameters
N	coef_bending	coef_stretching	coef_vol_expansion sp_curve
1024  2.50          25.00            1000000.00        2.0
radius	pos_bot_wall	sigma	epsilon   theta_attractive
1.00     -1.05          0.17     4.00     0.52
## Montecarlo parameters
Dfac	kBT mc_total_iters	mc_dump_iter
4       1.00  50000            100
## Afm Tip parameters
tip_radius	tip_pos_z	afm_sigma	afm_epsilon
0.20         1.05       0.17         4.00
`````

Values of the parameter are stored directly below the name. For instance the number
of points used to represent the membrane above is 1024 (The number below N). The other parameters which can be varied are:

* **Membrane specific parameter**
    +  coef_bending : 
    +  coef_stretching : the Young's modulus doc/paper.pdf
    +  coef_vol_expansion : the bulk's modulus in doc/paper.pdf
    +  sp_curve :: spontaneous curvature. 
    +  radius :: radius (We should probably remove this ) 

* **Parameters for the bottom stick wall**
    + pos_bot_wall :: z-position of the bottom wall. The value must be smaller than the radius
    + sigma :: The $\sigma$ of the bottom LJ potential
    + epsilon :: Relative strength LJ potential
    + theta_attractive :: All the points for which $\theta$ (see fig 11) is less this value will be affected by the attractive surface.

* **Monte Carlo Parameters**
    + Dfac :: Each  
    + kBT ::
    + mc_total_iters ::
    + mc_dump_iters ::
* **AFM parameters**
    + tip_radius ::
    + tip_pos_z ::
    + afm_sigma ::
    + afm_epsilon ::


stored and the folder to write the simulation
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

Before equilibration      |  After equilibration
:-------------------------:|:-------------------------:
![](./doc/figs/surf_mc_random.png)   |  ![](./doc/figs/surf_mc_lattice.png)


