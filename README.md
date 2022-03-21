# MeMC
A Monte-Carlo code for the simulation of fluctuating enclosed membranes. Such system can be
of relevance for the study of viruses, exosomes, e.t.c.

## Installation

If all the prequisites are satisfied, the installation is easy. The following
commands should do the job:

```bash
git clone https://github.com/vipinagrawal25/MeMC
cd MeMC
make
```
If successful, one should find an extra "bin" directory in the folder. The bin
directory will contain the binaries "exe_start" and "exe_memc" inside

## Prequisites

The code requires the following:

1) A c++ compiler. We have tested the code against gnu [g++](https://gcc.gnu.org/) version 11.2.0 on x86_64
CPU.
2) [Hdf5](https://www.hdfgroup.org/solutions/hdf5) libraries for reading and writing
data.
3) [Python](https://www.python.org/) version 3.8 (NOTE:: Cross check versions we use) with [scipy](https://www.scipy.org), [numpy](https://www.numpy.org), [h5py](https://www.h5py.org) and
[numpy-quartenion](https://https://pypi.org/project/numpy-quaternion/) installed.

The installation instruction is provided on the pages linked with each package. In
Ubuntu,
```bash
sudo apt install g++ hdf5
```
should install both. The different python libraries can be installed using pip. For
example
```bash
pip install numpy-quartenion
```
installs numpy-quartenion library.


# Using the Code:

We shall now dive deeper and explain the different part of the code. As stated
previously, the main purpose of the MeMC is the Monte-Carlo simulation of enclosed
Membranes. For the details check the document on doc/paper.pdf. The usage can be
divided into three distinct part. Since the number is three, we shall borrow the
quotes from the movie "The Prestige" and call them "the pledge", "the turn" and "the
prestige"

##  The pledge (Starting the simulation)

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

