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

```bash
./bin/exe_start 1024 sph data_sph

```

## The turn (Generating the connections)

For MeMC, one needs to triangulate the points generated as described in para above.
We rely on the python libraries to do so. The code to do the trick is  utils/gen_memc_conf.py, which take the data to be triangulated as an input.  Executing the following

```bash
python utils/gen_memc_conf.py data_sph/part_pos0003.bin 
```
will write the input file with all the connections in conf/dmemc_conf.h5.  


The surface config shown below::

![plot](./doc/figs/surf_mc_random.png)

After 60000 monte-carlo steps become::

![plot](./doc/figs/surf_mc_lattice.png)


## The prestige  
