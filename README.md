# MeMC

MeMC is an open-source Monte Carlo package for simulating elastic membranes. It is designed to study the mechanics of biological nano-vesicles (e.g. exosomes) where thermal fluctuations renormalize elastic constants. The package supports both spherical shells and flat membranes, fluid and tethered topology, and AFM indentation experiments.

For the physical background see `paper/paper.pdf`.

---

## Repository layout

```
MeMC/
├── start/       — generates equilibrated initial positions (LJ randomisation)
├── memc/        — main Monte Carlo simulation
├── utils/       — pre/post-processing Python scripts
├── docs/        — API documentation (Doxygen)
└── paper/       — manuscript and figures
```

---

## Prerequisites

| Dependency | Notes |
|---|---|
| `mpic++` | Tested with GCC 5.4+ on x86-64 |
| `gfortran` | Required for Fortran namelist reader |
| HDF5 | System package or manual install |
| Python 3.8+ | `numpy`, `scipy`, `h5py`, `numpy-quaternion` |

On Ubuntu/Debian:

```bash
apt install g++ gfortran libhdf5-dev mpich
pip install numpy scipy h5py numpy-quaternion
```

---

## Building

Each sub-package has its own Makefile.

```bash
# Build the start executable
cd start && make
# Produces: start/exe_start

# Build the main MC executable
cd memc && make
# Produces: memc/bin/exe_memc
```

---

## Workflow

### 1. Generate an initial configuration

`start/exe_start` places N points on a sphere or flat plane using a Lennard-Jones
Monte Carlo and outputs HDF5 snapshots.

```bash
cd start
./exe_start <N> <metric> <outfolder> <mc_steps>
```

| Argument | Description |
|---|---|
| `N` | Number of mesh nodes |
| `metric` | `sph` for sphere, `cart` for flat plane |
| `outfolder` | Directory for output snapshots |
| `mc_steps` | Number of MC steps for randomisation |

Example — 1024 nodes on a flat plane:

```bash
./exe_start 1024 cart run_flat 50000
```

Snapshots are written as `run_flat/snap_00000.h5`, `snap_00001.h5`, ...

---

### 2. Build the mesh connectivity

`utils/makeinput.py` reads positions from a start snapshot, triangulates the
surface and writes `input.h5` (positions + neighbour list) into the same folder.
It auto-detects spherical vs flat coordinates.

```bash
python utils/makeinput.py <snapshot.h5> <outfolder>
```

Example:

```bash
python utils/makeinput.py start/run_flat/snap_00099.h5 start/run_flat
# writes: start/run_flat/input.h5
#         start/run_flat/input.vtk
```

Use the last (most equilibrated) snapshot as input.

---

### 3. Set up a simulation folder

`exe_memc` expects a zero-padded folder (e.g. `00000/`) containing:
- `input.h5` — mesh from step 2
- `para_file.in` — simulation parameters (Fortran namelist format)

```bash
mkdir -p sim/00000
cp start/run_flat/input.h5 sim/00000/
cp <para_file.in>           sim/00000/
```

#### Parameter file format

Parameters are read as Fortran namelists. A minimal `para_file.in`:

```fortran
&Celllist
isselfrepulsive = F
boxsize = 10.0
cutoff = 1.5
minL = 0.0
sig = 1.0
eps = 0.0
/

&Meshpara
N = 1024
nghst = 12
bdry_cdt = 1        ! 0=channel  1=frame  other=periodic
radius = 1.0        ! set to 1.0 for flat membranes
/

&Membrane
N = 1024
coef_bend = 4.0
YY = 1e4
radius = 1.0
bdry_type = 1
/

&Bendpara
coef_bend = 4.0
minC = 0.0
maxC = 0.0
theta = 0.0
spcurv = 0.0
/

&spcurvpara
spcurv_which = 'none'
minC = 0.0
maxC = 0.0
theta = 0.0
/

&StickPara
pos_bot_wall = -1.0
sigma = 1.0
eps1 = 0.0
eps2 = 0.0
/

&mcpara
Mcalgo = 'mpolis'
dfac = 16.0
kbt = 1.0
is_restart = F
tot_mc_iter = 5000
dump_skip = 100
is_fluid = F
min_allowed_nbr = 4
fac_len_vertices = 1.3
fluidize_every = 10
/

&springpara
do_spring = F
icompute = 0
nPole_eq_z = 0.0
sPole_eq_z = 0.0
/

&afmpara
do_afm = F
tip_rad = 1.0
tip_pos_z = 0.0
sigma = 1.0
epsilon = 0.0
/

&Actpara
act_which = 'none'
doactivity = F
maxA = 0.0
minA = 0.0
/

&fluidpara
is_fluid = F
min_allowed_nbr = 4
fluidize_every = 10
fac_len_vertices = 1.3
/

&Volpara
do_volume = F
is_pressurized = F
coef_vol_exp = 0.0
pressure = 0.0
/

&Areapara
do_area = F
coef_area_exp = 0.0
/

&Stretchpara
YY = 1e4
do_volume = F
is_pressurized = F
is_pressure_ideal = F
coef_vol_expansion = 0.0
pext = 0.0
pint = 0.0
coef_area_expansion = 0.0
do_area = F
/
```

Key parameters:

| Namelist | Parameter | Description |
|---|---|---|
| `Meshpara` | `bdry_cdt` | 0=channel, 1=frame (flat), other=periodic |
| `Meshpara` | `radius` | Set to 1.0 for flat membranes |
| `Bendpara` | `coef_bend` | Bending rigidity κ |
| `Stretchpara` | `YY` | 2D Young's modulus |
| `mcpara` | `dfac` | Step-size divisor (larger = smaller steps) |
| `mcpara` | `is_fluid` | Enable bond-flip moves for fluid membrane |
| `mcpara` | `fluidize_every` | Bond-flip sweep every N MC iterations |
| `Volpara` | `do_volume` | Enable volume constraint |
| `afmpara` | `do_afm` | Enable AFM tip potential |

---

### 4. Run the simulation

`exe_memc` takes a single integer argument — the folder index. With MPI:

```bash
cd sim
mpirun -n <nproc> ../memc/bin/exe_memc <start_index>
```

With `start_index=0` and 1 MPI rank the code runs in `sim/00000/`. With 4 ranks
it runs simultaneously in `00000/`, `00001/`, `00002/`, `00003/`.

Example:

```bash
cd sim
mpirun -n 1 ../memc/bin/exe_memc 0
```

---

## Output

All output files are written inside each numbered folder (e.g. `00000/`).

| File | Description |
|---|---|
| `mc_log` | Per-iteration energies (see below) |
| `terminal.out` | Acceptance rate and total energy summary |
| `snap_NNNNN.h5` | Position + connectivity snapshots |
| `restartindex.txt` | Last completed iteration (used for restart) |
| `*.out` | Parameter echo files (`bendpara.out`, `mcpara.out`, …) |

#### `mc_log` columns

```
iter  accepted_frac  bend_e  stretch_e  stick_e  [pressure_e]  [selfrep_e]  total_e  volume  area  volt0
```

Optional columns `pressure_e` and `selfrep_e` appear only when the respective
features are enabled.

---

## Restarting

Set `is_restart = T` in `mcpara`. The code reads the last snapshot recorded in
`restartindex.txt` and continues from that point.

---

## Visualization

Convert HDF5 snapshots to VTK for [ParaView](https://www.paraview.org) or
[VisIt](https://visit-dav.github.io/visit-website/):

```bash
python utils/viz_memc.py <snap.h5> <output.vtk>
```

---

## Utilities

| Script | Purpose |
|---|---|
| `utils/makeinput.py` | Build `input.h5` from a start snapshot (auto-detects flat/sphere) |
| `utils/viz_memc.py` | Convert snapshot to VTK |
| `utils/check_status.py` | Diagnose bad runs |
| `utils/check_plot.py` | Plot energy traces |
| `utils/paramio.py` | Read/write parameter files |

---

## Developers

API documentation: [https://vipinagrawal25.github.io/MeMC/](https://vipinagrawal25.github.io/MeMC/)
