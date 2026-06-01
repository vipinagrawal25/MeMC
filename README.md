# MeMC — semisolid branch

MeMC is an open-source Monte Carlo package for simulating elastic membranes. It is designed to study the mechanics of biological nano-vesicles (e.g. exosomes) where thermal fluctuations renormalize elastic constants. The package supports both spherical shells and flat membranes, fluid and tethered topology, and AFM indentation experiments.

This branch (`semisolid`) extends the base code with:
- **Semisolid membranes** — a subset of nodes is pinned to a fixed triangulation (no bond flips), while the rest remains fluid. The solid node set is loaded from `solid_index.h5`.
- **Shear deformation** — an affine shear strain is applied to the initial configuration before the MC run begins. The frame is held fixed at the sheared positions throughout the simulation.

For the physical background see `paper/paper.pdf`.

---

## Repository layout

```
MeMC/
├── start/        — generates equilibrated initial positions (LJ randomisation)
├── memc/         — main Monte Carlo simulation (semisolid + shear)
│   └── old/      — reference: previous flat-struct implementation
├── utils/        — pre/post-processing Python scripts
├── docs/         — API documentation (Doxygen)
└── paper/        — manuscript and figures
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

## Quick-start workflow

**1. Compile `start`**
```bash
cd start && make
# produces: start/exe_start
```

**2. Generate initial positions**

Run the LJ Monte Carlo to place N points on a flat plane or sphere:
```bash
cd start
./exe_start <N> cart <outfolder> <mc_steps>
# example: ./exe_start 1024 cart run_flat 50000
# writes:  run_flat/snap_00000.h5 ... snap_000NN.h5
```
Use the last (most equilibrated) snapshot in the next step.

**3. Build the mesh connectivity**

Triangulate the point cloud and write `input.h5` (positions + neighbour list).
Auto-detects flat vs spherical coordinates:
```bash
python utils/makeinput.py start/run_flat/snap_00499.h5 start/run_flat
# writes: start/run_flat/input.h5
```

**4. (Semisolid only) Generate solid node indices**

Pick `num_solid_points` random bulk nodes and write `solid_index.h5`:
```bash
python utils/makesolidpts.py start/run_flat/input.h5 <num_solid_points> --bdry_type 1 --seed 42
# writes: start/run_flat/solid_index.h5
```
Note `num_solid_points` is hardcoded as 2048 in the main cpp code. Thus choose 2048 or change in the code. 

**5. Compile `memc`**
```bash
cd memc && make
# produces: memc/bin/exe_memc
```

**6. Set up a run folder**

`exe_memc` reads from a zero-padded folder (e.g. `memc/00000/`):
```bash
mkdir -p memc/00000
cp start/run_flat/input.h5     memc/00000/
cp start/run_flat/solid_index.h5  memc/00000/   # semisolid only
cp <your_para_file.in>         memc/00000/para_file.in
```

**7. Run the simulation**
```bash
mpirun -n <nproc> memc/bin/exe_memc <start_index>
# example (single process, folder 00000/):
mpirun -n 1 memc/bin/exe_memc 0
```

With `nproc` MPI ranks and `start_index=K`, the ranks run in folders `K/`, `K+1/`, ..., `K+nproc-1/` simultaneously.

---

## Parameter file

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
bdry_cdt = 1        ! 0=channel  1=frame (flat)  other=periodic
radius = 1.0        ! set to 1.0 for flat membranes
/

&Bendpara
coef_bend = 4.0
minC = 0.0
maxC = 0.0
theta = 0.0
spcurv = 0.0
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
tot_mc_iter = 1000
dump_skip = 1000
is_fluid = F
is_semisolid = F
num_solid_points = 0
min_allowed_nbr = 4
fac_len_vertices = 1.3
fluidize_every = 10
/

&shearpara
do_shear = F
slope = 0.0
constant = 0.0
shear_every = 0
/

&Actpara
act_which = 'none'
doactivity = F
maxA = 0.0
minA = 0.0
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
| `mcpara` | `is_semisolid` | Enable semisolid mode (requires `solid_index.h5`) |
| `mcpara` | `num_solid_points` | Number of solid (non-flipping) nodes |
| `mcpara` | `fluidize_every` | Bond-flip sweep every N MC iterations |
| `shearpara` | `do_shear` | Apply affine shear to initial config |
| `shearpara` | `slope` | Shear strain γ (Δx per unit y) |
| `shearpara` | `constant` | Frame spring constant (unused in set-and-run mode) |

#### Semisolid namelist

```fortran
&mcpara
...
is_fluid        = T
is_semisolid    = T
num_solid_points = 200
fluidize_every  = 10
/
```

#### Shear namelist

```fortran
&shearpara
do_shear   = T
slope      = 0.1
constant   = 0.0
shear_every = 0
/
```

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
| `utils/makesolidpts.py` | Generate `solid_index.h5` for semisolid simulations |
| `utils/viz_memc.py` | Convert snapshot to VTK |
| `utils/check_status.py` | Diagnose bad runs |
| `utils/check_plot.py` | Plot energy traces |
| `utils/paramio.py` | Read/write parameter files |

#### Generating solid points

```bash
python utils/makesolidpts.py 00000/input.h5 200 --bdry_type 1 --seed 42
# writes: 00000/solid_index.h5
```

Nodes in the frame and their direct neighbours are excluded from selection.
The remaining `num_solid_points` nodes are chosen at random.

---

## Developers

API documentation: [https://vipinagrawal25.github.io/MeMC/](https://vipinagrawal25.github.io/MeMC/)
