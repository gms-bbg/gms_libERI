# libERI

GPU-accelerated two-electron repulsion integral (ERI) library for quantum chemistry, designed as an engine plugin for [GAMESS](https://www.msg.chem.iastate.edu/gamess/). libERI offloads the computationally dominant Fock-matrix build to GPUs using OpenMP target offload, achieving significant speedups over conventional CPU-based integral evaluation.

## Overview

In Hartree-Fock and DFT calculations the two-electron integral step is the bottleneck. libERI replaces this step by computing all required shell-quartet classes entirely on the GPU. It provides a handle-based Fortran API with zero global state, and can integrate with GAMESS either through direct linking or as an [MDI]([https://molssi-mdi.github.io/MDI_Library/html/index.html](https://github.com/MolSSI-MDI/MDI_Library/tree/master)) engine plugin.

### Supported integral classes

Integrals are organized by the angular momentum of the four shells in a quartet (s=0, p=1, d=2, f=3):

| Method | Shell types | Kernel count | Description |
|---|---|---|---|
| **Rotated-axis** | sp (ss through dd) | ~23 | McMurchie-Davidson rotated-axis method for low angular momentum |
| **ERIC** | f-containing (sf, pf, ...) | 7 | Obara-Saika style for f-function quartets |
| **Rys quadrature** | dd and higher (generated) | ~45 | Rys-polynomial based evaluation for higher angular momentum |

F-function integrals (ERIC + Rys) are optional and gated behind a compile flag since they significantly increase build time.

### Key features

- GPU offload via OpenMP target directives (portable across vendors)
- Schwarz screening for integral prescreening
- Shell-pair precomputation for efficient quartet formation
- Boys function evaluation using tabulated grid interpolation
- MPI-parallel with optional MDI plugin architecture
- Supports RHF 

## Project structure

```
libERI/
  CMakeLists.txt              # Top-level build configuration
  src/
    liberi_types.f90          # Kind parameters, basis_t, shell_pair_t, shell_pair_container_t, eri_resources_t
    liberi_parameters.F90     # Physical/numerical constants, Rys quadrature weights
    liberi_boys.F90           # Boys function tabulated grid data (~18k points)
    liberi_shell_pair.F90     # Shell-pair precomputation into shell_pair_container_t
    liberi_driver.F90         # compute_integrals(pairs, density, fock) — integral dispatch
    interface/
      liberi_interface.F90    # Handle-based API: create/setup/fock_build/cleanup/destroy
      gms_liberi_wrapper.F90  # GAMESS direct-link wrapper (liberi_init/fock/finalize)
      liberi_mdi_api.F90      # MDI engine plugin (local mdi_state_t, no globals)
    rhf/
      rot_axis/               # Rotated-axis integral kernels (s, p, d shells)
        rot_axis_kernels.F90  # Module interface (Fortran submodules)
        int0000.F90           # (ss|ss) kernel
        int0001.F90           # (ss|sp) kernel
        ...                   # 21 kernels total
      eric/                   # ERIC integral kernels (f-function, optional)
        eric_kernels.F90      # Module interface
        int0030_ericgen.F90   # (ss|sf) kernel
        ...                   # 7 kernels total
      rys/                    # Rys quadrature kernels (generated, optional)
        rys_kernels.F90       # Module interface
        int3333_rysgen.F90    # (ff|ff) kernel (~39k lines, machine-generated)
        ...                   # ~45 kernels total
  cmake/
    CMakeLists.txt            # Compiler flag auto-detection
    modules/                  # FindMDI, FindLIBERI, status summary
  tests/
    fortran/                  # Fortran regression tests (H2, H2O)
    test_plugin.py            # MDI plugin unit tests (pytest + mdi4py)
    general_tests/            # GAMESS integration tests (exam01-exam06)
```

## Building

### Requirements

- CMake >= 3.22
- MPI (any implementation)
- OpenMP-capable Fortran compiler (with GPU offload support for GPU acceleration)
- [MDI Library](https://github.com/MolSSI-MDI/MDI_Library) (optional, for MDI plugin mode)

### Supported compilers

| Compiler | GPU offload | Notes |
|---|---|---|
| **NVHPC** (`nvfortran`) | NVIDIA GPUs | Primary target; best performance |
| **Cray** (`ftn`) | AMD GPUs | Used on Frontier |
| **Intel** (`ifx`) | Intel GPUs | Used on Sunspot / Aurora |
| **GNU** (`gfortran`) | None | Compiles but runs serially; poor performance |

### Quick start (CMake)

```bash
mkdir build && cd build
cmake ..
make
```

This produces `libERI.so` (shared library when MDI is enabled) or `libERI.a` (static library otherwise).

### CMake options

| Option | Default | Description |
|---|---|---|
| `ERI_GPU_ARCH` | `cc70` | NVIDIA compute capability (e.g. `cc70`, `cc80`) |
| `ERI_ENABLE_GPU` | `ON`* | Enable GPU offloading flags (*OFF for gfortran) |
| `ERI_ENABLE_F` | `OFF` | Enable f-function integrals (significantly increases compile time) |
| `ERI_ENABLE_MDI` | `ON` | Enable MDI library support for GAMESS integration |
| `ERI_USE_MPI` | `ON` | MPI support in MDI (only used if `ERI_ENABLE_MDI=ON`) |
| `ERI_ENABLE_TESTING` | `ON` | Build regression tests |
| `ERI_ENABLE_BENCHMARKS` | `OFF` | Build performance benchmarks |

### Example configurations

```bash
# Basic build (no MDI, no GPU)
cmake ..

# With MDI support for GAMESS integration
cmake -DERI_ENABLE_MDI=ON ..

# NVIDIA GPU (V100)
cmake -DERI_GPU_ARCH=cc70 ..

# NVIDIA GPU (A100)
cmake -DERI_GPU_ARCH=cc80 ..

# With f-function integrals (increases compile time)
cmake -DERI_ENABLE_F=ON ..

# Full build with MDI, GPU, and benchmarks
cmake -DERI_ENABLE_MDI=ON -DERI_GPU_ARCH=cc80 -DERI_ENABLE_BENCHMARKS=ON ..
```

## Usage

### Handle-based API

The primary API uses an opaque handle that owns all state (basis data, shell pairs, GPU mappings). There are no global variables anywhere in the library.

```fortran
use liberi_types, only: dp
use liberi_interface, only: liberi_handle_t, liberi_create, liberi_setup, &
                            liberi_fock_build, liberi_cleanup, liberi_destroy

type(liberi_handle_t) :: handle
real(dp) :: density(n_tri), fock(n_tri)

! 1. Create handle
call liberi_create(handle)

! 2. Set up basis and compute shell pairs (once per geometry)
call liberi_setup(handle, nsh, natoms, num_bas, mxgtot, &
                  ang_mom, contr_num, sh_loc, atom_num, atom_loc, &
                  start_bas, end_bas, exponents, &
                  contr_coef_s, contr_coef_p, contr_coef_d, contr_coef_f, &
                  coords, schwrz_int, my_rank, num_procs)

! 3. Build Fock matrix (every SCF iteration)
call liberi_fock_build(handle, density, fock, n_tri)

! 4. Clean up
call liberi_cleanup(handle)
call liberi_destroy(handle)
```

See `tests/fortran/` for complete working examples.

### Direct GAMESS integration

The `gamess_liberi_wrapper` module provides a simplified three-call interface for GAMESS with a module-level `save`d handle:

```fortran
use gamess_liberi_wrapper, only: liberi_init, liberi_fock, liberi_finalize

call liberi_init(...)       ! Once per geometry (before SCF loop)
call liberi_fock(density, fock, n)  ! Every SCF iteration
call liberi_finalize()      ! After convergence
```

Build libERI as a static library (`-DERI_ENABLE_MDI=OFF`) and link with GAMESS. The GAMESS link step **must** include `-mp=gpu -gpu=<arch> -gpu=mem:separate`.

### As an MDI plugin

When built with `-DERI_ENABLE_MDI=ON`, libERI can run as an MDI engine plugin, loaded by GAMESS at runtime. The MDI API uses a local `mdi_state_t` and calls the same handle API internally.

**Inputs received from driver (GAMESS):**
- Basis set info: atom coordinates, shell angular momenta, exponents, contraction coefficients
- Schwarz screening integrals
- AO density matrix

**Outputs sent to driver:**
- Fock matrix

## Testing

### Fortran regression tests

```bash
# Via CTest (from build directory):
ctest

# Directly:
mpirun -np 1 ./test_h2_sonly
mpirun -np 1 ./test_h2_ccpvdz
mpirun -np 1 ./test_h2o_ccpvdz

# Regenerate reference values:
mpirun -np 1 ./test_h2_sonly --generate
```

Tests use the handle API with local arrays (no globals, no MDI, no GAMESS) and validate Fock matrix results against hardcoded reference values for small molecules (H2, H2O) with various basis sets.


## Benchmarks

Performance benchmarks are available for C2H6 and C60 molecules:

```bash
cmake -DERI_ENABLE_BENCHMARKS=ON ..
make
mpirun -np 1 ./bench_c2h6
mpirun -np 1 ./bench_c60
```


## Contributors

- **Melisa Alkan** (NVIDIA) — Original idea and implementation, developed through a [MolSSI](https://molssi.org/) fellowship
- **Taylor Barnes** (MolSSI) — MDI integration and plugin architecture
- **Jeff Hammond** (NVIDIA) — Technical guidance and support
- **Daniel Del Angel** (Iowa State University) — F-function integrals (PhD project)
- **Jorge Galvez** (Australian National University) — Software engineering and HPC consulting

With heavy relation and sponsorship through the GAMESS-ECP project under Professor Mark S. Gordon 
at Iowa State University / Ames Laboratory.
