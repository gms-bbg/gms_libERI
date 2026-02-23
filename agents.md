# libERI

GPU-accelerated two-electron repulsion integral (ERI) library for quantum chemistry, written in Fortran with OpenMP target offloading. Integrates with GAMESS as either an MDI plugin or a directly linked library.

## Quick Reference

### Build

```bash
# NVIDIA GPU (primary target):
cmake -B build -DERI_GPU_ARCH=cc70 -DERI_ENABLE_MDI=OFF
cmake --build build

# With MDI plugin support:
cmake -B build -DERI_GPU_ARCH=cc70 -DERI_ENABLE_MDI=ON
cmake --build build

# Run tests:
cd build && ctest
# Or directly:  mpirun -np 1 ./test_h2_sonly
```

### CMake Options

| Option | Default | Description |
|---|---|---|
| `ERI_GPU_ARCH` | `cc70` | NVIDIA compute capability |
| `ERI_ENABLE_GPU` | ON (NVHPC) | OpenMP target offloading |
| `ERI_ENABLE_MDI` | ON | MDI plugin support (builds shared lib) |
| `ERI_ENABLE_F` | OFF | F-function integrals (eric + rys kernels) |
| `ERI_ENABLE_TESTING` | ON | Build test executables |
| `ERI_ENABLE_BENCHMARKS` | OFF | Build benchmark executables |

## Architecture

**Zero global mutable state.** `liberi_globals` has been deleted. All data flows through explicit arguments and opaque handles.

### Public API (`src/interface/liberi_interface.F90`)

```
liberi_create(handle)       -- allocate handle
liberi_setup(handle, ...)   -- once per geometry: store basis, compute shell pairs, map to GPU
liberi_fock_build(handle, density, fock, size)  -- every SCF iteration
liberi_cleanup(handle)      -- unmap GPU, deallocate
liberi_destroy(handle)      -- destroy handle
```

### GAMESS Wrapper (`src/interface/gms_liberi_wrapper.F90`)

Thin module with `save`d handle. Three entry points for GAMESS: `liberi_init`, `liberi_fock`, `liberi_finalize`.

### Data Flow

```
Caller --> liberi_setup() --> shell_pair() --> target enter data (GPU)
                                |
Caller --> liberi_fock_build() --> compute_integrals(pairs, density, fock)
                                       |
                                  target enter data (density, fock)
                                  call int0000(pairs%ss, density, fock, res)
                                  call int0001(pairs%ss, pairs%sp, density, fock, res)
                                  ... (62 kernels total)
                                  target exit data (fock)
```

### Key Types (`src/liberi_types.f90`)

- **`basis_t`** — basis set data (exponents, coefficients, shell metadata, coordinates, Schwarz integrals)
- **`shell_pair_t`** — precomputed primitive pair data for one shell-pair class (exponents, coefficients, pair locations)
- **`shell_pair_container_t`** — holds all 10 shell-pair types (ss, sp, pp, sd, pd, dd, sf, pf, df, ff) + `eri_resources_t`
- **`eri_resources_t`** — helper arrays for kernels (coord_sh, ia, contr_num, atom_loc, shell classification arrays, MPI rank/size, num_bas)
- **`liberi_handle_t`** — opaque handle owning `basis_t` + `shell_pair_container_t`

### Module Map

| Module | File | Role |
|---|---|---|
| `liberi_types` | `src/liberi_types.f90` | Kind params, all derived type definitions |
| `liberi_parameters` | `src/liberi_parameters.F90` | Physical/numerical constants, Rys weights, `cutoff_schwarz` |
| `liberi_boys` | `src/liberi_boys.F90` | Boys function grid data (3 static arrays) |
| `liberi_shell_pair` | `src/liberi_shell_pair.F90` | Shell-pair precomputation + GPU mapping |
| `liberi_integral_driver` | `src/liberi_driver.F90` | `compute_integrals(pairs, density, fock)` — single driver entry point |
| `liberi_interface` | `src/interface/liberi_interface.F90` | Handle-based API |
| `gamess_liberi_wrapper` | `src/interface/gms_liberi_wrapper.F90` | GAMESS direct-link wrapper |
| `mdi_api` | `src/interface/liberi_mdi_api.F90` | MDI plugin API (local `mdi_state_t`, no globals) |
| `rot_axis_kernels` | `src/rhf/rot_axis/rot_axis_kernels.F90` | Module interface for 21 rot_axis kernel submodules |

### Kernel Files

- `src/rhf/rot_axis/` — 21 rotated-axis kernels (14 hand-written + 7 rysgen)
- `src/rhf/eric/` — eric kernels (f-function, requires `ERI_ENABLE_F`)
- `src/rhf/rys/` — rys kernels (f-function, requires `ERI_ENABLE_F`)

All 62 kernels have the same signature pattern: `(bra_pair[, ket_pair], density, fock, res)`. No kernel uses any globals.

### GPU Data Mapping Strategy

| Data | Mapped In | Lifetime |
|---|---|---|
| Shell-pair arrays | `shell_pair()` via `target enter data` | Setup → `shell_pair_deallocate()` |
| `eri_resources_t` | `shell_pair()` | Setup → `shell_pair_deallocate()` |
| `boys_grid_zero`, `exponent_grid` | `shell_pair()` | Setup → `shell_pair_deallocate()` |
| `density`, `fock` | `compute_integrals()` | Per SCF iteration (`enter`/`exit`) |

Kernels use `shared()` only — zero `map()` clauses. All data is pre-mapped.

## Compiler: nvfortran (NVIDIA HPC SDK)

### Critical Gotchas

1. **`%LOC` reserved** — nvfortran uses `%LOC` as built-in address-of. Never name a derived-type component `loc`. Use `pair_loc` or similar.

2. **OMP `shared()` must use whole types** — `shared(ss_pair)` works; `shared(ss_pair%field)` is a syntax error.

3. **OMP `map(to:)` deep copy** — Must list parent AND each allocatable member: `map(to: pairs%ss, pairs%ss%d_coeff, pairs%ss%pair_loc, ...)`.

4. **OMP target name collision** — Inside `!$omp target default(none)`, nvfortran may resolve a local variable name to a global module variable of the same name, even if the global is not imported. Avoid naming subroutine arguments the same as global variables.

5. **Large `parameter` arrays crash `fort2`** — Arrays above ~18000 elements as `parameter` cause the compiler backend to segfault. `boys_grid_zero` (18040) and `exponent_grid` (5005) must use `data` statements or runtime init, not `parameter`.

6. **`-gpu=mem:separate` is required** — Without it, CUDA unified memory masks broken mappings. The code is validated with explicit mappings only.

## GAMESS Integration

### Direct Link (preferred, no MDI)

GAMESS calls `liberi_init`/`liberi_fock`/`liberi_finalize` from `gamess_liberi_wrapper`. GAMESS compiles with `-i8` (64-bit default integers) — all cross-library integer args must use `integer(c_int)` (32-bit) via `iso_c_binding`.

**Critical link flags**: GAMESS link step MUST include `-mp=gpu -gpu=<arch> -gpu=mem:separate`. Without these, OMP target regions fall back to host emulation, breaking derived-type deep copy.

### MDI Plugin (legacy)

`liberi_mdi_api.F90` receives basis data via MDI commands into a local `mdi_state_t`, calls the handle API. Three plugin launches per SCF job (setup, per-iteration fock build, cleanup).

### Basis Set Convention

- `nosp_basis%cc` is a UNIFIED contraction coefficient array. SP shells are split into separate S and P with duplicated primitives.
- Contraction coefficients are split into 4 arrays: `contr_coef_s`, `contr_coef_p`, `contr_coef_d`, `contr_coef_f` (only one non-zero per shell).
- Angular momentum: 1=s, 2=p, 3=d, 4=f.

## Tests

Tests under `tests/fortran/` use the handle API with local arrays (no globals, no MDI, no GAMESS).

| Test | Molecule | Basis | Kernels |
|---|---|---|---|
| `test_h2_sonly` | H2 | s-only | int0000 |
| `test_h2_ccpvdz` | H2 | cc-pVDZ | 6 s/p kernels |
| `test_h2o_ccpvdz` | H2O | cc-pVDZ | 15+ s/p/d kernels |

Run with `ctest` from build dir, or `mpirun -np 1 ./test_h2_sonly`. Use `--generate` flag to regenerate reference Fock values.

## Open Issues

- **Hardcoded magic constants** — Boys function `f_increment` values inlined in all kernels instead of using the array
- **Axis rotation duplication** — ~60 identical lines in every kernel, should be a device subroutine
- **`liberi_setup()` 20-argument signature** — should accept a `basis_t` struct, but requires GAMESS-side changes too
- **Runtime GPU memory query** — chunk sizing still uses compile-time `#ifdef HPC_PM`
