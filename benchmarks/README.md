# libERI Performance Benchmarks

Standalone Fortran benchmarks for timing `integral_driver()` across compilers and code changes. These are **not** ctest tests — they're manual benchmarks for performance comparisons.

## Build

```bash
cmake -DERI_ENABLE_BENCHMARKS=ON ..
make bench_c2h6 bench_c60
```

The benchmark executables link against `libERI.so` (the shared library), same as the test programs. They do not compile any library sources themselves.

## Run

```bash
mpirun -np 1 ./benchmarks/bench_c2h6
mpirun -np 1 ./benchmarks/bench_c60
```

Set `OMP_NUM_THREADS` to control thread count:

```bash
OMP_NUM_THREADS=4 mpirun -np 1 ./benchmarks/bench_c2h6
```

## Benchmarks

### `bench_c2h6` — Ethane (C2H6, 8 atoms)

Coordinates from libcint `fortran_time_c2h6.F90` (Bohr). Runs two basis sets sequentially:

| Basis   | Shells | Basis Functions | Angular Momentum | mxgtot |
|---------|--------|-----------------|------------------|--------|
| 6-31G   | 22     | 30              | s, p             | 6      |
| cc-pVDZ | 30     | 60              | s, p, d          | 8      |

**Shell breakdown (6-31G):**
- Carbon (2 atoms): 5 shells each — s(6), s(3), s(1), p(3), p(1). L-shells split into separate S and P.
- Hydrogen (6 atoms): 2 shells each — s(3), s(1)

**Shell breakdown (cc-pVDZ):**
- Carbon (2 atoms): 6 shells each — s(8), s(8), s(1), p(3), p(1), d(1). General contraction split into two separate s-shells with duplicated exponents.
- Hydrogen (6 atoms): 3 shells each — s(3), s(1), p(1)

### `bench_c60` — Fullerene (C60, 60 atoms)

Coordinates from libcint `time_c60.c` (Bohr). Runs three basis sets sequentially:

| Basis   | Shells | Basis Functions | Angular Momentum | mxgtot |
|---------|--------|-----------------|------------------|--------|
| STO-3G  | 180    | 300             | s, p             | 3      |
| 6-31G   | 300    | 540             | s, p             | 6      |
| cc-pVDZ | 360    | 900             | s, p, d          | 8      |

**Shell breakdown per Carbon atom:**

| Basis   | Shells per C | BF per C | Shell types                          |
|---------|-------------|----------|--------------------------------------|
| STO-3G  | 3           | 5        | s(3), s(3), p(3)                     |
| 6-31G   | 5           | 9        | s(6), s(3), s(1), p(3), p(1)         |
| cc-pVDZ | 6           | 15       | s(8), s(8), s(1), p(3), p(1), d(1)   |

## Implementation Details

### How it works

Each benchmark:
1. Calls `MPI_Init` (single rank)
2. For each basis set:
   - Populates `liberi_globals` (coords, ang_mom, contr_num, sh_loc, atom_num, atom_loc, start_bas, end_bas, exponents, contr_coef_s/p/d/f)
   - Sets `density = 0.0` (zero density — pure timing, no correctness check)
   - Sets `schwrz_int = 1.0d10` (disables Schwarz screening to measure raw integral time)
   - Calls `shell_pair()` to precompute shell pair data
   - Times `integral_driver()` with `omp_get_wtime()`
   - Prints summary (shells, basis functions, shell pair counts, wall time)
   - Deallocates all globals and shell-pair arrays before the next basis set
3. Calls `MPI_Finalize`

### Normalization

Contraction coefficients use GAMESS normalization, pre-multiplied at initialization:

```
s: N = (2*alpha/pi)^(3/4)
p: N = 2*sqrt(alpha) * (2*alpha/pi)^(3/4)
d: N = 4*alpha * (2*alpha/pi)^(3/4)
```

### Conventions

- **Angular momentum**: 1=s, 2=p, 3=d (libERI convention, 1-indexed)
- **Cartesian d-functions**: 6 components (xx, yy, zz, xy, xz, yz)
- **General contractions**: Split into separate shells with duplicated exponents
- **L-shells (SP)**: Split into separate S and P shells
- **Exponent storage**: Flat array indexed via `sh_loc`, stride = `mxgtot`

### Basis set sources

- **STO-3G**: Hehre, Stewart, Pople (1969)
- **6-31G**: Hehre, Ditchfield, Pople (1972) — exponents/coefficients from libcint `time_c2h6.c`
- **cc-pVDZ**: Dunning (1989) — exponents/coefficients from libcint `time_c60.c` and BSE

## Example Output

```
=== libERI Benchmark: C60 ===
Threads:  4

--- STO-3G ---
  Atoms:           60
  Shells:          180
  Basis functions:  300
  Shell pairs (ss): ...
  Shell pairs (sp): ...
  Shell pairs (pp): ...
  Wall time:          X.XXX s

--- 6-31G ---
  Atoms:           60
  Shells:          300
  Basis functions:  540
  Shell pairs (ss): ...
  Shell pairs (sp): ...
  Shell pairs (pp): ...
  Wall time:         XX.XXX s

--- cc-pVDZ ---
  Atoms:           60
  Shells:          360
  Basis functions:  900
  Shell pairs (ss): ...
  Shell pairs (sp): ...
  Shell pairs (pp): ...
  Shell pairs (sd): ...
  Shell pairs (pd): ...
  Shell pairs (dd): ...
  Wall time:        XXX.XXX s
```
