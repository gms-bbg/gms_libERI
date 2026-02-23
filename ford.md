---
src_dir: src
         app
exclude_dir: src/rhf/rot_axis
             src/rhf/eric
             src/rhf/rys
exclude: src/liberi_boys.F90
output_dir: docs
project: libERI
summary: GPU-accelerated two-electron repulsion integral library for quantum chemistry
project_github: https://github.com/gamess-bbg/gamess-libERI
author: libERI developers
github: https://github.com/JorgeG94
graph: true
graph_maxnodes: 200
graph_maxdepth: 4
coloured_edges: true
display: public
         protected
source: true
proc_internals: false
sort: permission-alpha
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            mpi:https://www.mpi-forum.org/docs/
            mdi:https://molssi-mdi.github.io/MDI_Library/html/index.html
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
               markdown.extensions.tables
---

[TOC]

libERI API Documentation
========================

This is the API documentation for **libERI**, a GPU-accelerated two-electron repulsion integral (ERI) library for quantum chemistry.

## Overview

libERI computes two-electron integrals and builds Fock matrices on GPUs using OpenMP target offload. It can be used:

1. **As a standalone library** - directly call the integral driver from your Fortran code
2. **As an MDI plugin** - integrate with GAMESS via the MolSSI Driver Interface

## Module Organization

### Core Modules

| Module | Description |
|--------|-------------|
| [[liberi_types]] | Kind parameters for portable integer and real types |
| [[liberi_constants]] | Fundamental mathematical constants |
| [[parameters]] | Derived constants, screening thresholds |
| [[liberi_globals]] | Global arrays for basis set, density, and Fock matrices |

### Computation Modules

| Module | Description |
|--------|-------------|
| [Boys Function](#boys-function-module) | Boys function tabulated grid for integral evaluation |
| [[liberi_shell_pair]] | Shell-pair precomputation |
| [[liberi_integral_driver]] | Top-level integral dispatch |
| [Integral Kernels](#integral-kernels-reference) | Machine-generated integral kernels (rot_axis, eric, rys) |

### Interface Modules

| Module | Description |
|--------|-------------|
| [[mdi_api]] | MDI engine interface for GAMESS integration |

## Quick Start

### Standalone Usage

```fortran
use liberi_globals
use liberi_shell_pair, only: shell_pair
use liberi_integral_driver, only: integral_driver

! Set MPI parameters (0, 1 for serial)
n_rank = 0
n_size = 1

! Populate basis set arrays...
! (natoms, nsh, num_bas, ang_mom, exponents, etc.)

! Allocate and set density matrix
allocate(density(num_bas*(num_bas+1)/2))
allocate(fa(num_bas*(num_bas+1)/2))

! Compute integrals
call shell_pair()       ! Precompute shell pairs
call integral_driver()  ! Build Fock matrix

! Results are in fa(:)
```

### MDI Plugin Usage

When compiled with `-DUSE_MDI`, libERI can be loaded as an MDI engine plugin by GAMESS. The driver sends basis set data and density matrices via MDI commands, and libERI returns the computed Fock matrix.

## Integral Kernels

libERI uses three types of integral kernels:

### Rotated-Axis Kernels

For s, p, d shell combinations using the Head-Gordon & Pople method. These are the default kernels for basis sets without f-functions.

### ERIC Kernels

McMurchie-Davidson style kernels for f-function quartets. Requires `ENABLE_F` compile flag.

### Rys Quadrature Kernels

Machine-generated kernels using Rys polynomial quadrature for higher angular momentum. Requires `ENABLE_F` compile flag.

## Building Documentation

Generate this documentation with FORD:

```bash
ford ford.md
```

The output will be in the `docs/` directory.

## Boys Function Module

**Location:** `src/liberi_boys.F90`

The `boys` module provides pre-tabulated Boys function values for efficient evaluation during integral computation. The Boys function \( F_m(t) \) is defined as:

$$ F_m(t) = \int_0^1 u^{2m} e^{-tu^2} du $$

### Public Interface

| Name | Type | Description |
|------|------|-------------|
| `t_max` | parameter | Maximum t value (25.0) for tabulated grid; use asymptotic form above |
| `m_increment` | parameter | Grid spacing increment |
| `f_increment` | array(8) | Increment factors for different angular momentum orders |
| `boys_grid_zero` | array(18040) | Tabulated Boys function values |
| `exponent_grid` | array | Exponent grid for interpolation |

This module is ~7,000 lines (mostly tabulated data) and excluded from FORD source parsing.

---

## Integral Kernels Reference

The integral kernels are machine-generated Fortran subroutines excluded from the main API documentation due to their size.

### rot_axis_kernels (`src/rhf/rot_axis/`)

Rotated-axis kernels for s, p, d shells using Head-Gordon & Pople method:

| Kernel | Quartet | Kernel | Quartet | Kernel | Quartet |
|--------|---------|--------|---------|--------|---------|
| `int0000` | (ss\|ss) | `int0001` | (ss\|sp) | `int0011` | (ss\|pp) |
| `int0101` | (sp\|sp) | `int0111` | (sp\|pp) | `int1111` | (pp\|pp) |
| `int0002` | (ss\|sd) | `int0012` | (ss\|pd) | `int0102` | (sp\|sd) |
| `int0112` | (sp\|pd) | `int0022` | (ss\|dd) | `int0122` | (sp\|dd) |
| `int1102` | (pp\|sd) | `int0202` | (sd\|sd) | `int2111` | (pd\|pp) |
| `int2120` | (pd\|pd) | `int2121` | (pd\|pd) | `int2211` | (dd\|pp) |
| `int2220` | (dd\|pd) | `int2221` | (dd\|pd) | `int2222gen` | (dd\|dd) |

### eric_kernels (`src/rhf/eric/`) - requires ENABLE_F

ERIC (McMurchie-Davidson) kernels for f-function quartets:

| Kernel | Quartet |
|--------|---------|
| `int0030` | (ss\|sf) |
| `int0031` | (ss\|pf) |
| `int0032` | (ss\|df) |
| `int1030` | (sp\|sf) |

### rys_kernels (`src/rhf/rys/`) - requires ENABLE_F

Rys quadrature kernels for higher angular momentum. Machine-generated, very large files.

### Naming Convention

`intABCD` where A,B,C,D are angular momentum: 0=s, 1=p, 2=d, 3=f

---

## References

- Head-Gordon, M.; Pople, J.A. "A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations" *J. Chem. Phys.* **89**, 5777 (1988)
- Rys, J.; Dupuis, M.; King, H.F. "Computation of electron repulsion integrals using the Rys quadrature method" *J. Comput. Chem.* **4**, 154 (1983)
- Boys, S.F. "Electronic wave functions. I." *Proc. R. Soc. London, Ser. A* **200**, 542 (1950)

{!README.md!}
