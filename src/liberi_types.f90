!! Kind parameters for libERI
module liberi_types
  !! Defines portable kind parameters for integer and real types.
  !!
  !! This module provides consistent type definitions across the library:
  !!
  !! - `int32`, `int64`: 32-bit and 64-bit integers from iso_fortran_env
  !! - `sp`: Single precision real (6 decimal digits)
  !! - `dp`: Double precision real (15 decimal digits)
  !!
  !! All computational kernels use `dp` for floating-point operations.
  use, intrinsic :: iso_fortran_env, only: int32, int64
  implicit none

  private

  public :: int32, int64

  integer, parameter, public :: sp = SELECTED_REAL_KIND(6, 37)
    !! Single precision real kind (approximately 6 decimal digits)
  integer, parameter, public :: dp = SELECTED_REAL_KIND(15, 307)
    !! Double precision real kind (approximately 15 decimal digits)

  integer, parameter, public :: default_real = dp

  type, public :: basis_t
    !! Basis set data container.
    !!
    !! Holds all shell metadata, exponents, contraction coefficients,
    !! coordinates, and Schwarz integrals needed for integral evaluation.
    integer :: nsh = 0
      !! Number of shells
    integer :: natoms = 0
      !! Number of atoms
    integer :: num_bas = 0
      !! Number of basis functions
    integer :: mxgtot = 0
      !! Maximum primitives per shell
    integer, allocatable :: ang_mom(:)
      !! Angular momentum per shell (1=s, 2=p, 3=d, 4=f)
    integer, allocatable :: contr_num(:)
      !! Contraction length per shell
    integer, allocatable :: sh_loc(:)
      !! Starting primitive index per shell
    integer, allocatable :: atom_num(:)
      !! Atom center per shell
    integer, allocatable :: atom_loc(:)
      !! Basis function offset per shell
    integer, allocatable :: start_bas(:)
      !! First basis function per shell
    integer, allocatable :: end_bas(:)
      !! Last basis function per shell
    real(dp), allocatable :: exponents(:)
      !! Gaussian exponents (mxgtot * nsh)
    real(dp), allocatable :: contr_coef_s(:)
      !! s-type contraction coefficients
    real(dp), allocatable :: contr_coef_p(:)
      !! p-type contraction coefficients
    real(dp), allocatable :: contr_coef_d(:)
      !! d-type contraction coefficients
    real(dp), allocatable :: contr_coef_f(:)
      !! f-type contraction coefficients
    real(dp), allocatable :: coords(:)
      !! Atomic coordinates (3 * natoms)
    real(dp), allocatable :: schwrz_int(:)
      !! Schwarz integrals (nsh * (nsh+1) / 2)
  end type basis_t

  type, public :: eri_resources_t
    real(dp), allocatable    :: coord_sh(:, :)
    integer(int64), allocatable :: ia(:)
    integer, allocatable     :: contr_num(:)
    integer, allocatable     :: atom_loc(:)
    integer(int64), allocatable :: i_s_shl(:), i_p_shl(:), i_d_shl(:), i_f_shl(:)
    integer(int64) :: n_s_shl = 0, n_p_shl = 0, n_d_shl = 0, n_f_shl = 0
    integer :: n_rank = 0, n_size = 1
    integer :: num_bas = 0
  end type eri_resources_t

  type, public :: shell_pair_t
    integer(int64) :: n_pairs = 0
    integer(int64), allocatable :: pair_loc(:)
    integer(int64), allocatable :: ismlp(:)
    real(dp), allocatable :: xints(:)
    real(dp), allocatable :: expon_a(:), expon_b(:)
    real(dp), allocatable :: t_alpha(:), t_beta(:)
    real(dp), allocatable :: t_expon_ab(:)
    real(dp), allocatable :: t_inverse_expon_ab(:)
    real(dp), allocatable :: d_coeff(:)
    real(dp), allocatable :: d_coeff_alt(:)
    real(dp), allocatable :: sq(:)
  end type shell_pair_t

  type, public :: shell_pair_container_t
    !! Container holding all shell pair types for integral evaluation.
    !! This encapsulates all precomputed shell-pair data and is passed
    !! to the integral driver and kernels.
    type(shell_pair_t) :: ss
    type(shell_pair_t) :: sp
    type(shell_pair_t) :: pp
    type(shell_pair_t) :: sd
    type(shell_pair_t) :: pd
    type(shell_pair_t) :: dd
    type(shell_pair_t) :: sf
    type(shell_pair_t) :: pf
    type(shell_pair_t) :: df
    type(shell_pair_t) :: ff
    type(eri_resources_t) :: resources
  end type shell_pair_container_t

    !! Default real kind used throughout the library
end module liberi_types
