!! GAMESS wrapper for direct libERI integration (no MDI)
!!
!! This module provides a drop-in replacement for the MDI-based interface.
!! It calls liberi_interface directly, eliminating the plugin overhead.
!!
!! Usage in GAMESS (replace MDI calls):
!!   - Call LIBERI_INIT() once before SCF loop
!!   - Call LIBERI_FOCK() each SCF iteration
!!   - Call LIBERI_FINALIZE() after convergence
!!
!! Build libERI as a static library and link with GAMESS.
module gamess_liberi_wrapper
  use liberi_interface
  use liberi_types, only: dp
  implicit none
  private

  public :: liberi_init
  public :: liberi_fock
  public :: liberi_finalize

  ! Module-level handle persists across calls
  type(liberi_handle_t), save :: handle
  logical, save :: is_initialized = .false.
  integer, save :: matrix_size = 0

contains

  subroutine liberi_init(natoms, num_bas, nsh, mxgtot, &
                         ang_mom, contr_num, sh_loc, &
                         atom_num, atom_loc, start_bas, end_bas, &
                         exponents, contr_coef_s, contr_coef_p, &
                         contr_coef_d, contr_coef_f, &
                         coords, schwrz_int, my_rank, num_procs)
    !! Initialize libERI with basis set data.
    !!
    !! Call this once per geometry, before the SCF iteration loop.
    !! Equivalent to MDI_COMPUTE_SHELL_PAIR in the old interface.
    integer, intent(in) :: natoms
      !! Number of atoms
    integer, intent(in) :: num_bas
      !! Number of basis functions
    integer, intent(in) :: nsh
      !! Number of shells
    integer, intent(in) :: mxgtot
      !! Maximum primitives per shell
    integer, intent(in) :: ang_mom(nsh)
      !! Angular momentum per shell (1=s, 2=p, 3=d, 4=f)
    integer, intent(in) :: contr_num(nsh)
      !! Contraction length per shell
    integer, intent(in) :: sh_loc(nsh)
      !! Starting primitive index per shell
    integer, intent(in) :: atom_num(nsh)
      !! Atom center per shell
    integer, intent(in) :: atom_loc(nsh)
      !! Basis function offset per shell
    integer, intent(in) :: start_bas(nsh)
      !! First basis function per shell
    integer, intent(in) :: end_bas(nsh)
      !! Last basis function per shell
    real(dp), intent(in) :: exponents(mxgtot * nsh)
      !! Gaussian exponents
    real(dp), intent(in) :: contr_coef_s(mxgtot * nsh)
      !! s-type contraction coefficients
    real(dp), intent(in) :: contr_coef_p(mxgtot * nsh)
      !! p-type contraction coefficients
    real(dp), intent(in) :: contr_coef_d(mxgtot * nsh)
      !! d-type contraction coefficients
    real(dp), intent(in) :: contr_coef_f(mxgtot * nsh)
      !! f-type contraction coefficients
    real(dp), intent(in) :: coords(3 * natoms)
      !! Atomic coordinates (packed as x1,y1,z1,x2,...)
    real(dp), intent(in) :: schwrz_int((nsh * (nsh + 1)) / 2)
      !! Schwarz integrals
    integer, intent(in) :: my_rank
      !! MPI rank (0-based)
    integer, intent(in) :: num_procs
      !! Number of MPI processes

    ! Clean up if previously initialized
    if (is_initialized) then
      call liberi_cleanup(handle)
      call liberi_destroy(handle)
    end if

    ! Create handle and set up
    call liberi_create(handle)
    call liberi_setup(handle, nsh, natoms, num_bas, mxgtot, &
                      ang_mom, contr_num, sh_loc, &
                      atom_num, atom_loc, start_bas, end_bas, &
                      exponents, &
                      contr_coef_s, contr_coef_p, &
                      contr_coef_d, contr_coef_f, &
                      coords, schwrz_int, &
                      my_rank, num_procs)

    matrix_size = (num_bas * (num_bas + 1)) / 2
    is_initialized = .true.

  end subroutine liberi_init


  subroutine liberi_fock(density, fock, n)
    !! Compute Fock matrix contribution from density matrix.
    !!
    !! Call this every SCF iteration.
    !! Equivalent to MDI_DRIVE_PLUGIN_GEN in the old interface.
    !!
    !! After this call, caller (GAMESS) must reduce fock across MPI ranks.
    integer, intent(in) :: n
      !! Size of packed matrices: num_bas*(num_bas+1)/2
    real(dp), intent(in) :: density(n)
      !! Density matrix (packed lower triangular)
    real(dp), intent(inout) :: fock(n)
      !! Fock matrix (packed lower triangular), this rank's contribution

    if (.not. is_initialized) then
      write(*,*) "ERROR: liberi_fock called before liberi_init"
      error stop 1
    end if

    if (n /= matrix_size) then
      write(*,*) "ERROR: matrix size mismatch in liberi_fock"
      write(*,*) "  Expected:", matrix_size
      write(*,*) "  Got:", n
      error stop 1
    end if

    call liberi_fock_build(handle, density, fock, n)

  end subroutine liberi_fock


  subroutine liberi_finalize()
    !! Clean up libERI resources.
    !!
    !! Call this after SCF convergence (or on error).
    !! Equivalent to MDI_CLEANUP in the old interface.

    if (is_initialized) then
      call liberi_cleanup(handle)
      call liberi_destroy(handle)
      is_initialized = .false.
    end if

  end subroutine liberi_finalize

end module gamess_liberi_wrapper
