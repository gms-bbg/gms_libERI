!! Direct Fortran interface for GAMESS integration (no MDI)
!!
!! This module provides the primary API for libERI when linked directly
!! to GAMESS (or other Fortran drivers) without the MDI plugin layer.
!!
!! The interface follows the libcchem/GMSHPC pattern with an opaque handle:
!!
!! 1. liberi_create() - create a handle
!! 2. liberi_setup() - called once per geometry (before SCF loop)
!! 3. liberi_fock_build() - called every SCF iteration
!! 4. liberi_cleanup() - called once after convergence
!! 5. liberi_destroy() - destroy the handle
!!
!! MPI reduction of the Fock matrix is handled by the caller (GAMESS).
!! Each rank computes its portion; GAMESS calls DDI_GSUMF afterward.
module liberi_interface
  use liberi_types, only: int32, int64, dp, basis_t, shell_pair_container_t
  use liberi_shell_pair, only: shell_pair, shell_pair_deallocate
  use liberi_integral_driver, only: compute_integrals
  implicit none
  private

  public :: liberi_handle_t
  public :: liberi_create
  public :: liberi_destroy
  public :: liberi_setup
  public :: liberi_fock_build
  public :: liberi_cleanup

  type :: liberi_handle_t
    !! Opaque handle for libERI state.
    !!
    !! Holds all persistent state for a libERI calculation:
    !! basis set data, shell-pair data (GPU-mapped), rank/size for
    !! work distribution, and initialization flag.
    !!
    !! Create with liberi_create(), destroy with liberi_destroy().
    private
    type(basis_t) :: basis
    type(shell_pair_container_t) :: pairs
    logical :: is_initialized = .false.
  end type liberi_handle_t

contains

  subroutine liberi_create(handle)
    !! Create a new libERI handle.
    !!
    !! Initializes a handle for libERI operations.
    !! The handle must be destroyed with liberi_destroy() when done.
    type(liberi_handle_t), intent(out) :: handle
      !! The newly created handle

    handle%is_initialized = .false.
  end subroutine liberi_create

  subroutine liberi_destroy(handle)
    !! Destroy a libERI handle.
    !!
    !! Cleans up all resources associated with the handle.
    !! Calls liberi_cleanup() if still initialized.
    type(liberi_handle_t), intent(inout) :: handle
      !! The handle to destroy

    if (handle%is_initialized) then
      call liberi_cleanup(handle)
    end if
  end subroutine liberi_destroy

  ! TODO: this can be done better from gamess and other programs, with how things are passed in
  ! having these many arguments is not ideal
  subroutine liberi_setup(handle, nsh_in, natoms_in, num_bas_in, mxgtot_in, &
                          ang_mom_in, contr_num_in, sh_loc_in, &
                          atom_num_in, atom_loc_in, &
                          start_bas_in, end_bas_in, &
                          exponents_in, &
                          contr_coef_s_in, contr_coef_p_in, &
                          contr_coef_d_in, contr_coef_f_in, &
                          coords_in, schwrz_int_in, &
                          my_rank_in, num_procs_in)
    !! Initialize libERI with basis set data.
    !!
    !! Called once per molecular geometry, before the SCF iteration loop.
    !! This routine stores all basis set data in the handle, computes
    !! shell pairs, and maps data to GPU.
    !!
    !! After this call, GPU data persists until liberi_cleanup() is called.
    type(liberi_handle_t), intent(inout) :: handle
      !! The libERI handle (from liberi_create)
    integer, intent(in) :: nsh_in
      !! Number of shells
    integer, intent(in) :: natoms_in
      !! Number of atoms
    integer, intent(in) :: num_bas_in
      !! Number of basis functions
    integer, intent(in) :: mxgtot_in
      !! Maximum primitives per shell
    integer, intent(in) :: my_rank_in
      !! MPI rank (0-based, from GAMESS)
    integer, intent(in) :: num_procs_in
      !! Number of MPI processes (from GAMESS)
    integer, intent(in) :: ang_mom_in(nsh_in)
      !! Angular momentum per shell (1=s, 2=p, 3=d, 4=f)
    integer, intent(in) :: contr_num_in(nsh_in)
      !! Contraction length per shell
    integer, intent(in) :: sh_loc_in(nsh_in)
      !! Starting primitive index per shell
    integer, intent(in) :: atom_num_in(nsh_in)
      !! Atom center per shell
    integer, intent(in) :: atom_loc_in(nsh_in)
      !! Basis function offset per shell
    integer, intent(in) :: start_bas_in(nsh_in)
      !! First basis function per shell
    integer, intent(in) :: end_bas_in(nsh_in)
      !! Last basis function per shell
    real(dp), intent(in) :: exponents_in(mxgtot_in*nsh_in)
      !! Gaussian exponents
    real(dp), intent(in) :: contr_coef_s_in(mxgtot_in*nsh_in)
      !! s-type contraction coefficients (normalized)
    real(dp), intent(in) :: contr_coef_p_in(mxgtot_in*nsh_in)
      !! p-type contraction coefficients (normalized)
    real(dp), intent(in) :: contr_coef_d_in(mxgtot_in*nsh_in)
      !! d-type contraction coefficients (normalized)
    real(dp), intent(in) :: contr_coef_f_in(mxgtot_in*nsh_in)
      !! f-type contraction coefficients (normalized)
    real(dp), intent(in) :: coords_in(3*natoms_in)
      !! Atomic coordinates (packed as x1,y1,z1,x2,...)
    real(dp), intent(in) :: schwrz_int_in((nsh_in*(nsh_in + 1))/2)
      !! Schwarz integrals

    integer :: coef_size

    ! Safety check: cleanup if already initialized
    if (handle%is_initialized) then
      call liberi_cleanup(handle)
    end if

    ! Store scalar dimensions
    handle%basis%nsh = nsh_in
    handle%basis%natoms = natoms_in
    handle%basis%num_bas = num_bas_in
    handle%basis%mxgtot = mxgtot_in
    coef_size = mxgtot_in*nsh_in

    ! Allocate and copy shell metadata arrays
    allocate (handle%basis%ang_mom(nsh_in))
    handle%basis%ang_mom = ang_mom_in

    allocate (handle%basis%contr_num(nsh_in))
    handle%basis%contr_num = contr_num_in

    allocate (handle%basis%sh_loc(nsh_in))
    handle%basis%sh_loc = sh_loc_in

    allocate (handle%basis%atom_num(nsh_in))
    handle%basis%atom_num = atom_num_in

    allocate (handle%basis%atom_loc(nsh_in))
    handle%basis%atom_loc = atom_loc_in

    allocate (handle%basis%start_bas(nsh_in))
    handle%basis%start_bas = start_bas_in

    allocate (handle%basis%end_bas(nsh_in))
    handle%basis%end_bas = end_bas_in

    ! Allocate and copy exponents and contraction coefficients
    allocate (handle%basis%exponents(coef_size))
    handle%basis%exponents = exponents_in

    allocate (handle%basis%contr_coef_s(coef_size))
    handle%basis%contr_coef_s = contr_coef_s_in

    allocate (handle%basis%contr_coef_p(coef_size))
    handle%basis%contr_coef_p = contr_coef_p_in

    allocate (handle%basis%contr_coef_d(coef_size))
    handle%basis%contr_coef_d = contr_coef_d_in

    allocate (handle%basis%contr_coef_f(coef_size))
    handle%basis%contr_coef_f = contr_coef_f_in

    ! Allocate and copy coordinates
    allocate (handle%basis%coords(3*natoms_in))
    handle%basis%coords = coords_in

    ! Allocate and copy Schwarz integrals
    allocate (handle%basis%schwrz_int((nsh_in*(nsh_in + 1))/2))
    handle%basis%schwrz_int = schwrz_int_in

    ! Store rank/size in resources (will be copied during shell_pair)
    handle%pairs%resources%n_rank = my_rank_in
    handle%pairs%resources%n_size = num_procs_in

    ! Compute shell pairs and map to GPU
    call shell_pair(handle%basis, handle%pairs)

    handle%is_initialized = .true.

  end subroutine liberi_setup

  subroutine liberi_fock_build(handle, density_in, fock_out, matrix_size)
    !! Compute Fock matrix contribution from density matrix.
    !!
    !! Called every SCF iteration. Uses pre-computed shell-pair data
    !! that was mapped to GPU during liberi_setup().
    !!
    !! Returns this rank's partial Fock contribution. The caller (GAMESS)
    !! is responsible for MPI reduction (DDI_GSUMF).
    type(liberi_handle_t), intent(inout) :: handle
      !! The libERI handle (from liberi_setup)
    integer, intent(in) :: matrix_size
      !! Size of packed matrices: num_bas*(num_bas+1)/2
    real(dp), intent(in) :: density_in(matrix_size)
      !! Density matrix (packed lower triangular)
    real(dp), intent(inout) :: fock_out(matrix_size)
      !! Fock matrix (packed lower triangular), this rank's contribution

    if (.not. handle%is_initialized) then
      write (*, *) "ERROR: liberi_fock_build called before liberi_setup"
      error stop 1
    end if

    ! Compute integrals - writes directly to fock_out
    call compute_integrals(handle%pairs, density_in, fock_out)

  end subroutine liberi_fock_build

  subroutine liberi_cleanup(handle)
    !! Clean up libERI resources.
    !!
    !! Called once after SCF convergence (or on error).
    !! Unmaps GPU data and deallocates all arrays.
    type(liberi_handle_t), intent(inout) :: handle
      !! The libERI handle to clean up

    if (.not. handle%is_initialized) return

    ! Unmap shell-pair data from GPU
    call shell_pair_deallocate(handle%pairs)

    ! Deallocate basis data
    if (allocated(handle%basis%ang_mom)) deallocate (handle%basis%ang_mom)
    if (allocated(handle%basis%contr_num)) deallocate (handle%basis%contr_num)
    if (allocated(handle%basis%sh_loc)) deallocate (handle%basis%sh_loc)
    if (allocated(handle%basis%atom_num)) deallocate (handle%basis%atom_num)
    if (allocated(handle%basis%atom_loc)) deallocate (handle%basis%atom_loc)
    if (allocated(handle%basis%start_bas)) deallocate (handle%basis%start_bas)
    if (allocated(handle%basis%end_bas)) deallocate (handle%basis%end_bas)
    if (allocated(handle%basis%exponents)) deallocate (handle%basis%exponents)
    if (allocated(handle%basis%contr_coef_s)) deallocate (handle%basis%contr_coef_s)
    if (allocated(handle%basis%contr_coef_p)) deallocate (handle%basis%contr_coef_p)
    if (allocated(handle%basis%contr_coef_d)) deallocate (handle%basis%contr_coef_d)
    if (allocated(handle%basis%contr_coef_f)) deallocate (handle%basis%contr_coef_f)
    if (allocated(handle%basis%coords)) deallocate (handle%basis%coords)
    if (allocated(handle%basis%schwrz_int)) deallocate (handle%basis%schwrz_int)

    handle%is_initialized = .false.

  end subroutine liberi_cleanup

end module liberi_interface
