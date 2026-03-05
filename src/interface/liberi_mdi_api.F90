!! MDI engine interface for GAMESS integration
module mdi_api
  !! Provides the MolSSI Driver Interface (MDI) plugin API for libERI.
  !!
  !! This module implements the MDI engine protocol, allowing libERI to be
  !! loaded as a plugin by GAMESS (or other MDI-compatible drivers) for
  !! GPU-accelerated two-electron integral evaluation.
  !!
  !! ### MDI Commands
  !!
  !! The following MDI commands are registered:
  !!
  !! | Command | Description |
  !! |---------|-------------|
  !! | `>NATOMS` | Receive number of atoms |
  !! | `>COORDS` | Receive atomic coordinates |
  !! | `>NSH` | Receive number of shells |
  !! | `>MXGTOT` | Receive max primitives per shell |
  !! | `>NUM_BAS` | Receive number of basis functions |
  !! | `>DENSITY` | Receive density matrix |
  !! | `>DENSITY_B` | Receive beta density (UHF) |
  !! | `<FOCK` | Send Fock matrix back to driver |
  !! | `EXIT` | Clean up and exit |
  !!
  !! ### Compilation
  !!
  !! This module requires the `USE_MDI` preprocessor flag and the MDI library.
  !! Without `USE_MDI`, the module provides stub implementations.
  !!
  !! ### Usage
  !!
  !! The plugin is loaded by GAMESS at runtime. The driver sends basis set
  !! data and density matrices, then requests Fock matrix computation.
#ifdef USE_MDI
  use mpi, only: MPI_Comm_rank, MPI_Comm_size, MPI_Bcast, MPI_Reduce, MPI_CHAR, MPI_DOUBLE, MPI_INT, MPI_SUM
#endif
  use iso_c_binding, only: C_PTR, C_FUNLOC, C_LOC, C_F_POINTER
  use liberi_types, only: int32, int64, dp
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_setup, &
                              liberi_fock_build, liberi_cleanup, liberi_destroy
#ifdef USE_MDI
  use mdi, ONLY: &
    MDI_Init, MDI_Send, MDI_Send_Command, MDI_INT, MDI_CHAR, MDI_NAME_LENGTH, &
    MDI_Accept_communicator, MDI_Recv_command, MDI_Recv, &
    MDI_Set_execute_command_func, MDI_DOUBLE, MDI_BYTE, &
    MDI_ENGINE, MDI_Get_role, MDI_Register_command, MDI_Register_node, &
    MDI_Register_callback, MDI_COMMAND_LENGTH, &
    MDI_Plugin_get_argc, MDI_Plugin_get_arg, MDI_Get_communicator, &
    MDI_Get_method, MDI_Set_plugin_state, MDI_MPI_get_world_comm
#endif
  implicit none

  type :: mdi_state_t
    !! Bundles all state for an MDI plugin session.
    !! Replaces module-level variables and liberi_globals usage.
    type(liberi_handle_t) :: handle

    ! MDI / MPI state
    integer :: comm = 0
    integer :: world_comm = 0
    logical :: terminate_flag = .false.

    ! MPI rank/size
    integer :: n_rank = 0
    integer :: n_size = 1

    ! Basis set scalars (accumulated from MDI receives)
    integer :: natoms = 0
    integer :: num_bas = 0
    integer :: nsh = 0
    integer :: mxgtot = 0
    integer :: nsh2 = 0
    integer :: size_of_matrix = 0

    ! SCF type
    character(LEN=10) :: typscf = ' '

    ! Basis set arrays (accumulated from MDI receives, freed after setup)
    integer, allocatable :: ang_mom(:)
    integer, allocatable :: contr_num(:)
    integer, allocatable :: sh_loc(:)
    integer, allocatable :: atom_num(:)
    integer, allocatable :: atom_loc(:)
    integer, allocatable :: start_bas(:)
    integer, allocatable :: end_bas(:)
    real(dp), allocatable :: exponents(:)
    real(dp), allocatable :: contr_coef_s(:)
    real(dp), allocatable :: contr_coef_p(:)
    real(dp), allocatable :: contr_coef_d(:)
    real(dp), allocatable :: contr_coef_f(:)
    real(dp), allocatable :: coords(:)
    real(dp), allocatable :: schwrz_int(:)

    ! Per-iteration matrices
    real(dp), allocatable :: density(:)
    real(dp), allocatable :: density_b(:)
    real(dp), allocatable :: fa(:)
    real(dp), allocatable :: fa_final2(:)
    real(dp), allocatable :: fb(:)
    real(dp), allocatable :: fb_final(:)
  end type mdi_state_t

  ! Persistent state across plugin launches (GAMESS launches the plugin
  ! multiple times: once for setup, then once per SCF iteration for Fock builds).
  ! The library stays loaded (RTLD_NODELETE), so module-level state survives.
  type(mdi_state_t), save, target :: persistent_state
  logical, save :: state_initialized = .false.

contains
#ifdef USE_MDI
  function MDI_Plugin_init_ERI(plugin_state) bind(C, name="MDI_Plugin_init_ERI")
    TYPE(C_PTR), VALUE :: plugin_state
    integer :: MDI_Plugin_init_ERI
    integer :: ierr

    CALL MDI_Set_plugin_state(plugin_state, ierr)

    if (.not. state_initialized) then
      call liberi_create(persistent_state%handle)
      state_initialized = .true.
    end if

    ! Get the MPI intra-communicator over which this plugin will run
    CALL MDI_MPI_get_world_comm(persistent_state%world_comm, ierr)
    call MPI_Comm_rank(persistent_state%world_comm, persistent_state%n_rank, ierr)
    call MPI_COMM_SIZE(persistent_state%world_comm, persistent_state%n_size, ierr)

    ! Perform one-time operations required to establish a connection with the driver
    CALL initialize_mdi(persistent_state)

    ! Reset terminate flag for this launch
    persistent_state%terminate_flag = .false.

    ! Respond to commands from the driver
    CALL respond_to_commands(persistent_state)

    if (allocated(persistent_state%density)) deallocate (persistent_state%density)
    if (allocated(persistent_state%fa)) deallocate (persistent_state%fa)
    if (allocated(persistent_state%density_b)) deallocate (persistent_state%density_b)
    if (allocated(persistent_state%fb)) deallocate (persistent_state%fb)

    MDI_Plugin_init_ERI = 0
  END function MDI_Plugin_init_ERI

  function MDI_Plugin_open_ERI(plugin_state) bind(C, name="MDI_Plugin_open_ERI")
    TYPE(C_PTR), VALUE :: plugin_state
    integer :: MDI_Plugin_open_ERI
    integer :: ierr

    CALL MDI_Set_plugin_state(plugin_state, ierr)

    if (.not. state_initialized) then
      call liberi_create(persistent_state%handle)
      state_initialized = .true.
    end if

    ! Get the MPI intra-communicator over which this plugin will run
    CALL MDI_MPI_get_world_comm(persistent_state%world_comm, ierr)
    call MPI_Comm_rank(persistent_state%world_comm, persistent_state%n_rank, ierr)
    call MPI_COMM_SIZE(persistent_state%world_comm, persistent_state%n_size, ierr)

    ! Perform one-time operations required to establish a connection with the driver
    CALL initialize_mdi(persistent_state)

    persistent_state%terminate_flag = .false.

    CALL respond_to_commands(persistent_state)

    MDI_Plugin_open_ERI = 0
  END function MDI_Plugin_open_ERI

  function MDI_Plugin_close_ERI() bind(C, name="MDI_Plugin_close_ERI")

    integer :: MDI_Plugin_close_ERI

    MDI_Plugin_close_ERI = 0
  END function MDI_Plugin_close_ERI

  SUBROUTINE initialize_mdi(state)
    type(mdi_state_t), intent(inout), target :: state
    integer :: ierr, role

    PROCEDURE(execute_command_wrapper), POINTER :: generic_command => null()
    TYPE(C_PTR)                         :: class_obj
    generic_command => execute_command_wrapper

    ! Confirm that the code is being run as an ENGINE
    call MDI_Get_role(role, ierr)
    IF (role .ne. MDI_ENGINE) THEN
      WRITE (6, *) 'ERROR: Must run libERI as an ENGINE'
    END IF

    ! Register the commands
    CALL MDI_Register_node("@DEFAULT", ierr)
    CALL MDI_Register_command("@DEFAULT", "EXIT", ierr)
    CALL MDI_Register_command("@DEFAULT", "EXIT_CONVERGED", ierr)
    CALL MDI_Register_command("@DEFAULT", "<FOCK", ierr)
    CALL MDI_Register_command("@DEFAULT", ">FOCK", ierr)
    CALL MDI_Register_command("@DEFAULT", ">NATOMS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">RANK_ID", ierr)
    CALL MDI_Register_command("@DEFAULT", ">NUM_PROCS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">NUM_BAS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">NSH", ierr)
    CALL MDI_Register_command("@DEFAULT", ">mxgtot", ierr)
    CALL MDI_Register_command("@DEFAULT", ">ANGM", ierr)
    CALL MDI_Register_command("@DEFAULT", "<COORDS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">COORDS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CONTR", ierr)
    CALL MDI_Register_command("@DEFAULT", ">SHLLOC", ierr)
    CALL MDI_Register_command("@DEFAULT", ">ATMNUM", ierr)
    CALL MDI_Register_command("@DEFAULT", ">ATMLOC", ierr)
    CALL MDI_Register_command("@DEFAULT", ">SBAS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">EBAS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">EXPON", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CONTR_COEF_S", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CONTR_COEF_P", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CONTR_COEF_D", ierr)
    CALL MDI_Register_command("@DEFAULT", ">CONTR_COEF_F", ierr)
    CALL MDI_Register_command("@DEFAULT", ">SCHWZ_INT", ierr)
    CALL MDI_Register_command("@DEFAULT", ">AO_DENS", ierr)
    CALL MDI_Register_command("@DEFAULT", ">TYPSCF", ierr)
    CALL MDI_Register_command("@DEFAULT", "<FOCK_B", ierr)
    CALL MDI_Register_command("@DEFAULT", ">FOCK_B", ierr)
    CALL MDI_Register_command("@DEFAULT", ">AO_DENS_B", ierr)

    ! Connect to the driver
    CALL MDI_Accept_communicator(state%comm, ierr)

    ! Set the generic execute_command function, passing state as class_obj
    class_obj = C_LOC(state)
    CALL MDI_Set_execute_command_func(C_FUNLOC(generic_command), class_obj, ierr)

  END SUBROUTINE initialize_mdi

  SUBROUTINE respond_to_commands(state)
    type(mdi_state_t), intent(inout) :: state
    character(len=:), ALLOCATABLE :: command
    integer :: ierr

    ALLOCATE (character(MDI_COMMAND_LENGTH) :: command)
    ! Respond to the driver's commands
    response_loop: DO

      CALL MDI_Recv_command(command, state%comm, ierr)
      CALL MPI_Bcast(command, MDI_COMMAND_LENGTH, MPI_CHAR, 0, state%world_comm, ierr)
      CALL execute_command(command, state%comm, state, ierr)
      IF (state%terminate_flag) EXIT

    END DO response_loop

    DEALLOCATE (command)

  END SUBROUTINE respond_to_commands

  function execute_command_wrapper(command, comm, obj_ptr)

    character(LEN=*), intent(IN)  :: command
    integer, intent(IN)           :: comm
    TYPE(C_PTR), VALUE            :: obj_ptr
    integer                       :: execute_command_wrapper

    integer                       :: ierr
    type(mdi_state_t), pointer    :: state_ptr

    call C_F_POINTER(obj_ptr, state_ptr)
    CALL execute_command(command, comm, state_ptr, ierr)
    execute_command_wrapper = ierr

  END function execute_command_wrapper

  SUBROUTINE execute_command(command, comm, state, ierr)

    character(LEN=*), intent(IN)  :: command
    integer, intent(IN)           :: comm
    type(mdi_state_t), intent(inout) :: state
    integer, intent(OUT)          :: ierr

    SELECT CASE (TRIM(command))
    CASE ("EXIT")
      state%terminate_flag = .true.

    CASE ("EXIT_CONVERGED")
      call liberi_cleanup(state%handle)
      call liberi_destroy(state%handle)
      state_initialized = .false.
      state%terminate_flag = .true.

    CASE (">TYPSCF")
      CALL MDI_Recv(state%typscf, 10, MDI_CHAR, comm, ierr)
      CALL MPI_Bcast(state%typscf, 10, MPI_CHAR, 0, state%world_comm, ierr)
    CASE (">NATOMS")
      CALL MDI_Recv(state%natoms, 1, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%natoms, 1, MPI_INT, 0, state%world_comm, ierr)
    CASE (">NUM_BAS")
      CALL MDI_Recv(state%num_bas, 1, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%num_bas, 1, MPI_INT, 0, state%world_comm, ierr)
    CASE (">NSH")
      CALL MDI_Recv(state%nsh, 1, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%nsh, 1, MPI_INT, 0, state%world_comm, ierr)
    CASE (">mxgtot")
      CALL MDI_Recv(state%mxgtot, 1, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%mxgtot, 1, MPI_INT, 0, state%world_comm, ierr)

    CASE (">ANGM")
      if (.not. allocated(state%ang_mom)) allocate (state%ang_mom(state%nsh))
      CALL MDI_Recv(state%ang_mom, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%ang_mom, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">CONTR")
      if (.not. allocated(state%contr_num)) allocate (state%contr_num(state%nsh))
      CALL MDI_Recv(state%contr_num, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%contr_num, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">SHLLOC")
      if (.not. allocated(state%sh_loc)) allocate (state%sh_loc(state%nsh))
      CALL MDI_Recv(state%sh_loc, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%sh_loc, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">ATMNUM")
      if (.not. allocated(state%atom_num)) allocate (state%atom_num(state%nsh))
      CALL MDI_Recv(state%atom_num, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%atom_num, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">ATMLOC")
      if (.not. allocated(state%atom_loc)) allocate (state%atom_loc(state%nsh))
      CALL MDI_Recv(state%atom_loc, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%atom_loc, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">SBAS")
      if (.not. allocated(state%start_bas)) allocate (state%start_bas(state%nsh))
      CALL MDI_Recv(state%start_bas, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%start_bas, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">EBAS")
      if (.not. allocated(state%end_bas)) allocate (state%end_bas(state%nsh))
      CALL MDI_Recv(state%end_bas, state%nsh, MDI_INT, comm, ierr)
      CALL MPI_Bcast(state%end_bas, state%nsh, MPI_INT, 0, state%world_comm, ierr)
    CASE (">EXPON")
      if (.not. allocated(state%exponents)) allocate (state%exponents(state%mxgtot*state%nsh))
      CALL MDI_Recv(state%exponents, state%mxgtot*state%nsh, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%exponents, state%mxgtot*state%nsh, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">CONTR_COEF_S")
      if (.not. allocated(state%contr_coef_s)) allocate (state%contr_coef_s(state%mxgtot*state%nsh))
      CALL MDI_Recv(state%contr_coef_s, state%mxgtot*state%nsh, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%contr_coef_s, state%mxgtot*state%nsh, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">CONTR_COEF_P")
      if (.not. allocated(state%contr_coef_p)) allocate (state%contr_coef_p(state%mxgtot*state%nsh))
      CALL MDI_Recv(state%contr_coef_p, state%mxgtot*state%nsh, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%contr_coef_p, state%mxgtot*state%nsh, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">CONTR_COEF_D")
      if (.not. allocated(state%contr_coef_d)) allocate (state%contr_coef_d(state%mxgtot*state%nsh))
      CALL MDI_Recv(state%contr_coef_d, state%mxgtot*state%nsh, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%contr_coef_d, state%mxgtot*state%nsh, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">CONTR_COEF_F")
      if (.not. allocated(state%contr_coef_f)) allocate (state%contr_coef_f(state%mxgtot*state%nsh))
      CALL MDI_Recv(state%contr_coef_f, state%mxgtot*state%nsh, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%contr_coef_f, state%mxgtot*state%nsh, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">COORDS")
      if (.not. allocated(state%coords)) allocate (state%coords(3*state%natoms))
      CALL MDI_Recv(state%coords, 3*state%natoms, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%coords, 3*state%natoms, MPI_DOUBLE, 0, state%world_comm, ierr)

    CASE (">SCHWZ_INT")
      state%nsh2 = (state%nsh*state%nsh + state%nsh)/2
      if (.not. allocated(state%schwrz_int)) allocate (state%schwrz_int(state%nsh2))
      CALL MDI_Recv(state%schwrz_int, state%nsh2, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%schwrz_int, state%nsh2, MPI_DOUBLE, 0, state%world_comm, ierr)

      ! Set up shell pairs via the handle API
      call liberi_setup(state%handle, state%nsh, state%natoms, state%num_bas, state%mxgtot, &
                        state%ang_mom, state%contr_num, state%sh_loc, &
                        state%atom_num, state%atom_loc, &
                        state%start_bas, state%end_bas, &
                        state%exponents, &
                        state%contr_coef_s, state%contr_coef_p, &
                        state%contr_coef_d, state%contr_coef_f, &
                        state%coords, state%schwrz_int, &
                        state%n_rank, state%n_size)

      ! Free accumulation buffers — handle has its own copies now
      if (allocated(state%ang_mom)) deallocate (state%ang_mom)
      if (allocated(state%contr_num)) deallocate (state%contr_num)
      if (allocated(state%sh_loc)) deallocate (state%sh_loc)
      if (allocated(state%atom_num)) deallocate (state%atom_num)
      if (allocated(state%atom_loc)) deallocate (state%atom_loc)
      if (allocated(state%start_bas)) deallocate (state%start_bas)
      if (allocated(state%end_bas)) deallocate (state%end_bas)
      if (allocated(state%exponents)) deallocate (state%exponents)
      if (allocated(state%contr_coef_s)) deallocate (state%contr_coef_s)
      if (allocated(state%contr_coef_p)) deallocate (state%contr_coef_p)
      if (allocated(state%contr_coef_d)) deallocate (state%contr_coef_d)
      if (allocated(state%contr_coef_f)) deallocate (state%contr_coef_f)
      if (allocated(state%coords)) deallocate (state%coords)
      if (allocated(state%schwrz_int)) deallocate (state%schwrz_int)

    CASE (">AO_DENS")
      state%size_of_matrix = (state%num_bas*state%num_bas + state%num_bas)/2
      if (.not. allocated(state%density)) allocate (state%density(state%size_of_matrix))
      CALL MDI_Recv(state%density, state%size_of_matrix, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%density, state%size_of_matrix, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE (">AO_DENS_B")
      state%size_of_matrix = (state%num_bas*state%num_bas + state%num_bas)/2
      if (.not. allocated(state%density_b)) allocate (state%density_b(state%size_of_matrix))
      CALL MDI_Recv(state%density_b, state%size_of_matrix, MDI_DOUBLE, comm, ierr)
      CALL MPI_Bcast(state%density_b, state%size_of_matrix, MPI_DOUBLE, 0, state%world_comm, ierr)
    CASE ("<FOCK")
      if (state%typscf == 'RHF      ') then
        if (.not. allocated(state%fa)) allocate (state%fa(state%size_of_matrix))
        if (.not. allocated(state%fa_final2)) allocate (state%fa_final2(state%size_of_matrix))
      else if (state%typscf == 'ROHF      ') then
        if (.not. allocated(state%fa)) allocate (state%fa(state%size_of_matrix))
        if (.not. allocated(state%fa_final2)) allocate (state%fa_final2(state%size_of_matrix))
        if (.not. allocated(state%fb)) allocate (state%fb(state%size_of_matrix))
        if (.not. allocated(state%fb_final)) allocate (state%fb_final(state%size_of_matrix))
      end if

      call liberi_fock_build(state%handle, state%density, state%fa, state%size_of_matrix)

      call MPI_Reduce(state%fa, state%fa_final2, state%size_of_matrix, MPI_DOUBLE, MPI_SUM, 0, state%world_comm, ierr)
      CALL MDI_Send(state%fa_final2, state%size_of_matrix, MDI_DOUBLE, comm, ierr)

      if (allocated(state%fa_final2)) deallocate (state%fa_final2)
      if (allocated(state%fb_final)) deallocate (state%fb_final)

    CASE ("<FOCK_B")
      call MPI_Reduce(state%fb, state%fb_final, state%size_of_matrix, MPI_DOUBLE, MPI_SUM, 0, state%world_comm, ierr)
      CALL MDI_Send(state%fb_final, state%size_of_matrix, MDI_DOUBLE, comm, ierr)

    CASE DEFAULT
      WRITE (6, *) 'Error: command not recognized'
    END SELECT

    ierr = 0
  END SUBROUTINE execute_command
#endif
end module mdi_api
