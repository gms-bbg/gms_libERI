! Test for libERI direct interface: H2 with S-functions only
!
! Tests the handle-based direct interface (liberi_interface module)
! that bypasses MDI for direct GAMESS integration.
!
! Exercises: liberi_create, liberi_setup, liberi_fock_build, liberi_cleanup, liberi_destroy
!
! Usage:
!   ./test_interface_h2             # validate
!   ./test_interface_h2 --generate  # print Fock values
!
program test_interface_h2
  use liberi_types, only: dp, int64
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_destroy, &
                              liberi_setup, liberi_fock_build, liberi_cleanup
  implicit none

  type(liberi_handle_t) :: handle
  integer, parameter :: n_tri = 10
  integer :: i
  integer :: nargs
  character(len=256) :: arg
  logical :: generate
  real(dp) :: max_err, err
  real(dp), parameter :: tol = 1.0d-10
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp

  ! Basis set data
  integer :: nsh_val, natoms_val, num_bas_val, mxgtot_val, tri_size
  integer, allocatable :: ang_mom_arr(:), contr_num_arr(:), sh_loc_arr(:)
  integer, allocatable :: atom_num_arr(:), atom_loc_arr(:)
  integer, allocatable :: start_bas_arr(:), end_bas_arr(:)
  real(dp), allocatable :: exponents_arr(:)
  real(dp), allocatable :: contr_coef_s_arr(:), contr_coef_p_arr(:)
  real(dp), allocatable :: contr_coef_d_arr(:), contr_coef_f_arr(:)
  real(dp), allocatable :: coords_arr(:), schwrz_int_arr(:)
  real(dp), allocatable :: density_arr(:), fock_arr(:)

  ! Reference Fock matrix: H2/cc-pVDZ s-only, identity density, no screening
  real(dp), parameter :: fa_ref(n_tri) = [ &
                         3.54291268527184511d-01, &  !  1
                         6.08427042171697963d-01, &  !  2
                         6.49503616123108030d-01, &  !  3
                         3.39857129332252594d-01, &  !  4
                         4.65476191699381481d-01, &  !  5
                         3.54291268527184400d-01, &  !  6
                         4.65476191699381536d-01, &  !  7
                         1.11809229981003599d+00, &  !  8
                         6.08427042171697963d-01, &  !  9
                         6.49503616123108030d-01 &   ! 10
                         ]


  ! Parse --generate flag
  generate = .false.
  nargs = command_argument_count()
  do i = 1, nargs
    call get_command_argument(i, arg)
    if (trim(arg) == '--generate') generate = .true.
  end do

  ! ============================================
  ! Set up basis: H2, bond length 1.4 bohr
  ! ============================================
  natoms_val = 2
  nsh_val = 4          ! 2 s-shells per H
  num_bas_val = 4      ! 1 basis fn per s-shell
  mxgtot_val = 3       ! max primitives in any shell
  tri_size = num_bas_val * (num_bas_val + 1) / 2

  ! Allocate arrays
  allocate(ang_mom_arr(nsh_val))
  allocate(contr_num_arr(nsh_val))
  allocate(sh_loc_arr(nsh_val))
  allocate(atom_num_arr(nsh_val))
  allocate(atom_loc_arr(nsh_val))
  allocate(start_bas_arr(nsh_val))
  allocate(end_bas_arr(nsh_val))
  allocate(coords_arr(3 * natoms_val))
  allocate(exponents_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_s_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_p_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_d_arr(mxgtot_val * nsh_val))
  allocate(contr_coef_f_arr(mxgtot_val * nsh_val))
  allocate(schwrz_int_arr(nsh_val * (nsh_val + 1) / 2))
  allocate(density_arr(tri_size))
  allocate(fock_arr(tri_size))

  ! Initialize coefficient arrays
  exponents_arr = 0.0_dp
  contr_coef_s_arr = 0.0_dp
  contr_coef_p_arr = 0.0_dp
  contr_coef_d_arr = 0.0_dp
  contr_coef_f_arr = 0.0_dp

  ! Atomic coordinates (bohr): H at (0,0,-0.7) and (0,0,+0.7)
  coords_arr(1) = 0.0_dp; coords_arr(2) = 0.0_dp; coords_arr(3) = -0.7_dp
  coords_arr(4) = 0.0_dp; coords_arr(5) = 0.0_dp; coords_arr(6) = 0.7_dp

  ! Shell 1: s, 3 primitives (contracted) on atom 1
  ang_mom_arr(1) = 1
  contr_num_arr(1) = 3
  sh_loc_arr(1) = 1
  atom_num_arr(1) = 1
  atom_loc_arr(1) = 1
  start_bas_arr(1) = 1; end_bas_arr(1) = 1
  exponents_arr(1) = 13.0100000_dp
  exponents_arr(2) = 1.9620000_dp
  exponents_arr(3) = 0.4446000_dp
  contr_coef_s_arr(1) = 0.0196850_dp * gto_norm_s(13.0100000_dp)
  contr_coef_s_arr(2) = 0.1379770_dp * gto_norm_s(1.9620000_dp)
  contr_coef_s_arr(3) = 0.4781480_dp * gto_norm_s(0.4446000_dp)

  ! Shell 2: s, 1 primitive (diffuse) on atom 1
  ang_mom_arr(2) = 1
  contr_num_arr(2) = 1
  sh_loc_arr(2) = 4
  atom_num_arr(2) = 1
  atom_loc_arr(2) = 2
  start_bas_arr(2) = 2; end_bas_arr(2) = 2
  exponents_arr(4) = 0.1220000_dp
  contr_coef_s_arr(4) = 1.0000000_dp * gto_norm_s(0.1220000_dp)

  ! Shell 3: s, 3 primitives (contracted) on atom 2
  ang_mom_arr(3) = 1
  contr_num_arr(3) = 3
  sh_loc_arr(3) = 5
  atom_num_arr(3) = 2
  atom_loc_arr(3) = 3
  start_bas_arr(3) = 3; end_bas_arr(3) = 3
  exponents_arr(5) = 13.0100000_dp
  exponents_arr(6) = 1.9620000_dp
  exponents_arr(7) = 0.4446000_dp
  contr_coef_s_arr(5) = 0.0196850_dp * gto_norm_s(13.0100000_dp)
  contr_coef_s_arr(6) = 0.1379770_dp * gto_norm_s(1.9620000_dp)
  contr_coef_s_arr(7) = 0.4781480_dp * gto_norm_s(0.4446000_dp)

  ! Shell 4: s, 1 primitive (diffuse) on atom 2
  ang_mom_arr(4) = 1
  contr_num_arr(4) = 1
  sh_loc_arr(4) = 8
  atom_num_arr(4) = 2
  atom_loc_arr(4) = 4
  start_bas_arr(4) = 4; end_bas_arr(4) = 4
  exponents_arr(8) = 0.1220000_dp
  contr_coef_s_arr(8) = 1.0000000_dp * gto_norm_s(0.1220000_dp)

  ! Disable Schwarz screening (set large values)
  schwrz_int_arr = 1.0d10

  ! Identity density matrix (diagonal = 1)
  density_arr = 0.0_dp
  do i = 1, num_bas_val
    density_arr(i * (i + 1) / 2) = 1.0_dp
  end do
  fock_arr = 0.0_dp

  ! ============================================
  ! Test the new direct interface
  ! ============================================
  call liberi_create(handle)

  ! Step 2: Setup (equivalent to MDI basis transfer + shell_pair)
  call liberi_setup(handle, nsh_val, natoms_val, num_bas_val, mxgtot_val, &
                    ang_mom_arr, contr_num_arr, sh_loc_arr,               &
                    atom_num_arr, atom_loc_arr,                           &
                    start_bas_arr, end_bas_arr,                           &
                    exponents_arr,                                        &
                    contr_coef_s_arr, contr_coef_p_arr,                   &
                    contr_coef_d_arr, contr_coef_f_arr,                   &
                    coords_arr, schwrz_int_arr,                           &
                    0, 1)  ! my_rank=0, num_procs=1

  ! Step 3: Fock build (this is what gets called every SCF iteration)
  call liberi_fock_build(handle, density_arr, fock_arr, tri_size)

  ! Step 4: Output and verify results (only rank 0)
    print *, ""
    print *, "Fock matrix (triangular packed, 10 elements):"
    do i = 1, tri_size
      print '(A,I3,A,ES24.17)', "  fa(", i, ") = ", fock_arr(i)
    end do

    if (generate) then
      print *, ""
      print *, "Copy-paste for fa_ref:"
      do i = 1, tri_size
        if (i < tri_size) then
          print '(A,ES26.17,A,I0)', "                         ", fock_arr(i), ", &  !  ", i
        else
          print '(A,ES26.17,A,I0)', "                         ", fock_arr(i), " &   !  ", i
        end if
      end do
    end if

    if (.not. generate) then
      max_err = 0.0_dp
      do i = 1, tri_size
        err = abs(fock_arr(i) - fa_ref(i))
        if (err > max_err) max_err = err
        if (err > tol) then
          print '(A,I3,A,ES24.15,A,ES24.15,A,ES10.2)', &
            "  FAIL fa(", i, ") = ", fock_arr(i), "  ref = ", fa_ref(i), "  err = ", err
        end if
      end do

      print *, ""
      if (max_err <= tol) then
        print '(A,ES10.2)', " PASS - max error: ", max_err
      else
        print '(A,ES10.2)', " FAIL - max error: ", max_err
        call liberi_destroy(handle)
        error stop 1
      end if
    end if

  ! Step 5: Cleanup and destroy
    print *, ""
    print *, "4. Cleaning up..."
  call liberi_cleanup(handle)
  call liberi_destroy(handle)

  ! Deallocate local arrays
  deallocate(ang_mom_arr, contr_num_arr, sh_loc_arr)
  deallocate(atom_num_arr, atom_loc_arr, start_bas_arr, end_bas_arr)
  deallocate(coords_arr, exponents_arr)
  deallocate(contr_coef_s_arr, contr_coef_p_arr, contr_coef_d_arr, contr_coef_f_arr)
  deallocate(schwrz_int_arr, density_arr, fock_arr)



contains

  ! GAMESS normalization for s-type primitive: N = (2a/pi)^(3/4)
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp * alpha / pi_val)**0.75_dp
  end function gto_norm_s

end program test_interface_h2
