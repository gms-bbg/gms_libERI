! Standalone test for libERI: H2 / cc-pVDZ
!
! Exercises (ss|ss), (ss|sp), (ss|pp), (sp|sp), (sp|pp), (pp|pp) kernels.
! Uses the liberi_interface API (no globals, no MDI).
!
! Usage:
!   ./test_h2_ccpvdz
!   ./test_h2_ccpvdz --generate
!
program test_h2_ccpvdz
  use liberi_types, only: dp, int64
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_destroy, &
                              liberi_setup, liberi_fock_build, liberi_cleanup
  implicit none

  type(liberi_handle_t) :: handle
  integer, parameter :: n_tri = 55
  integer :: i, tri_size
  integer :: nargs, iunit
  character(len=256) :: arg
  logical :: generate
  real(dp) :: max_err, err
  real(dp), parameter :: tol = 1.0d-10
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp

  ! Basis set data
  integer :: nsh_val, natoms_val, num_bas_val, mxgtot_val
  integer, allocatable :: ang_mom_arr(:), contr_num_arr(:), sh_loc_arr(:)
  integer, allocatable :: atom_num_arr(:), atom_loc_arr(:)
  integer, allocatable :: start_bas_arr(:), end_bas_arr(:)
  real(dp), allocatable :: exponents_arr(:)
  real(dp), allocatable :: contr_coef_s_arr(:), contr_coef_p_arr(:)
  real(dp), allocatable :: contr_coef_d_arr(:), contr_coef_f_arr(:)
  real(dp), allocatable :: coords_arr(:), schwrz_int_arr(:)
  real(dp), allocatable :: density_arr(:), fock_arr(:)

  ! Reference Fock matrix: H2/cc-pVDZ, identity density, no screening
  real(dp), parameter :: fa_ref(n_tri) = [ &
                         1.62604241861396770d+00, &  !  1
                         3.25113803590695349d+00, &  !  2
                         3.24245511082653382d+00, &  !  3
                         0.00000000000000000d+00, &  !  4
                         0.00000000000000000d+00, &  !  5
                         4.55331209455179042d+00, &  !  6
                         0.00000000000000000d+00, &  !  7
                         0.00000000000000000d+00, &  !  8
                         0.00000000000000000d+00, &  !  9
                         4.55331209455178954d+00, &  ! 10
                         3.42480278809267513d-01, &  ! 11
                         4.53866185943757272d-01, &  ! 12
                         0.00000000000000000d+00, &  ! 13
                         0.00000000000000000d+00, &  ! 14
                         4.69669813975613071d+00, &  ! 15
                         1.82712553583075765d+00, &  ! 16
                         2.70464914998301476d+00, &  ! 17
                         0.00000000000000000d+00, &  ! 18
                         0.00000000000000000d+00, &  ! 19
                         3.20583208827434163d+00, &  ! 20
                         1.62604241861396770d+00, &  ! 21
                         2.70464914998301476d+00, &  ! 22
                         5.83443825838191721d+00, &  ! 23
                         0.00000000000000000d+00, &  ! 24
                         0.00000000000000000d+00, &  ! 25
                         1.72879053792554793d+00, &  ! 26
                         3.25113803590695349d+00, &  ! 27
                         3.24245511082653426d+00, &  ! 28
                         0.00000000000000000d+00, &  ! 29
                         0.00000000000000000d+00, &  ! 30
                         4.43985548181141265d+00, &  ! 31
                         0.00000000000000000d+00, &  ! 32
                         0.00000000000000000d+00, &  ! 33
                         0.00000000000000000d+00, &  ! 34
                         0.00000000000000000d+00, &  ! 35
                         4.55331209455179042d+00, &  ! 36
                         0.00000000000000000d+00, &  ! 37
                         0.00000000000000000d+00, &  ! 38
                         0.00000000000000000d+00, &  ! 39
                         4.43985548181141265d+00, &  ! 40
                         0.00000000000000000d+00, &  ! 41
                         0.00000000000000000d+00, &  ! 42
                         0.00000000000000000d+00, &  ! 43
                         0.00000000000000000d+00, &  ! 44
                         4.55331209455179042d+00, &  ! 45
                         -3.20583208827434163d+00, &  ! 46
                         -1.72879053792554793d+00, &  ! 47
                         0.00000000000000000d+00, &  ! 48
                         0.00000000000000000d+00, &  ! 49
                         -2.46963256595329916d+00, &  ! 50
                         -3.42480278809267513d-01, &  ! 51
                         -4.53866185943757272d-01, &  ! 52
                         0.00000000000000000d+00, &  ! 53
                         0.00000000000000000d+00, &  ! 54
                         4.69669813975613071d+00 &  ! 55
                         ]

  ! cc-pVDZ basis for hydrogen
  !   S  3  1.00
  !       13.0100000     0.0196850
  !        1.9620000     0.1379770
  !        0.4446000     0.4781480
  !   S  1  1.00
  !        0.1220000     1.0000000
  !   P  1  1.00
  !        0.7270000     1.0000000

  ! Parse --generate flag
  generate = .false.
  nargs = command_argument_count()
  do i = 1, nargs
    call get_command_argument(i, arg)
    if (trim(arg) == '--generate') generate = .true.
  end do

  ! ============================================
  ! Molecule: H2, bond length 1.4 bohr
  ! ============================================
  natoms_val = 2
  nsh_val = 6          ! 3 shells per H: s(3), s(1), p(1)
  num_bas_val = 10     ! (1+1+3) * 2 atoms
  mxgtot_val = 3       ! max primitives in any shell
  tri_size = num_bas_val * (num_bas_val + 1) / 2   ! 55

  ! ============================================
  ! Allocate arrays
  ! ============================================
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

  ! ============================================
  ! Atomic coordinates (bohr)
  ! H1 at (0, 0, -0.7), H2 at (0, 0, +0.7)
  ! ============================================
  coords_arr(1) = 0.0_dp; coords_arr(2) = 0.0_dp; coords_arr(3) = -0.7_dp
  coords_arr(4) = 0.0_dp; coords_arr(5) = 0.0_dp; coords_arr(6) = 0.7_dp

  ! ============================================
  ! Basis set: cc-pVDZ for H (atom 1)
  ! ============================================

  ! Shell 1: s, 3 primitives
  ang_mom_arr(1) = 1      ! 1=s in libERI convention
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

  ! Shell 2: s, 1 primitive
  ang_mom_arr(2) = 1
  contr_num_arr(2) = 1
  sh_loc_arr(2) = 4
  atom_num_arr(2) = 1
  atom_loc_arr(2) = 2
  start_bas_arr(2) = 2; end_bas_arr(2) = 2
  exponents_arr(4) = 0.1220000_dp
  contr_coef_s_arr(4) = 1.0000000_dp * gto_norm_s(0.1220000_dp)

  ! Shell 3: p, 1 primitive
  ang_mom_arr(3) = 2      ! 2=p in libERI convention
  contr_num_arr(3) = 1
  sh_loc_arr(3) = 5
  atom_num_arr(3) = 1
  atom_loc_arr(3) = 3
  start_bas_arr(3) = 3; end_bas_arr(3) = 5
  exponents_arr(5) = 0.7270000_dp
  contr_coef_p_arr(5) = 1.0000000_dp * gto_norm_p(0.7270000_dp)

  ! ============================================
  ! Basis set: cc-pVDZ for H (atom 2)
  ! Same exponents/coefficients, different center
  ! ============================================

  ! Shell 4: s, 3 primitives
  ang_mom_arr(4) = 1
  contr_num_arr(4) = 3
  sh_loc_arr(4) = 6
  atom_num_arr(4) = 2
  atom_loc_arr(4) = 6
  start_bas_arr(4) = 6; end_bas_arr(4) = 6
  exponents_arr(6) = 13.0100000_dp
  exponents_arr(7) = 1.9620000_dp
  exponents_arr(8) = 0.4446000_dp
  contr_coef_s_arr(6) = 0.0196850_dp * gto_norm_s(13.0100000_dp)
  contr_coef_s_arr(7) = 0.1379770_dp * gto_norm_s(1.9620000_dp)
  contr_coef_s_arr(8) = 0.4781480_dp * gto_norm_s(0.4446000_dp)

  ! Shell 5: s, 1 primitive
  ang_mom_arr(5) = 1
  contr_num_arr(5) = 1
  sh_loc_arr(5) = 9
  atom_num_arr(5) = 2
  atom_loc_arr(5) = 7
  start_bas_arr(5) = 7; end_bas_arr(5) = 7
  exponents_arr(9) = 0.1220000_dp
  contr_coef_s_arr(9) = 1.0000000_dp * gto_norm_s(0.1220000_dp)

  ! Shell 6: p, 1 primitive
  ang_mom_arr(6) = 2
  contr_num_arr(6) = 1
  sh_loc_arr(6) = 10
  atom_num_arr(6) = 2
  atom_loc_arr(6) = 8
  start_bas_arr(6) = 8; end_bas_arr(6) = 10
  exponents_arr(10) = 0.7270000_dp
  contr_coef_p_arr(10) = 1.0000000_dp * gto_norm_p(0.7270000_dp)

  ! ============================================
  ! Disable Schwarz screening (set large values)
  ! ============================================
  schwrz_int_arr = 1.0d10

  ! ============================================
  ! Density matrix — identity (diagonal = 1)
  ! ============================================
  density_arr = 0.0_dp
  do i = 1, num_bas_val
    density_arr(i * (i + 1) / 2) = 1.0_dp
  end do
  fock_arr = 0.0_dp

  ! ============================================
  ! Use liberi_interface API
  ! ============================================
  print *, "Creating handle..."
  call liberi_create(handle)

  print *, "Setting up basis and computing shell pairs..."
  call liberi_setup(handle, nsh_val, natoms_val, num_bas_val, mxgtot_val, &
                    ang_mom_arr, contr_num_arr, sh_loc_arr,               &
                    atom_num_arr, atom_loc_arr,                           &
                    start_bas_arr, end_bas_arr,                           &
                    exponents_arr,                                        &
                    contr_coef_s_arr, contr_coef_p_arr,                   &
                    contr_coef_d_arr, contr_coef_f_arr,                   &
                    coords_arr, schwrz_int_arr,                           &
                    0, 1)  ! my_rank=0, num_procs=1

  print *, "Computing Fock matrix..."
  call liberi_fock_build(handle, density_arr, fock_arr, tri_size)

  ! ============================================
  ! Generate or validate Fock matrix
  ! ============================================
  if (generate) then
    iunit = 20
    open(unit=iunit, file='fa_h2_reference.dat', form='formatted', status='replace')
    write(iunit, '(I6)') tri_size
    do i = 1, tri_size
      write(iunit, '(ES26.17)') fock_arr(i)
    end do
    close(iunit)
    print *, ""
    print *, "Wrote reference to: fa_h2_reference.dat"
    do i = 1, tri_size
      print '(A,I3,A,ES24.15)', "  fa(", i, ") = ", fock_arr(i)
    end do
  else
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

  ! ============================================
  ! Cleanup
  ! ============================================
  print *, "Cleaning up..."
  call liberi_cleanup(handle)
  call liberi_destroy(handle)

  deallocate(ang_mom_arr, contr_num_arr, sh_loc_arr)
  deallocate(atom_num_arr, atom_loc_arr, start_bas_arr, end_bas_arr)
  deallocate(coords_arr, exponents_arr)
  deallocate(contr_coef_s_arr, contr_coef_p_arr, contr_coef_d_arr, contr_coef_f_arr)
  deallocate(schwrz_int_arr, density_arr, fock_arr)

  print *, "Done!"

contains

  ! GAMESS normalization for s-type primitive: N = (2α/π)^(3/4)
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp * alpha / pi_val)**0.75_dp
  end function gto_norm_s

  ! GAMESS normalization for p-type primitive: N = 2√α × (2α/π)^(3/4)
  pure function gto_norm_p(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 2.0_dp * sqrt(alpha) * (2.0_dp * alpha / pi_val)**0.75_dp
  end function gto_norm_p

end program test_h2_ccpvdz
