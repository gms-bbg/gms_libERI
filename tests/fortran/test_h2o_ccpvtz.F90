! Standalone test for libERI: H2O / cc-pVTZ
!
! Exercises s, p, d, and f kernels (requires ENABLE_F).
! Uses the liberi_interface API (no globals, no MDI).
!
! Usage:
!   ./test_h2o_ccpvtz
!   ./test_h2o_ccpvtz --generate
!
program test_h2o_ccpvtz
  use liberi_types, only: dp, int64
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_destroy, &
                              liberi_setup, liberi_fock_build, liberi_cleanup
  implicit none

  type(liberi_handle_t) :: handle
  integer, parameter :: n_tri = 2145   ! 65*(65+1)/2
  integer :: i, tri_size, ref_size
  integer :: nargs, iunit
  character(len=256) :: arg, ref_file
  logical :: generate, file_exists
  real(dp) :: max_err, err, trace
  real(dp), parameter :: tol = 1.0d-10
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp
  real(dp), allocatable :: fa_ref(:)

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

  ! Parse --generate flag
  generate = .false.
  nargs = command_argument_count()
  do i = 1, nargs
    call get_command_argument(i, arg)
    if (trim(arg) == '--generate') generate = .true.
  end do

  ! ============================================
  ! Molecule: H2O
  !   O at origin, H atoms symmetric in yz-plane
  !   r(OH) = 1.8096 bohr, angle(HOH) = 104.51 deg
  ! ============================================
  natoms_val = 3
  nsh_val = 22         ! O: 10 shells, each H: 6 shells
  num_bas_val = 65     ! O: 35, H: 15 each, total 65
  mxgtot_val = 10      ! max primitives in any shell (O general-contracted s)
  tri_size = num_bas_val * (num_bas_val + 1) / 2   ! 2145

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
  !   O at origin, H1 and H2 in yz-plane
  ! ============================================
  coords_arr(1) = 0.0_dp; coords_arr(2) = 0.0_dp; coords_arr(3) = 0.0_dp
  coords_arr(4) = 0.0_dp; coords_arr(5) = 1.431390_dp; coords_arr(6) = 1.107160_dp
  coords_arr(7) = 0.0_dp; coords_arr(8) = -1.431390_dp; coords_arr(9) = 1.107160_dp

  ! ============================================
  ! Basis set: cc-pVTZ for O (atom 1)
  ! From Basis Set Exchange
  ! ============================================
  call setup_oxygen_basis()
  call setup_hydrogen_basis(2, 11)   ! H atom 2, starting at shell 11
  call setup_hydrogen_basis(3, 17)   ! H atom 3, starting at shell 17

  ! ============================================
  ! Disable Schwarz screening (set large values)
  ! ============================================
  schwrz_int_arr = 1.0d10

  ! ============================================
  ! Density matrix - identity (diagonal = 1)
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
    open(unit=iunit, file='fa_h2o_ccpvtz_reference.dat', form='formatted', status='replace')
    write(iunit, '(I6)') tri_size
    do i = 1, tri_size
      write(iunit, '(ES26.17)') fock_arr(i)
    end do
    close(iunit)
    print *, ""
    print *, "Wrote reference to: fa_h2o_ccpvtz_reference.dat"
    trace = 0.0_dp
    do i = 1, num_bas_val
      trace = trace + fock_arr(i * (i + 1) / 2)
    end do
    print '(A,ES24.15)', "  Trace = ", trace
  else
    ! Try to find reference file
    ref_file = 'fa_h2o_ccpvtz_reference.dat'
    inquire(file=ref_file, exist=file_exists)
    if (.not. file_exists) then
      ref_file = '../tests/fortran/fa_h2o_ccpvtz_reference.dat'
      inquire(file=ref_file, exist=file_exists)
    end if

    if (file_exists) then
      ! Read reference and validate
      iunit = 20
      open(unit=iunit, file=ref_file, form='formatted', status='old')
      read(iunit, '(I6)') ref_size
      if (ref_size /= tri_size) then
        print *, "ERROR: Reference size mismatch"
        call liberi_destroy(handle)
        error stop 1
      end if
      allocate(fa_ref(tri_size))
      do i = 1, tri_size
        read(iunit, '(ES26.17)') fa_ref(i)
      end do
      close(iunit)

      max_err = 0.0_dp
      do i = 1, tri_size
        err = abs(fock_arr(i) - fa_ref(i))
        if (err > max_err) max_err = err
        if (err > tol) then
          print '(A,I4,A,ES24.15,A,ES24.15,A,ES10.2)', &
            "  FAIL fa(", i, ") = ", fock_arr(i), "  ref = ", fa_ref(i), "  err = ", err
        end if
      end do

      print *, ""
      if (max_err <= tol) then
        print '(A,ES10.2)', " PASS - max error: ", max_err
      else
        print '(A,ES10.2)', " FAIL - max error: ", max_err
        deallocate(fa_ref)
        call liberi_destroy(handle)
        error stop 1
      end if
      deallocate(fa_ref)
    else
      ! No reference file, just print trace
      trace = 0.0_dp
      do i = 1, num_bas_val
        trace = trace + fock_arr(i * (i + 1) / 2)
      end do
      print *, ""
      print '(A,ES24.15)', " Fock matrix trace: ", trace
      print *, " PASS - H2O/cc-pVTZ calculation completed (no reference file for validation)"
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

  subroutine setup_oxygen_basis()
    integer :: iloc

    ! Shell 1: S 10 primitives (1st general contraction)
    ang_mom_arr(1) = 1
    contr_num_arr(1) = 10
    sh_loc_arr(1) = 1
    atom_num_arr(1) = 1
    atom_loc_arr(1) = 1
    start_bas_arr(1) = 1; end_bas_arr(1) = 1
    iloc = 1
    exponents_arr(iloc) = 15330.0_dp; contr_coef_s_arr(iloc) = 0.000508_dp*gto_norm_s(15330.0_dp)
    exponents_arr(iloc + 1) = 2299.0_dp; contr_coef_s_arr(iloc + 1) = 0.003929_dp*gto_norm_s(2299.0_dp)
    exponents_arr(iloc + 2) = 522.4_dp; contr_coef_s_arr(iloc + 2) = 0.020243_dp*gto_norm_s(522.4_dp)
    exponents_arr(iloc + 3) = 147.3_dp; contr_coef_s_arr(iloc + 3) = 0.079181_dp*gto_norm_s(147.3_dp)
    exponents_arr(iloc + 4) = 47.55_dp; contr_coef_s_arr(iloc + 4) = 0.230687_dp*gto_norm_s(47.55_dp)
    exponents_arr(iloc + 5) = 16.76_dp; contr_coef_s_arr(iloc + 5) = 0.433118_dp*gto_norm_s(16.76_dp)
    exponents_arr(iloc + 6) = 6.207_dp; contr_coef_s_arr(iloc + 6) = 0.350260_dp*gto_norm_s(6.207_dp)
    exponents_arr(iloc + 7) = 1.752_dp; contr_coef_s_arr(iloc + 7) = 0.042728_dp*gto_norm_s(1.752_dp)
    exponents_arr(iloc + 8) = 0.6882_dp; contr_coef_s_arr(iloc + 8) = -0.008154_dp*gto_norm_s(0.6882_dp)
    exponents_arr(iloc + 9) = 0.2384_dp; contr_coef_s_arr(iloc + 9) = 0.002381_dp*gto_norm_s(0.2384_dp)

    ! Shell 2: S 1 primitive (uncontracted from general)
    ang_mom_arr(2) = 1
    contr_num_arr(2) = 1
    sh_loc_arr(2) = 11
    atom_num_arr(2) = 1
    atom_loc_arr(2) = 2
    start_bas_arr(2) = 2; end_bas_arr(2) = 2
    iloc = 11
    exponents_arr(iloc) = 1.752_dp
    contr_coef_s_arr(iloc) = 1.0_dp*gto_norm_s(1.752_dp)

    ! Shell 3: S 10 primitives (2nd general contraction)
    ang_mom_arr(3) = 1
    contr_num_arr(3) = 10
    sh_loc_arr(3) = 12
    atom_num_arr(3) = 1
    atom_loc_arr(3) = 3
    start_bas_arr(3) = 3; end_bas_arr(3) = 3
    iloc = 12
    exponents_arr(iloc) = 15330.0_dp; contr_coef_s_arr(iloc) = -0.000115_dp*gto_norm_s(15330.0_dp)
    exponents_arr(iloc + 1) = 2299.0_dp; contr_coef_s_arr(iloc + 1) = -0.000895_dp*gto_norm_s(2299.0_dp)
    exponents_arr(iloc + 2) = 522.4_dp; contr_coef_s_arr(iloc + 2) = -0.004636_dp*gto_norm_s(522.4_dp)
    exponents_arr(iloc + 3) = 147.3_dp; contr_coef_s_arr(iloc + 3) = -0.018724_dp*gto_norm_s(147.3_dp)
    exponents_arr(iloc + 4) = 47.55_dp; contr_coef_s_arr(iloc + 4) = -0.058463_dp*gto_norm_s(47.55_dp)
    exponents_arr(iloc + 5) = 16.76_dp; contr_coef_s_arr(iloc + 5) = -0.136463_dp*gto_norm_s(16.76_dp)
    exponents_arr(iloc + 6) = 6.207_dp; contr_coef_s_arr(iloc + 6) = -0.175740_dp*gto_norm_s(6.207_dp)
    exponents_arr(iloc + 7) = 1.752_dp; contr_coef_s_arr(iloc + 7) = 0.160934_dp*gto_norm_s(1.752_dp)
    exponents_arr(iloc + 8) = 0.6882_dp; contr_coef_s_arr(iloc + 8) = 0.603418_dp*gto_norm_s(0.6882_dp)
    exponents_arr(iloc + 9) = 0.2384_dp; contr_coef_s_arr(iloc + 9) = 0.378765_dp*gto_norm_s(0.2384_dp)

    ! Shell 4: S 1 primitive (uncontracted)
    ang_mom_arr(4) = 1
    contr_num_arr(4) = 1
    sh_loc_arr(4) = 22
    atom_num_arr(4) = 1
    atom_loc_arr(4) = 4
    start_bas_arr(4) = 4; end_bas_arr(4) = 4
    iloc = 22
    exponents_arr(iloc) = 0.2384_dp
    contr_coef_s_arr(iloc) = 1.0_dp*gto_norm_s(0.2384_dp)

    ! Shell 5: P 1 primitive (uncontracted)
    ang_mom_arr(5) = 2
    contr_num_arr(5) = 1
    sh_loc_arr(5) = 23
    atom_num_arr(5) = 1
    atom_loc_arr(5) = 5
    start_bas_arr(5) = 5; end_bas_arr(5) = 7
    iloc = 23
    exponents_arr(iloc) = 0.7156_dp
    contr_coef_p_arr(iloc) = 1.0_dp*gto_norm_p(0.7156_dp)

    ! Shell 6: P 5 primitives (general contraction)
    ang_mom_arr(6) = 2
    contr_num_arr(6) = 5
    sh_loc_arr(6) = 24
    atom_num_arr(6) = 1
    atom_loc_arr(6) = 8
    start_bas_arr(6) = 8; end_bas_arr(6) = 10
    iloc = 24
    exponents_arr(iloc) = 34.46_dp; contr_coef_p_arr(iloc) = 0.015928_dp*gto_norm_p(34.46_dp)
    exponents_arr(iloc + 1) = 7.749_dp; contr_coef_p_arr(iloc + 1) = 0.09974_dp*gto_norm_p(7.749_dp)
    exponents_arr(iloc + 2) = 2.280_dp; contr_coef_p_arr(iloc + 2) = 0.310492_dp*gto_norm_p(2.280_dp)
    exponents_arr(iloc + 3) = 0.7156_dp; contr_coef_p_arr(iloc + 3) = 0.491026_dp*gto_norm_p(0.7156_dp)
    exponents_arr(iloc + 4) = 0.214_dp; contr_coef_p_arr(iloc + 4) = 0.336337_dp*gto_norm_p(0.214_dp)

    ! Shell 7: P 1 primitive (uncontracted)
    ang_mom_arr(7) = 2
    contr_num_arr(7) = 1
    sh_loc_arr(7) = 29
    atom_num_arr(7) = 1
    atom_loc_arr(7) = 11
    start_bas_arr(7) = 11; end_bas_arr(7) = 13
    iloc = 29
    exponents_arr(iloc) = 0.214_dp
    contr_coef_p_arr(iloc) = 1.0_dp*gto_norm_p(0.214_dp)

    ! Shell 8: D 1 primitive
    ang_mom_arr(8) = 3
    contr_num_arr(8) = 1
    sh_loc_arr(8) = 30
    atom_num_arr(8) = 1
    atom_loc_arr(8) = 14
    start_bas_arr(8) = 14; end_bas_arr(8) = 19
    iloc = 30
    exponents_arr(iloc) = 2.314_dp
    contr_coef_d_arr(iloc) = 1.0_dp*gto_norm_d(2.314_dp)

    ! Shell 9: D 1 primitive
    ang_mom_arr(9) = 3
    contr_num_arr(9) = 1
    sh_loc_arr(9) = 31
    atom_num_arr(9) = 1
    atom_loc_arr(9) = 20
    start_bas_arr(9) = 20; end_bas_arr(9) = 25
    iloc = 31
    exponents_arr(iloc) = 0.645_dp
    contr_coef_d_arr(iloc) = 1.0_dp*gto_norm_d(0.645_dp)

    ! Shell 10: F 1 primitive
    ang_mom_arr(10) = 4
    contr_num_arr(10) = 1
    sh_loc_arr(10) = 32
    atom_num_arr(10) = 1
    atom_loc_arr(10) = 26
    start_bas_arr(10) = 26; end_bas_arr(10) = 35
    iloc = 32
    exponents_arr(iloc) = 1.428_dp
    contr_coef_f_arr(iloc) = 1.0_dp*gto_norm_f(1.428_dp)

  end subroutine setup_oxygen_basis

  subroutine setup_hydrogen_basis(atom_id, shell_start)
    integer, intent(in) :: atom_id, shell_start
    integer :: sh, iloc, bas_start

    ! Calculate starting basis function index
    if (atom_id == 2) then
      bas_start = 36  ! After O's 35 basis functions
    else
      bas_start = 51  ! After O + H1
    end if

    sh = shell_start

    ! Shell: S 1 primitive (uncontracted)
    ang_mom_arr(sh) = 1
    contr_num_arr(sh) = 1
    sh_loc_arr(sh) = 33 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start
    start_bas_arr(sh) = bas_start; end_bas_arr(sh) = bas_start
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 0.3258_dp
    contr_coef_s_arr(iloc) = 1.0_dp*gto_norm_s(0.3258_dp)
    sh = sh + 1

    ! Shell: S 5 primitives (general contraction)
    ang_mom_arr(sh) = 1
    contr_num_arr(sh) = 5
    sh_loc_arr(sh) = 34 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start + 1
    start_bas_arr(sh) = bas_start + 1; end_bas_arr(sh) = bas_start + 1
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 33.87_dp; contr_coef_s_arr(iloc) = 0.006068_dp*gto_norm_s(33.87_dp)
    exponents_arr(iloc + 1) = 5.095_dp; contr_coef_s_arr(iloc + 1) = 0.045308_dp*gto_norm_s(5.095_dp)
    exponents_arr(iloc + 2) = 1.159_dp; contr_coef_s_arr(iloc + 2) = 0.202822_dp*gto_norm_s(1.159_dp)
    exponents_arr(iloc + 3) = 0.3258_dp; contr_coef_s_arr(iloc + 3) = 0.503903_dp*gto_norm_s(0.3258_dp)
    exponents_arr(iloc + 4) = 0.1027_dp; contr_coef_s_arr(iloc + 4) = 0.383421_dp*gto_norm_s(0.1027_dp)
    sh = sh + 1

    ! Shell: S 1 primitive (uncontracted)
    ang_mom_arr(sh) = 1
    contr_num_arr(sh) = 1
    sh_loc_arr(sh) = 39 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start + 2
    start_bas_arr(sh) = bas_start + 2; end_bas_arr(sh) = bas_start + 2
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 0.1027_dp
    contr_coef_s_arr(iloc) = 1.0_dp*gto_norm_s(0.1027_dp)
    sh = sh + 1

    ! Shell: P 1 primitive
    ang_mom_arr(sh) = 2
    contr_num_arr(sh) = 1
    sh_loc_arr(sh) = 40 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start + 3
    start_bas_arr(sh) = bas_start + 3; end_bas_arr(sh) = bas_start + 5
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 1.407_dp
    contr_coef_p_arr(iloc) = 1.0_dp*gto_norm_p(1.407_dp)
    sh = sh + 1

    ! Shell: P 1 primitive
    ang_mom_arr(sh) = 2
    contr_num_arr(sh) = 1
    sh_loc_arr(sh) = 41 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start + 6
    start_bas_arr(sh) = bas_start + 6; end_bas_arr(sh) = bas_start + 8
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 0.388_dp
    contr_coef_p_arr(iloc) = 1.0_dp*gto_norm_p(0.388_dp)
    sh = sh + 1

    ! Shell: D 1 primitive
    ang_mom_arr(sh) = 3
    contr_num_arr(sh) = 1
    sh_loc_arr(sh) = 42 + (atom_id - 2)*10
    atom_num_arr(sh) = atom_id
    atom_loc_arr(sh) = bas_start + 9
    start_bas_arr(sh) = bas_start + 9; end_bas_arr(sh) = bas_start + 14
    iloc = sh_loc_arr(sh)
    exponents_arr(iloc) = 1.057_dp
    contr_coef_d_arr(iloc) = 1.0_dp*gto_norm_d(1.057_dp)

  end subroutine setup_hydrogen_basis

  ! GAMESS normalization for s-type primitive: N = (2a/pi)^(3/4)
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp*alpha/pi_val)**0.75_dp
  end function

  ! GAMESS normalization for p-type primitive: N = 2*sqrt(a) * (2a/pi)^(3/4)
  pure function gto_norm_p(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 2.0_dp*sqrt(alpha)*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

  ! GAMESS normalization for d-type primitive: N = 4*a * (2a/pi)^(3/4)
  pure function gto_norm_d(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 4.0_dp*alpha*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

  ! GAMESS normalization for f-type primitive: N = 8*a^(3/2) * (2a/pi)^(3/4)
  pure function gto_norm_f(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 8.0_dp*alpha*sqrt(alpha)*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

end program test_h2o_ccpvtz
