! Standalone test for libERI: H2 with S-functions only
!
! Exercises ONLY the (ss|ss) kernel (int0000) and the Boys function.
! Uses cc-pVDZ s-shells (no p, no d, no f).
!
! Basis per H:
!   S  3  (contracted):  exp = 13.01, 1.962, 0.4446
!   S  1  (diffuse):     exp = 0.122
!
! 4 shells, 4 basis functions, tri_size = 10
!
! Usage:
!   mpirun -np 1 ./test_h2_sonly             # validate
!   mpirun -np 1 ./test_h2_sonly --generate  # print Fock values
!
program test_h2_sonly
  use liberi_types, only: dp, int64, int32
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_setup, &
                              liberi_fock_build, liberi_cleanup, liberi_destroy
  implicit none

  type(liberi_handle_t) :: handle

  integer, parameter :: nsh_val = 4
  integer, parameter :: natoms_val = 2
  integer, parameter :: num_bas_val = 4
  integer, parameter :: mxgtot_val = 3
  integer, parameter :: n_tri = 10

  integer :: ang_mom(nsh_val)
  integer :: contr_num_arr(nsh_val)
  integer :: sh_loc(nsh_val)
  integer :: atom_num(nsh_val)
  integer :: atom_loc(nsh_val)
  integer :: start_bas(nsh_val)
  integer :: end_bas(nsh_val)
  real(dp) :: coords(3*natoms_val)
  real(dp) :: exponents(mxgtot_val*nsh_val)
  real(dp) :: contr_coef_s(mxgtot_val*nsh_val)
  real(dp) :: contr_coef_p(mxgtot_val*nsh_val)
  real(dp) :: contr_coef_d(mxgtot_val*nsh_val)
  real(dp) :: contr_coef_f(mxgtot_val*nsh_val)
  real(dp) :: schwrz_int(nsh_val*(nsh_val + 1)/2)
  real(dp) :: density(n_tri)
  real(dp) :: fa(n_tri)

  integer(int64) :: i
  integer :: nargs
  character(len=256) :: arg
  logical :: generate
  real(dp) :: max_err, err
  real(dp), parameter :: tol = 1.0d-10
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp

  ! Reference Fock matrix: H2/cc-pVDZ s-only, identity density, no screening
  ! Generated with NVHPC (known-good compiler)
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
    call get_command_argument(int(i), arg)
    if (trim(arg) == '--generate') generate = .true.
  end do

  ! ============================================
  ! 0-initialize coefficient arrays
  ! ============================================
  exponents = 0.0_dp
  contr_coef_s = 0.0_dp
  contr_coef_p = 0.0_dp
  contr_coef_d = 0.0_dp
  contr_coef_f = 0.0_dp

  ! ============================================
  ! Atomic coordinates (bohr)
  ! H1 at (0, 0, -0.7), H2 at (0, 0, +0.7)
  ! ============================================
  coords(1) = 0.0_dp; coords(2) = 0.0_dp; coords(3) = -0.7_dp
  coords(4) = 0.0_dp; coords(5) = 0.0_dp; coords(6) = 0.7_dp

  ! ============================================
  ! Basis set: cc-pVDZ s-shells for H (atom 1)
  ! ============================================

  ! Shell 1: s, 3 primitives (contracted)
  ang_mom(1) = 1      ! 1=s in libERI convention
  contr_num_arr(1) = 3
  sh_loc(1) = 1
  atom_num(1) = 1
  atom_loc(1) = 1
  start_bas(1) = 1; end_bas(1) = 1
  exponents(1) = 13.0100000_dp
  exponents(2) = 1.9620000_dp
  exponents(3) = 0.4446000_dp
  contr_coef_s(1) = 0.0196850_dp*gto_norm_s(13.0100000_dp)
  contr_coef_s(2) = 0.1379770_dp*gto_norm_s(1.9620000_dp)
  contr_coef_s(3) = 0.4781480_dp*gto_norm_s(0.4446000_dp)

  ! Shell 2: s, 1 primitive (diffuse)
  ang_mom(2) = 1
  contr_num_arr(2) = 1
  sh_loc(2) = 4
  atom_num(2) = 1
  atom_loc(2) = 2
  start_bas(2) = 2; end_bas(2) = 2
  exponents(4) = 0.1220000_dp
  contr_coef_s(4) = 1.0000000_dp*gto_norm_s(0.1220000_dp)

  ! ============================================
  ! Basis set: cc-pVDZ s-shells for H (atom 2)
  ! ============================================

  ! Shell 3: s, 3 primitives (contracted)
  ang_mom(3) = 1
  contr_num_arr(3) = 3
  sh_loc(3) = 5
  atom_num(3) = 2
  atom_loc(3) = 3
  start_bas(3) = 3; end_bas(3) = 3
  exponents(5) = 13.0100000_dp
  exponents(6) = 1.9620000_dp
  exponents(7) = 0.4446000_dp
  contr_coef_s(5) = 0.0196850_dp*gto_norm_s(13.0100000_dp)
  contr_coef_s(6) = 0.1379770_dp*gto_norm_s(1.9620000_dp)
  contr_coef_s(7) = 0.4781480_dp*gto_norm_s(0.4446000_dp)

  ! Shell 4: s, 1 primitive (diffuse)
  ang_mom(4) = 1
  contr_num_arr(4) = 1
  sh_loc(4) = 8
  atom_num(4) = 2
  atom_loc(4) = 4
  start_bas(4) = 4; end_bas(4) = 4
  exponents(8) = 0.1220000_dp
  contr_coef_s(8) = 1.0000000_dp*gto_norm_s(0.1220000_dp)

  ! ============================================
  ! Disable Schwarz screening (set large values)
  ! ============================================
  schwrz_int = 1.0d10

  ! ============================================
  ! Set up handle and compute shell pairs
  ! ============================================
  print *, "Computing shell pairs..."
  call liberi_create(handle)
  call liberi_setup(handle, nsh_val, natoms_val, num_bas_val, mxgtot_val, &
                    ang_mom, contr_num_arr, sh_loc, &
                    atom_num, atom_loc, start_bas, end_bas, &
                    exponents, contr_coef_s, contr_coef_p, &
                    contr_coef_d, contr_coef_f, &
                    coords, schwrz_int, 0, 1)

  ! ============================================
  ! Density matrix — identity (diagonal = 1)
  ! ============================================
  density = 0.0_dp
  do i = 1, num_bas_val
    density(i*(i + 1)/2) = 1.0_dp
  end do
  fa = 0.0_dp

  ! ============================================
  ! Run integral driver (only int0000 will fire)
  ! ============================================
  print *, "Running integral driver (s-only, int0000 only)..."
  call liberi_fock_build(handle, density, fa, n_tri)

  ! ============================================
  ! Output results
  ! ============================================
  print *, ""
  print *, "Fock matrix (triangular packed, 10 elements):"
  do i = 1, n_tri
    print '(A,I3,A,ES24.17)', "  fa(", i, ") = ", fa(i)
  end do

  if (generate) then
    print *, ""
    print *, "Copy-paste for fa_ref:"
    do i = 1, n_tri
      if (i < n_tri) then
        print '(A,ES26.17,A,I0)', "                         ", fa(i), ", &  !  ", i
      else
        print '(A,ES26.17,A,I0)', "                         ", fa(i), " &   !  ", i
      end if
    end do
  end if

  if (.not. generate) then
    max_err = 0.0_dp
    do i = 1, n_tri
      err = abs(fa(i) - fa_ref(i))
      if (err > max_err) max_err = err
      if (err > tol) then
        print '(A,I3,A,ES24.15,A,ES24.15,A,ES10.2)', &
          "  FAIL fa(", i, ") = ", fa(i), "  ref = ", fa_ref(i), "  err = ", err
      end if
    end do

    print *, ""
    if (max_err <= tol) then
      print '(A,ES10.2)', " PASS - max error: ", max_err
    else
      print '(A,ES10.2)', " FAIL - max error: ", max_err
    error stop
    end if
  end if

  ! ============================================
  ! Cleanup
  ! ============================================
  call liberi_destroy(handle)

contains

  ! GAMESS normalization for s-type primitive: N = (2a/pi)^(3/4)
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp*alpha/pi_val)**0.75_dp
  end function

end program test_h2_sonly
