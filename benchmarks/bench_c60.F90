! Performance benchmark for libERI: C60 / Fullerene
!
! Three basis sets: STO-3G (s,p), 6-31G (s,p), cc-pVDZ (s,p,d).
! Coordinates from libcint time_c60.c (Bohr).
! This is a manual benchmark for timing — NOT a ctest test.
!
! Basis summary (60 Carbon atoms):
!   STO-3G:  180 shells,  300 bf (s,p only)
!   6-31G:   300 shells,  540 bf (s,p only)
!   cc-pVDZ: 360 shells,  900 bf (s,p,d)
!
! Usage:
!   mpirun -np 1 ./bench_c60              (runs all three)
!   mpirun -np 1 ./bench_c60 --sto        (STO-3G only)
!   mpirun -np 1 ./bench_c60 --631g       (6-31G only)
!   mpirun -np 1 ./bench_c60 --ccd        (cc-pVDZ only)
!   mpirun -np 1 ./bench_c60 --sto --ccd  (STO-3G + cc-pVDZ)
!
program bench_c60
  use omp_lib
  use liberi_types, only: dp, int64, int32
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_setup, &
                              liberi_fock_build, liberi_cleanup, liberi_destroy
  implicit none

  integer :: iat, nargs, iarg
  integer(int64) :: i, tri_size, bf_idx, sh_idx, exp_idx
  real(dp) :: t_start, t_end
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp
  logical :: run_sto, run_631g, run_ccd
  character(len=32) :: arg

  ! C60 coordinates (Bohr) — from libcint time_c60.c
  integer, parameter :: nat = 60
  real(dp), parameter :: xyz(3, nat) = reshape([ &
                                               -0.779_dp, 0.757_dp, 0.000_dp, &  !  1
                                               0.684_dp, 0.757_dp, 0.000_dp, &  !  2
                                               1.376_dp, 1.956_dp, 0.000_dp, &  !  3
                                               0.644_dp, 3.224_dp, -0.000_dp, &  !  4
                                               -0.739_dp, 3.224_dp, 0.000_dp, &  !  5
                                               -1.471_dp, 1.956_dp, 0.000_dp, &  !  6
                                               -1.231_dp, -0.348_dp, -0.845_dp, &  !  7
                                               -0.047_dp, -1.032_dp, -1.367_dp, &  !  8
                                               1.137_dp, -0.348_dp, -0.845_dp, &  !  9
                                               2.257_dp, -0.196_dp, -1.644_dp, &  ! 10
                                               2.561_dp, 2.118_dp, -0.845_dp, &  ! 11
                                               1.376_dp, 4.169_dp, -0.845_dp, &  ! 12
                                               0.684_dp, 5.063_dp, -1.644_dp, &  ! 13
                                               -0.779_dp, 5.063_dp, -1.644_dp, &  ! 14
                                               -1.471_dp, 4.169_dp, -0.845_dp, &  ! 15
                                               -2.656_dp, 3.485_dp, -1.367_dp, &  ! 16
                                               -2.655_dp, 2.118_dp, -0.845_dp, &  ! 17
                                               -3.083_dp, 1.071_dp, -1.644_dp, &  ! 18
                                               -2.351_dp, -0.196_dp, -1.644_dp, &  ! 19
                                               -0.047_dp, -1.526_dp, -2.661_dp, &  ! 20
                                               -1.231_dp, -1.365_dp, -3.506_dp, &  ! 21
                                               -2.351_dp, -0.718_dp, -3.012_dp, &  ! 22
                                               -3.083_dp, 0.226_dp, -3.857_dp, &  ! 23
                                               -3.536_dp, 1.332_dp, -3.012_dp, &  ! 24
                                               -3.536_dp, 2.626_dp, -3.506_dp, &  ! 25
                                               -3.083_dp, 3.732_dp, -2.660_dp, &  ! 26
                                               -2.352_dp, 4.677_dp, -3.506_dp, &  ! 27
                                               -1.231_dp, 5.324_dp, -3.012_dp, &  ! 28
                                               2.561_dp, 3.485_dp, -1.367_dp, &  ! 29
                                               1.136_dp, 4.307_dp, -5.673_dp, &  ! 30
                                               2.256_dp, 4.155_dp, -4.873_dp, &  ! 31
                                               2.988_dp, 2.887_dp, -4.873_dp, &  ! 32
                                               2.561_dp, 1.841_dp, -5.673_dp, &  ! 33
                                               1.376_dp, 2.002_dp, -6.518_dp, &  ! 34
                                               -0.779_dp, 3.201_dp, -6.518_dp, &  ! 35
                                               -1.231_dp, 4.307_dp, -5.673_dp, &  ! 36
                                               -0.047_dp, 4.991_dp, -5.150_dp, &  ! 37
                                               -0.047_dp, 5.485_dp, -3.857_dp, &  ! 38
                                               1.136_dp, 5.324_dp, -3.012_dp, &  ! 39
                                               2.257_dp, 4.677_dp, -3.506_dp, &  ! 40
                                               3.441_dp, 2.626_dp, -3.506_dp, &  ! 41
                                               3.441_dp, 1.332_dp, -3.012_dp, &  ! 42
                                               2.989_dp, 0.226_dp, -3.857_dp, &  ! 43
                                               2.561_dp, 0.473_dp, -5.150_dp, &  ! 44
                                               1.376_dp, -0.210_dp, -5.673_dp, &  ! 45
                                               0.644_dp, 0.734_dp, -6.518_dp, &  ! 46
                                               -0.739_dp, 0.734_dp, -6.518_dp, &  ! 47
                                               -1.471_dp, 2.002_dp, -6.518_dp, &  ! 48
                                               -2.352_dp, 4.155_dp, -4.873_dp, &  ! 49
                                               -3.084_dp, 2.887_dp, -4.873_dp, &  ! 50
                                               -2.656_dp, 1.840_dp, -5.673_dp, &  ! 51
                                               -2.656_dp, 0.473_dp, -5.150_dp, &  ! 52
                                               -1.471_dp, -0.210_dp, -5.673_dp, &  ! 53
                                               -0.779_dp, -1.104_dp, -4.873_dp, &  ! 54
                                               0.684_dp, -1.104_dp, -4.873_dp, &  ! 55
                                               1.136_dp, -1.365_dp, -3.506_dp, &  ! 56
                                               2.257_dp, -0.718_dp, -3.012_dp, &  ! 57
                                               2.988_dp, 3.732_dp, -2.661_dp, &  ! 58
                                               2.989_dp, 1.071_dp, -1.644_dp, &  ! 59
                                               0.684_dp, 3.201_dp, -6.518_dp &  ! 60
                                               ], [3, nat])

  ! Local basis arrays (sized per setup routine)
  integer :: natoms, nsh, num_bas, mxgtot
  integer, allocatable :: ang_mom(:), contr_num(:), sh_loc(:)
  integer, allocatable :: atom_num(:), atom_loc(:), start_bas(:), end_bas(:)
  real(dp), allocatable :: coords(:), exponents(:)
  real(dp), allocatable :: contr_coef_s(:), contr_coef_p(:), contr_coef_d(:), contr_coef_f(:)
  real(dp), allocatable :: schwrz_int(:), density(:), fa(:)

  ! ============================================
  ! Parse command-line arguments
  ! ============================================
  run_sto = .false.
  run_631g = .false.
  run_ccd = .false.

  nargs = command_argument_count()
  do iarg = 1, nargs
    call get_command_argument(iarg, arg)
    select case (trim(arg))
    case ('--sto')
      run_sto = .true.
    case ('--631g')
      run_631g = .true.
    case ('--ccd')
      run_ccd = .true.
    case default
      write (*, '(A,A)') "Unknown argument: ", trim(arg)
      stop 1
    end select
  end do

  ! No flags => run all
  if (.not. run_sto .and. .not. run_631g .and. .not. run_ccd) then
    run_sto = .true.
    run_631g = .true.
    run_ccd = .true.
  end if

  write (*, '(A)') "=== libERI Benchmark: C60 ==="
  write (*, '(A,I0)') "Threads:  ", omp_get_max_threads()
  write (*, '(A)') ""

  ! ============================================
  ! Benchmark 1: STO-3G
  ! ============================================
  if (run_sto) then
    call setup_sto3g()
    call run_benchmark("STO-3G")
    call cleanup_all()
  end if

  ! ============================================
  ! Benchmark 2: 6-31G
  ! ============================================
  if (run_631g) then
    call setup_631g()
    call run_benchmark("6-31G")
    call cleanup_all()
  end if

  ! ============================================
  ! Benchmark 3: cc-pVDZ
  ! ============================================
  if (run_ccd) then
    call setup_ccpvdz()
    call run_benchmark("cc-pVDZ")
    call cleanup_all()
  end if

contains

  ! ============================================
  ! Run benchmark: liberi_setup + liberi_fock_build
  ! ============================================
  subroutine run_benchmark(basis_name)
    character(len=*), intent(in) :: basis_name
    type(liberi_handle_t) :: handle

    write (*, '(A,A,A)') "Computing shell pairs (", basis_name, ")..."
    call liberi_create(handle)
    call liberi_setup(handle, nsh, natoms, num_bas, mxgtot, &
                      ang_mom, contr_num, sh_loc, &
                      atom_num, atom_loc, start_bas, end_bas, &
                      exponents, contr_coef_s, contr_coef_p, &
                      contr_coef_d, contr_coef_f, &
                      coords, schwrz_int, 0, 1)

    ! 0.0_dp density — pure timing
    density = 0.0_dp
    fa = 0.0_dp

    write (*, '(A,A,A)') "Running integral driver (", basis_name, ")..."
    t_start = omp_get_wtime()
    call liberi_fock_build(handle, density, fa, int(tri_size))
    t_end = omp_get_wtime()

    write (*, '(A)') ""
    write (*, '(A,A,A)') "--- ", basis_name, " ---"
    write (*, '(A,I0)') "  Atoms:           ", natoms
    write (*, '(A,I0)') "  Shells:          ", nsh
    write (*, '(A,I0)') "  Basis functions:  ", num_bas
    write (*, '(A,F12.3,A)') "  Wall time:       ", t_end - t_start, " s"
    write (*, '(A)') ""

    call liberi_destroy(handle)
  end subroutine run_benchmark

  ! ============================================
  ! Fill coordinates into coords array
  ! ============================================
  subroutine fill_coords()
    do i = 1, natoms
      coords(3*(i - 1) + 1) = xyz(1, int(i))
      coords(3*(i - 1) + 2) = xyz(2, int(i))
      coords(3*(i - 1) + 3) = xyz(3, int(i))
    end do
  end subroutine fill_coords

  ! ============================================
  ! Setup STO-3G basis for C60
  ! ============================================
  subroutine setup_sto3g()
    ! Carbon STO-3G: 3 shells (s3, s3, p3) => 5 bf per C
    ! Total: 60*3 = 180 shells, 60*5 = 300 bf, mxgtot = 3
    natoms = nat
    nsh = 180
    num_bas = 300
    mxgtot = 3

    call allocate_arrays()
    call fill_coords()

    sh_idx = 0
    bf_idx = 0
    exp_idx = 0

    do iat = 1, nat
      call add_c_sto3g(iat)
    end do

    schwrz_int = 1.0d10
  end subroutine setup_sto3g

  ! ============================================
  ! Add Carbon STO-3G basis
  ! ============================================
  subroutine add_c_sto3g(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 3 primitives (inner)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 71.6168370_dp
    exponents(base + 2) = 13.0450960_dp
    exponents(base + 3) = 3.5305122_dp
    contr_coef_s(base + 1) = 0.15432897_dp*gto_norm_s(71.6168370_dp)
    contr_coef_s(base + 2) = 0.53532814_dp*gto_norm_s(13.0450960_dp)
    contr_coef_s(base + 3) = 0.44463454_dp*gto_norm_s(3.5305122_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 2: s from SP, 3 primitives
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 2.9412494_dp
    exponents(base + 2) = 0.6834831_dp
    exponents(base + 3) = 0.2222899_dp
    contr_coef_s(base + 1) = -0.09996723_dp*gto_norm_s(2.9412494_dp)
    contr_coef_s(base + 2) = 0.39951283_dp*gto_norm_s(0.6834831_dp)
    contr_coef_s(base + 3) = 0.70011547_dp*gto_norm_s(0.2222899_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 3: p from SP, 3 primitives
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 2
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 2
    base = exp_idx
    exponents(base + 1) = 2.9412494_dp
    exponents(base + 2) = 0.6834831_dp
    exponents(base + 3) = 0.2222899_dp
    contr_coef_p(base + 1) = 0.15591627_dp*gto_norm_p(2.9412494_dp)
    contr_coef_p(base + 2) = 0.60768372_dp*gto_norm_p(0.6834831_dp)
    contr_coef_p(base + 3) = 0.39195739_dp*gto_norm_p(0.2222899_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot
  end subroutine add_c_sto3g

  ! ============================================
  ! Setup 6-31G basis for C60
  ! ============================================
  subroutine setup_631g()
    ! Carbon 6-31G: 5 shells (s6, s3, s1, p3, p1) => 9 bf per C
    ! Total: 60*5 = 300 shells, 60*9 = 540 bf, mxgtot = 6
    natoms = nat
    nsh = 300
    num_bas = 540
    mxgtot = 6

    call allocate_arrays()
    call fill_coords()

    sh_idx = 0
    bf_idx = 0
    exp_idx = 0

    do iat = 1, nat
      call add_c_631g(iat)
    end do

    schwrz_int = 1.0d10
  end subroutine setup_631g

  ! ============================================
  ! Add Carbon 6-31G basis
  ! ============================================
  subroutine add_c_631g(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 6 primitives (inner)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 6
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 3047.5249000_dp
    exponents(base + 2) = 457.3695100_dp
    exponents(base + 3) = 103.9486900_dp
    exponents(base + 4) = 29.2101550_dp
    exponents(base + 5) = 9.2866630_dp
    exponents(base + 6) = 3.1639270_dp
    contr_coef_s(base + 1) = 0.0018347_dp*gto_norm_s(3047.5249000_dp)
    contr_coef_s(base + 2) = 0.0140373_dp*gto_norm_s(457.3695100_dp)
    contr_coef_s(base + 3) = 0.0688426_dp*gto_norm_s(103.9486900_dp)
    contr_coef_s(base + 4) = 0.2321844_dp*gto_norm_s(29.2101550_dp)
    contr_coef_s(base + 5) = 0.4679413_dp*gto_norm_s(9.2866630_dp)
    contr_coef_s(base + 6) = 0.3623120_dp*gto_norm_s(3.1639270_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 2: s from SP, 3 primitives
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 7.8682724_dp
    exponents(base + 2) = 1.8812885_dp
    exponents(base + 3) = 0.5442493_dp
    contr_coef_s(base + 1) = -0.1193324_dp*gto_norm_s(7.8682724_dp)
    contr_coef_s(base + 2) = -0.1608542_dp*gto_norm_s(1.8812885_dp)
    contr_coef_s(base + 3) = 1.1434564_dp*gto_norm_s(0.5442493_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 3: s, 1 primitive (diffuse)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 1
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 0.1687144_dp
    contr_coef_s(base + 1) = 1.0000000_dp*gto_norm_s(0.1687144_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 4: p from SP, 3 primitives
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 2
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 2
    base = exp_idx
    exponents(base + 1) = 7.8682724_dp
    exponents(base + 2) = 1.8812885_dp
    exponents(base + 3) = 0.5442493_dp
    contr_coef_p(base + 1) = 0.0689991_dp*gto_norm_p(7.8682724_dp)
    contr_coef_p(base + 2) = 0.3164240_dp*gto_norm_p(1.8812885_dp)
    contr_coef_p(base + 3) = 0.7443083_dp*gto_norm_p(0.5442493_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot

    ! Shell 5: p, 1 primitive (diffuse)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 2
    contr_num(sh_idx) = 1
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 2
    base = exp_idx
    exponents(base + 1) = 0.1687144_dp
    contr_coef_p(base + 1) = 1.0000000_dp*gto_norm_p(0.1687144_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot
  end subroutine add_c_631g

  ! ============================================
  ! Setup cc-pVDZ basis for C60
  ! ============================================
  subroutine setup_ccpvdz()
    ! Carbon cc-pVDZ: 6 shells (s8, s8, s1, p3, p1, d1) => 15 bf per C
    ! Total: 60*6 = 360 shells, 60*15 = 900 bf, mxgtot = 8
    natoms = nat
    nsh = 360
    num_bas = 900
    mxgtot = 8

    call allocate_arrays()
    call fill_coords()

    sh_idx = 0
    bf_idx = 0
    exp_idx = 0

    do iat = 1, nat
      call add_c_ccpvdz(iat)
    end do

    schwrz_int = 1.0d10
  end subroutine setup_ccpvdz

  ! ============================================
  ! Add Carbon cc-pVDZ basis (general contraction split)
  ! ============================================
  subroutine add_c_ccpvdz(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 8 primitives (general contraction, 1st function)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 8
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 6665.0000000_dp
    exponents(base + 2) = 1000.0000000_dp
    exponents(base + 3) = 228.0000000_dp
    exponents(base + 4) = 64.7100000_dp
    exponents(base + 5) = 21.0600000_dp
    exponents(base + 6) = 7.4950000_dp
    exponents(base + 7) = 2.7970000_dp
    exponents(base + 8) = 0.5215000_dp
    contr_coef_s(base + 1) = 0.000692_dp*gto_norm_s(6665.0000000_dp)
    contr_coef_s(base + 2) = 0.005329_dp*gto_norm_s(1000.0000000_dp)
    contr_coef_s(base + 3) = 0.027077_dp*gto_norm_s(228.0000000_dp)
    contr_coef_s(base + 4) = 0.101718_dp*gto_norm_s(64.7100000_dp)
    contr_coef_s(base + 5) = 0.274740_dp*gto_norm_s(21.0600000_dp)
    contr_coef_s(base + 6) = 0.448564_dp*gto_norm_s(7.4950000_dp)
    contr_coef_s(base + 7) = 0.285074_dp*gto_norm_s(2.7970000_dp)
    contr_coef_s(base + 8) = 0.015204_dp*gto_norm_s(0.5215000_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 2: s, 8 primitives (general contraction, 2nd function)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 8
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 6665.0000000_dp
    exponents(base + 2) = 1000.0000000_dp
    exponents(base + 3) = 228.0000000_dp
    exponents(base + 4) = 64.7100000_dp
    exponents(base + 5) = 21.0600000_dp
    exponents(base + 6) = 7.4950000_dp
    exponents(base + 7) = 2.7970000_dp
    exponents(base + 8) = 0.5215000_dp
    contr_coef_s(base + 1) = -0.000146_dp*gto_norm_s(6665.0000000_dp)
    contr_coef_s(base + 2) = -0.001154_dp*gto_norm_s(1000.0000000_dp)
    contr_coef_s(base + 3) = -0.005725_dp*gto_norm_s(228.0000000_dp)
    contr_coef_s(base + 4) = -0.023312_dp*gto_norm_s(64.7100000_dp)
    contr_coef_s(base + 5) = -0.063955_dp*gto_norm_s(21.0600000_dp)
    contr_coef_s(base + 6) = -0.149981_dp*gto_norm_s(7.4950000_dp)
    contr_coef_s(base + 7) = -0.127262_dp*gto_norm_s(2.7970000_dp)
    contr_coef_s(base + 8) = 0.544529_dp*gto_norm_s(0.5215000_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 3: s, 1 primitive
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 1
    contr_num(sh_idx) = 1
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx
    base = exp_idx
    exponents(base + 1) = 0.1596000_dp
    contr_coef_s(base + 1) = 1.0000000_dp*gto_norm_s(0.1596000_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 4: p, 3 primitives
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 2
    contr_num(sh_idx) = 3
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 2
    base = exp_idx
    exponents(base + 1) = 9.4390000_dp
    exponents(base + 2) = 2.0020000_dp
    exponents(base + 3) = 0.5456000_dp
    contr_coef_p(base + 1) = 0.038109_dp*gto_norm_p(9.4390000_dp)
    contr_coef_p(base + 2) = 0.209480_dp*gto_norm_p(2.0020000_dp)
    contr_coef_p(base + 3) = 0.508557_dp*gto_norm_p(0.5456000_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot

    ! Shell 5: p, 1 primitive
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 2
    contr_num(sh_idx) = 1
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 2
    base = exp_idx
    exponents(base + 1) = 0.1517000_dp
    contr_coef_p(base + 1) = 1.0000000_dp*gto_norm_p(0.1517000_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot

    ! Shell 6: d, 1 primitive (6 Cartesian components)
    sh_idx = sh_idx + 1
    ang_mom(sh_idx) = 3
    contr_num(sh_idx) = 1
    sh_loc(sh_idx) = exp_idx + 1
    atom_num(sh_idx) = iatom
    bf_idx = bf_idx + 1
    atom_loc(sh_idx) = bf_idx
    start_bas(sh_idx) = bf_idx
    end_bas(sh_idx) = bf_idx + 5
    base = exp_idx
    exponents(base + 1) = 0.5500000_dp
    contr_coef_d(base + 1) = 1.0000000_dp*gto_norm_d(0.5500000_dp)
    bf_idx = bf_idx + 5
    exp_idx = exp_idx + mxgtot
  end subroutine add_c_ccpvdz

  ! ============================================
  ! Allocate local arrays
  ! ============================================
  subroutine allocate_arrays()
    allocate (ang_mom(nsh))
    allocate (contr_num(nsh))
    allocate (sh_loc(nsh))
    allocate (atom_num(nsh))
    allocate (atom_loc(nsh))
    allocate (start_bas(nsh))
    allocate (end_bas(nsh))
    allocate (coords(3*natoms))
    allocate (exponents(mxgtot*nsh))
    allocate (contr_coef_s(mxgtot*nsh))
    allocate (contr_coef_p(mxgtot*nsh))
    allocate (contr_coef_d(mxgtot*nsh))
    allocate (contr_coef_f(mxgtot*nsh))

    tri_size = num_bas*(num_bas + 1)/2
    allocate (schwrz_int(nsh*(nsh + 1)/2))
    allocate (density(tri_size))
    allocate (fa(tri_size))

    exponents = 0.0_dp
    contr_coef_s = 0.0_dp
    contr_coef_p = 0.0_dp
    contr_coef_d = 0.0_dp
    contr_coef_f = 0.0_dp
  end subroutine allocate_arrays

  ! ============================================
  ! Deallocate everything
  ! ============================================
  subroutine cleanup_all()
    if (allocated(ang_mom)) deallocate (ang_mom)
    if (allocated(contr_num)) deallocate (contr_num)
    if (allocated(sh_loc)) deallocate (sh_loc)
    if (allocated(atom_num)) deallocate (atom_num)
    if (allocated(atom_loc)) deallocate (atom_loc)
    if (allocated(start_bas)) deallocate (start_bas)
    if (allocated(end_bas)) deallocate (end_bas)
    if (allocated(coords)) deallocate (coords)
    if (allocated(exponents)) deallocate (exponents)
    if (allocated(contr_coef_s)) deallocate (contr_coef_s)
    if (allocated(contr_coef_p)) deallocate (contr_coef_p)
    if (allocated(contr_coef_d)) deallocate (contr_coef_d)
    if (allocated(contr_coef_f)) deallocate (contr_coef_f)
    if (allocated(schwrz_int)) deallocate (schwrz_int)
    if (allocated(density)) deallocate (density)
    if (allocated(fa)) deallocate (fa)
  end subroutine cleanup_all

  ! ============================================
  ! GAMESS normalization functions
  ! ============================================
  pure function gto_norm_s(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = (2.0_dp*alpha/pi_val)**0.75_dp
  end function

  pure function gto_norm_p(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 2.0_dp*sqrt(alpha)*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

  pure function gto_norm_d(alpha) result(norm)
    real(dp), intent(in) :: alpha
    real(dp) :: norm
    norm = 4.0_dp*alpha*(2.0_dp*alpha/pi_val)**0.75_dp
  end function

end program bench_c60
