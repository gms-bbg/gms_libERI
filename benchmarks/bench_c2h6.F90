! Performance benchmark for libERI: C2H6 / Ethane
!
! 2.0_dp basis sets: 6-31G (s,p) and cc-pVDZ (s,p,d).
! Coordinates from libcint fortran_time_c2h6.F90 (Bohr).
! This is a manual benchmark for timing — NOT a ctest test.
!
! Usage:
!   mpirun -np 1 ./bench_c2h6
!
program bench_c2h6
  use omp_lib
  use liberi_types, only: dp, int64, int32
  use liberi_interface, only: liberi_handle_t, liberi_create, liberi_setup, &
                              liberi_fock_build, liberi_cleanup, liberi_destroy
  implicit none

  integer(int64) :: i, tri_size, bf_idx, sh_idx, exp_idx
  real(dp) :: t_start, t_end
  real(dp), parameter :: pi_val = 3.14159265358979323846264338327950288_dp

  ! C2H6 coordinates (Bohr) — from libcint
  ! Atoms: C1, H1, H2, H3, C2, H4, H5, H6
  integer, parameter :: nat = 8
  real(dp), parameter :: xyz(3, nat) = reshape([ &
                                               0.000_dp, 0.000_dp, 0.769_dp, &  ! C1
                                               0.000_dp, 1.014_dp, 1.174_dp, &  ! H1
                                               -0.878_dp, -0.507_dp, 1.174_dp, &  ! H2
                                               0.878_dp, -0.507_dp, 1.174_dp, &  ! H3
                                               0.000_dp, 0.000_dp, -0.769_dp, &  ! C2
                                               0.000_dp, 1.014_dp, -1.174_dp, &  ! H4
                                               -0.878_dp, -0.507_dp, -1.174_dp, &  ! H5
                                               0.878_dp, -0.507_dp, -1.174_dp &  ! H6
                                               ], [3, nat])

  ! Local basis arrays (sized per setup routine)
  integer :: natoms, nsh, num_bas, mxgtot
  integer, allocatable :: ang_mom(:), contr_num(:), sh_loc(:)
  integer, allocatable :: atom_num(:), atom_loc(:), start_bas(:), end_bas(:)
  real(dp), allocatable :: coords(:), exponents(:)
  real(dp), allocatable :: contr_coef_s(:), contr_coef_p(:), contr_coef_d(:), contr_coef_f(:)
  real(dp), allocatable :: schwrz_int(:), density(:), fa(:)

  write (*, '(A)') "=== libERI Benchmark: C2H6 ==="
  write (*, '(A,I0)') "Threads:  ", omp_get_max_threads()
  write (*, '(A)') ""

  ! ============================================
  ! Benchmark 1: 6-31G
  ! ============================================
  call setup_631g()
  call run_benchmark("6-31G")
  call cleanup_all()

  ! ============================================
  ! Benchmark 2: cc-pVDZ
  ! ============================================
  call setup_ccpvdz()
  call run_benchmark("cc-pVDZ")
  call cleanup_all()

contains

  ! ============================================
  ! Run benchmark: liberi_setup + liberi_fock_build
  ! ============================================
  subroutine run_benchmark(basis_name)
    character(len=*), intent(in) :: basis_name
    type(liberi_handle_t) :: handle

    call liberi_create(handle)
    call liberi_setup(handle, nsh, natoms, num_bas, mxgtot, &
                      ang_mom, contr_num, sh_loc, &
                      atom_num, atom_loc, start_bas, end_bas, &
                      exponents, contr_coef_s, contr_coef_p, &
                      contr_coef_d, contr_coef_f, &
                      coords, schwrz_int, 0, 1)

    ! 0.0_dp density — pure timing, no correctness check
    density = 0.0_dp
    fa = 0.0_dp

    t_start = omp_get_wtime()
    call liberi_fock_build(handle, density, fa, int(tri_size))
    t_end = omp_get_wtime()

    write (*, '(A,A,A)') "--- ", basis_name, " ---"
    write (*, '(A,I0)') "  Shells:          ", nsh
    write (*, '(A,I0)') "  Basis functions:  ", num_bas
    write (*, '(A,F12.3,A)') "  Wall time:       ", t_end - t_start, " s"
    write (*, '(A)') ""

    call liberi_destroy(handle)
  end subroutine run_benchmark

  ! ============================================
  ! Setup 6-31G basis for C2H6
  ! ============================================
  subroutine setup_631g()
    ! Carbon 6-31G: 5 shells (s6, s3, s1, p3, p1) => 9 bf per C
    ! Hydrogen 6-31G: 2 shells (s3, s1) => 2 bf per H
    ! Total: 2*5 + 6*2 = 22 shells, 2*9 + 6*2 = 30 bf, mxgtot = 6

    natoms = nat
    nsh = 22
    num_bas = 30
    mxgtot = 6

    call allocate_arrays()

    ! Coordinates
    do i = 1, natoms
      coords(3*(i - 1) + 1) = xyz(1, i)
      coords(3*(i - 1) + 2) = xyz(2, i)
      coords(3*(i - 1) + 3) = xyz(3, i)
    end do

    ! Carbon 6-31G exponents and coefficients
    ! Shell layout per Carbon: s(6), s(3), s(1), p(3), p(1)
    sh_idx = 0
    bf_idx = 0
    exp_idx = 0

    ! --- Carbon 1 (atom 1) ---
    call add_c_631g(1)
    ! --- Carbon 2 (atom 5) ---
    call add_c_631g(5)

    ! --- Hydrogens (atoms 2,3,4,6,7,8) ---
    call add_h_631g(2)
    call add_h_631g(3)
    call add_h_631g(4)
    call add_h_631g(6)
    call add_h_631g(7)
    call add_h_631g(8)

    schwrz_int = 1.0d10
  end subroutine setup_631g

  ! ============================================
  ! Add Carbon 6-31G basis
  ! ============================================
  subroutine add_c_631g(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 6 primitives
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

    ! Shell 2: s, 3 primitives (from SP shell)
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

    ! Shell 4: p, 3 primitives (from SP shell)
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
    bf_idx = bf_idx + 2   ! p = 3 functions, already counted 1
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
    bf_idx = bf_idx + 2   ! p = 3 functions
    exp_idx = exp_idx + mxgtot
  end subroutine add_c_631g

  ! ============================================
  ! Add Hydrogen 6-31G basis
  ! ============================================
  subroutine add_h_631g(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 3 primitives
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
    exponents(base + 1) = 18.7311370_dp
    exponents(base + 2) = 2.8253937_dp
    exponents(base + 3) = 0.6401217_dp
    contr_coef_s(base + 1) = 0.0334946_dp*gto_norm_s(18.7311370_dp)
    contr_coef_s(base + 2) = 0.2347269_dp*gto_norm_s(2.8253937_dp)
    contr_coef_s(base + 3) = 0.8137573_dp*gto_norm_s(0.6401217_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 2: s, 1 primitive (diffuse)
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
    exponents(base + 1) = 0.1612778_dp
    contr_coef_s(base + 1) = 1.0000000_dp*gto_norm_s(0.1612778_dp)
    exp_idx = exp_idx + mxgtot
  end subroutine add_h_631g

  ! ============================================
  ! Setup cc-pVDZ basis for C2H6
  ! ============================================
  subroutine setup_ccpvdz()
    ! Carbon cc-pVDZ: 6 shells (s8+s8 gen contr, s1, p3, p1, d1) => 15 bf per C
    ! Hydrogen cc-pVDZ: 3 shells (s3, s1, p1) => 5 bf per H
    ! Total: 2*6 + 6*3 = 30 shells, 2*15 + 6*5 = 60 bf, mxgtot = 8

    natoms = nat
    nsh = 30
    num_bas = 60
    mxgtot = 8

    call allocate_arrays()

    ! Coordinates
    do i = 1, natoms
      coords(3*(i - 1) + 1) = xyz(1, i)
      coords(3*(i - 1) + 2) = xyz(2, i)
      coords(3*(i - 1) + 3) = xyz(3, i)
    end do

    sh_idx = 0
    bf_idx = 0
    exp_idx = 0

    ! --- Carbon 1 (atom 1) ---
    call add_c_ccpvdz(1)
    ! --- Carbon 2 (atom 5) ---
    call add_c_ccpvdz(5)

    ! --- Hydrogens ---
    call add_h_ccpvdz(2)
    call add_h_ccpvdz(3)
    call add_h_ccpvdz(4)
    call add_h_ccpvdz(6)
    call add_h_ccpvdz(7)
    call add_h_ccpvdz(8)

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
  ! Add Hydrogen cc-pVDZ basis
  ! ============================================
  subroutine add_h_ccpvdz(iatom)
    integer, intent(in) :: iatom
    integer :: base

    ! Shell 1: s, 3 primitives
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
    exponents(base + 1) = 13.0100000_dp
    exponents(base + 2) = 1.9620000_dp
    exponents(base + 3) = 0.4446000_dp
    contr_coef_s(base + 1) = 0.0196850_dp*gto_norm_s(13.0100000_dp)
    contr_coef_s(base + 2) = 0.1379770_dp*gto_norm_s(1.9620000_dp)
    contr_coef_s(base + 3) = 0.4781480_dp*gto_norm_s(0.4446000_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 2: s, 1 primitive
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
    exponents(base + 1) = 0.1220000_dp
    contr_coef_s(base + 1) = 1.0000000_dp*gto_norm_s(0.1220000_dp)
    exp_idx = exp_idx + mxgtot

    ! Shell 3: p, 1 primitive
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
    exponents(base + 1) = 0.7270000_dp
    contr_coef_p(base + 1) = 1.0000000_dp*gto_norm_p(0.7270000_dp)
    bf_idx = bf_idx + 2
    exp_idx = exp_idx + mxgtot
  end subroutine add_h_ccpvdz

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

end program bench_c2h6
