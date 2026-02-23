! The total angular momentum of this class is:          10
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3322_impl
contains
  module subroutine int3322(ff_pair, dd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, dd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n22ket(:)
    real(dp), allocatable :: xint33bra(:), xint22ket(:)
    integer(kind=int64) :: nffbra, nddket
    real(dp) :: scutffbra, scutddket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxj2, maxl, maxl2
    integer(kind=int64) :: n, i1, i3, i4, i5, k3, k4, nn, nm, km, iaa, ib, nj, ni, nl, nk
    real(dp) :: cp10, c10, cp01, c01
    integer(kind=int64) :: nx, ny, nz, mx, my, mz
    integer(kind=int64) :: bra_loop, ket_loop, ijtop, kltop
    real(dp) :: t_expon_ab, t_expon_cd, t_inverse_expon_ab, t_inverse_expon_cd
    real(dp) :: t_expon_abcd_inverse, rho, expe, dum, rab, rcd
    real(dp) :: brrk, akxk, akyk, akzk, t_expon_c, t_expon_d, t_expon_a, t_expon_b
    real(dp) :: xa, ya, za, axak, ayak, azak, axai, ayai, azai, bbrrk, xb, yb, zb, bxbk
    real(dp) :: bybk, bzbk, bxbi, bybi, bzbi, xx, c1x, c2x, c3x, c4x, c1y, c2y, c3y, c4y
    real(dp) :: c1z, c2z, c3z, c4z, f00, u2, duminv, dm2inv, bp01, b00, b10, xcp00, xc00
    real(dp) :: ycp00, zcp00, zc00, yc00, dij, dxij, dyij, dzij, dxkl, dykl, dzkl
    real(dp) :: buff(9)
    real(dp) :: roots(6), wghts(6)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(40), wgrid(40), p0(40), p1(40), p2(40)
    real(dp) :: rts(6), wts(6), alpha(6), beta(6), wrk(6)
    real(dp) :: xin(864), yin(864), zin(864)
    real(dp) :: eri_value(3600)
    real(dp) :: d33bra(100), d22ket(36)
    integer(kind=int64) :: ix(10), jx(10), kx(6), lx(6)
    integer(kind=int64) :: iy(10), jy(10), ky(6), ly(6)
    integer(kind=int64) :: iz(10), jz(10), kz(6), lz(6)
    integer(kind=int64) :: in(7), in1(7), kn(5)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(36), kly(36), klz(36)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj, kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 37
    in1(3) = 73
    in1(4) = 109
    in1(5) = 118
    in1(6) = 127
    in1(7) = 136

    kn(1) = 0
    kn(2) = 3
    kn(3) = 6
    kn(4) = 7
    kn(5) = 8

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 2
    lx(2) = 0
    lx(3) = 0
    lx(4) = 1
    lx(5) = 1
    lx(6) = 0

    kx(1) = 6
    kx(2) = 0
    kx(3) = 0
    kx(4) = 3
    kx(5) = 3
    kx(6) = 0

    jx(1) = 27
    jx(2) = 0
    jx(3) = 0
    jx(4) = 18
    jx(5) = 18
    jx(6) = 9
    jx(7) = 0
    jx(8) = 9
    jx(9) = 0
    jx(10) = 9

    ix(1) = 109
    ix(2) = 1
    ix(3) = 1
    ix(4) = 73
    ix(5) = 73
    ix(6) = 37
    ix(7) = 1
    ix(8) = 37
    ix(9) = 1
    ix(10) = 37

    ! y-arrays

    ly(1) = 0
    ly(2) = 2
    ly(3) = 0
    ly(4) = 1
    ly(5) = 0
    ly(6) = 1

    ky(1) = 0
    ky(2) = 6
    ky(3) = 0
    ky(4) = 3
    ky(5) = 0
    ky(6) = 3

    jy(1) = 0
    jy(2) = 27
    jy(3) = 0
    jy(4) = 9
    jy(5) = 0
    jy(6) = 18
    jy(7) = 18
    jy(8) = 0
    jy(9) = 9
    jy(10) = 9

    iy(1) = 1
    iy(2) = 109
    iy(3) = 1
    iy(4) = 37
    iy(5) = 1
    iy(6) = 73
    iy(7) = 73
    iy(8) = 1
    iy(9) = 37
    iy(10) = 37

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 2
    lz(4) = 0
    lz(5) = 1
    lz(6) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 6
    kz(4) = 0
    kz(5) = 3
    kz(6) = 3

    jz(1) = 0
    jz(2) = 0
    jz(3) = 27
    jz(4) = 0
    jz(5) = 9
    jz(6) = 0
    jz(7) = 9
    jz(8) = 18
    jz(9) = 18
    jz(10) = 9

    iz(1) = 1
    iz(2) = 1
    iz(3) = 109
    iz(4) = 1
    iz(5) = 37
    iz(6) = 1
    iz(7) = 37
    iz(8) = 73
    iz(9) = 73
    iz(10) = 37

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 136
    ijx(2) = 109
    ijx(3) = 109
    ijx(4) = 127
    ijx(5) = 127
    ijx(6) = 118
    ijx(7) = 109
    ijx(8) = 118
    ijx(9) = 109
    ijx(10) = 118
    ijx(11) = 28
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 19
    ijx(15) = 19
    ijx(16) = 10
    ijx(17) = 1
    ijx(18) = 10
    ijx(19) = 1
    ijx(20) = 10
    ijx(21) = 28
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 19
    ijx(25) = 19
    ijx(26) = 10
    ijx(27) = 1
    ijx(28) = 10
    ijx(29) = 1
    ijx(30) = 10
    ijx(31) = 100
    ijx(32) = 73
    ijx(33) = 73
    ijx(34) = 91
    ijx(35) = 91
    ijx(36) = 82
    ijx(37) = 73
    ijx(38) = 82
    ijx(39) = 73
    ijx(40) = 82
    ijx(41) = 100
    ijx(42) = 73
    ijx(43) = 73
    ijx(44) = 91
    ijx(45) = 91
    ijx(46) = 82
    ijx(47) = 73
    ijx(48) = 82
    ijx(49) = 73
    ijx(50) = 82
    ijx(51) = 64
    ijx(52) = 37
    ijx(53) = 37
    ijx(54) = 55
    ijx(55) = 55
    ijx(56) = 46
    ijx(57) = 37
    ijx(58) = 46
    ijx(59) = 37
    ijx(60) = 46
    ijx(61) = 28
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 19
    ijx(65) = 19
    ijx(66) = 10
    ijx(67) = 1
    ijx(68) = 10
    ijx(69) = 1
    ijx(70) = 10
    ijx(71) = 64
    ijx(72) = 37
    ijx(73) = 37
    ijx(74) = 55
    ijx(75) = 55
    ijx(76) = 46
    ijx(77) = 37
    ijx(78) = 46
    ijx(79) = 37
    ijx(80) = 46
    ijx(81) = 28
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 19
    ijx(85) = 19
    ijx(86) = 10
    ijx(87) = 1
    ijx(88) = 10
    ijx(89) = 1
    ijx(90) = 10
    ijx(91) = 64
    ijx(92) = 37
    ijx(93) = 37
    ijx(94) = 55
    ijx(95) = 55
    ijx(96) = 46
    ijx(97) = 37
    ijx(98) = 46
    ijx(99) = 37
    ijx(100) = 46

    ijy(1) = 1
    ijy(2) = 28
    ijy(3) = 1
    ijy(4) = 10
    ijy(5) = 1
    ijy(6) = 19
    ijy(7) = 19
    ijy(8) = 1
    ijy(9) = 10
    ijy(10) = 10
    ijy(11) = 109
    ijy(12) = 136
    ijy(13) = 109
    ijy(14) = 118
    ijy(15) = 109
    ijy(16) = 127
    ijy(17) = 127
    ijy(18) = 109
    ijy(19) = 118
    ijy(20) = 118
    ijy(21) = 1
    ijy(22) = 28
    ijy(23) = 1
    ijy(24) = 10
    ijy(25) = 1
    ijy(26) = 19
    ijy(27) = 19
    ijy(28) = 1
    ijy(29) = 10
    ijy(30) = 10
    ijy(31) = 37
    ijy(32) = 64
    ijy(33) = 37
    ijy(34) = 46
    ijy(35) = 37
    ijy(36) = 55
    ijy(37) = 55
    ijy(38) = 37
    ijy(39) = 46
    ijy(40) = 46
    ijy(41) = 1
    ijy(42) = 28
    ijy(43) = 1
    ijy(44) = 10
    ijy(45) = 1
    ijy(46) = 19
    ijy(47) = 19
    ijy(48) = 1
    ijy(49) = 10
    ijy(50) = 10
    ijy(51) = 73
    ijy(52) = 100
    ijy(53) = 73
    ijy(54) = 82
    ijy(55) = 73
    ijy(56) = 91
    ijy(57) = 91
    ijy(58) = 73
    ijy(59) = 82
    ijy(60) = 82
    ijy(61) = 73
    ijy(62) = 100
    ijy(63) = 73
    ijy(64) = 82
    ijy(65) = 73
    ijy(66) = 91
    ijy(67) = 91
    ijy(68) = 73
    ijy(69) = 82
    ijy(70) = 82
    ijy(71) = 1
    ijy(72) = 28
    ijy(73) = 1
    ijy(74) = 10
    ijy(75) = 1
    ijy(76) = 19
    ijy(77) = 19
    ijy(78) = 1
    ijy(79) = 10
    ijy(80) = 10
    ijy(81) = 37
    ijy(82) = 64
    ijy(83) = 37
    ijy(84) = 46
    ijy(85) = 37
    ijy(86) = 55
    ijy(87) = 55
    ijy(88) = 37
    ijy(89) = 46
    ijy(90) = 46
    ijy(91) = 37
    ijy(92) = 64
    ijy(93) = 37
    ijy(94) = 46
    ijy(95) = 37
    ijy(96) = 55
    ijy(97) = 55
    ijy(98) = 37
    ijy(99) = 46
    ijy(100) = 46

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 28
    ijz(4) = 1
    ijz(5) = 10
    ijz(6) = 1
    ijz(7) = 10
    ijz(8) = 19
    ijz(9) = 19
    ijz(10) = 10
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 28
    ijz(14) = 1
    ijz(15) = 10
    ijz(16) = 1
    ijz(17) = 10
    ijz(18) = 19
    ijz(19) = 19
    ijz(20) = 10
    ijz(21) = 109
    ijz(22) = 109
    ijz(23) = 136
    ijz(24) = 109
    ijz(25) = 118
    ijz(26) = 109
    ijz(27) = 118
    ijz(28) = 127
    ijz(29) = 127
    ijz(30) = 118
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 28
    ijz(34) = 1
    ijz(35) = 10
    ijz(36) = 1
    ijz(37) = 10
    ijz(38) = 19
    ijz(39) = 19
    ijz(40) = 10
    ijz(41) = 37
    ijz(42) = 37
    ijz(43) = 64
    ijz(44) = 37
    ijz(45) = 46
    ijz(46) = 37
    ijz(47) = 46
    ijz(48) = 55
    ijz(49) = 55
    ijz(50) = 46
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 28
    ijz(54) = 1
    ijz(55) = 10
    ijz(56) = 1
    ijz(57) = 10
    ijz(58) = 19
    ijz(59) = 19
    ijz(60) = 10
    ijz(61) = 37
    ijz(62) = 37
    ijz(63) = 64
    ijz(64) = 37
    ijz(65) = 46
    ijz(66) = 37
    ijz(67) = 46
    ijz(68) = 55
    ijz(69) = 55
    ijz(70) = 46
    ijz(71) = 73
    ijz(72) = 73
    ijz(73) = 100
    ijz(74) = 73
    ijz(75) = 82
    ijz(76) = 73
    ijz(77) = 82
    ijz(78) = 91
    ijz(79) = 91
    ijz(80) = 82
    ijz(81) = 73
    ijz(82) = 73
    ijz(83) = 100
    ijz(84) = 73
    ijz(85) = 82
    ijz(86) = 73
    ijz(87) = 82
    ijz(88) = 91
    ijz(89) = 91
    ijz(90) = 82
    ijz(91) = 37
    ijz(92) = 37
    ijz(93) = 64
    ijz(94) = 37
    ijz(95) = 46
    ijz(96) = 37
    ijz(97) = 46
    ijz(98) = 55
    ijz(99) = 55
    ijz(100) = 46

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 8
    klx(2) = 6
    klx(3) = 6
    klx(4) = 7
    klx(5) = 7
    klx(6) = 6
    klx(7) = 2
    klx(8) = 0
    klx(9) = 0
    klx(10) = 1
    klx(11) = 1
    klx(12) = 0
    klx(13) = 2
    klx(14) = 0
    klx(15) = 0
    klx(16) = 1
    klx(17) = 1
    klx(18) = 0
    klx(19) = 5
    klx(20) = 3
    klx(21) = 3
    klx(22) = 4
    klx(23) = 4
    klx(24) = 3
    klx(25) = 5
    klx(26) = 3
    klx(27) = 3
    klx(28) = 4
    klx(29) = 4
    klx(30) = 3
    klx(31) = 2
    klx(32) = 0
    klx(33) = 0
    klx(34) = 1
    klx(35) = 1
    klx(36) = 0

    kly(1) = 0
    kly(2) = 2
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 1
    kly(7) = 6
    kly(8) = 8
    kly(9) = 6
    kly(10) = 7
    kly(11) = 6
    kly(12) = 7
    kly(13) = 0
    kly(14) = 2
    kly(15) = 0
    kly(16) = 1
    kly(17) = 0
    kly(18) = 1
    kly(19) = 3
    kly(20) = 5
    kly(21) = 3
    kly(22) = 4
    kly(23) = 3
    kly(24) = 4
    kly(25) = 0
    kly(26) = 2
    kly(27) = 0
    kly(28) = 1
    kly(29) = 0
    kly(30) = 1
    kly(31) = 3
    kly(32) = 5
    kly(33) = 3
    kly(34) = 4
    kly(35) = 3
    kly(36) = 4

    klz(1) = 0
    klz(2) = 0
    klz(3) = 2
    klz(4) = 0
    klz(5) = 1
    klz(6) = 1
    klz(7) = 0
    klz(8) = 0
    klz(9) = 2
    klz(10) = 0
    klz(11) = 1
    klz(12) = 1
    klz(13) = 6
    klz(14) = 6
    klz(15) = 8
    klz(16) = 6
    klz(17) = 7
    klz(18) = 7
    klz(19) = 0
    klz(20) = 0
    klz(21) = 2
    klz(22) = 0
    klz(23) = 1
    klz(24) = 1
    klz(25) = 3
    klz(26) = 3
    klz(27) = 5
    klz(28) = 3
    klz(29) = 4
    klz(30) = 4
    klz(31) = 3
    klz(32) = 3
    klz(33) = 5
    klz(34) = 3
    klz(35) = 4
    klz(36) = 4

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n22ket(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22ket(res%n_d_shl*(res%n_d_shl + 1)/2))

    ! Start screening

    scutffbra = cutoff_schwarz/maxval(ff_pair%xints)
    nffbra = 0
    do ij = 1, res%n_f_shl*(res%n_f_shl + 1)/2
      if (ff_pair%xints(ij) .ge. scutffbra) then
        nffbra = nffbra + 1
        xint33bra(nffbra) = ff_pair%xints(ij)
        n33bra(nffbra) = ij
      end if
    end do

    scutddket = cutoff_schwarz/maxval(dd_pair%xints)
    nddket = 0
    do ij = 1, res%n_d_shl*(res%n_d_shl + 1)/2
      if (dd_pair%xints(ij) .ge. scutddket) then
        nddket = nddket + 1
        xint22ket(nddket) = dd_pair%xints(ij)
        n22ket(nddket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*nddket) .le. nchunksize_int64) nchunksize_int64 = nffbra*nddket
    ntile = int(nffbra*nddket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*nddket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint22ket, n22ket, dd_pair, ff_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d22ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d33bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,maxl,maxl2,iandj,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nffbra) + 1
              kl_tmp = (iquart - 1)/nffbra + 1

              test = xint33bra(ij_tmp)*xint22ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n22ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_d_shl(ksh_tmp)
                lsh = res%i_d_shl(lsh_tmp)

                ijtop = res%contr_num(ish)*res%contr_num(jsh)

                kltop = res%contr_num(ksh)*res%contr_num(lsh)

                ! --- Distances between shells ----

                rab = ((res%coord_sh(ish, 1) - res%coord_sh(jsh, 1))*(res%coord_sh(ish, 1) - res%coord_sh(jsh, 1)) + &
                       (res%coord_sh(ish, 2) - res%coord_sh(jsh, 2))*(res%coord_sh(ish, 2) - res%coord_sh(jsh, 2)) + &
                       (res%coord_sh(ish, 3) - res%coord_sh(jsh, 3))*(res%coord_sh(ish, 3) - res%coord_sh(jsh, 3)))

                rcd = ((res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1))*(res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1)) + &
                       (res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2))*(res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2)) + &
                       (res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3))*(res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3)))

                ! Distances used in the bra transfer equation

                dxij = res%coord_sh(ish, 1) - res%coord_sh(jsh, 1)
                dyij = res%coord_sh(ish, 2) - res%coord_sh(jsh, 2)
                dzij = res%coord_sh(ish, 3) - res%coord_sh(jsh, 3)

                ! Distances used in the ket transfer equation

                dxkl = res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1)
                dykl = res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2)
                dzkl = res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3)

                ! Prepare to begin looping

                eri_value = 0.0_dp

                ket_loop = 0

                ! --- Start looping over primitives ---

                do k = 1, kltop

                  ket_loop = ket_loop + 1

                  t_expon_cd = dd_pair%t_expon_ab(dd_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = dd_pair%expon_a(dd_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = dd_pair%expon_b(dd_pair%pair_loc(kl) + ket_loop)
                  t_inverse_expon_cd = 1.0_dp/t_expon_cd

                  brrk = t_expon_c*rcd                ! Same variable names
                  akxk = t_expon_c*res%coord_sh(ksh, 1)    ! as in the original
                  akyk = t_expon_c*res%coord_sh(ksh, 2)    ! GAMESS subroutine
                  akzk = t_expon_c*res%coord_sh(ksh, 3)

                  bbrrk = t_expon_d*brrk*t_inverse_expon_cd
                  xb = (akxk + t_expon_d*res%coord_sh(lsh, 1))*t_inverse_expon_cd
                  yb = (akyk + t_expon_d*res%coord_sh(lsh, 2))*t_inverse_expon_cd
                  zb = (akzk + t_expon_d*res%coord_sh(lsh, 3))*t_inverse_expon_cd

                  bxbk = t_expon_cd*(xb - res%coord_sh(ksh, 1))
                  bybk = t_expon_cd*(yb - res%coord_sh(ksh, 2))
                  bzbk = t_expon_cd*(zb - res%coord_sh(ksh, 3))

                  bxbi = t_expon_cd*(xb - res%coord_sh(ish, 1))
                  bybi = t_expon_cd*(yb - res%coord_sh(ish, 2))
                  bzbi = t_expon_cd*(zb - res%coord_sh(ish, 3))

                  ! Contraction coefficient factor + angular momentum normalization
                  ! This could be optimized by moving twopi_5_2 to the bra
                  ! Whether this would have much impact or be benfictial at all
                  ! is debatable

                  d22ket(1) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(2) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(3) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(4) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(5) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(6) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(7) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(8) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(9) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(10) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(11) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(12) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(13) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(14) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(15) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(16) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(17) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(18) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(19) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(20) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(21) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(22) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(23) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(24) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(25) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(26) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(27) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(28) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(29) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(30) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(31) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(32) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(33) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(34) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(35) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(36) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3

                  bra_loop = 0

                  do i = 1, ijtop

                    bra_loop = bra_loop + 1

                    t_expon_ab = ff_pair%t_expon_ab(ff_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = ff_pair%expon_a(ff_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = ff_pair%expon_b(ff_pair%pair_loc(ij) + bra_loop)
                    t_inverse_expon_ab = 1.0_dp/t_expon_ab

                    dum = bbrrk + t_expon_a*t_expon_b*rab*t_inverse_expon_ab
                    t_expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd)
                    expe = exp(-dum)*sqrt(t_expon_abcd_inverse)
                    rho = t_expon_ab*t_expon_cd*t_expon_abcd_inverse

                    xa = (t_expon_a*res%coord_sh(ish, 1) + t_expon_b*res%coord_sh(jsh, 1))
                    ya = (t_expon_a*res%coord_sh(ish, 2) + t_expon_b*res%coord_sh(jsh, 2))
                    za = (t_expon_a*res%coord_sh(ish, 3) + t_expon_b*res%coord_sh(jsh, 3))

                    xa = xa*t_inverse_expon_ab
                    ya = ya*t_inverse_expon_ab
                    za = za*t_inverse_expon_ab

                    axak = t_expon_ab*(xa - res%coord_sh(ksh, 1))
                    ayak = t_expon_ab*(ya - res%coord_sh(ksh, 2))
                    azak = t_expon_ab*(za - res%coord_sh(ksh, 3))

                    axai = t_expon_ab*(xa - res%coord_sh(ish, 1))
                    ayai = t_expon_ab*(ya - res%coord_sh(ish, 2))
                    azai = t_expon_ab*(za - res%coord_sh(ish, 3))

                    c1x = bxbk + axak
                    c2x = t_expon_ab*bxbk
                    c3x = bxbi + axai
                    c4x = t_expon_cd*axai

                    c1y = bybk + ayak
                    c2y = t_expon_ab*bybk
                    c3y = bybi + ayai
                    c4y = t_expon_cd*ayai

                    c1z = bzbk + azak
                    c2z = t_expon_ab*bzbk
                    c3z = bzbi + azai
                    c4z = t_expon_cd*azai

                    ! Contraction coefficient factor + angular momentum normalization
                    ! Multiply only when the array element is not 1.0_dp to begin with

                    d33bra(1) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(2) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(3) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(4) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(5) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(6) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(7) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(8) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(9) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(10) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(11) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(12) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(13) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(14) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(15) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(16) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(17) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(18) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(19) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(20) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(21) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(22) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(23) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)
                    d33bra(24) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(25) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(26) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(27) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(28) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(29) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(30) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(31) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(32) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(33) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(34) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(35) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(36) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(37) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(38) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(39) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(40) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(41) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(42) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(43) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(44) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(45) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(46) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(47) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(48) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(49) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(50) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(51) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(52) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(53) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(54) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(55) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(56) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(57) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(58) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(59) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(60) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(61) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(62) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(63) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(64) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(65) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(66) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(67) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(68) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(69) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(70) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(71) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(72) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(73) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(74) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(75) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(76) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(77) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(78) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(79) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(80) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(81) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(82) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(83) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d33bra(84) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(85) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(86) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(87) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(88) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(89) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5
                    d33bra(90) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt5*sqrt3
                    d33bra(91) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(92) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(93) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d33bra(94) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(95) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(96) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(97) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(98) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(99) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5
                    d33bra(100) = ff_pair%d_coeff_alt(ff_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt5*sqrt3

                    ! Compute xx and evaluate for the respective number of roots and weights

                    xx = rho*((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb) + (za - zb)*(za - zb))

                    ! Evaluate roots and weights

                    if (xx .ge. 60.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.9874701406848084D-01
                      rts(2) = factr*0.8983028345696176D+00
                      rts(3) = factr*0.2552589802668170D+01
                      rts(4) = factr*0.5196152530054469D+01
                      rts(5) = factr*0.9124248037531169D+01
                      rts(6) = factr*0.1512995978110809D+02

                      wts(1) = factw*0.5701352362624800D+00
                      wts(2) = factw*0.2604923102641603D+00
                      wts(3) = factw*0.5160798561588373D-01
                      wts(4) = factw*0.3905390584629056D-02
                      wts(5) = factw*0.8573687043587825D-04
                      wts(6) = factw*0.2658551684356290D-06

                    else ! "regular" evaluation

                      rgrid(1) = 0.7764167660645311D-06
                      rgrid(2) = 0.2150066216486457D-04
                      rgrid(3) = 0.1292774686851154D-03
                      rgrid(4) = 0.4427485262711727D-03
                      rgrid(5) = 0.1128529682849472D-02
                      rgrid(6) = 0.2396160899229440D-02
                      rgrid(7) = 0.4491713694776070D-02
                      rgrid(8) = 0.7690217393316887D-02
                      rgrid(9) = 0.1228709604735500D-01
                      rgrid(10) = 0.1858883348816642D-01
                      rgrid(11) = 0.2690310419233024D-01
                      rgrid(12) = 0.3752862210285295D-01
                      rgrid(13) = 0.5074496784251226D-01
                      rgrid(14) = 0.6680265670059675D-01
                      rgrid(15) = 0.8591370530679703D-01
                      rgrid(16) = 0.1082429441270552D+00
                      rgrid(17) = 0.1339003060774142D+00
                      rgrid(18) = 0.1629342990513548D+00
                      rgrid(19) = 0.1953268425285067D+00
                      rgrid(20) = 0.2309896163367905D+00
                      rgrid(21) = 0.2697620338428412D+00
                      rgrid(22) = 0.3114109132037622D+00
                      rgrid(23) = 0.3556318797527262D+00
                      rgrid(24) = 0.4020524910846677D+00
                      rgrid(25) = 0.4502370349528138D+00
                      rgrid(26) = 0.4996929096784019D+00
                      rgrid(27) = 0.5498784583867753D+00
                      rgrid(28) = 0.6002120929376404D+00
                      rgrid(29) = 0.6500825117708330D+00
                      rgrid(30) = 0.6988597888065096D+00
                      rgrid(31) = 0.7459070886780932D+00
                      rgrid(32) = 0.7905927474738742D+00
                      rgrid(33) = 0.8323024482266286D+00
                      rgrid(34) = 0.8704512169070354D+00
                      rgrid(35) = 0.9044949678681030D+00
                      rgrid(36) = 0.9339413379615261D+00
                      rgrid(37) = 0.9583595677400624D+00
                      rgrid(38) = 0.9773892274524588D+00
                      rgrid(39) = 0.9907477393616217D+00
                      rgrid(40) = 0.9982384861273251D+00

                      wgrid(1) = 0.2260638549266816D-02*exp(-xx*0.7764167660645311D-06)
                      wgrid(2) = 0.5249142265576177D-02*exp(-xx*0.2150066216486457D-04)
                      wgrid(3) = 0.8210529190954175D-02*exp(-xx*0.1292774686851154D-03)
                      wgrid(4) = 0.1112292459708341D-01*exp(-xx*0.4427485262711727D-03)
                      wgrid(5) = 0.1396850349001161D-01*exp(-xx*0.1128529682849472D-02)
                      wgrid(6) = 0.1673009764127369D-01*exp(-xx*0.2396160899229440D-02)
                      wgrid(7) = 0.1939108398723607D-01*exp(-xx*0.4491713694776070D-02)
                      wgrid(8) = 0.2193545409283397D-01*exp(-xx*0.7690217393316887D-02)
                      wgrid(9) = 0.2434790381753858D-01*exp(-xx*0.1228709604735500D-01)
                      wgrid(10) = 0.2661392349196827D-01*exp(-xx*0.1858883348816642D-01)
                      wgrid(11) = 0.2871988454969593D-01*exp(-xx*0.2690310419233024D-01)
                      wgrid(12) = 0.3065312124646436D-01*exp(-xx*0.3752862210285295D-01)
                      wgrid(13) = 0.3240200672830060D-01*exp(-xx*0.5074496784251226D-01)
                      wgrid(14) = 0.3395602290761694D-01*exp(-xx*0.6680265670059675D-01)
                      wgrid(15) = 0.3530582369564367D-01*exp(-xx*0.8591370530679703D-01)
                      wgrid(16) = 0.3644329119790184D-01*exp(-xx*0.1082429441270552D+00)
                      wgrid(17) = 0.3736158452898425D-01*exp(-xx*0.1339003060774142D+00)
                      wgrid(18) = 0.3805518095031275D-01*exp(-xx*0.1629342990513548D+00)
                      wgrid(19) = 0.3851990908212372D-01*exp(-xx*0.1953268425285067D+00)
                      wgrid(20) = 0.3875297398921261D-01*exp(-xx*0.2309896163367905D+00)
                      wgrid(21) = 0.3875297398921173D-01*exp(-xx*0.2697620338428412D+00)
                      wgrid(22) = 0.3851990908212474D-01*exp(-xx*0.3114109132037622D+00)
                      wgrid(23) = 0.3805518095031294D-01*exp(-xx*0.3556318797527262D+00)
                      wgrid(24) = 0.3736158452898516D-01*exp(-xx*0.4020524910846677D+00)
                      wgrid(25) = 0.3644329119790151D-01*exp(-xx*0.4502370349528138D+00)
                      wgrid(26) = 0.3530582369564279D-01*exp(-xx*0.4996929096784019D+00)
                      wgrid(27) = 0.3395602290761690D-01*exp(-xx*0.5498784583867753D+00)
                      wgrid(28) = 0.3240200672830062D-01*exp(-xx*0.6002120929376404D+00)
                      wgrid(29) = 0.3065312124646425D-01*exp(-xx*0.6500825117708330D+00)
                      wgrid(30) = 0.2871988454969628D-01*exp(-xx*0.6988597888065096D+00)
                      wgrid(31) = 0.2661392349196852D-01*exp(-xx*0.7459070886780932D+00)
                      wgrid(32) = 0.2434790381753621D-01*exp(-xx*0.7905927474738742D+00)
                      wgrid(33) = 0.2193545409283680D-01*exp(-xx*0.8323024482266286D+00)
                      wgrid(34) = 0.1939108398723554D-01*exp(-xx*0.8704512169070354D+00)
                      wgrid(35) = 0.1673009764127440D-01*exp(-xx*0.9044949678681030D+00)
                      wgrid(36) = 0.1396850349001189D-01*exp(-xx*0.9339413379615261D+00)
                      wgrid(37) = 0.1112292459708395D-01*exp(-xx*0.9583595677400624D+00)
                      wgrid(38) = 0.8210529190953284D-02*exp(-xx*0.9773892274524588D+00)
                      wgrid(39) = 0.5249142265576468D-02*exp(-xx*0.9907477393616217D+00)
                      wgrid(40) = 0.2260638549266104D-02*exp(-xx*0.9982384861273251D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 40
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 40
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 5

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 40

                          if (wgrid(m) .eq. 0.0d+00) goto 30

                          p0(m) = p1(m)
                          p1(m) = p2(m)
                          p2(m) = (rgrid(m) - alpha(kk))*p1(m) - beta(kk)*p0(m)
                          t = wgrid(m)*p2(m)*p2(m)
                          sum1 = sum1 + t
                          sum2 = sum2 + t*rgrid(m)

30                        continue

                          alpha(kk + 1) = sum2/sum1
                          beta(kk + 1) = sum1/sum0
                          sum0 = sum1

                        end do

                        ! End of RYSDS

                        ! Call to RYSGW

                        rts(1) = alpha(1)
                        wts(1) = 1.0D+00
                        wrk(6) = 0.0D+00
                        do 100 kk = 2, 6

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 6

                            jj = 0

105                         do 110 m = l, 6
                              if (m .eq. 6) go to 120
                              if (abs(wrk(m)) .le. (1.0D-14)*(abs(rts(m)) + abs(rts(m + 1)))) go to 120
110                           continue

120                           dpp = rts(l)
                              if (m .eq. l) go to 240
                              jj = jj + 1

                              dg = (rts(l + 1) - dpp)/(2.0D+00*wrk(l))
                              dr = sqrt(dg*dg + 1.0D+00)
                              dg = rts(m) - dpp + wrk(l)/(dg + sign(dr, dg))
                              ds = 1.0D+00
                              dc = 1.0D+00
                              dpp = 0.0D+00
                              mml = m - l

                              do 200 ii = 1, mml

                                mmii = m - ii
                                df = ds*wrk(mmii)
                                db = dc*wrk(mmii)
                                if (abs(df) .lt. abs(dg)) go to 150
                                dc = dg/df
                                dr = sqrt(dc*dc + 1.0D+00)
                                wrk(mmii + 1) = df*dr
                                ds = 1.0D+00/dr
                                dc = dc*ds
                                go to 160
150                             ds = df/dg
                                dr = sqrt(ds*ds + 1.0D+00)
                                wrk(mmii + 1) = dg*dr
                                dc = 1.0D+00/dr
                                ds = ds*dc
160                             dg = rts(mmii + 1) - dpp
                                dr = (rts(mmii) - dg)*ds + 2.0D+00*dc*db
                                dpp = ds*dr
                                rts(mmii + 1) = dg + dpp
                                dg = dc*dr - db
                                df = wts(mmii + 1)
                                wts(mmii + 1) = ds*wts(mmii) + dc*df
                                wts(mmii) = dc*wts(mmii) - ds*df

200                             continue

                                rts(l) = rts(l) - dpp
                                wrk(l) = dg
                                wrk(m) = 0.0D+00
                                go to 105

240                             continue

                                do 300 ii = 2, 6

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 6
                                    if (rts(jj) .ge. dpp) go to 260
                                    kk = jj
                                    dpp = rts(jj)
260                                 continue

                                    if (kk .eq. iim1) go to 300

                                    rts(kk) = rts(iim1)
                                    rts(iim1) = dpp
                                    dpp = wts(iim1)
                                    wts(iim1) = wts(kk)
                                    wts(kk) = dpp

300                                 continue

                                    do 310 kk = 1, 6
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 6
                                        roots(kk) = rts(kk)/(1.0_dp - rts(kk))
                                        wghts(kk) = wts(kk)
                                      end do

                                      ! Start computation of 2.0_dp-electron integrals for each root

                                      ! mm = 0
                                      ! do m = 1, nroots

                                      u2 = roots(1)*rho
                                      f00 = expe*wghts(1)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =    1

                                      xin(1) = 1.0_dp
                                      yin(1) = 1.0_dp
                                      zin(1) = f00

                                      ! i2 = in(2) =   37
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(37) = xc00
                                      yin(37) = yc00
                                      zin(37) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    4

                                      xin(4) = xcp00
                                      yin(4) = ycp00
                                      zin(4) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   40
                                      ! i2 =   37

                                      xin(40) = xcp00*xin(37) + cp10
                                      yin(40) = ycp00*yin(37) + cp10
                                      zin(40) = zcp00*zin(37) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   37

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   73
                                      ! i3 =    1
                                      ! i4 =   37

                                      xin(73) = c10*xin(1) + xc00*xin(37)
                                      yin(73) = c10*yin(1) + yc00*yin(37)
                                      zin(73) = c10*zin(1) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   76
                                      ! i5 =   73
                                      ! i4 =   37

                                      xin(76) = xcp00*xin(73) + cp10*xin(37)
                                      yin(76) = ycp00*yin(73) + cp10*yin(37)
                                      zin(76) = zcp00*zin(73) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   73

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  109
                                      ! i3 =   37
                                      ! i4 =   73

                                      xin(109) = c10*xin(37) + xc00*xin(73)
                                      yin(109) = c10*yin(37) + yc00*yin(73)
                                      zin(109) = c10*zin(37) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  112
                                      ! i5 =  109
                                      ! i4 =   73

                                      xin(112) = xcp00*xin(109) + cp10*xin(73)
                                      yin(112) = ycp00*yin(109) + cp10*yin(73)
                                      zin(112) = zcp00*zin(109) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =  109

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  118
                                      ! i3 =   73
                                      ! i4 =  109

                                      xin(118) = c10*xin(73) + xc00*xin(109)
                                      yin(118) = c10*yin(73) + yc00*yin(109)
                                      zin(118) = c10*zin(73) + zc00*zin(109)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  121
                                      ! i5 =  118
                                      ! i4 =  109

                                      xin(121) = xcp00*xin(118) + cp10*xin(109)
                                      yin(121) = ycp00*yin(118) + cp10*yin(109)
                                      zin(121) = zcp00*zin(118) + cp10*zin(109)

                                      ! ------------------

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  118

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  127
                                      ! i3 =  109
                                      ! i4 =  118

                                      xin(127) = c10*xin(109) + xc00*xin(118)
                                      yin(127) = c10*yin(109) + yc00*yin(118)
                                      zin(127) = c10*zin(109) + zc00*zin(118)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  130
                                      ! i5 =  127
                                      ! i4 =  118

                                      xin(130) = xcp00*xin(127) + cp10*xin(118)
                                      yin(130) = ycp00*yin(127) + cp10*yin(118)
                                      zin(130) = zcp00*zin(127) + cp10*zin(118)

                                      ! ------------------

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  136
                                      ! i3 =  118
                                      ! i4 =  127

                                      xin(136) = c10*xin(118) + xc00*xin(127)
                                      yin(136) = c10*yin(118) + yc00*yin(127)
                                      zin(136) = c10*zin(118) + zc00*zin(127)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  139
                                      ! i5 =  136
                                      ! i4 =  127

                                      xin(139) = xcp00*xin(136) + cp10*xin(127)
                                      yin(139) = ycp00*yin(136) + cp10*yin(127)
                                      zin(139) = zcp00*zin(136) + cp10*zin(127)

                                      ! ------------------

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    4

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    7
                                      ! i3 =    1
                                      ! i4 =    4

                                      xin(7) = cp01*xin(1) + xcp00*xin(4)
                                      yin(7) = cp01*yin(1) + ycp00*yin(4)
                                      zin(7) = cp01*zin(1) + zcp00*zin(4)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   43

                                      xin(43) = xc00*xin(7) + c01*xin(4)
                                      yin(43) = yc00*yin(7) + c01*yin(4)
                                      zin(43) = zc00*zin(7) + c01*zin(4)

                                      ! ------------------

                                      ! i3 = i4 =    4
                                      ! i4 = i5 =    7

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    8
                                      ! i3 =    4
                                      ! i4 =    7

                                      xin(8) = cp01*xin(4) + xcp00*xin(7)
                                      yin(8) = cp01*yin(4) + ycp00*yin(7)
                                      zin(8) = cp01*zin(4) + zcp00*zin(7)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   44

                                      xin(44) = xc00*xin(8) + c01*xin(7)
                                      yin(44) = yc00*yin(8) + c01*yin(7)
                                      zin(44) = zc00*zin(8) + c01*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    9
                                      ! i3 =    7
                                      ! i4 =    8

                                      xin(9) = cp01*xin(7) + xcp00*xin(8)
                                      yin(9) = cp01*yin(7) + ycp00*yin(8)
                                      zin(9) = cp01*zin(7) + zcp00*zin(8)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   45

                                      xin(45) = xc00*xin(9) + c01*xin(8)
                                      yin(45) = yc00*yin(9) + c01*yin(8)
                                      zin(45) = zc00*zin(9) + c01*zin(8)

                                      ! ------------------

                                      ! i3 = i4 =    8
                                      ! i4 = i5 =    9

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   37

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   73

                                      xin(79) = c10*xin(7) + xc00*xin(43) + c01*xin(40)
                                      yin(79) = c10*yin(7) + yc00*yin(43) + c01*yin(40)
                                      zin(79) = c10*zin(7) + zc00*zin(43) + c01*zin(40)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  109

                                      xin(115) = c10*xin(43) + xc00*xin(79) + c01*xin(76)
                                      yin(115) = c10*yin(43) + yc00*yin(79) + c01*yin(76)
                                      zin(115) = c10*zin(43) + zc00*zin(79) + c01*zin(76)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =  109

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  118

                                      xin(124) = c10*xin(79) + xc00*xin(115) + c01*xin(112)
                                      yin(124) = c10*yin(79) + yc00*yin(115) + c01*yin(112)
                                      zin(124) = c10*zin(79) + zc00*zin(115) + c01*zin(112)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  118

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  127

                                      xin(133) = c10*xin(115) + xc00*xin(124) + c01*xin(121)
                                      yin(133) = c10*yin(115) + yc00*yin(124) + c01*yin(121)
                                      zin(133) = c10*zin(115) + zc00*zin(124) + c01*zin(121)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  136

                                      xin(142) = c10*xin(124) + xc00*xin(133) + c01*xin(130)
                                      yin(142) = c10*yin(124) + yc00*yin(133) + c01*yin(130)
                                      zin(142) = c10*zin(124) + zc00*zin(133) + c01*zin(130)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   37

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   73

                                      xin(80) = c10*xin(8) + xc00*xin(44) + c01*xin(43)
                                      yin(80) = c10*yin(8) + yc00*yin(44) + c01*yin(43)
                                      zin(80) = c10*zin(8) + zc00*zin(44) + c01*zin(43)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  109

                                      xin(116) = c10*xin(44) + xc00*xin(80) + c01*xin(79)
                                      yin(116) = c10*yin(44) + yc00*yin(80) + c01*yin(79)
                                      zin(116) = c10*zin(44) + zc00*zin(80) + c01*zin(79)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =  109

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  118

                                      xin(125) = c10*xin(80) + xc00*xin(116) + c01*xin(115)
                                      yin(125) = c10*yin(80) + yc00*yin(116) + c01*yin(115)
                                      zin(125) = c10*zin(80) + zc00*zin(116) + c01*zin(115)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  118

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  127

                                      xin(134) = c10*xin(116) + xc00*xin(125) + c01*xin(124)
                                      yin(134) = c10*yin(116) + yc00*yin(125) + c01*yin(124)
                                      zin(134) = c10*zin(116) + zc00*zin(125) + c01*zin(124)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  136

                                      xin(143) = c10*xin(125) + xc00*xin(134) + c01*xin(133)
                                      yin(143) = c10*yin(125) + yc00*yin(134) + c01*yin(133)
                                      zin(143) = c10*zin(125) + zc00*zin(134) + c01*zin(133)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   37

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   73

                                      xin(81) = c10*xin(9) + xc00*xin(45) + c01*xin(44)
                                      yin(81) = c10*yin(9) + yc00*yin(45) + c01*yin(44)
                                      zin(81) = c10*zin(9) + zc00*zin(45) + c01*zin(44)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  109

                                      xin(117) = c10*xin(45) + xc00*xin(81) + c01*xin(80)
                                      yin(117) = c10*yin(45) + yc00*yin(81) + c01*yin(80)
                                      zin(117) = c10*zin(45) + zc00*zin(81) + c01*zin(80)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =  109

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  118

                                      xin(126) = c10*xin(81) + xc00*xin(117) + c01*xin(116)
                                      yin(126) = c10*yin(81) + yc00*yin(117) + c01*yin(116)
                                      zin(126) = c10*zin(81) + zc00*zin(117) + c01*zin(116)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  118

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  127

                                      xin(135) = c10*xin(117) + xc00*xin(126) + c01*xin(125)
                                      yin(135) = c10*yin(117) + yc00*yin(126) + c01*yin(125)
                                      zin(135) = c10*zin(117) + zc00*zin(126) + c01*zin(125)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  136

                                      xin(144) = c10*xin(126) + xc00*xin(135) + c01*xin(134)
                                      yin(144) = c10*yin(126) + yc00*yin(135) + c01*yin(134)
                                      zin(144) = c10*zin(126) + zc00*zin(135) + c01*zin(134)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  136

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  136

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  127

                                      xin(136) = xin(136) + dxij*xin(127)
                                      yin(136) = yin(136) + dyij*yin(127)
                                      zin(136) = zin(136) + dzij*zin(127)

                                      ! i3 = i4 =  127
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  118

                                      xin(127) = xin(127) + dxij*xin(118)
                                      yin(127) = yin(127) + dyij*yin(118)
                                      zin(127) = zin(127) + dzij*zin(118)

                                      ! i3 = i4 =  118
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  109

                                      xin(118) = xin(118) + dxij*xin(109)
                                      yin(118) = yin(118) + dyij*yin(109)
                                      zin(118) = zin(118) + dzij*zin(109)

                                      ! i3 = i4 =  109
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  136

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  127

                                      xin(136) = xin(136) + dxij*xin(127)
                                      yin(136) = yin(136) + dyij*yin(127)
                                      zin(136) = zin(136) + dzij*zin(127)

                                      ! i3 = i4 =  127
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  118

                                      xin(127) = xin(127) + dxij*xin(118)
                                      yin(127) = yin(127) + dyij*yin(118)
                                      zin(127) = zin(127) + dzij*zin(118)

                                      ! i3 = i4 =  118
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  136

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  127

                                      xin(136) = xin(136) + dxij*xin(127)
                                      yin(136) = yin(136) + dyij*yin(127)
                                      zin(136) = zin(136) + dzij*zin(127)

                                      ! i3 = i4 =  127
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   10

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(37) + dxij*xin(1)
                                      yin(10) = yin(37) + dyij*yin(1)
                                      zin(10) = zin(37) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    2

                                      xin(46) = xin(73) + dxij*xin(37)
                                      yin(46) = yin(73) + dyij*yin(37)
                                      zin(46) = zin(73) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   82

                                      ! ni =    3

                                      xin(82) = xin(109) + dxij*xin(73)
                                      yin(82) = yin(109) + dyij*yin(73)
                                      zin(82) = zin(109) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  118

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    3

                                      xin(19) = xin(46) + dxij*xin(10)
                                      yin(19) = yin(46) + dyij*yin(10)
                                      zin(19) = zin(46) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    2

                                      xin(55) = xin(82) + dxij*xin(46)
                                      yin(55) = yin(82) + dyij*yin(46)
                                      zin(55) = zin(82) + dzij*zin(46)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    3

                                      xin(91) = xin(118) + dxij*xin(82)
                                      yin(91) = yin(118) + dyij*yin(82)
                                      zin(91) = zin(118) + dzij*zin(82)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   28

                                      ! nj =    3

                                      ! i4 = i3 =   28

                                      ! do ni = 1,    3

                                      xin(28) = xin(55) + dxij*xin(19)
                                      yin(28) = yin(55) + dyij*yin(19)
                                      zin(28) = zin(55) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    2

                                      xin(64) = xin(91) + dxij*xin(55)
                                      yin(64) = yin(91) + dyij*yin(55)
                                      zin(64) = zin(91) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  100

                                      ! ni =    3

                                      xin(100) = xin(127) + dxij*xin(91)
                                      yin(100) = yin(127) + dyij*yin(91)
                                      zin(100) = zin(127) + dzij*zin(91)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  136

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   37

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  130

                                      xin(139) = xin(139) + dxij*xin(130)
                                      yin(139) = yin(139) + dyij*yin(130)
                                      zin(139) = zin(139) + dzij*zin(130)

                                      ! i3 = i4 =  130
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  121

                                      xin(130) = xin(130) + dxij*xin(121)
                                      yin(130) = yin(130) + dyij*yin(121)
                                      zin(130) = zin(130) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  112

                                      xin(121) = xin(121) + dxij*xin(112)
                                      yin(121) = yin(121) + dyij*yin(112)
                                      zin(121) = zin(121) + dzij*zin(112)

                                      ! i3 = i4 =  112
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  130

                                      xin(139) = xin(139) + dxij*xin(130)
                                      yin(139) = yin(139) + dyij*yin(130)
                                      zin(139) = zin(139) + dzij*zin(130)

                                      ! i3 = i4 =  130
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  121

                                      xin(130) = xin(130) + dxij*xin(121)
                                      yin(130) = yin(130) + dyij*yin(121)
                                      zin(130) = zin(130) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  130

                                      xin(139) = xin(139) + dxij*xin(130)
                                      yin(139) = yin(139) + dyij*yin(130)
                                      zin(139) = zin(139) + dzij*zin(130)

                                      ! i3 = i4 =  130
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   13

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(40) + dxij*xin(4)
                                      yin(13) = yin(40) + dyij*yin(4)
                                      zin(13) = zin(40) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! ni =    2

                                      xin(49) = xin(76) + dxij*xin(40)
                                      yin(49) = yin(76) + dyij*yin(40)
                                      zin(49) = zin(76) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    3

                                      xin(85) = xin(112) + dxij*xin(76)
                                      yin(85) = yin(112) + dyij*yin(76)
                                      zin(85) = zin(112) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   22

                                      ! nj =    2

                                      ! i4 = i3 =   22

                                      ! do ni = 1,    3

                                      xin(22) = xin(49) + dxij*xin(13)
                                      yin(22) = yin(49) + dyij*yin(13)
                                      zin(22) = zin(49) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    2

                                      xin(58) = xin(85) + dxij*xin(49)
                                      yin(58) = yin(85) + dyij*yin(49)
                                      zin(58) = zin(85) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    3

                                      xin(94) = xin(121) + dxij*xin(85)
                                      yin(94) = yin(121) + dyij*yin(85)
                                      zin(94) = zin(121) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   31

                                      ! nj =    3

                                      ! i4 = i3 =   31

                                      ! do ni = 1,    3

                                      xin(31) = xin(58) + dxij*xin(22)
                                      yin(31) = yin(58) + dyij*yin(22)
                                      zin(31) = zin(58) + dzij*zin(22)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    2

                                      xin(67) = xin(94) + dxij*xin(58)
                                      yin(67) = yin(94) + dyij*yin(58)
                                      zin(67) = zin(94) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    3

                                      xin(103) = xin(130) + dxij*xin(94)
                                      yin(103) = yin(130) + dyij*yin(94)
                                      zin(103) = zin(130) + dzij*zin(94)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   40

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(142) = xin(142) + dxij*xin(133)
                                      yin(142) = yin(142) + dyij*yin(133)
                                      zin(142) = zin(142) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  124

                                      xin(133) = xin(133) + dxij*xin(124)
                                      yin(133) = yin(133) + dyij*yin(124)
                                      zin(133) = zin(133) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  115

                                      xin(124) = xin(124) + dxij*xin(115)
                                      yin(124) = yin(124) + dyij*yin(115)
                                      zin(124) = zin(124) + dzij*zin(115)

                                      ! i3 = i4 =  115
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(142) = xin(142) + dxij*xin(133)
                                      yin(142) = yin(142) + dyij*yin(133)
                                      zin(142) = zin(142) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  124

                                      xin(133) = xin(133) + dxij*xin(124)
                                      yin(133) = yin(133) + dyij*yin(124)
                                      zin(133) = zin(133) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(142) = xin(142) + dxij*xin(133)
                                      yin(142) = yin(142) + dyij*yin(133)
                                      zin(142) = zin(142) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   16

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   16

                                      ! do ni = 1,    3

                                      xin(16) = xin(43) + dxij*xin(7)
                                      yin(16) = yin(43) + dyij*yin(7)
                                      zin(16) = zin(43) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   52

                                      ! ni =    2

                                      xin(52) = xin(79) + dxij*xin(43)
                                      yin(52) = yin(79) + dyij*yin(43)
                                      zin(52) = zin(79) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    3

                                      xin(88) = xin(115) + dxij*xin(79)
                                      yin(88) = yin(115) + dyij*yin(79)
                                      zin(88) = zin(115) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  124

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    2

                                      ! i4 = i3 =   25

                                      ! do ni = 1,    3

                                      xin(25) = xin(52) + dxij*xin(16)
                                      yin(25) = yin(52) + dyij*yin(16)
                                      zin(25) = zin(52) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    2

                                      xin(61) = xin(88) + dxij*xin(52)
                                      yin(61) = yin(88) + dyij*yin(52)
                                      zin(61) = zin(88) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! ni =    3

                                      xin(97) = xin(124) + dxij*xin(88)
                                      yin(97) = yin(124) + dyij*yin(88)
                                      zin(97) = zin(124) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   34

                                      ! nj =    3

                                      ! i4 = i3 =   34

                                      ! do ni = 1,    3

                                      xin(34) = xin(61) + dxij*xin(25)
                                      yin(34) = yin(61) + dyij*yin(25)
                                      zin(34) = zin(61) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    2

                                      xin(70) = xin(97) + dxij*xin(61)
                                      yin(70) = yin(97) + dyij*yin(61)
                                      zin(70) = zin(97) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                      ! ni =    3

                                      xin(106) = xin(133) + dxij*xin(97)
                                      yin(106) = yin(133) + dyij*yin(97)
                                      zin(106) = zin(133) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   43

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  134

                                      xin(143) = xin(143) + dxij*xin(134)
                                      yin(143) = yin(143) + dyij*yin(134)
                                      zin(143) = zin(143) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  125

                                      xin(134) = xin(134) + dxij*xin(125)
                                      yin(134) = yin(134) + dyij*yin(125)
                                      zin(134) = zin(134) + dzij*zin(125)

                                      ! i3 = i4 =  125
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  116

                                      xin(125) = xin(125) + dxij*xin(116)
                                      yin(125) = yin(125) + dyij*yin(116)
                                      zin(125) = zin(125) + dzij*zin(116)

                                      ! i3 = i4 =  116
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  134

                                      xin(143) = xin(143) + dxij*xin(134)
                                      yin(143) = yin(143) + dyij*yin(134)
                                      zin(143) = zin(143) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  125

                                      xin(134) = xin(134) + dxij*xin(125)
                                      yin(134) = yin(134) + dyij*yin(125)
                                      zin(134) = zin(134) + dzij*zin(125)

                                      ! i3 = i4 =  125
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  134

                                      xin(143) = xin(143) + dxij*xin(134)
                                      yin(143) = yin(143) + dyij*yin(134)
                                      zin(143) = zin(143) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   17

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    3

                                      xin(17) = xin(44) + dxij*xin(8)
                                      yin(17) = yin(44) + dyij*yin(8)
                                      zin(17) = zin(44) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    2

                                      xin(53) = xin(80) + dxij*xin(44)
                                      yin(53) = yin(80) + dyij*yin(44)
                                      zin(53) = zin(80) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    3

                                      xin(89) = xin(116) + dxij*xin(80)
                                      yin(89) = yin(116) + dyij*yin(80)
                                      zin(89) = zin(116) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   26

                                      ! nj =    2

                                      ! i4 = i3 =   26

                                      ! do ni = 1,    3

                                      xin(26) = xin(53) + dxij*xin(17)
                                      yin(26) = yin(53) + dyij*yin(17)
                                      zin(26) = zin(53) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   62

                                      ! ni =    2

                                      xin(62) = xin(89) + dxij*xin(53)
                                      yin(62) = yin(89) + dyij*yin(53)
                                      zin(62) = zin(89) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   98

                                      ! ni =    3

                                      xin(98) = xin(125) + dxij*xin(89)
                                      yin(98) = yin(125) + dyij*yin(89)
                                      zin(98) = zin(125) + dzij*zin(89)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  134

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   35

                                      ! nj =    3

                                      ! i4 = i3 =   35

                                      ! do ni = 1,    3

                                      xin(35) = xin(62) + dxij*xin(26)
                                      yin(35) = yin(62) + dyij*yin(26)
                                      zin(35) = zin(62) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    2

                                      xin(71) = xin(98) + dxij*xin(62)
                                      yin(71) = yin(98) + dyij*yin(62)
                                      zin(71) = zin(98) + dzij*zin(62)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    3

                                      xin(107) = xin(134) + dxij*xin(98)
                                      yin(107) = yin(134) + dyij*yin(98)
                                      zin(107) = zin(134) + dzij*zin(98)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   44

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(144) = xin(144) + dxij*xin(135)
                                      yin(144) = yin(144) + dyij*yin(135)
                                      zin(144) = zin(144) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  126

                                      xin(135) = xin(135) + dxij*xin(126)
                                      yin(135) = yin(135) + dyij*yin(126)
                                      zin(135) = zin(135) + dzij*zin(126)

                                      ! i3 = i4 =  126
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  117

                                      xin(126) = xin(126) + dxij*xin(117)
                                      yin(126) = yin(126) + dyij*yin(117)
                                      zin(126) = zin(126) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(144) = xin(144) + dxij*xin(135)
                                      yin(144) = yin(144) + dyij*yin(135)
                                      zin(144) = zin(144) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  126

                                      xin(135) = xin(135) + dxij*xin(126)
                                      yin(135) = yin(135) + dyij*yin(126)
                                      zin(135) = zin(135) + dzij*zin(126)

                                      ! i3 = i4 =  126
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(144) = xin(144) + dxij*xin(135)
                                      yin(144) = yin(144) + dyij*yin(135)
                                      zin(144) = zin(144) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   18

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(45) + dxij*xin(9)
                                      yin(18) = yin(45) + dyij*yin(9)
                                      zin(18) = zin(45) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    2

                                      xin(54) = xin(81) + dxij*xin(45)
                                      yin(54) = yin(81) + dyij*yin(45)
                                      zin(54) = zin(81) + dzij*zin(45)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    3

                                      xin(90) = xin(117) + dxij*xin(81)
                                      yin(90) = yin(117) + dyij*yin(81)
                                      zin(90) = zin(117) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    2

                                      ! i4 = i3 =   27

                                      ! do ni = 1,    3

                                      xin(27) = xin(54) + dxij*xin(18)
                                      yin(27) = yin(54) + dyij*yin(18)
                                      zin(27) = zin(54) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    2

                                      xin(63) = xin(90) + dxij*xin(54)
                                      yin(63) = yin(90) + dyij*yin(54)
                                      zin(63) = zin(90) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   99

                                      ! ni =    3

                                      xin(99) = xin(126) + dxij*xin(90)
                                      yin(99) = yin(126) + dyij*yin(90)
                                      zin(99) = zin(126) + dzij*zin(90)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   36

                                      ! nj =    3

                                      ! i4 = i3 =   36

                                      ! do ni = 1,    3

                                      xin(36) = xin(63) + dxij*xin(27)
                                      yin(36) = yin(63) + dyij*yin(27)
                                      zin(36) = zin(63) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    2

                                      xin(72) = xin(99) + dxij*xin(63)
                                      yin(72) = yin(99) + dyij*yin(63)
                                      zin(72) = zin(99) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    3

                                      xin(108) = xin(135) + dxij*xin(99)
                                      yin(108) = yin(135) + dyij*yin(99)
                                      zin(108) = zin(135) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   45

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    9

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    8

                                      xin(9) = xin(9) + dxkl*xin(8)
                                      yin(9) = yin(9) + dykl*yin(8)
                                      zin(9) = zin(9) + dzkl*zin(8)

                                      ! i3 = i4 =    8
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    9

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    8

                                      xin(9) = xin(9) + dxkl*xin(8)
                                      yin(9) = yin(9) + dykl*yin(8)
                                      zin(9) = zin(9) + dzkl*zin(8)

                                      ! i3 = i4 =    8
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    2

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    2

                                      xin(2) = xin(4) + dxkl*xin(1)
                                      yin(2) = yin(4) + dykl*yin(1)
                                      zin(2) = zin(4) + dzkl*zin(1)
                                      ! i4 = i4 + lang+1 =    5

                                      ! nk =    2

                                      xin(5) = xin(7) + dxkl*xin(4)
                                      yin(5) = yin(7) + dykl*yin(4)
                                      zin(5) = zin(7) + dzkl*zin(4)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! i4 = i3 =    3

                                      ! do nk = 1,    2

                                      xin(3) = xin(5) + dxkl*xin(2)
                                      yin(3) = yin(5) + dykl*yin(2)
                                      zin(3) = zin(5) + dzkl*zin(2)
                                      ! i4 = i4 + lang+1 =    6

                                      ! nk =    2

                                      xin(6) = xin(8) + dxkl*xin(5)
                                      yin(6) = yin(8) + dykl*yin(5)
                                      zin(6) = zin(8) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    9

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =    4

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   10

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   18

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   17

                                      xin(18) = xin(18) + dxkl*xin(17)
                                      yin(18) = yin(18) + dykl*yin(17)
                                      zin(18) = zin(18) + dzkl*zin(17)

                                      ! i3 = i4 =   17
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   16

                                      xin(17) = xin(17) + dxkl*xin(16)
                                      yin(17) = yin(17) + dykl*yin(16)
                                      zin(17) = zin(17) + dzkl*zin(16)

                                      ! i3 = i4 =   16
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   18

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   17

                                      xin(18) = xin(18) + dxkl*xin(17)
                                      yin(18) = yin(18) + dykl*yin(17)
                                      zin(18) = zin(18) + dzkl*zin(17)

                                      ! i3 = i4 =   17
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   11

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   11

                                      ! do nk = 1,    2

                                      xin(11) = xin(13) + dxkl*xin(10)
                                      yin(11) = yin(13) + dykl*yin(10)
                                      zin(11) = zin(13) + dzkl*zin(10)
                                      ! i4 = i4 + lang+1 =   14

                                      ! nk =    2

                                      xin(14) = xin(16) + dxkl*xin(13)
                                      yin(14) = yin(16) + dykl*yin(13)
                                      zin(14) = zin(16) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   17

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   12

                                      ! nl =    2

                                      ! i4 = i3 =   12

                                      ! do nk = 1,    2

                                      xin(12) = xin(14) + dxkl*xin(11)
                                      yin(12) = yin(14) + dykl*yin(11)
                                      zin(12) = zin(14) + dzkl*zin(11)
                                      ! i4 = i4 + lang+1 =   15

                                      ! nk =    2

                                      xin(15) = xin(17) + dxkl*xin(14)
                                      yin(15) = yin(17) + dykl*yin(14)
                                      zin(15) = zin(17) + dzkl*zin(14)
                                      ! i4 = i4 + lang+1 =   18

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   13

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   19

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   27

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   26

                                      xin(27) = xin(27) + dxkl*xin(26)
                                      yin(27) = yin(27) + dykl*yin(26)
                                      zin(27) = zin(27) + dzkl*zin(26)

                                      ! i3 = i4 =   26
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   25

                                      xin(26) = xin(26) + dxkl*xin(25)
                                      yin(26) = yin(26) + dykl*yin(25)
                                      zin(26) = zin(26) + dzkl*zin(25)

                                      ! i3 = i4 =   25
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   27

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   26

                                      xin(27) = xin(27) + dxkl*xin(26)
                                      yin(27) = yin(27) + dykl*yin(26)
                                      zin(27) = zin(27) + dzkl*zin(26)

                                      ! i3 = i4 =   26
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   20

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   20

                                      ! do nk = 1,    2

                                      xin(20) = xin(22) + dxkl*xin(19)
                                      yin(20) = yin(22) + dykl*yin(19)
                                      zin(20) = zin(22) + dzkl*zin(19)
                                      ! i4 = i4 + lang+1 =   23

                                      ! nk =    2

                                      xin(23) = xin(25) + dxkl*xin(22)
                                      yin(23) = yin(25) + dykl*yin(22)
                                      zin(23) = zin(25) + dzkl*zin(22)
                                      ! i4 = i4 + lang+1 =   26

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   21

                                      ! nl =    2

                                      ! i4 = i3 =   21

                                      ! do nk = 1,    2

                                      xin(21) = xin(23) + dxkl*xin(20)
                                      yin(21) = yin(23) + dykl*yin(20)
                                      zin(21) = zin(23) + dzkl*zin(20)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    2

                                      xin(24) = xin(26) + dxkl*xin(23)
                                      yin(24) = yin(26) + dykl*yin(23)
                                      zin(24) = zin(26) + dzkl*zin(23)
                                      ! i4 = i4 + lang+1 =   27

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   22

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   28

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   34

                                      xin(35) = xin(35) + dxkl*xin(34)
                                      yin(35) = yin(35) + dykl*yin(34)
                                      zin(35) = zin(35) + dzkl*zin(34)

                                      ! i3 = i4 =   34
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   29

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   29

                                      ! do nk = 1,    2

                                      xin(29) = xin(31) + dxkl*xin(28)
                                      yin(29) = yin(31) + dykl*yin(28)
                                      zin(29) = zin(31) + dzkl*zin(28)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    2

                                      xin(32) = xin(34) + dxkl*xin(31)
                                      yin(32) = yin(34) + dykl*yin(31)
                                      zin(32) = zin(34) + dzkl*zin(31)
                                      ! i4 = i4 + lang+1 =   35

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   30

                                      ! nl =    2

                                      ! i4 = i3 =   30

                                      ! do nk = 1,    2

                                      xin(30) = xin(32) + dxkl*xin(29)
                                      yin(30) = yin(32) + dykl*yin(29)
                                      zin(30) = zin(32) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   33

                                      ! nk =    2

                                      xin(33) = xin(35) + dxkl*xin(32)
                                      yin(33) = yin(35) + dykl*yin(32)
                                      zin(33) = zin(35) + dzkl*zin(32)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   31

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   37

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   45

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   44

                                      xin(45) = xin(45) + dxkl*xin(44)
                                      yin(45) = yin(45) + dykl*yin(44)
                                      zin(45) = zin(45) + dzkl*zin(44)

                                      ! i3 = i4 =   44
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   43

                                      xin(44) = xin(44) + dxkl*xin(43)
                                      yin(44) = yin(44) + dykl*yin(43)
                                      zin(44) = zin(44) + dzkl*zin(43)

                                      ! i3 = i4 =   43
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   45

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   44

                                      xin(45) = xin(45) + dxkl*xin(44)
                                      yin(45) = yin(45) + dykl*yin(44)
                                      zin(45) = zin(45) + dzkl*zin(44)

                                      ! i3 = i4 =   44
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   38

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   38

                                      ! do nk = 1,    2

                                      xin(38) = xin(40) + dxkl*xin(37)
                                      yin(38) = yin(40) + dykl*yin(37)
                                      zin(38) = zin(40) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   41

                                      ! nk =    2

                                      xin(41) = xin(43) + dxkl*xin(40)
                                      yin(41) = yin(43) + dykl*yin(40)
                                      zin(41) = zin(43) + dzkl*zin(40)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   39

                                      ! nl =    2

                                      ! i4 = i3 =   39

                                      ! do nk = 1,    2

                                      xin(39) = xin(41) + dxkl*xin(38)
                                      yin(39) = yin(41) + dykl*yin(38)
                                      zin(39) = zin(41) + dzkl*zin(38)
                                      ! i4 = i4 + lang+1 =   42

                                      ! nk =    2

                                      xin(42) = xin(44) + dxkl*xin(41)
                                      yin(42) = yin(44) + dykl*yin(41)
                                      zin(42) = zin(44) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   45

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   40

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   46

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   54

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   53

                                      xin(54) = xin(54) + dxkl*xin(53)
                                      yin(54) = yin(54) + dykl*yin(53)
                                      zin(54) = zin(54) + dzkl*zin(53)

                                      ! i3 = i4 =   53
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   52

                                      xin(53) = xin(53) + dxkl*xin(52)
                                      yin(53) = yin(53) + dykl*yin(52)
                                      zin(53) = zin(53) + dzkl*zin(52)

                                      ! i3 = i4 =   52
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   54

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   53

                                      xin(54) = xin(54) + dxkl*xin(53)
                                      yin(54) = yin(54) + dykl*yin(53)
                                      zin(54) = zin(54) + dzkl*zin(53)

                                      ! i3 = i4 =   53
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   47

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   47

                                      ! do nk = 1,    2

                                      xin(47) = xin(49) + dxkl*xin(46)
                                      yin(47) = yin(49) + dykl*yin(46)
                                      zin(47) = zin(49) + dzkl*zin(46)
                                      ! i4 = i4 + lang+1 =   50

                                      ! nk =    2

                                      xin(50) = xin(52) + dxkl*xin(49)
                                      yin(50) = yin(52) + dykl*yin(49)
                                      zin(50) = zin(52) + dzkl*zin(49)
                                      ! i4 = i4 + lang+1 =   53

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   48

                                      ! nl =    2

                                      ! i4 = i3 =   48

                                      ! do nk = 1,    2

                                      xin(48) = xin(50) + dxkl*xin(47)
                                      yin(48) = yin(50) + dykl*yin(47)
                                      zin(48) = zin(50) + dzkl*zin(47)
                                      ! i4 = i4 + lang+1 =   51

                                      ! nk =    2

                                      xin(51) = xin(53) + dxkl*xin(50)
                                      yin(51) = yin(53) + dykl*yin(50)
                                      zin(51) = zin(53) + dzkl*zin(50)
                                      ! i4 = i4 + lang+1 =   54

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   49

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   55

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   63

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   62

                                      xin(63) = xin(63) + dxkl*xin(62)
                                      yin(63) = yin(63) + dykl*yin(62)
                                      zin(63) = zin(63) + dzkl*zin(62)

                                      ! i3 = i4 =   62
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   61

                                      xin(62) = xin(62) + dxkl*xin(61)
                                      yin(62) = yin(62) + dykl*yin(61)
                                      zin(62) = zin(62) + dzkl*zin(61)

                                      ! i3 = i4 =   61
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   63

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   62

                                      xin(63) = xin(63) + dxkl*xin(62)
                                      yin(63) = yin(63) + dykl*yin(62)
                                      zin(63) = zin(63) + dzkl*zin(62)

                                      ! i3 = i4 =   62
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   56

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   56

                                      ! do nk = 1,    2

                                      xin(56) = xin(58) + dxkl*xin(55)
                                      yin(56) = yin(58) + dykl*yin(55)
                                      zin(56) = zin(58) + dzkl*zin(55)
                                      ! i4 = i4 + lang+1 =   59

                                      ! nk =    2

                                      xin(59) = xin(61) + dxkl*xin(58)
                                      yin(59) = yin(61) + dykl*yin(58)
                                      zin(59) = zin(61) + dzkl*zin(58)
                                      ! i4 = i4 + lang+1 =   62

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   57

                                      ! nl =    2

                                      ! i4 = i3 =   57

                                      ! do nk = 1,    2

                                      xin(57) = xin(59) + dxkl*xin(56)
                                      yin(57) = yin(59) + dykl*yin(56)
                                      zin(57) = zin(59) + dzkl*zin(56)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      xin(60) = xin(62) + dxkl*xin(59)
                                      yin(60) = yin(62) + dykl*yin(59)
                                      zin(60) = zin(62) + dzkl*zin(59)
                                      ! i4 = i4 + lang+1 =   63

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   58

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   64

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   70

                                      xin(71) = xin(71) + dxkl*xin(70)
                                      yin(71) = yin(71) + dykl*yin(70)
                                      zin(71) = zin(71) + dzkl*zin(70)

                                      ! i3 = i4 =   70
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   65

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   65

                                      ! do nk = 1,    2

                                      xin(65) = xin(67) + dxkl*xin(64)
                                      yin(65) = yin(67) + dykl*yin(64)
                                      zin(65) = zin(67) + dzkl*zin(64)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      xin(68) = xin(70) + dxkl*xin(67)
                                      yin(68) = yin(70) + dykl*yin(67)
                                      zin(68) = zin(70) + dzkl*zin(67)
                                      ! i4 = i4 + lang+1 =   71

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   66

                                      ! nl =    2

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    2

                                      xin(66) = xin(68) + dxkl*xin(65)
                                      yin(66) = yin(68) + dykl*yin(65)
                                      zin(66) = zin(68) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   69

                                      ! nk =    2

                                      xin(69) = xin(71) + dxkl*xin(68)
                                      yin(69) = yin(71) + dykl*yin(68)
                                      zin(69) = zin(71) + dzkl*zin(68)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   81

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   80

                                      xin(81) = xin(81) + dxkl*xin(80)
                                      yin(81) = yin(81) + dykl*yin(80)
                                      zin(81) = zin(81) + dzkl*zin(80)

                                      ! i3 = i4 =   80
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   81

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   80

                                      xin(81) = xin(81) + dxkl*xin(80)
                                      yin(81) = yin(81) + dykl*yin(80)
                                      zin(81) = zin(81) + dzkl*zin(80)

                                      ! i3 = i4 =   80
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    2

                                      xin(74) = xin(76) + dxkl*xin(73)
                                      yin(74) = yin(76) + dykl*yin(73)
                                      zin(74) = zin(76) + dzkl*zin(73)
                                      ! i4 = i4 + lang+1 =   77

                                      ! nk =    2

                                      xin(77) = xin(79) + dxkl*xin(76)
                                      yin(77) = yin(79) + dykl*yin(76)
                                      zin(77) = zin(79) + dzkl*zin(76)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! i4 = i3 =   75

                                      ! do nk = 1,    2

                                      xin(75) = xin(77) + dxkl*xin(74)
                                      yin(75) = yin(77) + dykl*yin(74)
                                      zin(75) = zin(77) + dzkl*zin(74)
                                      ! i4 = i4 + lang+1 =   78

                                      ! nk =    2

                                      xin(78) = xin(80) + dxkl*xin(77)
                                      yin(78) = yin(80) + dykl*yin(77)
                                      zin(78) = zin(80) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   81

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   76

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   82

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   90

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   89

                                      xin(90) = xin(90) + dxkl*xin(89)
                                      yin(90) = yin(90) + dykl*yin(89)
                                      zin(90) = zin(90) + dzkl*zin(89)

                                      ! i3 = i4 =   89
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   88

                                      xin(89) = xin(89) + dxkl*xin(88)
                                      yin(89) = yin(89) + dykl*yin(88)
                                      zin(89) = zin(89) + dzkl*zin(88)

                                      ! i3 = i4 =   88
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   90

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   89

                                      xin(90) = xin(90) + dxkl*xin(89)
                                      yin(90) = yin(90) + dykl*yin(89)
                                      zin(90) = zin(90) + dzkl*zin(89)

                                      ! i3 = i4 =   89
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   83

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   83

                                      ! do nk = 1,    2

                                      xin(83) = xin(85) + dxkl*xin(82)
                                      yin(83) = yin(85) + dykl*yin(82)
                                      zin(83) = zin(85) + dzkl*zin(82)
                                      ! i4 = i4 + lang+1 =   86

                                      ! nk =    2

                                      xin(86) = xin(88) + dxkl*xin(85)
                                      yin(86) = yin(88) + dykl*yin(85)
                                      zin(86) = zin(88) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   89

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   84

                                      ! nl =    2

                                      ! i4 = i3 =   84

                                      ! do nk = 1,    2

                                      xin(84) = xin(86) + dxkl*xin(83)
                                      yin(84) = yin(86) + dykl*yin(83)
                                      zin(84) = zin(86) + dzkl*zin(83)
                                      ! i4 = i4 + lang+1 =   87

                                      ! nk =    2

                                      xin(87) = xin(89) + dxkl*xin(86)
                                      yin(87) = yin(89) + dykl*yin(86)
                                      zin(87) = zin(89) + dzkl*zin(86)
                                      ! i4 = i4 + lang+1 =   90

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   85

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   91

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   99

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   98

                                      xin(99) = xin(99) + dxkl*xin(98)
                                      yin(99) = yin(99) + dykl*yin(98)
                                      zin(99) = zin(99) + dzkl*zin(98)

                                      ! i3 = i4 =   98
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   97

                                      xin(98) = xin(98) + dxkl*xin(97)
                                      yin(98) = yin(98) + dykl*yin(97)
                                      zin(98) = zin(98) + dzkl*zin(97)

                                      ! i3 = i4 =   97
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   99

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   98

                                      xin(99) = xin(99) + dxkl*xin(98)
                                      yin(99) = yin(99) + dykl*yin(98)
                                      zin(99) = zin(99) + dzkl*zin(98)

                                      ! i3 = i4 =   98
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   92

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   92

                                      ! do nk = 1,    2

                                      xin(92) = xin(94) + dxkl*xin(91)
                                      yin(92) = yin(94) + dykl*yin(91)
                                      zin(92) = zin(94) + dzkl*zin(91)
                                      ! i4 = i4 + lang+1 =   95

                                      ! nk =    2

                                      xin(95) = xin(97) + dxkl*xin(94)
                                      yin(95) = yin(97) + dykl*yin(94)
                                      zin(95) = zin(97) + dzkl*zin(94)
                                      ! i4 = i4 + lang+1 =   98

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   93

                                      ! nl =    2

                                      ! i4 = i3 =   93

                                      ! do nk = 1,    2

                                      xin(93) = xin(95) + dxkl*xin(92)
                                      yin(93) = yin(95) + dykl*yin(92)
                                      zin(93) = zin(95) + dzkl*zin(92)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    2

                                      xin(96) = xin(98) + dxkl*xin(95)
                                      yin(96) = yin(98) + dykl*yin(95)
                                      zin(96) = zin(98) + dzkl*zin(95)
                                      ! i4 = i4 + lang+1 =   99

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   94

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  100

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  106

                                      xin(107) = xin(107) + dxkl*xin(106)
                                      yin(107) = yin(107) + dykl*yin(106)
                                      zin(107) = zin(107) + dzkl*zin(106)

                                      ! i3 = i4 =  106
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  101

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  101

                                      ! do nk = 1,    2

                                      xin(101) = xin(103) + dxkl*xin(100)
                                      yin(101) = yin(103) + dykl*yin(100)
                                      zin(101) = zin(103) + dzkl*zin(100)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    2

                                      xin(104) = xin(106) + dxkl*xin(103)
                                      yin(104) = yin(106) + dykl*yin(103)
                                      zin(104) = zin(106) + dzkl*zin(103)
                                      ! i4 = i4 + lang+1 =  107

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  102

                                      ! nl =    2

                                      ! i4 = i3 =  102

                                      ! do nk = 1,    2

                                      xin(102) = xin(104) + dxkl*xin(101)
                                      yin(102) = yin(104) + dykl*yin(101)
                                      zin(102) = zin(104) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  105

                                      ! nk =    2

                                      xin(105) = xin(107) + dxkl*xin(104)
                                      yin(105) = yin(107) + dykl*yin(104)
                                      zin(105) = zin(107) + dzkl*zin(104)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  103

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  109

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  117

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  116

                                      xin(117) = xin(117) + dxkl*xin(116)
                                      yin(117) = yin(117) + dykl*yin(116)
                                      zin(117) = zin(117) + dzkl*zin(116)

                                      ! i3 = i4 =  116
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  115

                                      xin(116) = xin(116) + dxkl*xin(115)
                                      yin(116) = yin(116) + dykl*yin(115)
                                      zin(116) = zin(116) + dzkl*zin(115)

                                      ! i3 = i4 =  115
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  117

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  116

                                      xin(117) = xin(117) + dxkl*xin(116)
                                      yin(117) = yin(117) + dykl*yin(116)
                                      zin(117) = zin(117) + dzkl*zin(116)

                                      ! i3 = i4 =  116
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  110

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  110

                                      ! do nk = 1,    2

                                      xin(110) = xin(112) + dxkl*xin(109)
                                      yin(110) = yin(112) + dykl*yin(109)
                                      zin(110) = zin(112) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  113

                                      ! nk =    2

                                      xin(113) = xin(115) + dxkl*xin(112)
                                      yin(113) = yin(115) + dykl*yin(112)
                                      zin(113) = zin(115) + dzkl*zin(112)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  111

                                      ! nl =    2

                                      ! i4 = i3 =  111

                                      ! do nk = 1,    2

                                      xin(111) = xin(113) + dxkl*xin(110)
                                      yin(111) = yin(113) + dykl*yin(110)
                                      zin(111) = zin(113) + dzkl*zin(110)
                                      ! i4 = i4 + lang+1 =  114

                                      ! nk =    2

                                      xin(114) = xin(116) + dxkl*xin(113)
                                      yin(114) = yin(116) + dykl*yin(113)
                                      zin(114) = zin(116) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  117

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  112

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  118

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  126

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  125

                                      xin(126) = xin(126) + dxkl*xin(125)
                                      yin(126) = yin(126) + dykl*yin(125)
                                      zin(126) = zin(126) + dzkl*zin(125)

                                      ! i3 = i4 =  125
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  124

                                      xin(125) = xin(125) + dxkl*xin(124)
                                      yin(125) = yin(125) + dykl*yin(124)
                                      zin(125) = zin(125) + dzkl*zin(124)

                                      ! i3 = i4 =  124
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  126

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  125

                                      xin(126) = xin(126) + dxkl*xin(125)
                                      yin(126) = yin(126) + dykl*yin(125)
                                      zin(126) = zin(126) + dzkl*zin(125)

                                      ! i3 = i4 =  125
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  119

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  119

                                      ! do nk = 1,    2

                                      xin(119) = xin(121) + dxkl*xin(118)
                                      yin(119) = yin(121) + dykl*yin(118)
                                      zin(119) = zin(121) + dzkl*zin(118)
                                      ! i4 = i4 + lang+1 =  122

                                      ! nk =    2

                                      xin(122) = xin(124) + dxkl*xin(121)
                                      yin(122) = yin(124) + dykl*yin(121)
                                      zin(122) = zin(124) + dzkl*zin(121)
                                      ! i4 = i4 + lang+1 =  125

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  120

                                      ! nl =    2

                                      ! i4 = i3 =  120

                                      ! do nk = 1,    2

                                      xin(120) = xin(122) + dxkl*xin(119)
                                      yin(120) = yin(122) + dykl*yin(119)
                                      zin(120) = zin(122) + dzkl*zin(119)
                                      ! i4 = i4 + lang+1 =  123

                                      ! nk =    2

                                      xin(123) = xin(125) + dxkl*xin(122)
                                      yin(123) = yin(125) + dykl*yin(122)
                                      zin(123) = zin(125) + dzkl*zin(122)
                                      ! i4 = i4 + lang+1 =  126

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  121

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  127

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  135

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  134

                                      xin(135) = xin(135) + dxkl*xin(134)
                                      yin(135) = yin(135) + dykl*yin(134)
                                      zin(135) = zin(135) + dzkl*zin(134)

                                      ! i3 = i4 =  134
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  133

                                      xin(134) = xin(134) + dxkl*xin(133)
                                      yin(134) = yin(134) + dykl*yin(133)
                                      zin(134) = zin(134) + dzkl*zin(133)

                                      ! i3 = i4 =  133
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  135

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  134

                                      xin(135) = xin(135) + dxkl*xin(134)
                                      yin(135) = yin(135) + dykl*yin(134)
                                      zin(135) = zin(135) + dzkl*zin(134)

                                      ! i3 = i4 =  134
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  128

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  128

                                      ! do nk = 1,    2

                                      xin(128) = xin(130) + dxkl*xin(127)
                                      yin(128) = yin(130) + dykl*yin(127)
                                      zin(128) = zin(130) + dzkl*zin(127)
                                      ! i4 = i4 + lang+1 =  131

                                      ! nk =    2

                                      xin(131) = xin(133) + dxkl*xin(130)
                                      yin(131) = yin(133) + dykl*yin(130)
                                      zin(131) = zin(133) + dzkl*zin(130)
                                      ! i4 = i4 + lang+1 =  134

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  129

                                      ! nl =    2

                                      ! i4 = i3 =  129

                                      ! do nk = 1,    2

                                      xin(129) = xin(131) + dxkl*xin(128)
                                      yin(129) = yin(131) + dykl*yin(128)
                                      zin(129) = zin(131) + dzkl*zin(128)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    2

                                      xin(132) = xin(134) + dxkl*xin(131)
                                      yin(132) = yin(134) + dykl*yin(131)
                                      zin(132) = zin(134) + dzkl*zin(131)
                                      ! i4 = i4 + lang+1 =  135

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  130

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  136

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  142

                                      xin(143) = xin(143) + dxkl*xin(142)
                                      yin(143) = yin(143) + dykl*yin(142)
                                      zin(143) = zin(143) + dzkl*zin(142)

                                      ! i3 = i4 =  142
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  137

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  137

                                      ! do nk = 1,    2

                                      xin(137) = xin(139) + dxkl*xin(136)
                                      yin(137) = yin(139) + dykl*yin(136)
                                      zin(137) = zin(139) + dzkl*zin(136)
                                      ! i4 = i4 + lang+1 =  140

                                      ! nk =    2

                                      xin(140) = xin(142) + dxkl*xin(139)
                                      yin(140) = yin(142) + dykl*yin(139)
                                      zin(140) = zin(142) + dzkl*zin(139)
                                      ! i4 = i4 + lang+1 =  143

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  138

                                      ! nl =    2

                                      ! i4 = i3 =  138

                                      ! do nk = 1,    2

                                      xin(138) = xin(140) + dxkl*xin(137)
                                      yin(138) = yin(140) + dykl*yin(137)
                                      zin(138) = zin(140) + dzkl*zin(137)
                                      ! i4 = i4 + lang+1 =  141

                                      ! nk =    2

                                      xin(141) = xin(143) + dxkl*xin(140)
                                      yin(141) = yin(143) + dykl*yin(140)
                                      zin(141) = zin(143) + dzkl*zin(140)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  139

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

                                      u2 = roots(2)*rho
                                      f00 = expe*wghts(2)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  145

                                      xin(145) = 1.0_dp
                                      yin(145) = 1.0_dp
                                      zin(145) = f00

                                      ! i2 = in(2) =  181
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(181) = xc00
                                      yin(181) = yc00
                                      zin(181) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  148

                                      xin(148) = xcp00
                                      yin(148) = ycp00
                                      zin(148) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  184
                                      ! i2 =  181

                                      xin(184) = xcp00*xin(181) + cp10
                                      yin(184) = ycp00*yin(181) + cp10
                                      zin(184) = zcp00*zin(181) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  181

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  217
                                      ! i3 =  145
                                      ! i4 =  181

                                      xin(217) = c10*xin(145) + xc00*xin(181)
                                      yin(217) = c10*yin(145) + yc00*yin(181)
                                      zin(217) = c10*zin(145) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  220
                                      ! i5 =  217
                                      ! i4 =  181

                                      xin(220) = xcp00*xin(217) + cp10*xin(181)
                                      yin(220) = ycp00*yin(217) + cp10*yin(181)
                                      zin(220) = zcp00*zin(217) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  217

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  253
                                      ! i3 =  181
                                      ! i4 =  217

                                      xin(253) = c10*xin(181) + xc00*xin(217)
                                      yin(253) = c10*yin(181) + yc00*yin(217)
                                      zin(253) = c10*zin(181) + zc00*zin(217)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  256
                                      ! i5 =  253
                                      ! i4 =  217

                                      xin(256) = xcp00*xin(253) + cp10*xin(217)
                                      yin(256) = ycp00*yin(253) + cp10*yin(217)
                                      zin(256) = zcp00*zin(253) + cp10*zin(217)

                                      ! ------------------

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  253

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  262
                                      ! i3 =  217
                                      ! i4 =  253

                                      xin(262) = c10*xin(217) + xc00*xin(253)
                                      yin(262) = c10*yin(217) + yc00*yin(253)
                                      zin(262) = c10*zin(217) + zc00*zin(253)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  265
                                      ! i5 =  262
                                      ! i4 =  253

                                      xin(265) = xcp00*xin(262) + cp10*xin(253)
                                      yin(265) = ycp00*yin(262) + cp10*yin(253)
                                      zin(265) = zcp00*zin(262) + cp10*zin(253)

                                      ! ------------------

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  262

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  271
                                      ! i3 =  253
                                      ! i4 =  262

                                      xin(271) = c10*xin(253) + xc00*xin(262)
                                      yin(271) = c10*yin(253) + yc00*yin(262)
                                      zin(271) = c10*zin(253) + zc00*zin(262)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  274
                                      ! i5 =  271
                                      ! i4 =  262

                                      xin(274) = xcp00*xin(271) + cp10*xin(262)
                                      yin(274) = ycp00*yin(271) + cp10*yin(262)
                                      zin(274) = zcp00*zin(271) + cp10*zin(262)

                                      ! ------------------

                                      ! i3 = i4 =  262
                                      ! i4 = i5 =  271

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  280
                                      ! i3 =  262
                                      ! i4 =  271

                                      xin(280) = c10*xin(262) + xc00*xin(271)
                                      yin(280) = c10*yin(262) + yc00*yin(271)
                                      zin(280) = c10*zin(262) + zc00*zin(271)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  283
                                      ! i5 =  280
                                      ! i4 =  271

                                      xin(283) = xcp00*xin(280) + cp10*xin(271)
                                      yin(283) = ycp00*yin(280) + cp10*yin(271)
                                      zin(283) = zcp00*zin(280) + cp10*zin(271)

                                      ! ------------------

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  280

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  148

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  151
                                      ! i3 =  145
                                      ! i4 =  148

                                      xin(151) = cp01*xin(145) + xcp00*xin(148)
                                      yin(151) = cp01*yin(145) + ycp00*yin(148)
                                      zin(151) = cp01*zin(145) + zcp00*zin(148)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  187

                                      xin(187) = xc00*xin(151) + c01*xin(148)
                                      yin(187) = yc00*yin(151) + c01*yin(148)
                                      zin(187) = zc00*zin(151) + c01*zin(148)

                                      ! ------------------

                                      ! i3 = i4 =  148
                                      ! i4 = i5 =  151

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  152
                                      ! i3 =  148
                                      ! i4 =  151

                                      xin(152) = cp01*xin(148) + xcp00*xin(151)
                                      yin(152) = cp01*yin(148) + ycp00*yin(151)
                                      zin(152) = cp01*zin(148) + zcp00*zin(151)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  188

                                      xin(188) = xc00*xin(152) + c01*xin(151)
                                      yin(188) = yc00*yin(152) + c01*yin(151)
                                      zin(188) = zc00*zin(152) + c01*zin(151)

                                      ! ------------------

                                      ! i3 = i4 =  151
                                      ! i4 = i5 =  152

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  153
                                      ! i3 =  151
                                      ! i4 =  152

                                      xin(153) = cp01*xin(151) + xcp00*xin(152)
                                      yin(153) = cp01*yin(151) + ycp00*yin(152)
                                      zin(153) = cp01*zin(151) + zcp00*zin(152)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  189

                                      xin(189) = xc00*xin(153) + c01*xin(152)
                                      yin(189) = yc00*yin(153) + c01*yin(152)
                                      zin(189) = zc00*zin(153) + c01*zin(152)

                                      ! ------------------

                                      ! i3 = i4 =  152
                                      ! i4 = i5 =  153

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  181

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  217

                                      xin(223) = c10*xin(151) + xc00*xin(187) + c01*xin(184)
                                      yin(223) = c10*yin(151) + yc00*yin(187) + c01*yin(184)
                                      zin(223) = c10*zin(151) + zc00*zin(187) + c01*zin(184)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  253

                                      xin(259) = c10*xin(187) + xc00*xin(223) + c01*xin(220)
                                      yin(259) = c10*yin(187) + yc00*yin(223) + c01*yin(220)
                                      zin(259) = c10*zin(187) + zc00*zin(223) + c01*zin(220)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  253

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  262

                                      xin(268) = c10*xin(223) + xc00*xin(259) + c01*xin(256)
                                      yin(268) = c10*yin(223) + yc00*yin(259) + c01*yin(256)
                                      zin(268) = c10*zin(223) + zc00*zin(259) + c01*zin(256)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  262

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  271

                                      xin(277) = c10*xin(259) + xc00*xin(268) + c01*xin(265)
                                      yin(277) = c10*yin(259) + yc00*yin(268) + c01*yin(265)
                                      zin(277) = c10*zin(259) + zc00*zin(268) + c01*zin(265)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  262
                                      ! i4 = i5 =  271

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  280

                                      xin(286) = c10*xin(268) + xc00*xin(277) + c01*xin(274)
                                      yin(286) = c10*yin(268) + yc00*yin(277) + c01*yin(274)
                                      zin(286) = c10*zin(268) + zc00*zin(277) + c01*zin(274)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  280

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  181

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  217

                                      xin(224) = c10*xin(152) + xc00*xin(188) + c01*xin(187)
                                      yin(224) = c10*yin(152) + yc00*yin(188) + c01*yin(187)
                                      zin(224) = c10*zin(152) + zc00*zin(188) + c01*zin(187)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  253

                                      xin(260) = c10*xin(188) + xc00*xin(224) + c01*xin(223)
                                      yin(260) = c10*yin(188) + yc00*yin(224) + c01*yin(223)
                                      zin(260) = c10*zin(188) + zc00*zin(224) + c01*zin(223)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  253

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  262

                                      xin(269) = c10*xin(224) + xc00*xin(260) + c01*xin(259)
                                      yin(269) = c10*yin(224) + yc00*yin(260) + c01*yin(259)
                                      zin(269) = c10*zin(224) + zc00*zin(260) + c01*zin(259)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  262

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  271

                                      xin(278) = c10*xin(260) + xc00*xin(269) + c01*xin(268)
                                      yin(278) = c10*yin(260) + yc00*yin(269) + c01*yin(268)
                                      zin(278) = c10*zin(260) + zc00*zin(269) + c01*zin(268)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  262
                                      ! i4 = i5 =  271

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  280

                                      xin(287) = c10*xin(269) + xc00*xin(278) + c01*xin(277)
                                      yin(287) = c10*yin(269) + yc00*yin(278) + c01*yin(277)
                                      zin(287) = c10*zin(269) + zc00*zin(278) + c01*zin(277)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  280

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  181

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  217

                                      xin(225) = c10*xin(153) + xc00*xin(189) + c01*xin(188)
                                      yin(225) = c10*yin(153) + yc00*yin(189) + c01*yin(188)
                                      zin(225) = c10*zin(153) + zc00*zin(189) + c01*zin(188)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  253

                                      xin(261) = c10*xin(189) + xc00*xin(225) + c01*xin(224)
                                      yin(261) = c10*yin(189) + yc00*yin(225) + c01*yin(224)
                                      zin(261) = c10*zin(189) + zc00*zin(225) + c01*zin(224)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  253

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  262

                                      xin(270) = c10*xin(225) + xc00*xin(261) + c01*xin(260)
                                      yin(270) = c10*yin(225) + yc00*yin(261) + c01*yin(260)
                                      zin(270) = c10*zin(225) + zc00*zin(261) + c01*zin(260)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  262

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  271

                                      xin(279) = c10*xin(261) + xc00*xin(270) + c01*xin(269)
                                      yin(279) = c10*yin(261) + yc00*yin(270) + c01*yin(269)
                                      zin(279) = c10*zin(261) + zc00*zin(270) + c01*zin(269)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  262
                                      ! i4 = i5 =  271

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  280

                                      xin(288) = c10*xin(270) + xc00*xin(279) + c01*xin(278)
                                      yin(288) = c10*yin(270) + yc00*yin(279) + c01*yin(278)
                                      zin(288) = c10*zin(270) + zc00*zin(279) + c01*zin(278)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  280

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  280

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  280

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  271

                                      xin(280) = xin(280) + dxij*xin(271)
                                      yin(280) = yin(280) + dyij*yin(271)
                                      zin(280) = zin(280) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  262

                                      xin(271) = xin(271) + dxij*xin(262)
                                      yin(271) = yin(271) + dyij*yin(262)
                                      zin(271) = zin(271) + dzij*zin(262)

                                      ! i3 = i4 =  262
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  253

                                      xin(262) = xin(262) + dxij*xin(253)
                                      yin(262) = yin(262) + dyij*yin(253)
                                      zin(262) = zin(262) + dzij*zin(253)

                                      ! i3 = i4 =  253
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  280

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  271

                                      xin(280) = xin(280) + dxij*xin(271)
                                      yin(280) = yin(280) + dyij*yin(271)
                                      zin(280) = zin(280) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  262

                                      xin(271) = xin(271) + dxij*xin(262)
                                      yin(271) = yin(271) + dyij*yin(262)
                                      zin(271) = zin(271) + dzij*zin(262)

                                      ! i3 = i4 =  262
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  280

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  271

                                      xin(280) = xin(280) + dxij*xin(271)
                                      yin(280) = yin(280) + dyij*yin(271)
                                      zin(280) = zin(280) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  154

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  154

                                      ! do ni = 1,    3

                                      xin(154) = xin(181) + dxij*xin(145)
                                      yin(154) = yin(181) + dyij*yin(145)
                                      zin(154) = zin(181) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  190

                                      ! ni =    2

                                      xin(190) = xin(217) + dxij*xin(181)
                                      yin(190) = yin(217) + dyij*yin(181)
                                      zin(190) = zin(217) + dzij*zin(181)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  226

                                      ! ni =    3

                                      xin(226) = xin(253) + dxij*xin(217)
                                      yin(226) = yin(253) + dyij*yin(217)
                                      zin(226) = zin(253) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  262

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  163

                                      ! nj =    2

                                      ! i4 = i3 =  163

                                      ! do ni = 1,    3

                                      xin(163) = xin(190) + dxij*xin(154)
                                      yin(163) = yin(190) + dyij*yin(154)
                                      zin(163) = zin(190) + dzij*zin(154)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  199

                                      ! ni =    2

                                      xin(199) = xin(226) + dxij*xin(190)
                                      yin(199) = yin(226) + dyij*yin(190)
                                      zin(199) = zin(226) + dzij*zin(190)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    3

                                      xin(235) = xin(262) + dxij*xin(226)
                                      yin(235) = yin(262) + dyij*yin(226)
                                      zin(235) = zin(262) + dzij*zin(226)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  271

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  172

                                      ! nj =    3

                                      ! i4 = i3 =  172

                                      ! do ni = 1,    3

                                      xin(172) = xin(199) + dxij*xin(163)
                                      yin(172) = yin(199) + dyij*yin(163)
                                      zin(172) = zin(199) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                      ! ni =    2

                                      xin(208) = xin(235) + dxij*xin(199)
                                      yin(208) = yin(235) + dyij*yin(199)
                                      zin(208) = zin(235) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  244

                                      ! ni =    3

                                      xin(244) = xin(271) + dxij*xin(235)
                                      yin(244) = yin(271) + dyij*yin(235)
                                      zin(244) = zin(271) + dzij*zin(235)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  280

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  181

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  274

                                      xin(283) = xin(283) + dxij*xin(274)
                                      yin(283) = yin(283) + dyij*yin(274)
                                      zin(283) = zin(283) + dzij*zin(274)

                                      ! i3 = i4 =  274
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  265

                                      xin(274) = xin(274) + dxij*xin(265)
                                      yin(274) = yin(274) + dyij*yin(265)
                                      zin(274) = zin(274) + dzij*zin(265)

                                      ! i3 = i4 =  265
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  256

                                      xin(265) = xin(265) + dxij*xin(256)
                                      yin(265) = yin(265) + dyij*yin(256)
                                      zin(265) = zin(265) + dzij*zin(256)

                                      ! i3 = i4 =  256
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  274

                                      xin(283) = xin(283) + dxij*xin(274)
                                      yin(283) = yin(283) + dyij*yin(274)
                                      zin(283) = zin(283) + dzij*zin(274)

                                      ! i3 = i4 =  274
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  265

                                      xin(274) = xin(274) + dxij*xin(265)
                                      yin(274) = yin(274) + dyij*yin(265)
                                      zin(274) = zin(274) + dzij*zin(265)

                                      ! i3 = i4 =  265
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  274

                                      xin(283) = xin(283) + dxij*xin(274)
                                      yin(283) = yin(283) + dyij*yin(274)
                                      zin(283) = zin(283) + dzij*zin(274)

                                      ! i3 = i4 =  274
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  157

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  157

                                      ! do ni = 1,    3

                                      xin(157) = xin(184) + dxij*xin(148)
                                      yin(157) = yin(184) + dyij*yin(148)
                                      zin(157) = zin(184) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! ni =    2

                                      xin(193) = xin(220) + dxij*xin(184)
                                      yin(193) = yin(220) + dyij*yin(184)
                                      zin(193) = zin(220) + dzij*zin(184)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  229

                                      ! ni =    3

                                      xin(229) = xin(256) + dxij*xin(220)
                                      yin(229) = yin(256) + dyij*yin(220)
                                      zin(229) = zin(256) + dzij*zin(220)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  265

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  166

                                      ! nj =    2

                                      ! i4 = i3 =  166

                                      ! do ni = 1,    3

                                      xin(166) = xin(193) + dxij*xin(157)
                                      yin(166) = yin(193) + dyij*yin(157)
                                      zin(166) = zin(193) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  202

                                      ! ni =    2

                                      xin(202) = xin(229) + dxij*xin(193)
                                      yin(202) = yin(229) + dyij*yin(193)
                                      zin(202) = zin(229) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  238

                                      ! ni =    3

                                      xin(238) = xin(265) + dxij*xin(229)
                                      yin(238) = yin(265) + dyij*yin(229)
                                      zin(238) = zin(265) + dzij*zin(229)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  274

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  175

                                      ! nj =    3

                                      ! i4 = i3 =  175

                                      ! do ni = 1,    3

                                      xin(175) = xin(202) + dxij*xin(166)
                                      yin(175) = yin(202) + dyij*yin(166)
                                      zin(175) = zin(202) + dzij*zin(166)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    2

                                      xin(211) = xin(238) + dxij*xin(202)
                                      yin(211) = yin(238) + dyij*yin(202)
                                      zin(211) = zin(238) + dzij*zin(202)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                      ! ni =    3

                                      xin(247) = xin(274) + dxij*xin(238)
                                      yin(247) = yin(274) + dyij*yin(238)
                                      zin(247) = zin(274) + dzij*zin(238)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  184

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  286

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(286) = xin(286) + dxij*xin(277)
                                      yin(286) = yin(286) + dyij*yin(277)
                                      zin(286) = zin(286) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  268

                                      xin(277) = xin(277) + dxij*xin(268)
                                      yin(277) = yin(277) + dyij*yin(268)
                                      zin(277) = zin(277) + dzij*zin(268)

                                      ! i3 = i4 =  268
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  259

                                      xin(268) = xin(268) + dxij*xin(259)
                                      yin(268) = yin(268) + dyij*yin(259)
                                      zin(268) = zin(268) + dzij*zin(259)

                                      ! i3 = i4 =  259
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  286

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(286) = xin(286) + dxij*xin(277)
                                      yin(286) = yin(286) + dyij*yin(277)
                                      zin(286) = zin(286) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  268

                                      xin(277) = xin(277) + dxij*xin(268)
                                      yin(277) = yin(277) + dyij*yin(268)
                                      zin(277) = zin(277) + dzij*zin(268)

                                      ! i3 = i4 =  268
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  286

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(286) = xin(286) + dxij*xin(277)
                                      yin(286) = yin(286) + dyij*yin(277)
                                      zin(286) = zin(286) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  160

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  160

                                      ! do ni = 1,    3

                                      xin(160) = xin(187) + dxij*xin(151)
                                      yin(160) = yin(187) + dyij*yin(151)
                                      zin(160) = zin(187) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  196

                                      ! ni =    2

                                      xin(196) = xin(223) + dxij*xin(187)
                                      yin(196) = yin(223) + dyij*yin(187)
                                      zin(196) = zin(223) + dzij*zin(187)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  232

                                      ! ni =    3

                                      xin(232) = xin(259) + dxij*xin(223)
                                      yin(232) = yin(259) + dyij*yin(223)
                                      zin(232) = zin(259) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  268

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  169

                                      ! nj =    2

                                      ! i4 = i3 =  169

                                      ! do ni = 1,    3

                                      xin(169) = xin(196) + dxij*xin(160)
                                      yin(169) = yin(196) + dyij*yin(160)
                                      zin(169) = zin(196) + dzij*zin(160)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    2

                                      xin(205) = xin(232) + dxij*xin(196)
                                      yin(205) = yin(232) + dyij*yin(196)
                                      zin(205) = zin(232) + dzij*zin(196)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! ni =    3

                                      xin(241) = xin(268) + dxij*xin(232)
                                      yin(241) = yin(268) + dyij*yin(232)
                                      zin(241) = zin(268) + dzij*zin(232)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  178

                                      ! nj =    3

                                      ! i4 = i3 =  178

                                      ! do ni = 1,    3

                                      xin(178) = xin(205) + dxij*xin(169)
                                      yin(178) = yin(205) + dyij*yin(169)
                                      zin(178) = zin(205) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                      ! ni =    2

                                      xin(214) = xin(241) + dxij*xin(205)
                                      yin(214) = yin(241) + dyij*yin(205)
                                      zin(214) = zin(241) + dzij*zin(205)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  250

                                      ! ni =    3

                                      xin(250) = xin(277) + dxij*xin(241)
                                      yin(250) = yin(277) + dyij*yin(241)
                                      zin(250) = zin(277) + dzij*zin(241)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  286

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  187

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  278

                                      xin(287) = xin(287) + dxij*xin(278)
                                      yin(287) = yin(287) + dyij*yin(278)
                                      zin(287) = zin(287) + dzij*zin(278)

                                      ! i3 = i4 =  278
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  269

                                      xin(278) = xin(278) + dxij*xin(269)
                                      yin(278) = yin(278) + dyij*yin(269)
                                      zin(278) = zin(278) + dzij*zin(269)

                                      ! i3 = i4 =  269
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  260

                                      xin(269) = xin(269) + dxij*xin(260)
                                      yin(269) = yin(269) + dyij*yin(260)
                                      zin(269) = zin(269) + dzij*zin(260)

                                      ! i3 = i4 =  260
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  278

                                      xin(287) = xin(287) + dxij*xin(278)
                                      yin(287) = yin(287) + dyij*yin(278)
                                      zin(287) = zin(287) + dzij*zin(278)

                                      ! i3 = i4 =  278
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  269

                                      xin(278) = xin(278) + dxij*xin(269)
                                      yin(278) = yin(278) + dyij*yin(269)
                                      zin(278) = zin(278) + dzij*zin(269)

                                      ! i3 = i4 =  269
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  278

                                      xin(287) = xin(287) + dxij*xin(278)
                                      yin(287) = yin(287) + dyij*yin(278)
                                      zin(287) = zin(287) + dzij*zin(278)

                                      ! i3 = i4 =  278
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  161

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  161

                                      ! do ni = 1,    3

                                      xin(161) = xin(188) + dxij*xin(152)
                                      yin(161) = yin(188) + dyij*yin(152)
                                      zin(161) = zin(188) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  197

                                      ! ni =    2

                                      xin(197) = xin(224) + dxij*xin(188)
                                      yin(197) = yin(224) + dyij*yin(188)
                                      zin(197) = zin(224) + dzij*zin(188)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    3

                                      xin(233) = xin(260) + dxij*xin(224)
                                      yin(233) = yin(260) + dyij*yin(224)
                                      zin(233) = zin(260) + dzij*zin(224)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  269

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  170

                                      ! nj =    2

                                      ! i4 = i3 =  170

                                      ! do ni = 1,    3

                                      xin(170) = xin(197) + dxij*xin(161)
                                      yin(170) = yin(197) + dyij*yin(161)
                                      zin(170) = zin(197) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  206

                                      ! ni =    2

                                      xin(206) = xin(233) + dxij*xin(197)
                                      yin(206) = yin(233) + dyij*yin(197)
                                      zin(206) = zin(233) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  242

                                      ! ni =    3

                                      xin(242) = xin(269) + dxij*xin(233)
                                      yin(242) = yin(269) + dyij*yin(233)
                                      zin(242) = zin(269) + dzij*zin(233)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  278

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  179

                                      ! nj =    3

                                      ! i4 = i3 =  179

                                      ! do ni = 1,    3

                                      xin(179) = xin(206) + dxij*xin(170)
                                      yin(179) = yin(206) + dyij*yin(170)
                                      zin(179) = zin(206) + dzij*zin(170)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    2

                                      xin(215) = xin(242) + dxij*xin(206)
                                      yin(215) = yin(242) + dyij*yin(206)
                                      zin(215) = zin(242) + dzij*zin(206)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    3

                                      xin(251) = xin(278) + dxij*xin(242)
                                      yin(251) = yin(278) + dyij*yin(242)
                                      zin(251) = zin(278) + dzij*zin(242)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  188

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(288) = xin(288) + dxij*xin(279)
                                      yin(288) = yin(288) + dyij*yin(279)
                                      zin(288) = zin(288) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  270

                                      xin(279) = xin(279) + dxij*xin(270)
                                      yin(279) = yin(279) + dyij*yin(270)
                                      zin(279) = zin(279) + dzij*zin(270)

                                      ! i3 = i4 =  270
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  261

                                      xin(270) = xin(270) + dxij*xin(261)
                                      yin(270) = yin(270) + dyij*yin(261)
                                      zin(270) = zin(270) + dzij*zin(261)

                                      ! i3 = i4 =  261
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(288) = xin(288) + dxij*xin(279)
                                      yin(288) = yin(288) + dyij*yin(279)
                                      zin(288) = zin(288) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  270

                                      xin(279) = xin(279) + dxij*xin(270)
                                      yin(279) = yin(279) + dyij*yin(270)
                                      zin(279) = zin(279) + dzij*zin(270)

                                      ! i3 = i4 =  270
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(288) = xin(288) + dxij*xin(279)
                                      yin(288) = yin(288) + dyij*yin(279)
                                      zin(288) = zin(288) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  162

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  162

                                      ! do ni = 1,    3

                                      xin(162) = xin(189) + dxij*xin(153)
                                      yin(162) = yin(189) + dyij*yin(153)
                                      zin(162) = zin(189) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  198

                                      ! ni =    2

                                      xin(198) = xin(225) + dxij*xin(189)
                                      yin(198) = yin(225) + dyij*yin(189)
                                      zin(198) = zin(225) + dzij*zin(189)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  234

                                      ! ni =    3

                                      xin(234) = xin(261) + dxij*xin(225)
                                      yin(234) = yin(261) + dyij*yin(225)
                                      zin(234) = zin(261) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  270

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  171

                                      ! nj =    2

                                      ! i4 = i3 =  171

                                      ! do ni = 1,    3

                                      xin(171) = xin(198) + dxij*xin(162)
                                      yin(171) = yin(198) + dyij*yin(162)
                                      zin(171) = zin(198) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    2

                                      xin(207) = xin(234) + dxij*xin(198)
                                      yin(207) = yin(234) + dyij*yin(198)
                                      zin(207) = zin(234) + dzij*zin(198)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  243

                                      ! ni =    3

                                      xin(243) = xin(270) + dxij*xin(234)
                                      yin(243) = yin(270) + dyij*yin(234)
                                      zin(243) = zin(270) + dzij*zin(234)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  279

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  180

                                      ! nj =    3

                                      ! i4 = i3 =  180

                                      ! do ni = 1,    3

                                      xin(180) = xin(207) + dxij*xin(171)
                                      yin(180) = yin(207) + dyij*yin(171)
                                      zin(180) = zin(207) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    2

                                      xin(216) = xin(243) + dxij*xin(207)
                                      yin(216) = yin(243) + dyij*yin(207)
                                      zin(216) = zin(243) + dzij*zin(207)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  252

                                      ! ni =    3

                                      xin(252) = xin(279) + dxij*xin(243)
                                      yin(252) = yin(279) + dyij*yin(243)
                                      zin(252) = zin(279) + dzij*zin(243)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  189

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  145

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  153

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  152

                                      xin(153) = xin(153) + dxkl*xin(152)
                                      yin(153) = yin(153) + dykl*yin(152)
                                      zin(153) = zin(153) + dzkl*zin(152)

                                      ! i3 = i4 =  152
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  151

                                      xin(152) = xin(152) + dxkl*xin(151)
                                      yin(152) = yin(152) + dykl*yin(151)
                                      zin(152) = zin(152) + dzkl*zin(151)

                                      ! i3 = i4 =  151
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  153

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  152

                                      xin(153) = xin(153) + dxkl*xin(152)
                                      yin(153) = yin(153) + dykl*yin(152)
                                      zin(153) = zin(153) + dzkl*zin(152)

                                      ! i3 = i4 =  152
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    2

                                      xin(146) = xin(148) + dxkl*xin(145)
                                      yin(146) = yin(148) + dykl*yin(145)
                                      zin(146) = zin(148) + dzkl*zin(145)
                                      ! i4 = i4 + lang+1 =  149

                                      ! nk =    2

                                      xin(149) = xin(151) + dxkl*xin(148)
                                      yin(149) = yin(151) + dykl*yin(148)
                                      zin(149) = zin(151) + dzkl*zin(148)
                                      ! i4 = i4 + lang+1 =  152

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! i4 = i3 =  147

                                      ! do nk = 1,    2

                                      xin(147) = xin(149) + dxkl*xin(146)
                                      yin(147) = yin(149) + dykl*yin(146)
                                      zin(147) = zin(149) + dzkl*zin(146)
                                      ! i4 = i4 + lang+1 =  150

                                      ! nk =    2

                                      xin(150) = xin(152) + dxkl*xin(149)
                                      yin(150) = yin(152) + dykl*yin(149)
                                      zin(150) = zin(152) + dzkl*zin(149)
                                      ! i4 = i4 + lang+1 =  153

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  148

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  154

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  162

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  161

                                      xin(162) = xin(162) + dxkl*xin(161)
                                      yin(162) = yin(162) + dykl*yin(161)
                                      zin(162) = zin(162) + dzkl*zin(161)

                                      ! i3 = i4 =  161
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  160

                                      xin(161) = xin(161) + dxkl*xin(160)
                                      yin(161) = yin(161) + dykl*yin(160)
                                      zin(161) = zin(161) + dzkl*zin(160)

                                      ! i3 = i4 =  160
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  162

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  161

                                      xin(162) = xin(162) + dxkl*xin(161)
                                      yin(162) = yin(162) + dykl*yin(161)
                                      zin(162) = zin(162) + dzkl*zin(161)

                                      ! i3 = i4 =  161
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  155

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  155

                                      ! do nk = 1,    2

                                      xin(155) = xin(157) + dxkl*xin(154)
                                      yin(155) = yin(157) + dykl*yin(154)
                                      zin(155) = zin(157) + dzkl*zin(154)
                                      ! i4 = i4 + lang+1 =  158

                                      ! nk =    2

                                      xin(158) = xin(160) + dxkl*xin(157)
                                      yin(158) = yin(160) + dykl*yin(157)
                                      zin(158) = zin(160) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  161

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  156

                                      ! nl =    2

                                      ! i4 = i3 =  156

                                      ! do nk = 1,    2

                                      xin(156) = xin(158) + dxkl*xin(155)
                                      yin(156) = yin(158) + dykl*yin(155)
                                      zin(156) = zin(158) + dzkl*zin(155)
                                      ! i4 = i4 + lang+1 =  159

                                      ! nk =    2

                                      xin(159) = xin(161) + dxkl*xin(158)
                                      yin(159) = yin(161) + dykl*yin(158)
                                      zin(159) = zin(161) + dzkl*zin(158)
                                      ! i4 = i4 + lang+1 =  162

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  157

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  163

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  171

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  170

                                      xin(171) = xin(171) + dxkl*xin(170)
                                      yin(171) = yin(171) + dykl*yin(170)
                                      zin(171) = zin(171) + dzkl*zin(170)

                                      ! i3 = i4 =  170
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  169

                                      xin(170) = xin(170) + dxkl*xin(169)
                                      yin(170) = yin(170) + dykl*yin(169)
                                      zin(170) = zin(170) + dzkl*zin(169)

                                      ! i3 = i4 =  169
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  171

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  170

                                      xin(171) = xin(171) + dxkl*xin(170)
                                      yin(171) = yin(171) + dykl*yin(170)
                                      zin(171) = zin(171) + dzkl*zin(170)

                                      ! i3 = i4 =  170
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  164

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  164

                                      ! do nk = 1,    2

                                      xin(164) = xin(166) + dxkl*xin(163)
                                      yin(164) = yin(166) + dykl*yin(163)
                                      zin(164) = zin(166) + dzkl*zin(163)
                                      ! i4 = i4 + lang+1 =  167

                                      ! nk =    2

                                      xin(167) = xin(169) + dxkl*xin(166)
                                      yin(167) = yin(169) + dykl*yin(166)
                                      zin(167) = zin(169) + dzkl*zin(166)
                                      ! i4 = i4 + lang+1 =  170

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  165

                                      ! nl =    2

                                      ! i4 = i3 =  165

                                      ! do nk = 1,    2

                                      xin(165) = xin(167) + dxkl*xin(164)
                                      yin(165) = yin(167) + dykl*yin(164)
                                      zin(165) = zin(167) + dzkl*zin(164)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    2

                                      xin(168) = xin(170) + dxkl*xin(167)
                                      yin(168) = yin(170) + dykl*yin(167)
                                      zin(168) = zin(170) + dzkl*zin(167)
                                      ! i4 = i4 + lang+1 =  171

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  166

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  172

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  178

                                      xin(179) = xin(179) + dxkl*xin(178)
                                      yin(179) = yin(179) + dykl*yin(178)
                                      zin(179) = zin(179) + dzkl*zin(178)

                                      ! i3 = i4 =  178
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  173

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  173

                                      ! do nk = 1,    2

                                      xin(173) = xin(175) + dxkl*xin(172)
                                      yin(173) = yin(175) + dykl*yin(172)
                                      zin(173) = zin(175) + dzkl*zin(172)
                                      ! i4 = i4 + lang+1 =  176

                                      ! nk =    2

                                      xin(176) = xin(178) + dxkl*xin(175)
                                      yin(176) = yin(178) + dykl*yin(175)
                                      zin(176) = zin(178) + dzkl*zin(175)
                                      ! i4 = i4 + lang+1 =  179

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  174

                                      ! nl =    2

                                      ! i4 = i3 =  174

                                      ! do nk = 1,    2

                                      xin(174) = xin(176) + dxkl*xin(173)
                                      yin(174) = yin(176) + dykl*yin(173)
                                      zin(174) = zin(176) + dzkl*zin(173)
                                      ! i4 = i4 + lang+1 =  177

                                      ! nk =    2

                                      xin(177) = xin(179) + dxkl*xin(176)
                                      yin(177) = yin(179) + dykl*yin(176)
                                      zin(177) = zin(179) + dzkl*zin(176)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  175

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  181

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  189

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  188

                                      xin(189) = xin(189) + dxkl*xin(188)
                                      yin(189) = yin(189) + dykl*yin(188)
                                      zin(189) = zin(189) + dzkl*zin(188)

                                      ! i3 = i4 =  188
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  187

                                      xin(188) = xin(188) + dxkl*xin(187)
                                      yin(188) = yin(188) + dykl*yin(187)
                                      zin(188) = zin(188) + dzkl*zin(187)

                                      ! i3 = i4 =  187
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  189

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  188

                                      xin(189) = xin(189) + dxkl*xin(188)
                                      yin(189) = yin(189) + dykl*yin(188)
                                      zin(189) = zin(189) + dzkl*zin(188)

                                      ! i3 = i4 =  188
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  182

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  182

                                      ! do nk = 1,    2

                                      xin(182) = xin(184) + dxkl*xin(181)
                                      yin(182) = yin(184) + dykl*yin(181)
                                      zin(182) = zin(184) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  185

                                      ! nk =    2

                                      xin(185) = xin(187) + dxkl*xin(184)
                                      yin(185) = yin(187) + dykl*yin(184)
                                      zin(185) = zin(187) + dzkl*zin(184)
                                      ! i4 = i4 + lang+1 =  188

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  183

                                      ! nl =    2

                                      ! i4 = i3 =  183

                                      ! do nk = 1,    2

                                      xin(183) = xin(185) + dxkl*xin(182)
                                      yin(183) = yin(185) + dykl*yin(182)
                                      zin(183) = zin(185) + dzkl*zin(182)
                                      ! i4 = i4 + lang+1 =  186

                                      ! nk =    2

                                      xin(186) = xin(188) + dxkl*xin(185)
                                      yin(186) = yin(188) + dykl*yin(185)
                                      zin(186) = zin(188) + dzkl*zin(185)
                                      ! i4 = i4 + lang+1 =  189

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  184

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  190

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  198

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  197

                                      xin(198) = xin(198) + dxkl*xin(197)
                                      yin(198) = yin(198) + dykl*yin(197)
                                      zin(198) = zin(198) + dzkl*zin(197)

                                      ! i3 = i4 =  197
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  196

                                      xin(197) = xin(197) + dxkl*xin(196)
                                      yin(197) = yin(197) + dykl*yin(196)
                                      zin(197) = zin(197) + dzkl*zin(196)

                                      ! i3 = i4 =  196
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  198

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  197

                                      xin(198) = xin(198) + dxkl*xin(197)
                                      yin(198) = yin(198) + dykl*yin(197)
                                      zin(198) = zin(198) + dzkl*zin(197)

                                      ! i3 = i4 =  197
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  191

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  191

                                      ! do nk = 1,    2

                                      xin(191) = xin(193) + dxkl*xin(190)
                                      yin(191) = yin(193) + dykl*yin(190)
                                      zin(191) = zin(193) + dzkl*zin(190)
                                      ! i4 = i4 + lang+1 =  194

                                      ! nk =    2

                                      xin(194) = xin(196) + dxkl*xin(193)
                                      yin(194) = yin(196) + dykl*yin(193)
                                      zin(194) = zin(196) + dzkl*zin(193)
                                      ! i4 = i4 + lang+1 =  197

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  192

                                      ! nl =    2

                                      ! i4 = i3 =  192

                                      ! do nk = 1,    2

                                      xin(192) = xin(194) + dxkl*xin(191)
                                      yin(192) = yin(194) + dykl*yin(191)
                                      zin(192) = zin(194) + dzkl*zin(191)
                                      ! i4 = i4 + lang+1 =  195

                                      ! nk =    2

                                      xin(195) = xin(197) + dxkl*xin(194)
                                      yin(195) = yin(197) + dykl*yin(194)
                                      zin(195) = zin(197) + dzkl*zin(194)
                                      ! i4 = i4 + lang+1 =  198

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  193

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  199

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  207

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  206

                                      xin(207) = xin(207) + dxkl*xin(206)
                                      yin(207) = yin(207) + dykl*yin(206)
                                      zin(207) = zin(207) + dzkl*zin(206)

                                      ! i3 = i4 =  206
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  205

                                      xin(206) = xin(206) + dxkl*xin(205)
                                      yin(206) = yin(206) + dykl*yin(205)
                                      zin(206) = zin(206) + dzkl*zin(205)

                                      ! i3 = i4 =  205
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  207

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  206

                                      xin(207) = xin(207) + dxkl*xin(206)
                                      yin(207) = yin(207) + dykl*yin(206)
                                      zin(207) = zin(207) + dzkl*zin(206)

                                      ! i3 = i4 =  206
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  200

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  200

                                      ! do nk = 1,    2

                                      xin(200) = xin(202) + dxkl*xin(199)
                                      yin(200) = yin(202) + dykl*yin(199)
                                      zin(200) = zin(202) + dzkl*zin(199)
                                      ! i4 = i4 + lang+1 =  203

                                      ! nk =    2

                                      xin(203) = xin(205) + dxkl*xin(202)
                                      yin(203) = yin(205) + dykl*yin(202)
                                      zin(203) = zin(205) + dzkl*zin(202)
                                      ! i4 = i4 + lang+1 =  206

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  201

                                      ! nl =    2

                                      ! i4 = i3 =  201

                                      ! do nk = 1,    2

                                      xin(201) = xin(203) + dxkl*xin(200)
                                      yin(201) = yin(203) + dykl*yin(200)
                                      zin(201) = zin(203) + dzkl*zin(200)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    2

                                      xin(204) = xin(206) + dxkl*xin(203)
                                      yin(204) = yin(206) + dykl*yin(203)
                                      zin(204) = zin(206) + dzkl*zin(203)
                                      ! i4 = i4 + lang+1 =  207

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  202

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  208

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  214

                                      xin(215) = xin(215) + dxkl*xin(214)
                                      yin(215) = yin(215) + dykl*yin(214)
                                      zin(215) = zin(215) + dzkl*zin(214)

                                      ! i3 = i4 =  214
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  209

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  209

                                      ! do nk = 1,    2

                                      xin(209) = xin(211) + dxkl*xin(208)
                                      yin(209) = yin(211) + dykl*yin(208)
                                      zin(209) = zin(211) + dzkl*zin(208)
                                      ! i4 = i4 + lang+1 =  212

                                      ! nk =    2

                                      xin(212) = xin(214) + dxkl*xin(211)
                                      yin(212) = yin(214) + dykl*yin(211)
                                      zin(212) = zin(214) + dzkl*zin(211)
                                      ! i4 = i4 + lang+1 =  215

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  210

                                      ! nl =    2

                                      ! i4 = i3 =  210

                                      ! do nk = 1,    2

                                      xin(210) = xin(212) + dxkl*xin(209)
                                      yin(210) = yin(212) + dykl*yin(209)
                                      zin(210) = zin(212) + dzkl*zin(209)
                                      ! i4 = i4 + lang+1 =  213

                                      ! nk =    2

                                      xin(213) = xin(215) + dxkl*xin(212)
                                      yin(213) = yin(215) + dykl*yin(212)
                                      zin(213) = zin(215) + dzkl*zin(212)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  211

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  225

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  224

                                      xin(225) = xin(225) + dxkl*xin(224)
                                      yin(225) = yin(225) + dykl*yin(224)
                                      zin(225) = zin(225) + dzkl*zin(224)

                                      ! i3 = i4 =  224
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  223

                                      xin(224) = xin(224) + dxkl*xin(223)
                                      yin(224) = yin(224) + dykl*yin(223)
                                      zin(224) = zin(224) + dzkl*zin(223)

                                      ! i3 = i4 =  223
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  225

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  224

                                      xin(225) = xin(225) + dxkl*xin(224)
                                      yin(225) = yin(225) + dykl*yin(224)
                                      zin(225) = zin(225) + dzkl*zin(224)

                                      ! i3 = i4 =  224
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    2

                                      xin(218) = xin(220) + dxkl*xin(217)
                                      yin(218) = yin(220) + dykl*yin(217)
                                      zin(218) = zin(220) + dzkl*zin(217)
                                      ! i4 = i4 + lang+1 =  221

                                      ! nk =    2

                                      xin(221) = xin(223) + dxkl*xin(220)
                                      yin(221) = yin(223) + dykl*yin(220)
                                      zin(221) = zin(223) + dzkl*zin(220)
                                      ! i4 = i4 + lang+1 =  224

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! i4 = i3 =  219

                                      ! do nk = 1,    2

                                      xin(219) = xin(221) + dxkl*xin(218)
                                      yin(219) = yin(221) + dykl*yin(218)
                                      zin(219) = zin(221) + dzkl*zin(218)
                                      ! i4 = i4 + lang+1 =  222

                                      ! nk =    2

                                      xin(222) = xin(224) + dxkl*xin(221)
                                      yin(222) = yin(224) + dykl*yin(221)
                                      zin(222) = zin(224) + dzkl*zin(221)
                                      ! i4 = i4 + lang+1 =  225

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  220

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  226

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  234

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  233

                                      xin(234) = xin(234) + dxkl*xin(233)
                                      yin(234) = yin(234) + dykl*yin(233)
                                      zin(234) = zin(234) + dzkl*zin(233)

                                      ! i3 = i4 =  233
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  232

                                      xin(233) = xin(233) + dxkl*xin(232)
                                      yin(233) = yin(233) + dykl*yin(232)
                                      zin(233) = zin(233) + dzkl*zin(232)

                                      ! i3 = i4 =  232
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  234

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  233

                                      xin(234) = xin(234) + dxkl*xin(233)
                                      yin(234) = yin(234) + dykl*yin(233)
                                      zin(234) = zin(234) + dzkl*zin(233)

                                      ! i3 = i4 =  233
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  227

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  227

                                      ! do nk = 1,    2

                                      xin(227) = xin(229) + dxkl*xin(226)
                                      yin(227) = yin(229) + dykl*yin(226)
                                      zin(227) = zin(229) + dzkl*zin(226)
                                      ! i4 = i4 + lang+1 =  230

                                      ! nk =    2

                                      xin(230) = xin(232) + dxkl*xin(229)
                                      yin(230) = yin(232) + dykl*yin(229)
                                      zin(230) = zin(232) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  233

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  228

                                      ! nl =    2

                                      ! i4 = i3 =  228

                                      ! do nk = 1,    2

                                      xin(228) = xin(230) + dxkl*xin(227)
                                      yin(228) = yin(230) + dykl*yin(227)
                                      zin(228) = zin(230) + dzkl*zin(227)
                                      ! i4 = i4 + lang+1 =  231

                                      ! nk =    2

                                      xin(231) = xin(233) + dxkl*xin(230)
                                      yin(231) = yin(233) + dykl*yin(230)
                                      zin(231) = zin(233) + dzkl*zin(230)
                                      ! i4 = i4 + lang+1 =  234

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  229

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  235

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  243

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  242

                                      xin(243) = xin(243) + dxkl*xin(242)
                                      yin(243) = yin(243) + dykl*yin(242)
                                      zin(243) = zin(243) + dzkl*zin(242)

                                      ! i3 = i4 =  242
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  241

                                      xin(242) = xin(242) + dxkl*xin(241)
                                      yin(242) = yin(242) + dykl*yin(241)
                                      zin(242) = zin(242) + dzkl*zin(241)

                                      ! i3 = i4 =  241
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  243

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  242

                                      xin(243) = xin(243) + dxkl*xin(242)
                                      yin(243) = yin(243) + dykl*yin(242)
                                      zin(243) = zin(243) + dzkl*zin(242)

                                      ! i3 = i4 =  242
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  236

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  236

                                      ! do nk = 1,    2

                                      xin(236) = xin(238) + dxkl*xin(235)
                                      yin(236) = yin(238) + dykl*yin(235)
                                      zin(236) = zin(238) + dzkl*zin(235)
                                      ! i4 = i4 + lang+1 =  239

                                      ! nk =    2

                                      xin(239) = xin(241) + dxkl*xin(238)
                                      yin(239) = yin(241) + dykl*yin(238)
                                      zin(239) = zin(241) + dzkl*zin(238)
                                      ! i4 = i4 + lang+1 =  242

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  237

                                      ! nl =    2

                                      ! i4 = i3 =  237

                                      ! do nk = 1,    2

                                      xin(237) = xin(239) + dxkl*xin(236)
                                      yin(237) = yin(239) + dykl*yin(236)
                                      zin(237) = zin(239) + dzkl*zin(236)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    2

                                      xin(240) = xin(242) + dxkl*xin(239)
                                      yin(240) = yin(242) + dykl*yin(239)
                                      zin(240) = zin(242) + dzkl*zin(239)
                                      ! i4 = i4 + lang+1 =  243

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  238

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  244

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  250

                                      xin(251) = xin(251) + dxkl*xin(250)
                                      yin(251) = yin(251) + dykl*yin(250)
                                      zin(251) = zin(251) + dzkl*zin(250)

                                      ! i3 = i4 =  250
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  245

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  245

                                      ! do nk = 1,    2

                                      xin(245) = xin(247) + dxkl*xin(244)
                                      yin(245) = yin(247) + dykl*yin(244)
                                      zin(245) = zin(247) + dzkl*zin(244)
                                      ! i4 = i4 + lang+1 =  248

                                      ! nk =    2

                                      xin(248) = xin(250) + dxkl*xin(247)
                                      yin(248) = yin(250) + dykl*yin(247)
                                      zin(248) = zin(250) + dzkl*zin(247)
                                      ! i4 = i4 + lang+1 =  251

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  246

                                      ! nl =    2

                                      ! i4 = i3 =  246

                                      ! do nk = 1,    2

                                      xin(246) = xin(248) + dxkl*xin(245)
                                      yin(246) = yin(248) + dykl*yin(245)
                                      zin(246) = zin(248) + dzkl*zin(245)
                                      ! i4 = i4 + lang+1 =  249

                                      ! nk =    2

                                      xin(249) = xin(251) + dxkl*xin(248)
                                      yin(249) = yin(251) + dykl*yin(248)
                                      zin(249) = zin(251) + dzkl*zin(248)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  247

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  253

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  261

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  260

                                      xin(261) = xin(261) + dxkl*xin(260)
                                      yin(261) = yin(261) + dykl*yin(260)
                                      zin(261) = zin(261) + dzkl*zin(260)

                                      ! i3 = i4 =  260
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  259

                                      xin(260) = xin(260) + dxkl*xin(259)
                                      yin(260) = yin(260) + dykl*yin(259)
                                      zin(260) = zin(260) + dzkl*zin(259)

                                      ! i3 = i4 =  259
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  261

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  260

                                      xin(261) = xin(261) + dxkl*xin(260)
                                      yin(261) = yin(261) + dykl*yin(260)
                                      zin(261) = zin(261) + dzkl*zin(260)

                                      ! i3 = i4 =  260
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  254

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  254

                                      ! do nk = 1,    2

                                      xin(254) = xin(256) + dxkl*xin(253)
                                      yin(254) = yin(256) + dykl*yin(253)
                                      zin(254) = zin(256) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  257

                                      ! nk =    2

                                      xin(257) = xin(259) + dxkl*xin(256)
                                      yin(257) = yin(259) + dykl*yin(256)
                                      zin(257) = zin(259) + dzkl*zin(256)
                                      ! i4 = i4 + lang+1 =  260

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  255

                                      ! nl =    2

                                      ! i4 = i3 =  255

                                      ! do nk = 1,    2

                                      xin(255) = xin(257) + dxkl*xin(254)
                                      yin(255) = yin(257) + dykl*yin(254)
                                      zin(255) = zin(257) + dzkl*zin(254)
                                      ! i4 = i4 + lang+1 =  258

                                      ! nk =    2

                                      xin(258) = xin(260) + dxkl*xin(257)
                                      yin(258) = yin(260) + dykl*yin(257)
                                      zin(258) = zin(260) + dzkl*zin(257)
                                      ! i4 = i4 + lang+1 =  261

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  256

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  262

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  270

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  269

                                      xin(270) = xin(270) + dxkl*xin(269)
                                      yin(270) = yin(270) + dykl*yin(269)
                                      zin(270) = zin(270) + dzkl*zin(269)

                                      ! i3 = i4 =  269
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  268

                                      xin(269) = xin(269) + dxkl*xin(268)
                                      yin(269) = yin(269) + dykl*yin(268)
                                      zin(269) = zin(269) + dzkl*zin(268)

                                      ! i3 = i4 =  268
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  270

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  269

                                      xin(270) = xin(270) + dxkl*xin(269)
                                      yin(270) = yin(270) + dykl*yin(269)
                                      zin(270) = zin(270) + dzkl*zin(269)

                                      ! i3 = i4 =  269
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  263

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  263

                                      ! do nk = 1,    2

                                      xin(263) = xin(265) + dxkl*xin(262)
                                      yin(263) = yin(265) + dykl*yin(262)
                                      zin(263) = zin(265) + dzkl*zin(262)
                                      ! i4 = i4 + lang+1 =  266

                                      ! nk =    2

                                      xin(266) = xin(268) + dxkl*xin(265)
                                      yin(266) = yin(268) + dykl*yin(265)
                                      zin(266) = zin(268) + dzkl*zin(265)
                                      ! i4 = i4 + lang+1 =  269

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  264

                                      ! nl =    2

                                      ! i4 = i3 =  264

                                      ! do nk = 1,    2

                                      xin(264) = xin(266) + dxkl*xin(263)
                                      yin(264) = yin(266) + dykl*yin(263)
                                      zin(264) = zin(266) + dzkl*zin(263)
                                      ! i4 = i4 + lang+1 =  267

                                      ! nk =    2

                                      xin(267) = xin(269) + dxkl*xin(266)
                                      yin(267) = yin(269) + dykl*yin(266)
                                      zin(267) = zin(269) + dzkl*zin(266)
                                      ! i4 = i4 + lang+1 =  270

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  265

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  271

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  279

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  278

                                      xin(279) = xin(279) + dxkl*xin(278)
                                      yin(279) = yin(279) + dykl*yin(278)
                                      zin(279) = zin(279) + dzkl*zin(278)

                                      ! i3 = i4 =  278
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  277

                                      xin(278) = xin(278) + dxkl*xin(277)
                                      yin(278) = yin(278) + dykl*yin(277)
                                      zin(278) = zin(278) + dzkl*zin(277)

                                      ! i3 = i4 =  277
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  279

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  278

                                      xin(279) = xin(279) + dxkl*xin(278)
                                      yin(279) = yin(279) + dykl*yin(278)
                                      zin(279) = zin(279) + dzkl*zin(278)

                                      ! i3 = i4 =  278
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  272

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  272

                                      ! do nk = 1,    2

                                      xin(272) = xin(274) + dxkl*xin(271)
                                      yin(272) = yin(274) + dykl*yin(271)
                                      zin(272) = zin(274) + dzkl*zin(271)
                                      ! i4 = i4 + lang+1 =  275

                                      ! nk =    2

                                      xin(275) = xin(277) + dxkl*xin(274)
                                      yin(275) = yin(277) + dykl*yin(274)
                                      zin(275) = zin(277) + dzkl*zin(274)
                                      ! i4 = i4 + lang+1 =  278

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  273

                                      ! nl =    2

                                      ! i4 = i3 =  273

                                      ! do nk = 1,    2

                                      xin(273) = xin(275) + dxkl*xin(272)
                                      yin(273) = yin(275) + dykl*yin(272)
                                      zin(273) = zin(275) + dzkl*zin(272)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    2

                                      xin(276) = xin(278) + dxkl*xin(275)
                                      yin(276) = yin(278) + dykl*yin(275)
                                      zin(276) = zin(278) + dzkl*zin(275)
                                      ! i4 = i4 + lang+1 =  279

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  274

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  280

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  286

                                      xin(287) = xin(287) + dxkl*xin(286)
                                      yin(287) = yin(287) + dykl*yin(286)
                                      zin(287) = zin(287) + dzkl*zin(286)

                                      ! i3 = i4 =  286
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  281

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  281

                                      ! do nk = 1,    2

                                      xin(281) = xin(283) + dxkl*xin(280)
                                      yin(281) = yin(283) + dykl*yin(280)
                                      zin(281) = zin(283) + dzkl*zin(280)
                                      ! i4 = i4 + lang+1 =  284

                                      ! nk =    2

                                      xin(284) = xin(286) + dxkl*xin(283)
                                      yin(284) = yin(286) + dykl*yin(283)
                                      zin(284) = zin(286) + dzkl*zin(283)
                                      ! i4 = i4 + lang+1 =  287

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  282

                                      ! nl =    2

                                      ! i4 = i3 =  282

                                      ! do nk = 1,    2

                                      xin(282) = xin(284) + dxkl*xin(281)
                                      yin(282) = yin(284) + dykl*yin(281)
                                      zin(282) = zin(284) + dzkl*zin(281)
                                      ! i4 = i4 + lang+1 =  285

                                      ! nk =    2

                                      xin(285) = xin(287) + dxkl*xin(284)
                                      yin(285) = yin(287) + dykl*yin(284)
                                      zin(285) = zin(287) + dzkl*zin(284)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  283

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  288

                                      u2 = roots(3)*rho
                                      f00 = expe*wghts(3)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  289

                                      xin(289) = 1.0_dp
                                      yin(289) = 1.0_dp
                                      zin(289) = f00

                                      ! i2 = in(2) =  325
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(325) = xc00
                                      yin(325) = yc00
                                      zin(325) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  292

                                      xin(292) = xcp00
                                      yin(292) = ycp00
                                      zin(292) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  328
                                      ! i2 =  325

                                      xin(328) = xcp00*xin(325) + cp10
                                      yin(328) = ycp00*yin(325) + cp10
                                      zin(328) = zcp00*zin(325) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  325

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  361
                                      ! i3 =  289
                                      ! i4 =  325

                                      xin(361) = c10*xin(289) + xc00*xin(325)
                                      yin(361) = c10*yin(289) + yc00*yin(325)
                                      zin(361) = c10*zin(289) + zc00*zin(325)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  364
                                      ! i5 =  361
                                      ! i4 =  325

                                      xin(364) = xcp00*xin(361) + cp10*xin(325)
                                      yin(364) = ycp00*yin(361) + cp10*yin(325)
                                      zin(364) = zcp00*zin(361) + cp10*zin(325)

                                      ! ------------------

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  361

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  397
                                      ! i3 =  325
                                      ! i4 =  361

                                      xin(397) = c10*xin(325) + xc00*xin(361)
                                      yin(397) = c10*yin(325) + yc00*yin(361)
                                      zin(397) = c10*zin(325) + zc00*zin(361)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  400
                                      ! i5 =  397
                                      ! i4 =  361

                                      xin(400) = xcp00*xin(397) + cp10*xin(361)
                                      yin(400) = ycp00*yin(397) + cp10*yin(361)
                                      zin(400) = zcp00*zin(397) + cp10*zin(361)

                                      ! ------------------

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  397

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  406
                                      ! i3 =  361
                                      ! i4 =  397

                                      xin(406) = c10*xin(361) + xc00*xin(397)
                                      yin(406) = c10*yin(361) + yc00*yin(397)
                                      zin(406) = c10*zin(361) + zc00*zin(397)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  409
                                      ! i5 =  406
                                      ! i4 =  397

                                      xin(409) = xcp00*xin(406) + cp10*xin(397)
                                      yin(409) = ycp00*yin(406) + cp10*yin(397)
                                      zin(409) = zcp00*zin(406) + cp10*zin(397)

                                      ! ------------------

                                      ! i3 = i4 =  397
                                      ! i4 = i5 =  406

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  415
                                      ! i3 =  397
                                      ! i4 =  406

                                      xin(415) = c10*xin(397) + xc00*xin(406)
                                      yin(415) = c10*yin(397) + yc00*yin(406)
                                      zin(415) = c10*zin(397) + zc00*zin(406)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  418
                                      ! i5 =  415
                                      ! i4 =  406

                                      xin(418) = xcp00*xin(415) + cp10*xin(406)
                                      yin(418) = ycp00*yin(415) + cp10*yin(406)
                                      zin(418) = zcp00*zin(415) + cp10*zin(406)

                                      ! ------------------

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  424
                                      ! i3 =  406
                                      ! i4 =  415

                                      xin(424) = c10*xin(406) + xc00*xin(415)
                                      yin(424) = c10*yin(406) + yc00*yin(415)
                                      zin(424) = c10*zin(406) + zc00*zin(415)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  427
                                      ! i5 =  424
                                      ! i4 =  415

                                      xin(427) = xcp00*xin(424) + cp10*xin(415)
                                      yin(427) = ycp00*yin(424) + cp10*yin(415)
                                      zin(427) = zcp00*zin(424) + cp10*zin(415)

                                      ! ------------------

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  289
                                      ! i4 = i1+k2 =  292

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  295
                                      ! i3 =  289
                                      ! i4 =  292

                                      xin(295) = cp01*xin(289) + xcp00*xin(292)
                                      yin(295) = cp01*yin(289) + ycp00*yin(292)
                                      zin(295) = cp01*zin(289) + zcp00*zin(292)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  331

                                      xin(331) = xc00*xin(295) + c01*xin(292)
                                      yin(331) = yc00*yin(295) + c01*yin(292)
                                      zin(331) = zc00*zin(295) + c01*zin(292)

                                      ! ------------------

                                      ! i3 = i4 =  292
                                      ! i4 = i5 =  295

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  296
                                      ! i3 =  292
                                      ! i4 =  295

                                      xin(296) = cp01*xin(292) + xcp00*xin(295)
                                      yin(296) = cp01*yin(292) + ycp00*yin(295)
                                      zin(296) = cp01*zin(292) + zcp00*zin(295)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  332

                                      xin(332) = xc00*xin(296) + c01*xin(295)
                                      yin(332) = yc00*yin(296) + c01*yin(295)
                                      zin(332) = zc00*zin(296) + c01*zin(295)

                                      ! ------------------

                                      ! i3 = i4 =  295
                                      ! i4 = i5 =  296

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  297
                                      ! i3 =  295
                                      ! i4 =  296

                                      xin(297) = cp01*xin(295) + xcp00*xin(296)
                                      yin(297) = cp01*yin(295) + ycp00*yin(296)
                                      zin(297) = cp01*zin(295) + zcp00*zin(296)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  333

                                      xin(333) = xc00*xin(297) + c01*xin(296)
                                      yin(333) = yc00*yin(297) + c01*yin(296)
                                      zin(333) = zc00*zin(297) + c01*zin(296)

                                      ! ------------------

                                      ! i3 = i4 =  296
                                      ! i4 = i5 =  297

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  325

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  361

                                      xin(367) = c10*xin(295) + xc00*xin(331) + c01*xin(328)
                                      yin(367) = c10*yin(295) + yc00*yin(331) + c01*yin(328)
                                      zin(367) = c10*zin(295) + zc00*zin(331) + c01*zin(328)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  361

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  397

                                      xin(403) = c10*xin(331) + xc00*xin(367) + c01*xin(364)
                                      yin(403) = c10*yin(331) + yc00*yin(367) + c01*yin(364)
                                      zin(403) = c10*zin(331) + zc00*zin(367) + c01*zin(364)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  397

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  406

                                      xin(412) = c10*xin(367) + xc00*xin(403) + c01*xin(400)
                                      yin(412) = c10*yin(367) + yc00*yin(403) + c01*yin(400)
                                      zin(412) = c10*zin(367) + zc00*zin(403) + c01*zin(400)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  397
                                      ! i4 = i5 =  406

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  415

                                      xin(421) = c10*xin(403) + xc00*xin(412) + c01*xin(409)
                                      yin(421) = c10*yin(403) + yc00*yin(412) + c01*yin(409)
                                      zin(421) = c10*zin(403) + zc00*zin(412) + c01*zin(409)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  424

                                      xin(430) = c10*xin(412) + xc00*xin(421) + c01*xin(418)
                                      yin(430) = c10*yin(412) + yc00*yin(421) + c01*yin(418)
                                      zin(430) = c10*zin(412) + zc00*zin(421) + c01*zin(418)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  325

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  361

                                      xin(368) = c10*xin(296) + xc00*xin(332) + c01*xin(331)
                                      yin(368) = c10*yin(296) + yc00*yin(332) + c01*yin(331)
                                      zin(368) = c10*zin(296) + zc00*zin(332) + c01*zin(331)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  361

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  397

                                      xin(404) = c10*xin(332) + xc00*xin(368) + c01*xin(367)
                                      yin(404) = c10*yin(332) + yc00*yin(368) + c01*yin(367)
                                      zin(404) = c10*zin(332) + zc00*zin(368) + c01*zin(367)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  397

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  406

                                      xin(413) = c10*xin(368) + xc00*xin(404) + c01*xin(403)
                                      yin(413) = c10*yin(368) + yc00*yin(404) + c01*yin(403)
                                      zin(413) = c10*zin(368) + zc00*zin(404) + c01*zin(403)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  397
                                      ! i4 = i5 =  406

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  415

                                      xin(422) = c10*xin(404) + xc00*xin(413) + c01*xin(412)
                                      yin(422) = c10*yin(404) + yc00*yin(413) + c01*yin(412)
                                      zin(422) = c10*zin(404) + zc00*zin(413) + c01*zin(412)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  424

                                      xin(431) = c10*xin(413) + xc00*xin(422) + c01*xin(421)
                                      yin(431) = c10*yin(413) + yc00*yin(422) + c01*yin(421)
                                      zin(431) = c10*zin(413) + zc00*zin(422) + c01*zin(421)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  325

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  361

                                      xin(369) = c10*xin(297) + xc00*xin(333) + c01*xin(332)
                                      yin(369) = c10*yin(297) + yc00*yin(333) + c01*yin(332)
                                      zin(369) = c10*zin(297) + zc00*zin(333) + c01*zin(332)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  361

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  397

                                      xin(405) = c10*xin(333) + xc00*xin(369) + c01*xin(368)
                                      yin(405) = c10*yin(333) + yc00*yin(369) + c01*yin(368)
                                      zin(405) = c10*zin(333) + zc00*zin(369) + c01*zin(368)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  397

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  406

                                      xin(414) = c10*xin(369) + xc00*xin(405) + c01*xin(404)
                                      yin(414) = c10*yin(369) + yc00*yin(405) + c01*yin(404)
                                      zin(414) = c10*zin(369) + zc00*zin(405) + c01*zin(404)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  397
                                      ! i4 = i5 =  406

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  415

                                      xin(423) = c10*xin(405) + xc00*xin(414) + c01*xin(413)
                                      yin(423) = c10*yin(405) + yc00*yin(414) + c01*yin(413)
                                      zin(423) = c10*zin(405) + zc00*zin(414) + c01*zin(413)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  424

                                      xin(432) = c10*xin(414) + xc00*xin(423) + c01*xin(422)
                                      yin(432) = c10*yin(414) + yc00*yin(423) + c01*yin(422)
                                      zin(432) = c10*zin(414) + zc00*zin(423) + c01*zin(422)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  424

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  424

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  415

                                      xin(424) = xin(424) + dxij*xin(415)
                                      yin(424) = yin(424) + dyij*yin(415)
                                      zin(424) = zin(424) + dzij*zin(415)

                                      ! i3 = i4 =  415
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  406

                                      xin(415) = xin(415) + dxij*xin(406)
                                      yin(415) = yin(415) + dyij*yin(406)
                                      zin(415) = zin(415) + dzij*zin(406)

                                      ! i3 = i4 =  406
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  397

                                      xin(406) = xin(406) + dxij*xin(397)
                                      yin(406) = yin(406) + dyij*yin(397)
                                      zin(406) = zin(406) + dzij*zin(397)

                                      ! i3 = i4 =  397
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  424

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  415

                                      xin(424) = xin(424) + dxij*xin(415)
                                      yin(424) = yin(424) + dyij*yin(415)
                                      zin(424) = zin(424) + dzij*zin(415)

                                      ! i3 = i4 =  415
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  406

                                      xin(415) = xin(415) + dxij*xin(406)
                                      yin(415) = yin(415) + dyij*yin(406)
                                      zin(415) = zin(415) + dzij*zin(406)

                                      ! i3 = i4 =  406
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  424

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  415

                                      xin(424) = xin(424) + dxij*xin(415)
                                      yin(424) = yin(424) + dyij*yin(415)
                                      zin(424) = zin(424) + dzij*zin(415)

                                      ! i3 = i4 =  415
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  298

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  298

                                      ! do ni = 1,    3

                                      xin(298) = xin(325) + dxij*xin(289)
                                      yin(298) = yin(325) + dyij*yin(289)
                                      zin(298) = zin(325) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  334

                                      ! ni =    2

                                      xin(334) = xin(361) + dxij*xin(325)
                                      yin(334) = yin(361) + dyij*yin(325)
                                      zin(334) = zin(361) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  370

                                      ! ni =    3

                                      xin(370) = xin(397) + dxij*xin(361)
                                      yin(370) = yin(397) + dyij*yin(361)
                                      zin(370) = zin(397) + dzij*zin(361)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  406

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  307

                                      ! nj =    2

                                      ! i4 = i3 =  307

                                      ! do ni = 1,    3

                                      xin(307) = xin(334) + dxij*xin(298)
                                      yin(307) = yin(334) + dyij*yin(298)
                                      zin(307) = zin(334) + dzij*zin(298)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  343

                                      ! ni =    2

                                      xin(343) = xin(370) + dxij*xin(334)
                                      yin(343) = yin(370) + dyij*yin(334)
                                      zin(343) = zin(370) + dzij*zin(334)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  379

                                      ! ni =    3

                                      xin(379) = xin(406) + dxij*xin(370)
                                      yin(379) = yin(406) + dyij*yin(370)
                                      zin(379) = zin(406) + dzij*zin(370)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  415

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  316

                                      ! nj =    3

                                      ! i4 = i3 =  316

                                      ! do ni = 1,    3

                                      xin(316) = xin(343) + dxij*xin(307)
                                      yin(316) = yin(343) + dyij*yin(307)
                                      zin(316) = zin(343) + dzij*zin(307)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  352

                                      ! ni =    2

                                      xin(352) = xin(379) + dxij*xin(343)
                                      yin(352) = yin(379) + dyij*yin(343)
                                      zin(352) = zin(379) + dzij*zin(343)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  388

                                      ! ni =    3

                                      xin(388) = xin(415) + dxij*xin(379)
                                      yin(388) = yin(415) + dyij*yin(379)
                                      zin(388) = zin(415) + dzij*zin(379)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  424

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  325

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  427

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  418

                                      xin(427) = xin(427) + dxij*xin(418)
                                      yin(427) = yin(427) + dyij*yin(418)
                                      zin(427) = zin(427) + dzij*zin(418)

                                      ! i3 = i4 =  418
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  409

                                      xin(418) = xin(418) + dxij*xin(409)
                                      yin(418) = yin(418) + dyij*yin(409)
                                      zin(418) = zin(418) + dzij*zin(409)

                                      ! i3 = i4 =  409
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  400

                                      xin(409) = xin(409) + dxij*xin(400)
                                      yin(409) = yin(409) + dyij*yin(400)
                                      zin(409) = zin(409) + dzij*zin(400)

                                      ! i3 = i4 =  400
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  427

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  418

                                      xin(427) = xin(427) + dxij*xin(418)
                                      yin(427) = yin(427) + dyij*yin(418)
                                      zin(427) = zin(427) + dzij*zin(418)

                                      ! i3 = i4 =  418
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  409

                                      xin(418) = xin(418) + dxij*xin(409)
                                      yin(418) = yin(418) + dyij*yin(409)
                                      zin(418) = zin(418) + dzij*zin(409)

                                      ! i3 = i4 =  409
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  427

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  418

                                      xin(427) = xin(427) + dxij*xin(418)
                                      yin(427) = yin(427) + dyij*yin(418)
                                      zin(427) = zin(427) + dzij*zin(418)

                                      ! i3 = i4 =  418
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  301

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  301

                                      ! do ni = 1,    3

                                      xin(301) = xin(328) + dxij*xin(292)
                                      yin(301) = yin(328) + dyij*yin(292)
                                      zin(301) = zin(328) + dzij*zin(292)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  337

                                      ! ni =    2

                                      xin(337) = xin(364) + dxij*xin(328)
                                      yin(337) = yin(364) + dyij*yin(328)
                                      zin(337) = zin(364) + dzij*zin(328)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  373

                                      ! ni =    3

                                      xin(373) = xin(400) + dxij*xin(364)
                                      yin(373) = yin(400) + dyij*yin(364)
                                      zin(373) = zin(400) + dzij*zin(364)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  409

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  310

                                      ! nj =    2

                                      ! i4 = i3 =  310

                                      ! do ni = 1,    3

                                      xin(310) = xin(337) + dxij*xin(301)
                                      yin(310) = yin(337) + dyij*yin(301)
                                      zin(310) = zin(337) + dzij*zin(301)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  346

                                      ! ni =    2

                                      xin(346) = xin(373) + dxij*xin(337)
                                      yin(346) = yin(373) + dyij*yin(337)
                                      zin(346) = zin(373) + dzij*zin(337)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  382

                                      ! ni =    3

                                      xin(382) = xin(409) + dxij*xin(373)
                                      yin(382) = yin(409) + dyij*yin(373)
                                      zin(382) = zin(409) + dzij*zin(373)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  418

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  319

                                      ! nj =    3

                                      ! i4 = i3 =  319

                                      ! do ni = 1,    3

                                      xin(319) = xin(346) + dxij*xin(310)
                                      yin(319) = yin(346) + dyij*yin(310)
                                      zin(319) = zin(346) + dzij*zin(310)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  355

                                      ! ni =    2

                                      xin(355) = xin(382) + dxij*xin(346)
                                      yin(355) = yin(382) + dyij*yin(346)
                                      zin(355) = zin(382) + dzij*zin(346)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  391

                                      ! ni =    3

                                      xin(391) = xin(418) + dxij*xin(382)
                                      yin(391) = yin(418) + dyij*yin(382)
                                      zin(391) = zin(418) + dzij*zin(382)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  427

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  328

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  430

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  421

                                      xin(430) = xin(430) + dxij*xin(421)
                                      yin(430) = yin(430) + dyij*yin(421)
                                      zin(430) = zin(430) + dzij*zin(421)

                                      ! i3 = i4 =  421
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  412

                                      xin(421) = xin(421) + dxij*xin(412)
                                      yin(421) = yin(421) + dyij*yin(412)
                                      zin(421) = zin(421) + dzij*zin(412)

                                      ! i3 = i4 =  412
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  403

                                      xin(412) = xin(412) + dxij*xin(403)
                                      yin(412) = yin(412) + dyij*yin(403)
                                      zin(412) = zin(412) + dzij*zin(403)

                                      ! i3 = i4 =  403
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  430

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  421

                                      xin(430) = xin(430) + dxij*xin(421)
                                      yin(430) = yin(430) + dyij*yin(421)
                                      zin(430) = zin(430) + dzij*zin(421)

                                      ! i3 = i4 =  421
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  412

                                      xin(421) = xin(421) + dxij*xin(412)
                                      yin(421) = yin(421) + dyij*yin(412)
                                      zin(421) = zin(421) + dzij*zin(412)

                                      ! i3 = i4 =  412
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  430

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  421

                                      xin(430) = xin(430) + dxij*xin(421)
                                      yin(430) = yin(430) + dyij*yin(421)
                                      zin(430) = zin(430) + dzij*zin(421)

                                      ! i3 = i4 =  421
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  304

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  304

                                      ! do ni = 1,    3

                                      xin(304) = xin(331) + dxij*xin(295)
                                      yin(304) = yin(331) + dyij*yin(295)
                                      zin(304) = zin(331) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  340

                                      ! ni =    2

                                      xin(340) = xin(367) + dxij*xin(331)
                                      yin(340) = yin(367) + dyij*yin(331)
                                      zin(340) = zin(367) + dzij*zin(331)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  376

                                      ! ni =    3

                                      xin(376) = xin(403) + dxij*xin(367)
                                      yin(376) = yin(403) + dyij*yin(367)
                                      zin(376) = zin(403) + dzij*zin(367)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  412

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  313

                                      ! nj =    2

                                      ! i4 = i3 =  313

                                      ! do ni = 1,    3

                                      xin(313) = xin(340) + dxij*xin(304)
                                      yin(313) = yin(340) + dyij*yin(304)
                                      zin(313) = zin(340) + dzij*zin(304)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                      ! ni =    2

                                      xin(349) = xin(376) + dxij*xin(340)
                                      yin(349) = yin(376) + dyij*yin(340)
                                      zin(349) = zin(376) + dzij*zin(340)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  385

                                      ! ni =    3

                                      xin(385) = xin(412) + dxij*xin(376)
                                      yin(385) = yin(412) + dyij*yin(376)
                                      zin(385) = zin(412) + dzij*zin(376)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  421

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  322

                                      ! nj =    3

                                      ! i4 = i3 =  322

                                      ! do ni = 1,    3

                                      xin(322) = xin(349) + dxij*xin(313)
                                      yin(322) = yin(349) + dyij*yin(313)
                                      zin(322) = zin(349) + dzij*zin(313)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  358

                                      ! ni =    2

                                      xin(358) = xin(385) + dxij*xin(349)
                                      yin(358) = yin(385) + dyij*yin(349)
                                      zin(358) = zin(385) + dzij*zin(349)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  394

                                      ! ni =    3

                                      xin(394) = xin(421) + dxij*xin(385)
                                      yin(394) = yin(421) + dyij*yin(385)
                                      zin(394) = zin(421) + dzij*zin(385)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  430

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  331

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  431

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  422

                                      xin(431) = xin(431) + dxij*xin(422)
                                      yin(431) = yin(431) + dyij*yin(422)
                                      zin(431) = zin(431) + dzij*zin(422)

                                      ! i3 = i4 =  422
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  413

                                      xin(422) = xin(422) + dxij*xin(413)
                                      yin(422) = yin(422) + dyij*yin(413)
                                      zin(422) = zin(422) + dzij*zin(413)

                                      ! i3 = i4 =  413
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  404

                                      xin(413) = xin(413) + dxij*xin(404)
                                      yin(413) = yin(413) + dyij*yin(404)
                                      zin(413) = zin(413) + dzij*zin(404)

                                      ! i3 = i4 =  404
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  431

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  422

                                      xin(431) = xin(431) + dxij*xin(422)
                                      yin(431) = yin(431) + dyij*yin(422)
                                      zin(431) = zin(431) + dzij*zin(422)

                                      ! i3 = i4 =  422
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  413

                                      xin(422) = xin(422) + dxij*xin(413)
                                      yin(422) = yin(422) + dyij*yin(413)
                                      zin(422) = zin(422) + dzij*zin(413)

                                      ! i3 = i4 =  413
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  431

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  422

                                      xin(431) = xin(431) + dxij*xin(422)
                                      yin(431) = yin(431) + dyij*yin(422)
                                      zin(431) = zin(431) + dzij*zin(422)

                                      ! i3 = i4 =  422
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  305

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  305

                                      ! do ni = 1,    3

                                      xin(305) = xin(332) + dxij*xin(296)
                                      yin(305) = yin(332) + dyij*yin(296)
                                      zin(305) = zin(332) + dzij*zin(296)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  341

                                      ! ni =    2

                                      xin(341) = xin(368) + dxij*xin(332)
                                      yin(341) = yin(368) + dyij*yin(332)
                                      zin(341) = zin(368) + dzij*zin(332)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  377

                                      ! ni =    3

                                      xin(377) = xin(404) + dxij*xin(368)
                                      yin(377) = yin(404) + dyij*yin(368)
                                      zin(377) = zin(404) + dzij*zin(368)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  413

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  314

                                      ! nj =    2

                                      ! i4 = i3 =  314

                                      ! do ni = 1,    3

                                      xin(314) = xin(341) + dxij*xin(305)
                                      yin(314) = yin(341) + dyij*yin(305)
                                      zin(314) = zin(341) + dzij*zin(305)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  350

                                      ! ni =    2

                                      xin(350) = xin(377) + dxij*xin(341)
                                      yin(350) = yin(377) + dyij*yin(341)
                                      zin(350) = zin(377) + dzij*zin(341)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  386

                                      ! ni =    3

                                      xin(386) = xin(413) + dxij*xin(377)
                                      yin(386) = yin(413) + dyij*yin(377)
                                      zin(386) = zin(413) + dzij*zin(377)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  422

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  323

                                      ! nj =    3

                                      ! i4 = i3 =  323

                                      ! do ni = 1,    3

                                      xin(323) = xin(350) + dxij*xin(314)
                                      yin(323) = yin(350) + dyij*yin(314)
                                      zin(323) = zin(350) + dzij*zin(314)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  359

                                      ! ni =    2

                                      xin(359) = xin(386) + dxij*xin(350)
                                      yin(359) = yin(386) + dyij*yin(350)
                                      zin(359) = zin(386) + dzij*zin(350)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  395

                                      ! ni =    3

                                      xin(395) = xin(422) + dxij*xin(386)
                                      yin(395) = yin(422) + dyij*yin(386)
                                      zin(395) = zin(422) + dzij*zin(386)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  431

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  332

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  432

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  423

                                      xin(432) = xin(432) + dxij*xin(423)
                                      yin(432) = yin(432) + dyij*yin(423)
                                      zin(432) = zin(432) + dzij*zin(423)

                                      ! i3 = i4 =  423
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  414

                                      xin(423) = xin(423) + dxij*xin(414)
                                      yin(423) = yin(423) + dyij*yin(414)
                                      zin(423) = zin(423) + dzij*zin(414)

                                      ! i3 = i4 =  414
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  405

                                      xin(414) = xin(414) + dxij*xin(405)
                                      yin(414) = yin(414) + dyij*yin(405)
                                      zin(414) = zin(414) + dzij*zin(405)

                                      ! i3 = i4 =  405
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  432

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  423

                                      xin(432) = xin(432) + dxij*xin(423)
                                      yin(432) = yin(432) + dyij*yin(423)
                                      zin(432) = zin(432) + dzij*zin(423)

                                      ! i3 = i4 =  423
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  414

                                      xin(423) = xin(423) + dxij*xin(414)
                                      yin(423) = yin(423) + dyij*yin(414)
                                      zin(423) = zin(423) + dzij*zin(414)

                                      ! i3 = i4 =  414
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  432

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  423

                                      xin(432) = xin(432) + dxij*xin(423)
                                      yin(432) = yin(432) + dyij*yin(423)
                                      zin(432) = zin(432) + dzij*zin(423)

                                      ! i3 = i4 =  423
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  306

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  306

                                      ! do ni = 1,    3

                                      xin(306) = xin(333) + dxij*xin(297)
                                      yin(306) = yin(333) + dyij*yin(297)
                                      zin(306) = zin(333) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  342

                                      ! ni =    2

                                      xin(342) = xin(369) + dxij*xin(333)
                                      yin(342) = yin(369) + dyij*yin(333)
                                      zin(342) = zin(369) + dzij*zin(333)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  378

                                      ! ni =    3

                                      xin(378) = xin(405) + dxij*xin(369)
                                      yin(378) = yin(405) + dyij*yin(369)
                                      zin(378) = zin(405) + dzij*zin(369)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  414

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  315

                                      ! nj =    2

                                      ! i4 = i3 =  315

                                      ! do ni = 1,    3

                                      xin(315) = xin(342) + dxij*xin(306)
                                      yin(315) = yin(342) + dyij*yin(306)
                                      zin(315) = zin(342) + dzij*zin(306)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  351

                                      ! ni =    2

                                      xin(351) = xin(378) + dxij*xin(342)
                                      yin(351) = yin(378) + dyij*yin(342)
                                      zin(351) = zin(378) + dzij*zin(342)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  387

                                      ! ni =    3

                                      xin(387) = xin(414) + dxij*xin(378)
                                      yin(387) = yin(414) + dyij*yin(378)
                                      zin(387) = zin(414) + dzij*zin(378)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  423

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  324

                                      ! nj =    3

                                      ! i4 = i3 =  324

                                      ! do ni = 1,    3

                                      xin(324) = xin(351) + dxij*xin(315)
                                      yin(324) = yin(351) + dyij*yin(315)
                                      zin(324) = zin(351) + dzij*zin(315)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  360

                                      ! ni =    2

                                      xin(360) = xin(387) + dxij*xin(351)
                                      yin(360) = yin(387) + dyij*yin(351)
                                      zin(360) = zin(387) + dzij*zin(351)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  396

                                      ! ni =    3

                                      xin(396) = xin(423) + dxij*xin(387)
                                      yin(396) = yin(423) + dyij*yin(387)
                                      zin(396) = zin(423) + dzij*zin(387)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  432

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  333

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  289

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  297

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  296

                                      xin(297) = xin(297) + dxkl*xin(296)
                                      yin(297) = yin(297) + dykl*yin(296)
                                      zin(297) = zin(297) + dzkl*zin(296)

                                      ! i3 = i4 =  296
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  295

                                      xin(296) = xin(296) + dxkl*xin(295)
                                      yin(296) = yin(296) + dykl*yin(295)
                                      zin(296) = zin(296) + dzkl*zin(295)

                                      ! i3 = i4 =  295
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  297

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  296

                                      xin(297) = xin(297) + dxkl*xin(296)
                                      yin(297) = yin(297) + dykl*yin(296)
                                      zin(297) = zin(297) + dzkl*zin(296)

                                      ! i3 = i4 =  296
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    2

                                      xin(290) = xin(292) + dxkl*xin(289)
                                      yin(290) = yin(292) + dykl*yin(289)
                                      zin(290) = zin(292) + dzkl*zin(289)
                                      ! i4 = i4 + lang+1 =  293

                                      ! nk =    2

                                      xin(293) = xin(295) + dxkl*xin(292)
                                      yin(293) = yin(295) + dykl*yin(292)
                                      zin(293) = zin(295) + dzkl*zin(292)
                                      ! i4 = i4 + lang+1 =  296

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! i4 = i3 =  291

                                      ! do nk = 1,    2

                                      xin(291) = xin(293) + dxkl*xin(290)
                                      yin(291) = yin(293) + dykl*yin(290)
                                      zin(291) = zin(293) + dzkl*zin(290)
                                      ! i4 = i4 + lang+1 =  294

                                      ! nk =    2

                                      xin(294) = xin(296) + dxkl*xin(293)
                                      yin(294) = yin(296) + dykl*yin(293)
                                      zin(294) = zin(296) + dzkl*zin(293)
                                      ! i4 = i4 + lang+1 =  297

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  292

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  298

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  306

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  305

                                      xin(306) = xin(306) + dxkl*xin(305)
                                      yin(306) = yin(306) + dykl*yin(305)
                                      zin(306) = zin(306) + dzkl*zin(305)

                                      ! i3 = i4 =  305
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  304

                                      xin(305) = xin(305) + dxkl*xin(304)
                                      yin(305) = yin(305) + dykl*yin(304)
                                      zin(305) = zin(305) + dzkl*zin(304)

                                      ! i3 = i4 =  304
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  306

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  305

                                      xin(306) = xin(306) + dxkl*xin(305)
                                      yin(306) = yin(306) + dykl*yin(305)
                                      zin(306) = zin(306) + dzkl*zin(305)

                                      ! i3 = i4 =  305
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  299

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  299

                                      ! do nk = 1,    2

                                      xin(299) = xin(301) + dxkl*xin(298)
                                      yin(299) = yin(301) + dykl*yin(298)
                                      zin(299) = zin(301) + dzkl*zin(298)
                                      ! i4 = i4 + lang+1 =  302

                                      ! nk =    2

                                      xin(302) = xin(304) + dxkl*xin(301)
                                      yin(302) = yin(304) + dykl*yin(301)
                                      zin(302) = zin(304) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  305

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  300

                                      ! nl =    2

                                      ! i4 = i3 =  300

                                      ! do nk = 1,    2

                                      xin(300) = xin(302) + dxkl*xin(299)
                                      yin(300) = yin(302) + dykl*yin(299)
                                      zin(300) = zin(302) + dzkl*zin(299)
                                      ! i4 = i4 + lang+1 =  303

                                      ! nk =    2

                                      xin(303) = xin(305) + dxkl*xin(302)
                                      yin(303) = yin(305) + dykl*yin(302)
                                      zin(303) = zin(305) + dzkl*zin(302)
                                      ! i4 = i4 + lang+1 =  306

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  301

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  307

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  315

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  314

                                      xin(315) = xin(315) + dxkl*xin(314)
                                      yin(315) = yin(315) + dykl*yin(314)
                                      zin(315) = zin(315) + dzkl*zin(314)

                                      ! i3 = i4 =  314
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  313

                                      xin(314) = xin(314) + dxkl*xin(313)
                                      yin(314) = yin(314) + dykl*yin(313)
                                      zin(314) = zin(314) + dzkl*zin(313)

                                      ! i3 = i4 =  313
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  315

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  314

                                      xin(315) = xin(315) + dxkl*xin(314)
                                      yin(315) = yin(315) + dykl*yin(314)
                                      zin(315) = zin(315) + dzkl*zin(314)

                                      ! i3 = i4 =  314
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  308

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  308

                                      ! do nk = 1,    2

                                      xin(308) = xin(310) + dxkl*xin(307)
                                      yin(308) = yin(310) + dykl*yin(307)
                                      zin(308) = zin(310) + dzkl*zin(307)
                                      ! i4 = i4 + lang+1 =  311

                                      ! nk =    2

                                      xin(311) = xin(313) + dxkl*xin(310)
                                      yin(311) = yin(313) + dykl*yin(310)
                                      zin(311) = zin(313) + dzkl*zin(310)
                                      ! i4 = i4 + lang+1 =  314

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  309

                                      ! nl =    2

                                      ! i4 = i3 =  309

                                      ! do nk = 1,    2

                                      xin(309) = xin(311) + dxkl*xin(308)
                                      yin(309) = yin(311) + dykl*yin(308)
                                      zin(309) = zin(311) + dzkl*zin(308)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    2

                                      xin(312) = xin(314) + dxkl*xin(311)
                                      yin(312) = yin(314) + dykl*yin(311)
                                      zin(312) = zin(314) + dzkl*zin(311)
                                      ! i4 = i4 + lang+1 =  315

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  310

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  316

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  324

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  323

                                      xin(324) = xin(324) + dxkl*xin(323)
                                      yin(324) = yin(324) + dykl*yin(323)
                                      zin(324) = zin(324) + dzkl*zin(323)

                                      ! i3 = i4 =  323
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  322

                                      xin(323) = xin(323) + dxkl*xin(322)
                                      yin(323) = yin(323) + dykl*yin(322)
                                      zin(323) = zin(323) + dzkl*zin(322)

                                      ! i3 = i4 =  322
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  324

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  323

                                      xin(324) = xin(324) + dxkl*xin(323)
                                      yin(324) = yin(324) + dykl*yin(323)
                                      zin(324) = zin(324) + dzkl*zin(323)

                                      ! i3 = i4 =  323
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  317

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  317

                                      ! do nk = 1,    2

                                      xin(317) = xin(319) + dxkl*xin(316)
                                      yin(317) = yin(319) + dykl*yin(316)
                                      zin(317) = zin(319) + dzkl*zin(316)
                                      ! i4 = i4 + lang+1 =  320

                                      ! nk =    2

                                      xin(320) = xin(322) + dxkl*xin(319)
                                      yin(320) = yin(322) + dykl*yin(319)
                                      zin(320) = zin(322) + dzkl*zin(319)
                                      ! i4 = i4 + lang+1 =  323

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  318

                                      ! nl =    2

                                      ! i4 = i3 =  318

                                      ! do nk = 1,    2

                                      xin(318) = xin(320) + dxkl*xin(317)
                                      yin(318) = yin(320) + dykl*yin(317)
                                      zin(318) = zin(320) + dzkl*zin(317)
                                      ! i4 = i4 + lang+1 =  321

                                      ! nk =    2

                                      xin(321) = xin(323) + dxkl*xin(320)
                                      yin(321) = yin(323) + dykl*yin(320)
                                      zin(321) = zin(323) + dzkl*zin(320)
                                      ! i4 = i4 + lang+1 =  324

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  319

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  325

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  325

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  333

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  332

                                      xin(333) = xin(333) + dxkl*xin(332)
                                      yin(333) = yin(333) + dykl*yin(332)
                                      zin(333) = zin(333) + dzkl*zin(332)

                                      ! i3 = i4 =  332
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  331

                                      xin(332) = xin(332) + dxkl*xin(331)
                                      yin(332) = yin(332) + dykl*yin(331)
                                      zin(332) = zin(332) + dzkl*zin(331)

                                      ! i3 = i4 =  331
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  333

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  332

                                      xin(333) = xin(333) + dxkl*xin(332)
                                      yin(333) = yin(333) + dykl*yin(332)
                                      zin(333) = zin(333) + dzkl*zin(332)

                                      ! i3 = i4 =  332
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  326

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  326

                                      ! do nk = 1,    2

                                      xin(326) = xin(328) + dxkl*xin(325)
                                      yin(326) = yin(328) + dykl*yin(325)
                                      zin(326) = zin(328) + dzkl*zin(325)
                                      ! i4 = i4 + lang+1 =  329

                                      ! nk =    2

                                      xin(329) = xin(331) + dxkl*xin(328)
                                      yin(329) = yin(331) + dykl*yin(328)
                                      zin(329) = zin(331) + dzkl*zin(328)
                                      ! i4 = i4 + lang+1 =  332

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  327

                                      ! nl =    2

                                      ! i4 = i3 =  327

                                      ! do nk = 1,    2

                                      xin(327) = xin(329) + dxkl*xin(326)
                                      yin(327) = yin(329) + dykl*yin(326)
                                      zin(327) = zin(329) + dzkl*zin(326)
                                      ! i4 = i4 + lang+1 =  330

                                      ! nk =    2

                                      xin(330) = xin(332) + dxkl*xin(329)
                                      yin(330) = yin(332) + dykl*yin(329)
                                      zin(330) = zin(332) + dzkl*zin(329)
                                      ! i4 = i4 + lang+1 =  333

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  328

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  334

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  342

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  341

                                      xin(342) = xin(342) + dxkl*xin(341)
                                      yin(342) = yin(342) + dykl*yin(341)
                                      zin(342) = zin(342) + dzkl*zin(341)

                                      ! i3 = i4 =  341
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  340

                                      xin(341) = xin(341) + dxkl*xin(340)
                                      yin(341) = yin(341) + dykl*yin(340)
                                      zin(341) = zin(341) + dzkl*zin(340)

                                      ! i3 = i4 =  340
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  342

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  341

                                      xin(342) = xin(342) + dxkl*xin(341)
                                      yin(342) = yin(342) + dykl*yin(341)
                                      zin(342) = zin(342) + dzkl*zin(341)

                                      ! i3 = i4 =  341
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  335

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  335

                                      ! do nk = 1,    2

                                      xin(335) = xin(337) + dxkl*xin(334)
                                      yin(335) = yin(337) + dykl*yin(334)
                                      zin(335) = zin(337) + dzkl*zin(334)
                                      ! i4 = i4 + lang+1 =  338

                                      ! nk =    2

                                      xin(338) = xin(340) + dxkl*xin(337)
                                      yin(338) = yin(340) + dykl*yin(337)
                                      zin(338) = zin(340) + dzkl*zin(337)
                                      ! i4 = i4 + lang+1 =  341

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  336

                                      ! nl =    2

                                      ! i4 = i3 =  336

                                      ! do nk = 1,    2

                                      xin(336) = xin(338) + dxkl*xin(335)
                                      yin(336) = yin(338) + dykl*yin(335)
                                      zin(336) = zin(338) + dzkl*zin(335)
                                      ! i4 = i4 + lang+1 =  339

                                      ! nk =    2

                                      xin(339) = xin(341) + dxkl*xin(338)
                                      yin(339) = yin(341) + dykl*yin(338)
                                      zin(339) = zin(341) + dzkl*zin(338)
                                      ! i4 = i4 + lang+1 =  342

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  337

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  343

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  351

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  350

                                      xin(351) = xin(351) + dxkl*xin(350)
                                      yin(351) = yin(351) + dykl*yin(350)
                                      zin(351) = zin(351) + dzkl*zin(350)

                                      ! i3 = i4 =  350
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  349

                                      xin(350) = xin(350) + dxkl*xin(349)
                                      yin(350) = yin(350) + dykl*yin(349)
                                      zin(350) = zin(350) + dzkl*zin(349)

                                      ! i3 = i4 =  349
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  351

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  350

                                      xin(351) = xin(351) + dxkl*xin(350)
                                      yin(351) = yin(351) + dykl*yin(350)
                                      zin(351) = zin(351) + dzkl*zin(350)

                                      ! i3 = i4 =  350
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  344

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  344

                                      ! do nk = 1,    2

                                      xin(344) = xin(346) + dxkl*xin(343)
                                      yin(344) = yin(346) + dykl*yin(343)
                                      zin(344) = zin(346) + dzkl*zin(343)
                                      ! i4 = i4 + lang+1 =  347

                                      ! nk =    2

                                      xin(347) = xin(349) + dxkl*xin(346)
                                      yin(347) = yin(349) + dykl*yin(346)
                                      zin(347) = zin(349) + dzkl*zin(346)
                                      ! i4 = i4 + lang+1 =  350

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  345

                                      ! nl =    2

                                      ! i4 = i3 =  345

                                      ! do nk = 1,    2

                                      xin(345) = xin(347) + dxkl*xin(344)
                                      yin(345) = yin(347) + dykl*yin(344)
                                      zin(345) = zin(347) + dzkl*zin(344)
                                      ! i4 = i4 + lang+1 =  348

                                      ! nk =    2

                                      xin(348) = xin(350) + dxkl*xin(347)
                                      yin(348) = yin(350) + dykl*yin(347)
                                      zin(348) = zin(350) + dzkl*zin(347)
                                      ! i4 = i4 + lang+1 =  351

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  346

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  352

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  358

                                      xin(359) = xin(359) + dxkl*xin(358)
                                      yin(359) = yin(359) + dykl*yin(358)
                                      zin(359) = zin(359) + dzkl*zin(358)

                                      ! i3 = i4 =  358
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  353

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  353

                                      ! do nk = 1,    2

                                      xin(353) = xin(355) + dxkl*xin(352)
                                      yin(353) = yin(355) + dykl*yin(352)
                                      zin(353) = zin(355) + dzkl*zin(352)
                                      ! i4 = i4 + lang+1 =  356

                                      ! nk =    2

                                      xin(356) = xin(358) + dxkl*xin(355)
                                      yin(356) = yin(358) + dykl*yin(355)
                                      zin(356) = zin(358) + dzkl*zin(355)
                                      ! i4 = i4 + lang+1 =  359

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  354

                                      ! nl =    2

                                      ! i4 = i3 =  354

                                      ! do nk = 1,    2

                                      xin(354) = xin(356) + dxkl*xin(353)
                                      yin(354) = yin(356) + dykl*yin(353)
                                      zin(354) = zin(356) + dzkl*zin(353)
                                      ! i4 = i4 + lang+1 =  357

                                      ! nk =    2

                                      xin(357) = xin(359) + dxkl*xin(356)
                                      yin(357) = yin(359) + dykl*yin(356)
                                      zin(357) = zin(359) + dzkl*zin(356)
                                      ! i4 = i4 + lang+1 =  360

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  355

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  361

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  369

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  368

                                      xin(369) = xin(369) + dxkl*xin(368)
                                      yin(369) = yin(369) + dykl*yin(368)
                                      zin(369) = zin(369) + dzkl*zin(368)

                                      ! i3 = i4 =  368
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  367

                                      xin(368) = xin(368) + dxkl*xin(367)
                                      yin(368) = yin(368) + dykl*yin(367)
                                      zin(368) = zin(368) + dzkl*zin(367)

                                      ! i3 = i4 =  367
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  369

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  368

                                      xin(369) = xin(369) + dxkl*xin(368)
                                      yin(369) = yin(369) + dykl*yin(368)
                                      zin(369) = zin(369) + dzkl*zin(368)

                                      ! i3 = i4 =  368
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  362

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  362

                                      ! do nk = 1,    2

                                      xin(362) = xin(364) + dxkl*xin(361)
                                      yin(362) = yin(364) + dykl*yin(361)
                                      zin(362) = zin(364) + dzkl*zin(361)
                                      ! i4 = i4 + lang+1 =  365

                                      ! nk =    2

                                      xin(365) = xin(367) + dxkl*xin(364)
                                      yin(365) = yin(367) + dykl*yin(364)
                                      zin(365) = zin(367) + dzkl*zin(364)
                                      ! i4 = i4 + lang+1 =  368

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  363

                                      ! nl =    2

                                      ! i4 = i3 =  363

                                      ! do nk = 1,    2

                                      xin(363) = xin(365) + dxkl*xin(362)
                                      yin(363) = yin(365) + dykl*yin(362)
                                      zin(363) = zin(365) + dzkl*zin(362)
                                      ! i4 = i4 + lang+1 =  366

                                      ! nk =    2

                                      xin(366) = xin(368) + dxkl*xin(365)
                                      yin(366) = yin(368) + dykl*yin(365)
                                      zin(366) = zin(368) + dzkl*zin(365)
                                      ! i4 = i4 + lang+1 =  369

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  364

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  370

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  378

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  377

                                      xin(378) = xin(378) + dxkl*xin(377)
                                      yin(378) = yin(378) + dykl*yin(377)
                                      zin(378) = zin(378) + dzkl*zin(377)

                                      ! i3 = i4 =  377
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  376

                                      xin(377) = xin(377) + dxkl*xin(376)
                                      yin(377) = yin(377) + dykl*yin(376)
                                      zin(377) = zin(377) + dzkl*zin(376)

                                      ! i3 = i4 =  376
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  378

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  377

                                      xin(378) = xin(378) + dxkl*xin(377)
                                      yin(378) = yin(378) + dykl*yin(377)
                                      zin(378) = zin(378) + dzkl*zin(377)

                                      ! i3 = i4 =  377
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  371

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  371

                                      ! do nk = 1,    2

                                      xin(371) = xin(373) + dxkl*xin(370)
                                      yin(371) = yin(373) + dykl*yin(370)
                                      zin(371) = zin(373) + dzkl*zin(370)
                                      ! i4 = i4 + lang+1 =  374

                                      ! nk =    2

                                      xin(374) = xin(376) + dxkl*xin(373)
                                      yin(374) = yin(376) + dykl*yin(373)
                                      zin(374) = zin(376) + dzkl*zin(373)
                                      ! i4 = i4 + lang+1 =  377

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  372

                                      ! nl =    2

                                      ! i4 = i3 =  372

                                      ! do nk = 1,    2

                                      xin(372) = xin(374) + dxkl*xin(371)
                                      yin(372) = yin(374) + dykl*yin(371)
                                      zin(372) = zin(374) + dzkl*zin(371)
                                      ! i4 = i4 + lang+1 =  375

                                      ! nk =    2

                                      xin(375) = xin(377) + dxkl*xin(374)
                                      yin(375) = yin(377) + dykl*yin(374)
                                      zin(375) = zin(377) + dzkl*zin(374)
                                      ! i4 = i4 + lang+1 =  378

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  373

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  379

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  387

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  386

                                      xin(387) = xin(387) + dxkl*xin(386)
                                      yin(387) = yin(387) + dykl*yin(386)
                                      zin(387) = zin(387) + dzkl*zin(386)

                                      ! i3 = i4 =  386
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  385

                                      xin(386) = xin(386) + dxkl*xin(385)
                                      yin(386) = yin(386) + dykl*yin(385)
                                      zin(386) = zin(386) + dzkl*zin(385)

                                      ! i3 = i4 =  385
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  387

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  386

                                      xin(387) = xin(387) + dxkl*xin(386)
                                      yin(387) = yin(387) + dykl*yin(386)
                                      zin(387) = zin(387) + dzkl*zin(386)

                                      ! i3 = i4 =  386
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  380

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  380

                                      ! do nk = 1,    2

                                      xin(380) = xin(382) + dxkl*xin(379)
                                      yin(380) = yin(382) + dykl*yin(379)
                                      zin(380) = zin(382) + dzkl*zin(379)
                                      ! i4 = i4 + lang+1 =  383

                                      ! nk =    2

                                      xin(383) = xin(385) + dxkl*xin(382)
                                      yin(383) = yin(385) + dykl*yin(382)
                                      zin(383) = zin(385) + dzkl*zin(382)
                                      ! i4 = i4 + lang+1 =  386

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  381

                                      ! nl =    2

                                      ! i4 = i3 =  381

                                      ! do nk = 1,    2

                                      xin(381) = xin(383) + dxkl*xin(380)
                                      yin(381) = yin(383) + dykl*yin(380)
                                      zin(381) = zin(383) + dzkl*zin(380)
                                      ! i4 = i4 + lang+1 =  384

                                      ! nk =    2

                                      xin(384) = xin(386) + dxkl*xin(383)
                                      yin(384) = yin(386) + dykl*yin(383)
                                      zin(384) = zin(386) + dzkl*zin(383)
                                      ! i4 = i4 + lang+1 =  387

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  382

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  388

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  396

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  395

                                      xin(396) = xin(396) + dxkl*xin(395)
                                      yin(396) = yin(396) + dykl*yin(395)
                                      zin(396) = zin(396) + dzkl*zin(395)

                                      ! i3 = i4 =  395
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  394

                                      xin(395) = xin(395) + dxkl*xin(394)
                                      yin(395) = yin(395) + dykl*yin(394)
                                      zin(395) = zin(395) + dzkl*zin(394)

                                      ! i3 = i4 =  394
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  396

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  395

                                      xin(396) = xin(396) + dxkl*xin(395)
                                      yin(396) = yin(396) + dykl*yin(395)
                                      zin(396) = zin(396) + dzkl*zin(395)

                                      ! i3 = i4 =  395
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  389

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  389

                                      ! do nk = 1,    2

                                      xin(389) = xin(391) + dxkl*xin(388)
                                      yin(389) = yin(391) + dykl*yin(388)
                                      zin(389) = zin(391) + dzkl*zin(388)
                                      ! i4 = i4 + lang+1 =  392

                                      ! nk =    2

                                      xin(392) = xin(394) + dxkl*xin(391)
                                      yin(392) = yin(394) + dykl*yin(391)
                                      zin(392) = zin(394) + dzkl*zin(391)
                                      ! i4 = i4 + lang+1 =  395

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  390

                                      ! nl =    2

                                      ! i4 = i3 =  390

                                      ! do nk = 1,    2

                                      xin(390) = xin(392) + dxkl*xin(389)
                                      yin(390) = yin(392) + dykl*yin(389)
                                      zin(390) = zin(392) + dzkl*zin(389)
                                      ! i4 = i4 + lang+1 =  393

                                      ! nk =    2

                                      xin(393) = xin(395) + dxkl*xin(392)
                                      yin(393) = yin(395) + dykl*yin(392)
                                      zin(393) = zin(395) + dzkl*zin(392)
                                      ! i4 = i4 + lang+1 =  396

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  391

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  397

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  397

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  405

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  404

                                      xin(405) = xin(405) + dxkl*xin(404)
                                      yin(405) = yin(405) + dykl*yin(404)
                                      zin(405) = zin(405) + dzkl*zin(404)

                                      ! i3 = i4 =  404
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  403

                                      xin(404) = xin(404) + dxkl*xin(403)
                                      yin(404) = yin(404) + dykl*yin(403)
                                      zin(404) = zin(404) + dzkl*zin(403)

                                      ! i3 = i4 =  403
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  405

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  404

                                      xin(405) = xin(405) + dxkl*xin(404)
                                      yin(405) = yin(405) + dykl*yin(404)
                                      zin(405) = zin(405) + dzkl*zin(404)

                                      ! i3 = i4 =  404
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  398

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  398

                                      ! do nk = 1,    2

                                      xin(398) = xin(400) + dxkl*xin(397)
                                      yin(398) = yin(400) + dykl*yin(397)
                                      zin(398) = zin(400) + dzkl*zin(397)
                                      ! i4 = i4 + lang+1 =  401

                                      ! nk =    2

                                      xin(401) = xin(403) + dxkl*xin(400)
                                      yin(401) = yin(403) + dykl*yin(400)
                                      zin(401) = zin(403) + dzkl*zin(400)
                                      ! i4 = i4 + lang+1 =  404

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  399

                                      ! nl =    2

                                      ! i4 = i3 =  399

                                      ! do nk = 1,    2

                                      xin(399) = xin(401) + dxkl*xin(398)
                                      yin(399) = yin(401) + dykl*yin(398)
                                      zin(399) = zin(401) + dzkl*zin(398)
                                      ! i4 = i4 + lang+1 =  402

                                      ! nk =    2

                                      xin(402) = xin(404) + dxkl*xin(401)
                                      yin(402) = yin(404) + dykl*yin(401)
                                      zin(402) = zin(404) + dzkl*zin(401)
                                      ! i4 = i4 + lang+1 =  405

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  400

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  406

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  414

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  413

                                      xin(414) = xin(414) + dxkl*xin(413)
                                      yin(414) = yin(414) + dykl*yin(413)
                                      zin(414) = zin(414) + dzkl*zin(413)

                                      ! i3 = i4 =  413
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  412

                                      xin(413) = xin(413) + dxkl*xin(412)
                                      yin(413) = yin(413) + dykl*yin(412)
                                      zin(413) = zin(413) + dzkl*zin(412)

                                      ! i3 = i4 =  412
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  414

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  413

                                      xin(414) = xin(414) + dxkl*xin(413)
                                      yin(414) = yin(414) + dykl*yin(413)
                                      zin(414) = zin(414) + dzkl*zin(413)

                                      ! i3 = i4 =  413
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  407

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  407

                                      ! do nk = 1,    2

                                      xin(407) = xin(409) + dxkl*xin(406)
                                      yin(407) = yin(409) + dykl*yin(406)
                                      zin(407) = zin(409) + dzkl*zin(406)
                                      ! i4 = i4 + lang+1 =  410

                                      ! nk =    2

                                      xin(410) = xin(412) + dxkl*xin(409)
                                      yin(410) = yin(412) + dykl*yin(409)
                                      zin(410) = zin(412) + dzkl*zin(409)
                                      ! i4 = i4 + lang+1 =  413

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  408

                                      ! nl =    2

                                      ! i4 = i3 =  408

                                      ! do nk = 1,    2

                                      xin(408) = xin(410) + dxkl*xin(407)
                                      yin(408) = yin(410) + dykl*yin(407)
                                      zin(408) = zin(410) + dzkl*zin(407)
                                      ! i4 = i4 + lang+1 =  411

                                      ! nk =    2

                                      xin(411) = xin(413) + dxkl*xin(410)
                                      yin(411) = yin(413) + dykl*yin(410)
                                      zin(411) = zin(413) + dzkl*zin(410)
                                      ! i4 = i4 + lang+1 =  414

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  409

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  415

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  423

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  422

                                      xin(423) = xin(423) + dxkl*xin(422)
                                      yin(423) = yin(423) + dykl*yin(422)
                                      zin(423) = zin(423) + dzkl*zin(422)

                                      ! i3 = i4 =  422
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  421

                                      xin(422) = xin(422) + dxkl*xin(421)
                                      yin(422) = yin(422) + dykl*yin(421)
                                      zin(422) = zin(422) + dzkl*zin(421)

                                      ! i3 = i4 =  421
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  423

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  422

                                      xin(423) = xin(423) + dxkl*xin(422)
                                      yin(423) = yin(423) + dykl*yin(422)
                                      zin(423) = zin(423) + dzkl*zin(422)

                                      ! i3 = i4 =  422
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  416

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  416

                                      ! do nk = 1,    2

                                      xin(416) = xin(418) + dxkl*xin(415)
                                      yin(416) = yin(418) + dykl*yin(415)
                                      zin(416) = zin(418) + dzkl*zin(415)
                                      ! i4 = i4 + lang+1 =  419

                                      ! nk =    2

                                      xin(419) = xin(421) + dxkl*xin(418)
                                      yin(419) = yin(421) + dykl*yin(418)
                                      zin(419) = zin(421) + dzkl*zin(418)
                                      ! i4 = i4 + lang+1 =  422

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  417

                                      ! nl =    2

                                      ! i4 = i3 =  417

                                      ! do nk = 1,    2

                                      xin(417) = xin(419) + dxkl*xin(416)
                                      yin(417) = yin(419) + dykl*yin(416)
                                      zin(417) = zin(419) + dzkl*zin(416)
                                      ! i4 = i4 + lang+1 =  420

                                      ! nk =    2

                                      xin(420) = xin(422) + dxkl*xin(419)
                                      yin(420) = yin(422) + dykl*yin(419)
                                      zin(420) = zin(422) + dzkl*zin(419)
                                      ! i4 = i4 + lang+1 =  423

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  418

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  424

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  432

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  431

                                      xin(432) = xin(432) + dxkl*xin(431)
                                      yin(432) = yin(432) + dykl*yin(431)
                                      zin(432) = zin(432) + dzkl*zin(431)

                                      ! i3 = i4 =  431
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  430

                                      xin(431) = xin(431) + dxkl*xin(430)
                                      yin(431) = yin(431) + dykl*yin(430)
                                      zin(431) = zin(431) + dzkl*zin(430)

                                      ! i3 = i4 =  430
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  432

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  431

                                      xin(432) = xin(432) + dxkl*xin(431)
                                      yin(432) = yin(432) + dykl*yin(431)
                                      zin(432) = zin(432) + dzkl*zin(431)

                                      ! i3 = i4 =  431
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  425

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  425

                                      ! do nk = 1,    2

                                      xin(425) = xin(427) + dxkl*xin(424)
                                      yin(425) = yin(427) + dykl*yin(424)
                                      zin(425) = zin(427) + dzkl*zin(424)
                                      ! i4 = i4 + lang+1 =  428

                                      ! nk =    2

                                      xin(428) = xin(430) + dxkl*xin(427)
                                      yin(428) = yin(430) + dykl*yin(427)
                                      zin(428) = zin(430) + dzkl*zin(427)
                                      ! i4 = i4 + lang+1 =  431

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  426

                                      ! nl =    2

                                      ! i4 = i3 =  426

                                      ! do nk = 1,    2

                                      xin(426) = xin(428) + dxkl*xin(425)
                                      yin(426) = yin(428) + dykl*yin(425)
                                      zin(426) = zin(428) + dzkl*zin(425)
                                      ! i4 = i4 + lang+1 =  429

                                      ! nk =    2

                                      xin(429) = xin(431) + dxkl*xin(428)
                                      yin(429) = yin(431) + dykl*yin(428)
                                      zin(429) = zin(431) + dzkl*zin(428)
                                      ! i4 = i4 + lang+1 =  432

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  427

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  433

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  433

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  432

                                      u2 = roots(4)*rho
                                      f00 = expe*wghts(4)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  433

                                      xin(433) = 1.0_dp
                                      yin(433) = 1.0_dp
                                      zin(433) = f00

                                      ! i2 = in(2) =  469
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(469) = xc00
                                      yin(469) = yc00
                                      zin(469) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  436

                                      xin(436) = xcp00
                                      yin(436) = ycp00
                                      zin(436) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  472
                                      ! i2 =  469

                                      xin(472) = xcp00*xin(469) + cp10
                                      yin(472) = ycp00*yin(469) + cp10
                                      zin(472) = zcp00*zin(469) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  469

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  505
                                      ! i3 =  433
                                      ! i4 =  469

                                      xin(505) = c10*xin(433) + xc00*xin(469)
                                      yin(505) = c10*yin(433) + yc00*yin(469)
                                      zin(505) = c10*zin(433) + zc00*zin(469)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  508
                                      ! i5 =  505
                                      ! i4 =  469

                                      xin(508) = xcp00*xin(505) + cp10*xin(469)
                                      yin(508) = ycp00*yin(505) + cp10*yin(469)
                                      zin(508) = zcp00*zin(505) + cp10*zin(469)

                                      ! ------------------

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  505

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  541
                                      ! i3 =  469
                                      ! i4 =  505

                                      xin(541) = c10*xin(469) + xc00*xin(505)
                                      yin(541) = c10*yin(469) + yc00*yin(505)
                                      zin(541) = c10*zin(469) + zc00*zin(505)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  544
                                      ! i5 =  541
                                      ! i4 =  505

                                      xin(544) = xcp00*xin(541) + cp10*xin(505)
                                      yin(544) = ycp00*yin(541) + cp10*yin(505)
                                      zin(544) = zcp00*zin(541) + cp10*zin(505)

                                      ! ------------------

                                      ! i3 = i4 =  505
                                      ! i4 = i5 =  541

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  550
                                      ! i3 =  505
                                      ! i4 =  541

                                      xin(550) = c10*xin(505) + xc00*xin(541)
                                      yin(550) = c10*yin(505) + yc00*yin(541)
                                      zin(550) = c10*zin(505) + zc00*zin(541)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  553
                                      ! i5 =  550
                                      ! i4 =  541

                                      xin(553) = xcp00*xin(550) + cp10*xin(541)
                                      yin(553) = ycp00*yin(550) + cp10*yin(541)
                                      zin(553) = zcp00*zin(550) + cp10*zin(541)

                                      ! ------------------

                                      ! i3 = i4 =  541
                                      ! i4 = i5 =  550

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  559
                                      ! i3 =  541
                                      ! i4 =  550

                                      xin(559) = c10*xin(541) + xc00*xin(550)
                                      yin(559) = c10*yin(541) + yc00*yin(550)
                                      zin(559) = c10*zin(541) + zc00*zin(550)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  562
                                      ! i5 =  559
                                      ! i4 =  550

                                      xin(562) = xcp00*xin(559) + cp10*xin(550)
                                      yin(562) = ycp00*yin(559) + cp10*yin(550)
                                      zin(562) = zcp00*zin(559) + cp10*zin(550)

                                      ! ------------------

                                      ! i3 = i4 =  550
                                      ! i4 = i5 =  559

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  568
                                      ! i3 =  550
                                      ! i4 =  559

                                      xin(568) = c10*xin(550) + xc00*xin(559)
                                      yin(568) = c10*yin(550) + yc00*yin(559)
                                      zin(568) = c10*zin(550) + zc00*zin(559)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  571
                                      ! i5 =  568
                                      ! i4 =  559

                                      xin(571) = xcp00*xin(568) + cp10*xin(559)
                                      yin(571) = ycp00*yin(568) + cp10*yin(559)
                                      zin(571) = zcp00*zin(568) + cp10*zin(559)

                                      ! ------------------

                                      ! i3 = i4 =  559
                                      ! i4 = i5 =  568

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  433
                                      ! i4 = i1+k2 =  436

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  439
                                      ! i3 =  433
                                      ! i4 =  436

                                      xin(439) = cp01*xin(433) + xcp00*xin(436)
                                      yin(439) = cp01*yin(433) + ycp00*yin(436)
                                      zin(439) = cp01*zin(433) + zcp00*zin(436)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  475

                                      xin(475) = xc00*xin(439) + c01*xin(436)
                                      yin(475) = yc00*yin(439) + c01*yin(436)
                                      zin(475) = zc00*zin(439) + c01*zin(436)

                                      ! ------------------

                                      ! i3 = i4 =  436
                                      ! i4 = i5 =  439

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  440
                                      ! i3 =  436
                                      ! i4 =  439

                                      xin(440) = cp01*xin(436) + xcp00*xin(439)
                                      yin(440) = cp01*yin(436) + ycp00*yin(439)
                                      zin(440) = cp01*zin(436) + zcp00*zin(439)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  476

                                      xin(476) = xc00*xin(440) + c01*xin(439)
                                      yin(476) = yc00*yin(440) + c01*yin(439)
                                      zin(476) = zc00*zin(440) + c01*zin(439)

                                      ! ------------------

                                      ! i3 = i4 =  439
                                      ! i4 = i5 =  440

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  441
                                      ! i3 =  439
                                      ! i4 =  440

                                      xin(441) = cp01*xin(439) + xcp00*xin(440)
                                      yin(441) = cp01*yin(439) + ycp00*yin(440)
                                      zin(441) = cp01*zin(439) + zcp00*zin(440)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  477

                                      xin(477) = xc00*xin(441) + c01*xin(440)
                                      yin(477) = yc00*yin(441) + c01*yin(440)
                                      zin(477) = zc00*zin(441) + c01*zin(440)

                                      ! ------------------

                                      ! i3 = i4 =  440
                                      ! i4 = i5 =  441

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  469

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  505

                                      xin(511) = c10*xin(439) + xc00*xin(475) + c01*xin(472)
                                      yin(511) = c10*yin(439) + yc00*yin(475) + c01*yin(472)
                                      zin(511) = c10*zin(439) + zc00*zin(475) + c01*zin(472)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  505

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  541

                                      xin(547) = c10*xin(475) + xc00*xin(511) + c01*xin(508)
                                      yin(547) = c10*yin(475) + yc00*yin(511) + c01*yin(508)
                                      zin(547) = c10*zin(475) + zc00*zin(511) + c01*zin(508)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  505
                                      ! i4 = i5 =  541

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  550

                                      xin(556) = c10*xin(511) + xc00*xin(547) + c01*xin(544)
                                      yin(556) = c10*yin(511) + yc00*yin(547) + c01*yin(544)
                                      zin(556) = c10*zin(511) + zc00*zin(547) + c01*zin(544)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  541
                                      ! i4 = i5 =  550

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  559

                                      xin(565) = c10*xin(547) + xc00*xin(556) + c01*xin(553)
                                      yin(565) = c10*yin(547) + yc00*yin(556) + c01*yin(553)
                                      zin(565) = c10*zin(547) + zc00*zin(556) + c01*zin(553)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  550
                                      ! i4 = i5 =  559

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  568

                                      xin(574) = c10*xin(556) + xc00*xin(565) + c01*xin(562)
                                      yin(574) = c10*yin(556) + yc00*yin(565) + c01*yin(562)
                                      zin(574) = c10*zin(556) + zc00*zin(565) + c01*zin(562)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  559
                                      ! i4 = i5 =  568

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  469

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  505

                                      xin(512) = c10*xin(440) + xc00*xin(476) + c01*xin(475)
                                      yin(512) = c10*yin(440) + yc00*yin(476) + c01*yin(475)
                                      zin(512) = c10*zin(440) + zc00*zin(476) + c01*zin(475)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  505

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  541

                                      xin(548) = c10*xin(476) + xc00*xin(512) + c01*xin(511)
                                      yin(548) = c10*yin(476) + yc00*yin(512) + c01*yin(511)
                                      zin(548) = c10*zin(476) + zc00*zin(512) + c01*zin(511)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  505
                                      ! i4 = i5 =  541

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  550

                                      xin(557) = c10*xin(512) + xc00*xin(548) + c01*xin(547)
                                      yin(557) = c10*yin(512) + yc00*yin(548) + c01*yin(547)
                                      zin(557) = c10*zin(512) + zc00*zin(548) + c01*zin(547)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  541
                                      ! i4 = i5 =  550

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  559

                                      xin(566) = c10*xin(548) + xc00*xin(557) + c01*xin(556)
                                      yin(566) = c10*yin(548) + yc00*yin(557) + c01*yin(556)
                                      zin(566) = c10*zin(548) + zc00*zin(557) + c01*zin(556)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  550
                                      ! i4 = i5 =  559

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  568

                                      xin(575) = c10*xin(557) + xc00*xin(566) + c01*xin(565)
                                      yin(575) = c10*yin(557) + yc00*yin(566) + c01*yin(565)
                                      zin(575) = c10*zin(557) + zc00*zin(566) + c01*zin(565)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  559
                                      ! i4 = i5 =  568

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  469

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  505

                                      xin(513) = c10*xin(441) + xc00*xin(477) + c01*xin(476)
                                      yin(513) = c10*yin(441) + yc00*yin(477) + c01*yin(476)
                                      zin(513) = c10*zin(441) + zc00*zin(477) + c01*zin(476)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  505

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  541

                                      xin(549) = c10*xin(477) + xc00*xin(513) + c01*xin(512)
                                      yin(549) = c10*yin(477) + yc00*yin(513) + c01*yin(512)
                                      zin(549) = c10*zin(477) + zc00*zin(513) + c01*zin(512)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  505
                                      ! i4 = i5 =  541

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  550

                                      xin(558) = c10*xin(513) + xc00*xin(549) + c01*xin(548)
                                      yin(558) = c10*yin(513) + yc00*yin(549) + c01*yin(548)
                                      zin(558) = c10*zin(513) + zc00*zin(549) + c01*zin(548)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  541
                                      ! i4 = i5 =  550

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  559

                                      xin(567) = c10*xin(549) + xc00*xin(558) + c01*xin(557)
                                      yin(567) = c10*yin(549) + yc00*yin(558) + c01*yin(557)
                                      zin(567) = c10*zin(549) + zc00*zin(558) + c01*zin(557)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  550
                                      ! i4 = i5 =  559

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  568

                                      xin(576) = c10*xin(558) + xc00*xin(567) + c01*xin(566)
                                      yin(576) = c10*yin(558) + yc00*yin(567) + c01*yin(566)
                                      zin(576) = c10*zin(558) + zc00*zin(567) + c01*zin(566)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  559
                                      ! i4 = i5 =  568

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  568

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  568

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  559

                                      xin(568) = xin(568) + dxij*xin(559)
                                      yin(568) = yin(568) + dyij*yin(559)
                                      zin(568) = zin(568) + dzij*zin(559)

                                      ! i3 = i4 =  559
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  550

                                      xin(559) = xin(559) + dxij*xin(550)
                                      yin(559) = yin(559) + dyij*yin(550)
                                      zin(559) = zin(559) + dzij*zin(550)

                                      ! i3 = i4 =  550
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  541

                                      xin(550) = xin(550) + dxij*xin(541)
                                      yin(550) = yin(550) + dyij*yin(541)
                                      zin(550) = zin(550) + dzij*zin(541)

                                      ! i3 = i4 =  541
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  568

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  559

                                      xin(568) = xin(568) + dxij*xin(559)
                                      yin(568) = yin(568) + dyij*yin(559)
                                      zin(568) = zin(568) + dzij*zin(559)

                                      ! i3 = i4 =  559
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  550

                                      xin(559) = xin(559) + dxij*xin(550)
                                      yin(559) = yin(559) + dyij*yin(550)
                                      zin(559) = zin(559) + dzij*zin(550)

                                      ! i3 = i4 =  550
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  568

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  559

                                      xin(568) = xin(568) + dxij*xin(559)
                                      yin(568) = yin(568) + dyij*yin(559)
                                      zin(568) = zin(568) + dzij*zin(559)

                                      ! i3 = i4 =  559
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  442

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  442

                                      ! do ni = 1,    3

                                      xin(442) = xin(469) + dxij*xin(433)
                                      yin(442) = yin(469) + dyij*yin(433)
                                      zin(442) = zin(469) + dzij*zin(433)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  478

                                      ! ni =    2

                                      xin(478) = xin(505) + dxij*xin(469)
                                      yin(478) = yin(505) + dyij*yin(469)
                                      zin(478) = zin(505) + dzij*zin(469)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  514

                                      ! ni =    3

                                      xin(514) = xin(541) + dxij*xin(505)
                                      yin(514) = yin(541) + dyij*yin(505)
                                      zin(514) = zin(541) + dzij*zin(505)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  550

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  451

                                      ! nj =    2

                                      ! i4 = i3 =  451

                                      ! do ni = 1,    3

                                      xin(451) = xin(478) + dxij*xin(442)
                                      yin(451) = yin(478) + dyij*yin(442)
                                      zin(451) = zin(478) + dzij*zin(442)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  487

                                      ! ni =    2

                                      xin(487) = xin(514) + dxij*xin(478)
                                      yin(487) = yin(514) + dyij*yin(478)
                                      zin(487) = zin(514) + dzij*zin(478)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  523

                                      ! ni =    3

                                      xin(523) = xin(550) + dxij*xin(514)
                                      yin(523) = yin(550) + dyij*yin(514)
                                      zin(523) = zin(550) + dzij*zin(514)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  559

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  460

                                      ! nj =    3

                                      ! i4 = i3 =  460

                                      ! do ni = 1,    3

                                      xin(460) = xin(487) + dxij*xin(451)
                                      yin(460) = yin(487) + dyij*yin(451)
                                      zin(460) = zin(487) + dzij*zin(451)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  496

                                      ! ni =    2

                                      xin(496) = xin(523) + dxij*xin(487)
                                      yin(496) = yin(523) + dyij*yin(487)
                                      zin(496) = zin(523) + dzij*zin(487)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  532

                                      ! ni =    3

                                      xin(532) = xin(559) + dxij*xin(523)
                                      yin(532) = yin(559) + dyij*yin(523)
                                      zin(532) = zin(559) + dzij*zin(523)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  568

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  469

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  571

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  562

                                      xin(571) = xin(571) + dxij*xin(562)
                                      yin(571) = yin(571) + dyij*yin(562)
                                      zin(571) = zin(571) + dzij*zin(562)

                                      ! i3 = i4 =  562
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  553

                                      xin(562) = xin(562) + dxij*xin(553)
                                      yin(562) = yin(562) + dyij*yin(553)
                                      zin(562) = zin(562) + dzij*zin(553)

                                      ! i3 = i4 =  553
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  544

                                      xin(553) = xin(553) + dxij*xin(544)
                                      yin(553) = yin(553) + dyij*yin(544)
                                      zin(553) = zin(553) + dzij*zin(544)

                                      ! i3 = i4 =  544
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  571

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  562

                                      xin(571) = xin(571) + dxij*xin(562)
                                      yin(571) = yin(571) + dyij*yin(562)
                                      zin(571) = zin(571) + dzij*zin(562)

                                      ! i3 = i4 =  562
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  553

                                      xin(562) = xin(562) + dxij*xin(553)
                                      yin(562) = yin(562) + dyij*yin(553)
                                      zin(562) = zin(562) + dzij*zin(553)

                                      ! i3 = i4 =  553
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  571

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  562

                                      xin(571) = xin(571) + dxij*xin(562)
                                      yin(571) = yin(571) + dyij*yin(562)
                                      zin(571) = zin(571) + dzij*zin(562)

                                      ! i3 = i4 =  562
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  445

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  445

                                      ! do ni = 1,    3

                                      xin(445) = xin(472) + dxij*xin(436)
                                      yin(445) = yin(472) + dyij*yin(436)
                                      zin(445) = zin(472) + dzij*zin(436)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  481

                                      ! ni =    2

                                      xin(481) = xin(508) + dxij*xin(472)
                                      yin(481) = yin(508) + dyij*yin(472)
                                      zin(481) = zin(508) + dzij*zin(472)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  517

                                      ! ni =    3

                                      xin(517) = xin(544) + dxij*xin(508)
                                      yin(517) = yin(544) + dyij*yin(508)
                                      zin(517) = zin(544) + dzij*zin(508)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  553

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  454

                                      ! nj =    2

                                      ! i4 = i3 =  454

                                      ! do ni = 1,    3

                                      xin(454) = xin(481) + dxij*xin(445)
                                      yin(454) = yin(481) + dyij*yin(445)
                                      zin(454) = zin(481) + dzij*zin(445)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  490

                                      ! ni =    2

                                      xin(490) = xin(517) + dxij*xin(481)
                                      yin(490) = yin(517) + dyij*yin(481)
                                      zin(490) = zin(517) + dzij*zin(481)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  526

                                      ! ni =    3

                                      xin(526) = xin(553) + dxij*xin(517)
                                      yin(526) = yin(553) + dyij*yin(517)
                                      zin(526) = zin(553) + dzij*zin(517)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  562

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  463

                                      ! nj =    3

                                      ! i4 = i3 =  463

                                      ! do ni = 1,    3

                                      xin(463) = xin(490) + dxij*xin(454)
                                      yin(463) = yin(490) + dyij*yin(454)
                                      zin(463) = zin(490) + dzij*zin(454)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  499

                                      ! ni =    2

                                      xin(499) = xin(526) + dxij*xin(490)
                                      yin(499) = yin(526) + dyij*yin(490)
                                      zin(499) = zin(526) + dzij*zin(490)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  535

                                      ! ni =    3

                                      xin(535) = xin(562) + dxij*xin(526)
                                      yin(535) = yin(562) + dyij*yin(526)
                                      zin(535) = zin(562) + dzij*zin(526)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  571

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  472

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  574

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  565

                                      xin(574) = xin(574) + dxij*xin(565)
                                      yin(574) = yin(574) + dyij*yin(565)
                                      zin(574) = zin(574) + dzij*zin(565)

                                      ! i3 = i4 =  565
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  556

                                      xin(565) = xin(565) + dxij*xin(556)
                                      yin(565) = yin(565) + dyij*yin(556)
                                      zin(565) = zin(565) + dzij*zin(556)

                                      ! i3 = i4 =  556
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  547

                                      xin(556) = xin(556) + dxij*xin(547)
                                      yin(556) = yin(556) + dyij*yin(547)
                                      zin(556) = zin(556) + dzij*zin(547)

                                      ! i3 = i4 =  547
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  574

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  565

                                      xin(574) = xin(574) + dxij*xin(565)
                                      yin(574) = yin(574) + dyij*yin(565)
                                      zin(574) = zin(574) + dzij*zin(565)

                                      ! i3 = i4 =  565
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  556

                                      xin(565) = xin(565) + dxij*xin(556)
                                      yin(565) = yin(565) + dyij*yin(556)
                                      zin(565) = zin(565) + dzij*zin(556)

                                      ! i3 = i4 =  556
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  574

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  565

                                      xin(574) = xin(574) + dxij*xin(565)
                                      yin(574) = yin(574) + dyij*yin(565)
                                      zin(574) = zin(574) + dzij*zin(565)

                                      ! i3 = i4 =  565
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  448

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  448

                                      ! do ni = 1,    3

                                      xin(448) = xin(475) + dxij*xin(439)
                                      yin(448) = yin(475) + dyij*yin(439)
                                      zin(448) = zin(475) + dzij*zin(439)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  484

                                      ! ni =    2

                                      xin(484) = xin(511) + dxij*xin(475)
                                      yin(484) = yin(511) + dyij*yin(475)
                                      zin(484) = zin(511) + dzij*zin(475)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  520

                                      ! ni =    3

                                      xin(520) = xin(547) + dxij*xin(511)
                                      yin(520) = yin(547) + dyij*yin(511)
                                      zin(520) = zin(547) + dzij*zin(511)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  556

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  457

                                      ! nj =    2

                                      ! i4 = i3 =  457

                                      ! do ni = 1,    3

                                      xin(457) = xin(484) + dxij*xin(448)
                                      yin(457) = yin(484) + dyij*yin(448)
                                      zin(457) = zin(484) + dzij*zin(448)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  493

                                      ! ni =    2

                                      xin(493) = xin(520) + dxij*xin(484)
                                      yin(493) = yin(520) + dyij*yin(484)
                                      zin(493) = zin(520) + dzij*zin(484)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  529

                                      ! ni =    3

                                      xin(529) = xin(556) + dxij*xin(520)
                                      yin(529) = yin(556) + dyij*yin(520)
                                      zin(529) = zin(556) + dzij*zin(520)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  565

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  466

                                      ! nj =    3

                                      ! i4 = i3 =  466

                                      ! do ni = 1,    3

                                      xin(466) = xin(493) + dxij*xin(457)
                                      yin(466) = yin(493) + dyij*yin(457)
                                      zin(466) = zin(493) + dzij*zin(457)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  502

                                      ! ni =    2

                                      xin(502) = xin(529) + dxij*xin(493)
                                      yin(502) = yin(529) + dyij*yin(493)
                                      zin(502) = zin(529) + dzij*zin(493)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  538

                                      ! ni =    3

                                      xin(538) = xin(565) + dxij*xin(529)
                                      yin(538) = yin(565) + dyij*yin(529)
                                      zin(538) = zin(565) + dzij*zin(529)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  574

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  475

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  575

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  566

                                      xin(575) = xin(575) + dxij*xin(566)
                                      yin(575) = yin(575) + dyij*yin(566)
                                      zin(575) = zin(575) + dzij*zin(566)

                                      ! i3 = i4 =  566
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  557

                                      xin(566) = xin(566) + dxij*xin(557)
                                      yin(566) = yin(566) + dyij*yin(557)
                                      zin(566) = zin(566) + dzij*zin(557)

                                      ! i3 = i4 =  557
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  548

                                      xin(557) = xin(557) + dxij*xin(548)
                                      yin(557) = yin(557) + dyij*yin(548)
                                      zin(557) = zin(557) + dzij*zin(548)

                                      ! i3 = i4 =  548
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  575

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  566

                                      xin(575) = xin(575) + dxij*xin(566)
                                      yin(575) = yin(575) + dyij*yin(566)
                                      zin(575) = zin(575) + dzij*zin(566)

                                      ! i3 = i4 =  566
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  557

                                      xin(566) = xin(566) + dxij*xin(557)
                                      yin(566) = yin(566) + dyij*yin(557)
                                      zin(566) = zin(566) + dzij*zin(557)

                                      ! i3 = i4 =  557
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  575

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  566

                                      xin(575) = xin(575) + dxij*xin(566)
                                      yin(575) = yin(575) + dyij*yin(566)
                                      zin(575) = zin(575) + dzij*zin(566)

                                      ! i3 = i4 =  566
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  449

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  449

                                      ! do ni = 1,    3

                                      xin(449) = xin(476) + dxij*xin(440)
                                      yin(449) = yin(476) + dyij*yin(440)
                                      zin(449) = zin(476) + dzij*zin(440)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  485

                                      ! ni =    2

                                      xin(485) = xin(512) + dxij*xin(476)
                                      yin(485) = yin(512) + dyij*yin(476)
                                      zin(485) = zin(512) + dzij*zin(476)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  521

                                      ! ni =    3

                                      xin(521) = xin(548) + dxij*xin(512)
                                      yin(521) = yin(548) + dyij*yin(512)
                                      zin(521) = zin(548) + dzij*zin(512)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  557

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  458

                                      ! nj =    2

                                      ! i4 = i3 =  458

                                      ! do ni = 1,    3

                                      xin(458) = xin(485) + dxij*xin(449)
                                      yin(458) = yin(485) + dyij*yin(449)
                                      zin(458) = zin(485) + dzij*zin(449)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  494

                                      ! ni =    2

                                      xin(494) = xin(521) + dxij*xin(485)
                                      yin(494) = yin(521) + dyij*yin(485)
                                      zin(494) = zin(521) + dzij*zin(485)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  530

                                      ! ni =    3

                                      xin(530) = xin(557) + dxij*xin(521)
                                      yin(530) = yin(557) + dyij*yin(521)
                                      zin(530) = zin(557) + dzij*zin(521)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  566

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  467

                                      ! nj =    3

                                      ! i4 = i3 =  467

                                      ! do ni = 1,    3

                                      xin(467) = xin(494) + dxij*xin(458)
                                      yin(467) = yin(494) + dyij*yin(458)
                                      zin(467) = zin(494) + dzij*zin(458)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  503

                                      ! ni =    2

                                      xin(503) = xin(530) + dxij*xin(494)
                                      yin(503) = yin(530) + dyij*yin(494)
                                      zin(503) = zin(530) + dzij*zin(494)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  539

                                      ! ni =    3

                                      xin(539) = xin(566) + dxij*xin(530)
                                      yin(539) = yin(566) + dyij*yin(530)
                                      zin(539) = zin(566) + dzij*zin(530)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  575

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  476

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  576

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  567

                                      xin(576) = xin(576) + dxij*xin(567)
                                      yin(576) = yin(576) + dyij*yin(567)
                                      zin(576) = zin(576) + dzij*zin(567)

                                      ! i3 = i4 =  567
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  558

                                      xin(567) = xin(567) + dxij*xin(558)
                                      yin(567) = yin(567) + dyij*yin(558)
                                      zin(567) = zin(567) + dzij*zin(558)

                                      ! i3 = i4 =  558
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  549

                                      xin(558) = xin(558) + dxij*xin(549)
                                      yin(558) = yin(558) + dyij*yin(549)
                                      zin(558) = zin(558) + dzij*zin(549)

                                      ! i3 = i4 =  549
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  576

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  567

                                      xin(576) = xin(576) + dxij*xin(567)
                                      yin(576) = yin(576) + dyij*yin(567)
                                      zin(576) = zin(576) + dzij*zin(567)

                                      ! i3 = i4 =  567
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  558

                                      xin(567) = xin(567) + dxij*xin(558)
                                      yin(567) = yin(567) + dyij*yin(558)
                                      zin(567) = zin(567) + dzij*zin(558)

                                      ! i3 = i4 =  558
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  576

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  567

                                      xin(576) = xin(576) + dxij*xin(567)
                                      yin(576) = yin(576) + dyij*yin(567)
                                      zin(576) = zin(576) + dzij*zin(567)

                                      ! i3 = i4 =  567
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  450

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  450

                                      ! do ni = 1,    3

                                      xin(450) = xin(477) + dxij*xin(441)
                                      yin(450) = yin(477) + dyij*yin(441)
                                      zin(450) = zin(477) + dzij*zin(441)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  486

                                      ! ni =    2

                                      xin(486) = xin(513) + dxij*xin(477)
                                      yin(486) = yin(513) + dyij*yin(477)
                                      zin(486) = zin(513) + dzij*zin(477)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  522

                                      ! ni =    3

                                      xin(522) = xin(549) + dxij*xin(513)
                                      yin(522) = yin(549) + dyij*yin(513)
                                      zin(522) = zin(549) + dzij*zin(513)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  558

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  459

                                      ! nj =    2

                                      ! i4 = i3 =  459

                                      ! do ni = 1,    3

                                      xin(459) = xin(486) + dxij*xin(450)
                                      yin(459) = yin(486) + dyij*yin(450)
                                      zin(459) = zin(486) + dzij*zin(450)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  495

                                      ! ni =    2

                                      xin(495) = xin(522) + dxij*xin(486)
                                      yin(495) = yin(522) + dyij*yin(486)
                                      zin(495) = zin(522) + dzij*zin(486)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  531

                                      ! ni =    3

                                      xin(531) = xin(558) + dxij*xin(522)
                                      yin(531) = yin(558) + dyij*yin(522)
                                      zin(531) = zin(558) + dzij*zin(522)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  567

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  468

                                      ! nj =    3

                                      ! i4 = i3 =  468

                                      ! do ni = 1,    3

                                      xin(468) = xin(495) + dxij*xin(459)
                                      yin(468) = yin(495) + dyij*yin(459)
                                      zin(468) = zin(495) + dzij*zin(459)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  504

                                      ! ni =    2

                                      xin(504) = xin(531) + dxij*xin(495)
                                      yin(504) = yin(531) + dyij*yin(495)
                                      zin(504) = zin(531) + dzij*zin(495)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  540

                                      ! ni =    3

                                      xin(540) = xin(567) + dxij*xin(531)
                                      yin(540) = yin(567) + dyij*yin(531)
                                      zin(540) = zin(567) + dzij*zin(531)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  576

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  477

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  433

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  441

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  440

                                      xin(441) = xin(441) + dxkl*xin(440)
                                      yin(441) = yin(441) + dykl*yin(440)
                                      zin(441) = zin(441) + dzkl*zin(440)

                                      ! i3 = i4 =  440
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  439

                                      xin(440) = xin(440) + dxkl*xin(439)
                                      yin(440) = yin(440) + dykl*yin(439)
                                      zin(440) = zin(440) + dzkl*zin(439)

                                      ! i3 = i4 =  439
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  441

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  440

                                      xin(441) = xin(441) + dxkl*xin(440)
                                      yin(441) = yin(441) + dykl*yin(440)
                                      zin(441) = zin(441) + dzkl*zin(440)

                                      ! i3 = i4 =  440
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  434

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  434

                                      ! do nk = 1,    2

                                      xin(434) = xin(436) + dxkl*xin(433)
                                      yin(434) = yin(436) + dykl*yin(433)
                                      zin(434) = zin(436) + dzkl*zin(433)
                                      ! i4 = i4 + lang+1 =  437

                                      ! nk =    2

                                      xin(437) = xin(439) + dxkl*xin(436)
                                      yin(437) = yin(439) + dykl*yin(436)
                                      zin(437) = zin(439) + dzkl*zin(436)
                                      ! i4 = i4 + lang+1 =  440

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  435

                                      ! nl =    2

                                      ! i4 = i3 =  435

                                      ! do nk = 1,    2

                                      xin(435) = xin(437) + dxkl*xin(434)
                                      yin(435) = yin(437) + dykl*yin(434)
                                      zin(435) = zin(437) + dzkl*zin(434)
                                      ! i4 = i4 + lang+1 =  438

                                      ! nk =    2

                                      xin(438) = xin(440) + dxkl*xin(437)
                                      yin(438) = yin(440) + dykl*yin(437)
                                      zin(438) = zin(440) + dzkl*zin(437)
                                      ! i4 = i4 + lang+1 =  441

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  436

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  442

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  450

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  449

                                      xin(450) = xin(450) + dxkl*xin(449)
                                      yin(450) = yin(450) + dykl*yin(449)
                                      zin(450) = zin(450) + dzkl*zin(449)

                                      ! i3 = i4 =  449
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  448

                                      xin(449) = xin(449) + dxkl*xin(448)
                                      yin(449) = yin(449) + dykl*yin(448)
                                      zin(449) = zin(449) + dzkl*zin(448)

                                      ! i3 = i4 =  448
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  450

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  449

                                      xin(450) = xin(450) + dxkl*xin(449)
                                      yin(450) = yin(450) + dykl*yin(449)
                                      zin(450) = zin(450) + dzkl*zin(449)

                                      ! i3 = i4 =  449
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  443

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  443

                                      ! do nk = 1,    2

                                      xin(443) = xin(445) + dxkl*xin(442)
                                      yin(443) = yin(445) + dykl*yin(442)
                                      zin(443) = zin(445) + dzkl*zin(442)
                                      ! i4 = i4 + lang+1 =  446

                                      ! nk =    2

                                      xin(446) = xin(448) + dxkl*xin(445)
                                      yin(446) = yin(448) + dykl*yin(445)
                                      zin(446) = zin(448) + dzkl*zin(445)
                                      ! i4 = i4 + lang+1 =  449

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  444

                                      ! nl =    2

                                      ! i4 = i3 =  444

                                      ! do nk = 1,    2

                                      xin(444) = xin(446) + dxkl*xin(443)
                                      yin(444) = yin(446) + dykl*yin(443)
                                      zin(444) = zin(446) + dzkl*zin(443)
                                      ! i4 = i4 + lang+1 =  447

                                      ! nk =    2

                                      xin(447) = xin(449) + dxkl*xin(446)
                                      yin(447) = yin(449) + dykl*yin(446)
                                      zin(447) = zin(449) + dzkl*zin(446)
                                      ! i4 = i4 + lang+1 =  450

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  445

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  451

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  459

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  458

                                      xin(459) = xin(459) + dxkl*xin(458)
                                      yin(459) = yin(459) + dykl*yin(458)
                                      zin(459) = zin(459) + dzkl*zin(458)

                                      ! i3 = i4 =  458
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  457

                                      xin(458) = xin(458) + dxkl*xin(457)
                                      yin(458) = yin(458) + dykl*yin(457)
                                      zin(458) = zin(458) + dzkl*zin(457)

                                      ! i3 = i4 =  457
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  459

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  458

                                      xin(459) = xin(459) + dxkl*xin(458)
                                      yin(459) = yin(459) + dykl*yin(458)
                                      zin(459) = zin(459) + dzkl*zin(458)

                                      ! i3 = i4 =  458
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  452

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  452

                                      ! do nk = 1,    2

                                      xin(452) = xin(454) + dxkl*xin(451)
                                      yin(452) = yin(454) + dykl*yin(451)
                                      zin(452) = zin(454) + dzkl*zin(451)
                                      ! i4 = i4 + lang+1 =  455

                                      ! nk =    2

                                      xin(455) = xin(457) + dxkl*xin(454)
                                      yin(455) = yin(457) + dykl*yin(454)
                                      zin(455) = zin(457) + dzkl*zin(454)
                                      ! i4 = i4 + lang+1 =  458

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  453

                                      ! nl =    2

                                      ! i4 = i3 =  453

                                      ! do nk = 1,    2

                                      xin(453) = xin(455) + dxkl*xin(452)
                                      yin(453) = yin(455) + dykl*yin(452)
                                      zin(453) = zin(455) + dzkl*zin(452)
                                      ! i4 = i4 + lang+1 =  456

                                      ! nk =    2

                                      xin(456) = xin(458) + dxkl*xin(455)
                                      yin(456) = yin(458) + dykl*yin(455)
                                      zin(456) = zin(458) + dzkl*zin(455)
                                      ! i4 = i4 + lang+1 =  459

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  454

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  460

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  468

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  467

                                      xin(468) = xin(468) + dxkl*xin(467)
                                      yin(468) = yin(468) + dykl*yin(467)
                                      zin(468) = zin(468) + dzkl*zin(467)

                                      ! i3 = i4 =  467
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  466

                                      xin(467) = xin(467) + dxkl*xin(466)
                                      yin(467) = yin(467) + dykl*yin(466)
                                      zin(467) = zin(467) + dzkl*zin(466)

                                      ! i3 = i4 =  466
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  468

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  467

                                      xin(468) = xin(468) + dxkl*xin(467)
                                      yin(468) = yin(468) + dykl*yin(467)
                                      zin(468) = zin(468) + dzkl*zin(467)

                                      ! i3 = i4 =  467
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  461

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  461

                                      ! do nk = 1,    2

                                      xin(461) = xin(463) + dxkl*xin(460)
                                      yin(461) = yin(463) + dykl*yin(460)
                                      zin(461) = zin(463) + dzkl*zin(460)
                                      ! i4 = i4 + lang+1 =  464

                                      ! nk =    2

                                      xin(464) = xin(466) + dxkl*xin(463)
                                      yin(464) = yin(466) + dykl*yin(463)
                                      zin(464) = zin(466) + dzkl*zin(463)
                                      ! i4 = i4 + lang+1 =  467

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  462

                                      ! nl =    2

                                      ! i4 = i3 =  462

                                      ! do nk = 1,    2

                                      xin(462) = xin(464) + dxkl*xin(461)
                                      yin(462) = yin(464) + dykl*yin(461)
                                      zin(462) = zin(464) + dzkl*zin(461)
                                      ! i4 = i4 + lang+1 =  465

                                      ! nk =    2

                                      xin(465) = xin(467) + dxkl*xin(464)
                                      yin(465) = yin(467) + dykl*yin(464)
                                      zin(465) = zin(467) + dzkl*zin(464)
                                      ! i4 = i4 + lang+1 =  468

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  463

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  469

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  469

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  477

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  476

                                      xin(477) = xin(477) + dxkl*xin(476)
                                      yin(477) = yin(477) + dykl*yin(476)
                                      zin(477) = zin(477) + dzkl*zin(476)

                                      ! i3 = i4 =  476
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  475

                                      xin(476) = xin(476) + dxkl*xin(475)
                                      yin(476) = yin(476) + dykl*yin(475)
                                      zin(476) = zin(476) + dzkl*zin(475)

                                      ! i3 = i4 =  475
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  477

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  476

                                      xin(477) = xin(477) + dxkl*xin(476)
                                      yin(477) = yin(477) + dykl*yin(476)
                                      zin(477) = zin(477) + dzkl*zin(476)

                                      ! i3 = i4 =  476
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  470

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  470

                                      ! do nk = 1,    2

                                      xin(470) = xin(472) + dxkl*xin(469)
                                      yin(470) = yin(472) + dykl*yin(469)
                                      zin(470) = zin(472) + dzkl*zin(469)
                                      ! i4 = i4 + lang+1 =  473

                                      ! nk =    2

                                      xin(473) = xin(475) + dxkl*xin(472)
                                      yin(473) = yin(475) + dykl*yin(472)
                                      zin(473) = zin(475) + dzkl*zin(472)
                                      ! i4 = i4 + lang+1 =  476

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  471

                                      ! nl =    2

                                      ! i4 = i3 =  471

                                      ! do nk = 1,    2

                                      xin(471) = xin(473) + dxkl*xin(470)
                                      yin(471) = yin(473) + dykl*yin(470)
                                      zin(471) = zin(473) + dzkl*zin(470)
                                      ! i4 = i4 + lang+1 =  474

                                      ! nk =    2

                                      xin(474) = xin(476) + dxkl*xin(473)
                                      yin(474) = yin(476) + dykl*yin(473)
                                      zin(474) = zin(476) + dzkl*zin(473)
                                      ! i4 = i4 + lang+1 =  477

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  472

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  478

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  486

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  485

                                      xin(486) = xin(486) + dxkl*xin(485)
                                      yin(486) = yin(486) + dykl*yin(485)
                                      zin(486) = zin(486) + dzkl*zin(485)

                                      ! i3 = i4 =  485
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  484

                                      xin(485) = xin(485) + dxkl*xin(484)
                                      yin(485) = yin(485) + dykl*yin(484)
                                      zin(485) = zin(485) + dzkl*zin(484)

                                      ! i3 = i4 =  484
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  486

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  485

                                      xin(486) = xin(486) + dxkl*xin(485)
                                      yin(486) = yin(486) + dykl*yin(485)
                                      zin(486) = zin(486) + dzkl*zin(485)

                                      ! i3 = i4 =  485
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  479

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  479

                                      ! do nk = 1,    2

                                      xin(479) = xin(481) + dxkl*xin(478)
                                      yin(479) = yin(481) + dykl*yin(478)
                                      zin(479) = zin(481) + dzkl*zin(478)
                                      ! i4 = i4 + lang+1 =  482

                                      ! nk =    2

                                      xin(482) = xin(484) + dxkl*xin(481)
                                      yin(482) = yin(484) + dykl*yin(481)
                                      zin(482) = zin(484) + dzkl*zin(481)
                                      ! i4 = i4 + lang+1 =  485

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  480

                                      ! nl =    2

                                      ! i4 = i3 =  480

                                      ! do nk = 1,    2

                                      xin(480) = xin(482) + dxkl*xin(479)
                                      yin(480) = yin(482) + dykl*yin(479)
                                      zin(480) = zin(482) + dzkl*zin(479)
                                      ! i4 = i4 + lang+1 =  483

                                      ! nk =    2

                                      xin(483) = xin(485) + dxkl*xin(482)
                                      yin(483) = yin(485) + dykl*yin(482)
                                      zin(483) = zin(485) + dzkl*zin(482)
                                      ! i4 = i4 + lang+1 =  486

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  481

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  487

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  495

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  494

                                      xin(495) = xin(495) + dxkl*xin(494)
                                      yin(495) = yin(495) + dykl*yin(494)
                                      zin(495) = zin(495) + dzkl*zin(494)

                                      ! i3 = i4 =  494
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  493

                                      xin(494) = xin(494) + dxkl*xin(493)
                                      yin(494) = yin(494) + dykl*yin(493)
                                      zin(494) = zin(494) + dzkl*zin(493)

                                      ! i3 = i4 =  493
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  495

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  494

                                      xin(495) = xin(495) + dxkl*xin(494)
                                      yin(495) = yin(495) + dykl*yin(494)
                                      zin(495) = zin(495) + dzkl*zin(494)

                                      ! i3 = i4 =  494
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  488

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  488

                                      ! do nk = 1,    2

                                      xin(488) = xin(490) + dxkl*xin(487)
                                      yin(488) = yin(490) + dykl*yin(487)
                                      zin(488) = zin(490) + dzkl*zin(487)
                                      ! i4 = i4 + lang+1 =  491

                                      ! nk =    2

                                      xin(491) = xin(493) + dxkl*xin(490)
                                      yin(491) = yin(493) + dykl*yin(490)
                                      zin(491) = zin(493) + dzkl*zin(490)
                                      ! i4 = i4 + lang+1 =  494

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  489

                                      ! nl =    2

                                      ! i4 = i3 =  489

                                      ! do nk = 1,    2

                                      xin(489) = xin(491) + dxkl*xin(488)
                                      yin(489) = yin(491) + dykl*yin(488)
                                      zin(489) = zin(491) + dzkl*zin(488)
                                      ! i4 = i4 + lang+1 =  492

                                      ! nk =    2

                                      xin(492) = xin(494) + dxkl*xin(491)
                                      yin(492) = yin(494) + dykl*yin(491)
                                      zin(492) = zin(494) + dzkl*zin(491)
                                      ! i4 = i4 + lang+1 =  495

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  490

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  496

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  504

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  503

                                      xin(504) = xin(504) + dxkl*xin(503)
                                      yin(504) = yin(504) + dykl*yin(503)
                                      zin(504) = zin(504) + dzkl*zin(503)

                                      ! i3 = i4 =  503
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  502

                                      xin(503) = xin(503) + dxkl*xin(502)
                                      yin(503) = yin(503) + dykl*yin(502)
                                      zin(503) = zin(503) + dzkl*zin(502)

                                      ! i3 = i4 =  502
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  504

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  503

                                      xin(504) = xin(504) + dxkl*xin(503)
                                      yin(504) = yin(504) + dykl*yin(503)
                                      zin(504) = zin(504) + dzkl*zin(503)

                                      ! i3 = i4 =  503
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  497

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  497

                                      ! do nk = 1,    2

                                      xin(497) = xin(499) + dxkl*xin(496)
                                      yin(497) = yin(499) + dykl*yin(496)
                                      zin(497) = zin(499) + dzkl*zin(496)
                                      ! i4 = i4 + lang+1 =  500

                                      ! nk =    2

                                      xin(500) = xin(502) + dxkl*xin(499)
                                      yin(500) = yin(502) + dykl*yin(499)
                                      zin(500) = zin(502) + dzkl*zin(499)
                                      ! i4 = i4 + lang+1 =  503

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  498

                                      ! nl =    2

                                      ! i4 = i3 =  498

                                      ! do nk = 1,    2

                                      xin(498) = xin(500) + dxkl*xin(497)
                                      yin(498) = yin(500) + dykl*yin(497)
                                      zin(498) = zin(500) + dzkl*zin(497)
                                      ! i4 = i4 + lang+1 =  501

                                      ! nk =    2

                                      xin(501) = xin(503) + dxkl*xin(500)
                                      yin(501) = yin(503) + dykl*yin(500)
                                      zin(501) = zin(503) + dzkl*zin(500)
                                      ! i4 = i4 + lang+1 =  504

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  499

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  505

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  505

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  513

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  512

                                      xin(513) = xin(513) + dxkl*xin(512)
                                      yin(513) = yin(513) + dykl*yin(512)
                                      zin(513) = zin(513) + dzkl*zin(512)

                                      ! i3 = i4 =  512
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  511

                                      xin(512) = xin(512) + dxkl*xin(511)
                                      yin(512) = yin(512) + dykl*yin(511)
                                      zin(512) = zin(512) + dzkl*zin(511)

                                      ! i3 = i4 =  511
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  513

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  512

                                      xin(513) = xin(513) + dxkl*xin(512)
                                      yin(513) = yin(513) + dykl*yin(512)
                                      zin(513) = zin(513) + dzkl*zin(512)

                                      ! i3 = i4 =  512
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  506

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  506

                                      ! do nk = 1,    2

                                      xin(506) = xin(508) + dxkl*xin(505)
                                      yin(506) = yin(508) + dykl*yin(505)
                                      zin(506) = zin(508) + dzkl*zin(505)
                                      ! i4 = i4 + lang+1 =  509

                                      ! nk =    2

                                      xin(509) = xin(511) + dxkl*xin(508)
                                      yin(509) = yin(511) + dykl*yin(508)
                                      zin(509) = zin(511) + dzkl*zin(508)
                                      ! i4 = i4 + lang+1 =  512

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  507

                                      ! nl =    2

                                      ! i4 = i3 =  507

                                      ! do nk = 1,    2

                                      xin(507) = xin(509) + dxkl*xin(506)
                                      yin(507) = yin(509) + dykl*yin(506)
                                      zin(507) = zin(509) + dzkl*zin(506)
                                      ! i4 = i4 + lang+1 =  510

                                      ! nk =    2

                                      xin(510) = xin(512) + dxkl*xin(509)
                                      yin(510) = yin(512) + dykl*yin(509)
                                      zin(510) = zin(512) + dzkl*zin(509)
                                      ! i4 = i4 + lang+1 =  513

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  508

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  514

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  522

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  521

                                      xin(522) = xin(522) + dxkl*xin(521)
                                      yin(522) = yin(522) + dykl*yin(521)
                                      zin(522) = zin(522) + dzkl*zin(521)

                                      ! i3 = i4 =  521
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  520

                                      xin(521) = xin(521) + dxkl*xin(520)
                                      yin(521) = yin(521) + dykl*yin(520)
                                      zin(521) = zin(521) + dzkl*zin(520)

                                      ! i3 = i4 =  520
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  522

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  521

                                      xin(522) = xin(522) + dxkl*xin(521)
                                      yin(522) = yin(522) + dykl*yin(521)
                                      zin(522) = zin(522) + dzkl*zin(521)

                                      ! i3 = i4 =  521
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  515

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  515

                                      ! do nk = 1,    2

                                      xin(515) = xin(517) + dxkl*xin(514)
                                      yin(515) = yin(517) + dykl*yin(514)
                                      zin(515) = zin(517) + dzkl*zin(514)
                                      ! i4 = i4 + lang+1 =  518

                                      ! nk =    2

                                      xin(518) = xin(520) + dxkl*xin(517)
                                      yin(518) = yin(520) + dykl*yin(517)
                                      zin(518) = zin(520) + dzkl*zin(517)
                                      ! i4 = i4 + lang+1 =  521

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  516

                                      ! nl =    2

                                      ! i4 = i3 =  516

                                      ! do nk = 1,    2

                                      xin(516) = xin(518) + dxkl*xin(515)
                                      yin(516) = yin(518) + dykl*yin(515)
                                      zin(516) = zin(518) + dzkl*zin(515)
                                      ! i4 = i4 + lang+1 =  519

                                      ! nk =    2

                                      xin(519) = xin(521) + dxkl*xin(518)
                                      yin(519) = yin(521) + dykl*yin(518)
                                      zin(519) = zin(521) + dzkl*zin(518)
                                      ! i4 = i4 + lang+1 =  522

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  517

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  523

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  531

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  530

                                      xin(531) = xin(531) + dxkl*xin(530)
                                      yin(531) = yin(531) + dykl*yin(530)
                                      zin(531) = zin(531) + dzkl*zin(530)

                                      ! i3 = i4 =  530
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  529

                                      xin(530) = xin(530) + dxkl*xin(529)
                                      yin(530) = yin(530) + dykl*yin(529)
                                      zin(530) = zin(530) + dzkl*zin(529)

                                      ! i3 = i4 =  529
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  531

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  530

                                      xin(531) = xin(531) + dxkl*xin(530)
                                      yin(531) = yin(531) + dykl*yin(530)
                                      zin(531) = zin(531) + dzkl*zin(530)

                                      ! i3 = i4 =  530
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  524

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  524

                                      ! do nk = 1,    2

                                      xin(524) = xin(526) + dxkl*xin(523)
                                      yin(524) = yin(526) + dykl*yin(523)
                                      zin(524) = zin(526) + dzkl*zin(523)
                                      ! i4 = i4 + lang+1 =  527

                                      ! nk =    2

                                      xin(527) = xin(529) + dxkl*xin(526)
                                      yin(527) = yin(529) + dykl*yin(526)
                                      zin(527) = zin(529) + dzkl*zin(526)
                                      ! i4 = i4 + lang+1 =  530

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  525

                                      ! nl =    2

                                      ! i4 = i3 =  525

                                      ! do nk = 1,    2

                                      xin(525) = xin(527) + dxkl*xin(524)
                                      yin(525) = yin(527) + dykl*yin(524)
                                      zin(525) = zin(527) + dzkl*zin(524)
                                      ! i4 = i4 + lang+1 =  528

                                      ! nk =    2

                                      xin(528) = xin(530) + dxkl*xin(527)
                                      yin(528) = yin(530) + dykl*yin(527)
                                      zin(528) = zin(530) + dzkl*zin(527)
                                      ! i4 = i4 + lang+1 =  531

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  526

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  532

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  540

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  539

                                      xin(540) = xin(540) + dxkl*xin(539)
                                      yin(540) = yin(540) + dykl*yin(539)
                                      zin(540) = zin(540) + dzkl*zin(539)

                                      ! i3 = i4 =  539
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  538

                                      xin(539) = xin(539) + dxkl*xin(538)
                                      yin(539) = yin(539) + dykl*yin(538)
                                      zin(539) = zin(539) + dzkl*zin(538)

                                      ! i3 = i4 =  538
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  540

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  539

                                      xin(540) = xin(540) + dxkl*xin(539)
                                      yin(540) = yin(540) + dykl*yin(539)
                                      zin(540) = zin(540) + dzkl*zin(539)

                                      ! i3 = i4 =  539
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  533

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  533

                                      ! do nk = 1,    2

                                      xin(533) = xin(535) + dxkl*xin(532)
                                      yin(533) = yin(535) + dykl*yin(532)
                                      zin(533) = zin(535) + dzkl*zin(532)
                                      ! i4 = i4 + lang+1 =  536

                                      ! nk =    2

                                      xin(536) = xin(538) + dxkl*xin(535)
                                      yin(536) = yin(538) + dykl*yin(535)
                                      zin(536) = zin(538) + dzkl*zin(535)
                                      ! i4 = i4 + lang+1 =  539

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  534

                                      ! nl =    2

                                      ! i4 = i3 =  534

                                      ! do nk = 1,    2

                                      xin(534) = xin(536) + dxkl*xin(533)
                                      yin(534) = yin(536) + dykl*yin(533)
                                      zin(534) = zin(536) + dzkl*zin(533)
                                      ! i4 = i4 + lang+1 =  537

                                      ! nk =    2

                                      xin(537) = xin(539) + dxkl*xin(536)
                                      yin(537) = yin(539) + dykl*yin(536)
                                      zin(537) = zin(539) + dzkl*zin(536)
                                      ! i4 = i4 + lang+1 =  540

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  535

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  541

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  541

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  549

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  548

                                      xin(549) = xin(549) + dxkl*xin(548)
                                      yin(549) = yin(549) + dykl*yin(548)
                                      zin(549) = zin(549) + dzkl*zin(548)

                                      ! i3 = i4 =  548
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  547

                                      xin(548) = xin(548) + dxkl*xin(547)
                                      yin(548) = yin(548) + dykl*yin(547)
                                      zin(548) = zin(548) + dzkl*zin(547)

                                      ! i3 = i4 =  547
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  549

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  548

                                      xin(549) = xin(549) + dxkl*xin(548)
                                      yin(549) = yin(549) + dykl*yin(548)
                                      zin(549) = zin(549) + dzkl*zin(548)

                                      ! i3 = i4 =  548
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  542

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  542

                                      ! do nk = 1,    2

                                      xin(542) = xin(544) + dxkl*xin(541)
                                      yin(542) = yin(544) + dykl*yin(541)
                                      zin(542) = zin(544) + dzkl*zin(541)
                                      ! i4 = i4 + lang+1 =  545

                                      ! nk =    2

                                      xin(545) = xin(547) + dxkl*xin(544)
                                      yin(545) = yin(547) + dykl*yin(544)
                                      zin(545) = zin(547) + dzkl*zin(544)
                                      ! i4 = i4 + lang+1 =  548

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  543

                                      ! nl =    2

                                      ! i4 = i3 =  543

                                      ! do nk = 1,    2

                                      xin(543) = xin(545) + dxkl*xin(542)
                                      yin(543) = yin(545) + dykl*yin(542)
                                      zin(543) = zin(545) + dzkl*zin(542)
                                      ! i4 = i4 + lang+1 =  546

                                      ! nk =    2

                                      xin(546) = xin(548) + dxkl*xin(545)
                                      yin(546) = yin(548) + dykl*yin(545)
                                      zin(546) = zin(548) + dzkl*zin(545)
                                      ! i4 = i4 + lang+1 =  549

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  544

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  550

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  558

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  557

                                      xin(558) = xin(558) + dxkl*xin(557)
                                      yin(558) = yin(558) + dykl*yin(557)
                                      zin(558) = zin(558) + dzkl*zin(557)

                                      ! i3 = i4 =  557
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  556

                                      xin(557) = xin(557) + dxkl*xin(556)
                                      yin(557) = yin(557) + dykl*yin(556)
                                      zin(557) = zin(557) + dzkl*zin(556)

                                      ! i3 = i4 =  556
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  558

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  557

                                      xin(558) = xin(558) + dxkl*xin(557)
                                      yin(558) = yin(558) + dykl*yin(557)
                                      zin(558) = zin(558) + dzkl*zin(557)

                                      ! i3 = i4 =  557
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  551

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  551

                                      ! do nk = 1,    2

                                      xin(551) = xin(553) + dxkl*xin(550)
                                      yin(551) = yin(553) + dykl*yin(550)
                                      zin(551) = zin(553) + dzkl*zin(550)
                                      ! i4 = i4 + lang+1 =  554

                                      ! nk =    2

                                      xin(554) = xin(556) + dxkl*xin(553)
                                      yin(554) = yin(556) + dykl*yin(553)
                                      zin(554) = zin(556) + dzkl*zin(553)
                                      ! i4 = i4 + lang+1 =  557

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  552

                                      ! nl =    2

                                      ! i4 = i3 =  552

                                      ! do nk = 1,    2

                                      xin(552) = xin(554) + dxkl*xin(551)
                                      yin(552) = yin(554) + dykl*yin(551)
                                      zin(552) = zin(554) + dzkl*zin(551)
                                      ! i4 = i4 + lang+1 =  555

                                      ! nk =    2

                                      xin(555) = xin(557) + dxkl*xin(554)
                                      yin(555) = yin(557) + dykl*yin(554)
                                      zin(555) = zin(557) + dzkl*zin(554)
                                      ! i4 = i4 + lang+1 =  558

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  553

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  559

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  567

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  566

                                      xin(567) = xin(567) + dxkl*xin(566)
                                      yin(567) = yin(567) + dykl*yin(566)
                                      zin(567) = zin(567) + dzkl*zin(566)

                                      ! i3 = i4 =  566
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  565

                                      xin(566) = xin(566) + dxkl*xin(565)
                                      yin(566) = yin(566) + dykl*yin(565)
                                      zin(566) = zin(566) + dzkl*zin(565)

                                      ! i3 = i4 =  565
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  567

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  566

                                      xin(567) = xin(567) + dxkl*xin(566)
                                      yin(567) = yin(567) + dykl*yin(566)
                                      zin(567) = zin(567) + dzkl*zin(566)

                                      ! i3 = i4 =  566
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  560

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  560

                                      ! do nk = 1,    2

                                      xin(560) = xin(562) + dxkl*xin(559)
                                      yin(560) = yin(562) + dykl*yin(559)
                                      zin(560) = zin(562) + dzkl*zin(559)
                                      ! i4 = i4 + lang+1 =  563

                                      ! nk =    2

                                      xin(563) = xin(565) + dxkl*xin(562)
                                      yin(563) = yin(565) + dykl*yin(562)
                                      zin(563) = zin(565) + dzkl*zin(562)
                                      ! i4 = i4 + lang+1 =  566

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  561

                                      ! nl =    2

                                      ! i4 = i3 =  561

                                      ! do nk = 1,    2

                                      xin(561) = xin(563) + dxkl*xin(560)
                                      yin(561) = yin(563) + dykl*yin(560)
                                      zin(561) = zin(563) + dzkl*zin(560)
                                      ! i4 = i4 + lang+1 =  564

                                      ! nk =    2

                                      xin(564) = xin(566) + dxkl*xin(563)
                                      yin(564) = yin(566) + dykl*yin(563)
                                      zin(564) = zin(566) + dzkl*zin(563)
                                      ! i4 = i4 + lang+1 =  567

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  562

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  568

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  576

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  575

                                      xin(576) = xin(576) + dxkl*xin(575)
                                      yin(576) = yin(576) + dykl*yin(575)
                                      zin(576) = zin(576) + dzkl*zin(575)

                                      ! i3 = i4 =  575
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  574

                                      xin(575) = xin(575) + dxkl*xin(574)
                                      yin(575) = yin(575) + dykl*yin(574)
                                      zin(575) = zin(575) + dzkl*zin(574)

                                      ! i3 = i4 =  574
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  576

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  575

                                      xin(576) = xin(576) + dxkl*xin(575)
                                      yin(576) = yin(576) + dykl*yin(575)
                                      zin(576) = zin(576) + dzkl*zin(575)

                                      ! i3 = i4 =  575
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  569

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  569

                                      ! do nk = 1,    2

                                      xin(569) = xin(571) + dxkl*xin(568)
                                      yin(569) = yin(571) + dykl*yin(568)
                                      zin(569) = zin(571) + dzkl*zin(568)
                                      ! i4 = i4 + lang+1 =  572

                                      ! nk =    2

                                      xin(572) = xin(574) + dxkl*xin(571)
                                      yin(572) = yin(574) + dykl*yin(571)
                                      zin(572) = zin(574) + dzkl*zin(571)
                                      ! i4 = i4 + lang+1 =  575

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  570

                                      ! nl =    2

                                      ! i4 = i3 =  570

                                      ! do nk = 1,    2

                                      xin(570) = xin(572) + dxkl*xin(569)
                                      yin(570) = yin(572) + dykl*yin(569)
                                      zin(570) = zin(572) + dzkl*zin(569)
                                      ! i4 = i4 + lang+1 =  573

                                      ! nk =    2

                                      xin(573) = xin(575) + dxkl*xin(572)
                                      yin(573) = yin(575) + dykl*yin(572)
                                      zin(573) = zin(575) + dzkl*zin(572)
                                      ! i4 = i4 + lang+1 =  576

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  571

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  577

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  577

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  576

                                      u2 = roots(5)*rho
                                      f00 = expe*wghts(5)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  577

                                      xin(577) = 1.0_dp
                                      yin(577) = 1.0_dp
                                      zin(577) = f00

                                      ! i2 = in(2) =  613
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(613) = xc00
                                      yin(613) = yc00
                                      zin(613) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  580

                                      xin(580) = xcp00
                                      yin(580) = ycp00
                                      zin(580) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  616
                                      ! i2 =  613

                                      xin(616) = xcp00*xin(613) + cp10
                                      yin(616) = ycp00*yin(613) + cp10
                                      zin(616) = zcp00*zin(613) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  577
                                      ! i4 = i2 =  613

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  649
                                      ! i3 =  577
                                      ! i4 =  613

                                      xin(649) = c10*xin(577) + xc00*xin(613)
                                      yin(649) = c10*yin(577) + yc00*yin(613)
                                      zin(649) = c10*zin(577) + zc00*zin(613)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  652
                                      ! i5 =  649
                                      ! i4 =  613

                                      xin(652) = xcp00*xin(649) + cp10*xin(613)
                                      yin(652) = ycp00*yin(649) + cp10*yin(613)
                                      zin(652) = zcp00*zin(649) + cp10*zin(613)

                                      ! ------------------

                                      ! i3 = i4 =  613
                                      ! i4 = i5 =  649

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  685
                                      ! i3 =  613
                                      ! i4 =  649

                                      xin(685) = c10*xin(613) + xc00*xin(649)
                                      yin(685) = c10*yin(613) + yc00*yin(649)
                                      zin(685) = c10*zin(613) + zc00*zin(649)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  688
                                      ! i5 =  685
                                      ! i4 =  649

                                      xin(688) = xcp00*xin(685) + cp10*xin(649)
                                      yin(688) = ycp00*yin(685) + cp10*yin(649)
                                      zin(688) = zcp00*zin(685) + cp10*zin(649)

                                      ! ------------------

                                      ! i3 = i4 =  649
                                      ! i4 = i5 =  685

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  694
                                      ! i3 =  649
                                      ! i4 =  685

                                      xin(694) = c10*xin(649) + xc00*xin(685)
                                      yin(694) = c10*yin(649) + yc00*yin(685)
                                      zin(694) = c10*zin(649) + zc00*zin(685)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  697
                                      ! i5 =  694
                                      ! i4 =  685

                                      xin(697) = xcp00*xin(694) + cp10*xin(685)
                                      yin(697) = ycp00*yin(694) + cp10*yin(685)
                                      zin(697) = zcp00*zin(694) + cp10*zin(685)

                                      ! ------------------

                                      ! i3 = i4 =  685
                                      ! i4 = i5 =  694

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  703
                                      ! i3 =  685
                                      ! i4 =  694

                                      xin(703) = c10*xin(685) + xc00*xin(694)
                                      yin(703) = c10*yin(685) + yc00*yin(694)
                                      zin(703) = c10*zin(685) + zc00*zin(694)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  706
                                      ! i5 =  703
                                      ! i4 =  694

                                      xin(706) = xcp00*xin(703) + cp10*xin(694)
                                      yin(706) = ycp00*yin(703) + cp10*yin(694)
                                      zin(706) = zcp00*zin(703) + cp10*zin(694)

                                      ! ------------------

                                      ! i3 = i4 =  694
                                      ! i4 = i5 =  703

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  712
                                      ! i3 =  694
                                      ! i4 =  703

                                      xin(712) = c10*xin(694) + xc00*xin(703)
                                      yin(712) = c10*yin(694) + yc00*yin(703)
                                      zin(712) = c10*zin(694) + zc00*zin(703)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  715
                                      ! i5 =  712
                                      ! i4 =  703

                                      xin(715) = xcp00*xin(712) + cp10*xin(703)
                                      yin(715) = ycp00*yin(712) + cp10*yin(703)
                                      zin(715) = zcp00*zin(712) + cp10*zin(703)

                                      ! ------------------

                                      ! i3 = i4 =  703
                                      ! i4 = i5 =  712

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  577
                                      ! i4 = i1+k2 =  580

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  583
                                      ! i3 =  577
                                      ! i4 =  580

                                      xin(583) = cp01*xin(577) + xcp00*xin(580)
                                      yin(583) = cp01*yin(577) + ycp00*yin(580)
                                      zin(583) = cp01*zin(577) + zcp00*zin(580)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  619

                                      xin(619) = xc00*xin(583) + c01*xin(580)
                                      yin(619) = yc00*yin(583) + c01*yin(580)
                                      zin(619) = zc00*zin(583) + c01*zin(580)

                                      ! ------------------

                                      ! i3 = i4 =  580
                                      ! i4 = i5 =  583

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  584
                                      ! i3 =  580
                                      ! i4 =  583

                                      xin(584) = cp01*xin(580) + xcp00*xin(583)
                                      yin(584) = cp01*yin(580) + ycp00*yin(583)
                                      zin(584) = cp01*zin(580) + zcp00*zin(583)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  620

                                      xin(620) = xc00*xin(584) + c01*xin(583)
                                      yin(620) = yc00*yin(584) + c01*yin(583)
                                      zin(620) = zc00*zin(584) + c01*zin(583)

                                      ! ------------------

                                      ! i3 = i4 =  583
                                      ! i4 = i5 =  584

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  585
                                      ! i3 =  583
                                      ! i4 =  584

                                      xin(585) = cp01*xin(583) + xcp00*xin(584)
                                      yin(585) = cp01*yin(583) + ycp00*yin(584)
                                      zin(585) = cp01*zin(583) + zcp00*zin(584)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  621

                                      xin(621) = xc00*xin(585) + c01*xin(584)
                                      yin(621) = yc00*yin(585) + c01*yin(584)
                                      zin(621) = zc00*zin(585) + c01*zin(584)

                                      ! ------------------

                                      ! i3 = i4 =  584
                                      ! i4 = i5 =  585

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  577
                                      ! i4 = i2 =  613

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  649

                                      xin(655) = c10*xin(583) + xc00*xin(619) + c01*xin(616)
                                      yin(655) = c10*yin(583) + yc00*yin(619) + c01*yin(616)
                                      zin(655) = c10*zin(583) + zc00*zin(619) + c01*zin(616)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  613
                                      ! i4 = i5 =  649

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  685

                                      xin(691) = c10*xin(619) + xc00*xin(655) + c01*xin(652)
                                      yin(691) = c10*yin(619) + yc00*yin(655) + c01*yin(652)
                                      zin(691) = c10*zin(619) + zc00*zin(655) + c01*zin(652)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  649
                                      ! i4 = i5 =  685

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  694

                                      xin(700) = c10*xin(655) + xc00*xin(691) + c01*xin(688)
                                      yin(700) = c10*yin(655) + yc00*yin(691) + c01*yin(688)
                                      zin(700) = c10*zin(655) + zc00*zin(691) + c01*zin(688)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  685
                                      ! i4 = i5 =  694

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  703

                                      xin(709) = c10*xin(691) + xc00*xin(700) + c01*xin(697)
                                      yin(709) = c10*yin(691) + yc00*yin(700) + c01*yin(697)
                                      zin(709) = c10*zin(691) + zc00*zin(700) + c01*zin(697)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  694
                                      ! i4 = i5 =  703

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  712

                                      xin(718) = c10*xin(700) + xc00*xin(709) + c01*xin(706)
                                      yin(718) = c10*yin(700) + yc00*yin(709) + c01*yin(706)
                                      zin(718) = c10*zin(700) + zc00*zin(709) + c01*zin(706)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  703
                                      ! i4 = i5 =  712

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  577
                                      ! i4 = i2 =  613

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  649

                                      xin(656) = c10*xin(584) + xc00*xin(620) + c01*xin(619)
                                      yin(656) = c10*yin(584) + yc00*yin(620) + c01*yin(619)
                                      zin(656) = c10*zin(584) + zc00*zin(620) + c01*zin(619)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  613
                                      ! i4 = i5 =  649

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  685

                                      xin(692) = c10*xin(620) + xc00*xin(656) + c01*xin(655)
                                      yin(692) = c10*yin(620) + yc00*yin(656) + c01*yin(655)
                                      zin(692) = c10*zin(620) + zc00*zin(656) + c01*zin(655)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  649
                                      ! i4 = i5 =  685

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  694

                                      xin(701) = c10*xin(656) + xc00*xin(692) + c01*xin(691)
                                      yin(701) = c10*yin(656) + yc00*yin(692) + c01*yin(691)
                                      zin(701) = c10*zin(656) + zc00*zin(692) + c01*zin(691)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  685
                                      ! i4 = i5 =  694

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  703

                                      xin(710) = c10*xin(692) + xc00*xin(701) + c01*xin(700)
                                      yin(710) = c10*yin(692) + yc00*yin(701) + c01*yin(700)
                                      zin(710) = c10*zin(692) + zc00*zin(701) + c01*zin(700)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  694
                                      ! i4 = i5 =  703

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  712

                                      xin(719) = c10*xin(701) + xc00*xin(710) + c01*xin(709)
                                      yin(719) = c10*yin(701) + yc00*yin(710) + c01*yin(709)
                                      zin(719) = c10*zin(701) + zc00*zin(710) + c01*zin(709)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  703
                                      ! i4 = i5 =  712

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  577
                                      ! i4 = i2 =  613

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  649

                                      xin(657) = c10*xin(585) + xc00*xin(621) + c01*xin(620)
                                      yin(657) = c10*yin(585) + yc00*yin(621) + c01*yin(620)
                                      zin(657) = c10*zin(585) + zc00*zin(621) + c01*zin(620)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  613
                                      ! i4 = i5 =  649

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  685

                                      xin(693) = c10*xin(621) + xc00*xin(657) + c01*xin(656)
                                      yin(693) = c10*yin(621) + yc00*yin(657) + c01*yin(656)
                                      zin(693) = c10*zin(621) + zc00*zin(657) + c01*zin(656)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  649
                                      ! i4 = i5 =  685

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  694

                                      xin(702) = c10*xin(657) + xc00*xin(693) + c01*xin(692)
                                      yin(702) = c10*yin(657) + yc00*yin(693) + c01*yin(692)
                                      zin(702) = c10*zin(657) + zc00*zin(693) + c01*zin(692)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  685
                                      ! i4 = i5 =  694

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  703

                                      xin(711) = c10*xin(693) + xc00*xin(702) + c01*xin(701)
                                      yin(711) = c10*yin(693) + yc00*yin(702) + c01*yin(701)
                                      zin(711) = c10*zin(693) + zc00*zin(702) + c01*zin(701)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  694
                                      ! i4 = i5 =  703

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  712

                                      xin(720) = c10*xin(702) + xc00*xin(711) + c01*xin(710)
                                      yin(720) = c10*yin(702) + yc00*yin(711) + c01*yin(710)
                                      zin(720) = c10*zin(702) + zc00*zin(711) + c01*zin(710)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  703
                                      ! i4 = i5 =  712

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  712

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  712

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  703

                                      xin(712) = xin(712) + dxij*xin(703)
                                      yin(712) = yin(712) + dyij*yin(703)
                                      zin(712) = zin(712) + dzij*zin(703)

                                      ! i3 = i4 =  703
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  694

                                      xin(703) = xin(703) + dxij*xin(694)
                                      yin(703) = yin(703) + dyij*yin(694)
                                      zin(703) = zin(703) + dzij*zin(694)

                                      ! i3 = i4 =  694
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  685

                                      xin(694) = xin(694) + dxij*xin(685)
                                      yin(694) = yin(694) + dyij*yin(685)
                                      zin(694) = zin(694) + dzij*zin(685)

                                      ! i3 = i4 =  685
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  712

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  703

                                      xin(712) = xin(712) + dxij*xin(703)
                                      yin(712) = yin(712) + dyij*yin(703)
                                      zin(712) = zin(712) + dzij*zin(703)

                                      ! i3 = i4 =  703
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  694

                                      xin(703) = xin(703) + dxij*xin(694)
                                      yin(703) = yin(703) + dyij*yin(694)
                                      zin(703) = zin(703) + dzij*zin(694)

                                      ! i3 = i4 =  694
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  712

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  703

                                      xin(712) = xin(712) + dxij*xin(703)
                                      yin(712) = yin(712) + dyij*yin(703)
                                      zin(712) = zin(712) + dzij*zin(703)

                                      ! i3 = i4 =  703
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  586

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  586

                                      ! do ni = 1,    3

                                      xin(586) = xin(613) + dxij*xin(577)
                                      yin(586) = yin(613) + dyij*yin(577)
                                      zin(586) = zin(613) + dzij*zin(577)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  622

                                      ! ni =    2

                                      xin(622) = xin(649) + dxij*xin(613)
                                      yin(622) = yin(649) + dyij*yin(613)
                                      zin(622) = zin(649) + dzij*zin(613)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  658

                                      ! ni =    3

                                      xin(658) = xin(685) + dxij*xin(649)
                                      yin(658) = yin(685) + dyij*yin(649)
                                      zin(658) = zin(685) + dzij*zin(649)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  694

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  595

                                      ! nj =    2

                                      ! i4 = i3 =  595

                                      ! do ni = 1,    3

                                      xin(595) = xin(622) + dxij*xin(586)
                                      yin(595) = yin(622) + dyij*yin(586)
                                      zin(595) = zin(622) + dzij*zin(586)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  631

                                      ! ni =    2

                                      xin(631) = xin(658) + dxij*xin(622)
                                      yin(631) = yin(658) + dyij*yin(622)
                                      zin(631) = zin(658) + dzij*zin(622)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  667

                                      ! ni =    3

                                      xin(667) = xin(694) + dxij*xin(658)
                                      yin(667) = yin(694) + dyij*yin(658)
                                      zin(667) = zin(694) + dzij*zin(658)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  703

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  604

                                      ! nj =    3

                                      ! i4 = i3 =  604

                                      ! do ni = 1,    3

                                      xin(604) = xin(631) + dxij*xin(595)
                                      yin(604) = yin(631) + dyij*yin(595)
                                      zin(604) = zin(631) + dzij*zin(595)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  640

                                      ! ni =    2

                                      xin(640) = xin(667) + dxij*xin(631)
                                      yin(640) = yin(667) + dyij*yin(631)
                                      zin(640) = zin(667) + dzij*zin(631)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  676

                                      ! ni =    3

                                      xin(676) = xin(703) + dxij*xin(667)
                                      yin(676) = yin(703) + dyij*yin(667)
                                      zin(676) = zin(703) + dzij*zin(667)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  712

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  613

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  715

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  706

                                      xin(715) = xin(715) + dxij*xin(706)
                                      yin(715) = yin(715) + dyij*yin(706)
                                      zin(715) = zin(715) + dzij*zin(706)

                                      ! i3 = i4 =  706
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  697

                                      xin(706) = xin(706) + dxij*xin(697)
                                      yin(706) = yin(706) + dyij*yin(697)
                                      zin(706) = zin(706) + dzij*zin(697)

                                      ! i3 = i4 =  697
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  688

                                      xin(697) = xin(697) + dxij*xin(688)
                                      yin(697) = yin(697) + dyij*yin(688)
                                      zin(697) = zin(697) + dzij*zin(688)

                                      ! i3 = i4 =  688
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  715

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  706

                                      xin(715) = xin(715) + dxij*xin(706)
                                      yin(715) = yin(715) + dyij*yin(706)
                                      zin(715) = zin(715) + dzij*zin(706)

                                      ! i3 = i4 =  706
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  697

                                      xin(706) = xin(706) + dxij*xin(697)
                                      yin(706) = yin(706) + dyij*yin(697)
                                      zin(706) = zin(706) + dzij*zin(697)

                                      ! i3 = i4 =  697
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  715

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  706

                                      xin(715) = xin(715) + dxij*xin(706)
                                      yin(715) = yin(715) + dyij*yin(706)
                                      zin(715) = zin(715) + dzij*zin(706)

                                      ! i3 = i4 =  706
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  589

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  589

                                      ! do ni = 1,    3

                                      xin(589) = xin(616) + dxij*xin(580)
                                      yin(589) = yin(616) + dyij*yin(580)
                                      zin(589) = zin(616) + dzij*zin(580)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  625

                                      ! ni =    2

                                      xin(625) = xin(652) + dxij*xin(616)
                                      yin(625) = yin(652) + dyij*yin(616)
                                      zin(625) = zin(652) + dzij*zin(616)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  661

                                      ! ni =    3

                                      xin(661) = xin(688) + dxij*xin(652)
                                      yin(661) = yin(688) + dyij*yin(652)
                                      zin(661) = zin(688) + dzij*zin(652)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  697

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  598

                                      ! nj =    2

                                      ! i4 = i3 =  598

                                      ! do ni = 1,    3

                                      xin(598) = xin(625) + dxij*xin(589)
                                      yin(598) = yin(625) + dyij*yin(589)
                                      zin(598) = zin(625) + dzij*zin(589)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  634

                                      ! ni =    2

                                      xin(634) = xin(661) + dxij*xin(625)
                                      yin(634) = yin(661) + dyij*yin(625)
                                      zin(634) = zin(661) + dzij*zin(625)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  670

                                      ! ni =    3

                                      xin(670) = xin(697) + dxij*xin(661)
                                      yin(670) = yin(697) + dyij*yin(661)
                                      zin(670) = zin(697) + dzij*zin(661)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  706

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  607

                                      ! nj =    3

                                      ! i4 = i3 =  607

                                      ! do ni = 1,    3

                                      xin(607) = xin(634) + dxij*xin(598)
                                      yin(607) = yin(634) + dyij*yin(598)
                                      zin(607) = zin(634) + dzij*zin(598)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  643

                                      ! ni =    2

                                      xin(643) = xin(670) + dxij*xin(634)
                                      yin(643) = yin(670) + dyij*yin(634)
                                      zin(643) = zin(670) + dzij*zin(634)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  679

                                      ! ni =    3

                                      xin(679) = xin(706) + dxij*xin(670)
                                      yin(679) = yin(706) + dyij*yin(670)
                                      zin(679) = zin(706) + dzij*zin(670)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  715

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  616

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  718

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  709

                                      xin(718) = xin(718) + dxij*xin(709)
                                      yin(718) = yin(718) + dyij*yin(709)
                                      zin(718) = zin(718) + dzij*zin(709)

                                      ! i3 = i4 =  709
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  700

                                      xin(709) = xin(709) + dxij*xin(700)
                                      yin(709) = yin(709) + dyij*yin(700)
                                      zin(709) = zin(709) + dzij*zin(700)

                                      ! i3 = i4 =  700
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  691

                                      xin(700) = xin(700) + dxij*xin(691)
                                      yin(700) = yin(700) + dyij*yin(691)
                                      zin(700) = zin(700) + dzij*zin(691)

                                      ! i3 = i4 =  691
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  718

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  709

                                      xin(718) = xin(718) + dxij*xin(709)
                                      yin(718) = yin(718) + dyij*yin(709)
                                      zin(718) = zin(718) + dzij*zin(709)

                                      ! i3 = i4 =  709
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  700

                                      xin(709) = xin(709) + dxij*xin(700)
                                      yin(709) = yin(709) + dyij*yin(700)
                                      zin(709) = zin(709) + dzij*zin(700)

                                      ! i3 = i4 =  700
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  718

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  709

                                      xin(718) = xin(718) + dxij*xin(709)
                                      yin(718) = yin(718) + dyij*yin(709)
                                      zin(718) = zin(718) + dzij*zin(709)

                                      ! i3 = i4 =  709
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  592

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  592

                                      ! do ni = 1,    3

                                      xin(592) = xin(619) + dxij*xin(583)
                                      yin(592) = yin(619) + dyij*yin(583)
                                      zin(592) = zin(619) + dzij*zin(583)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  628

                                      ! ni =    2

                                      xin(628) = xin(655) + dxij*xin(619)
                                      yin(628) = yin(655) + dyij*yin(619)
                                      zin(628) = zin(655) + dzij*zin(619)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  664

                                      ! ni =    3

                                      xin(664) = xin(691) + dxij*xin(655)
                                      yin(664) = yin(691) + dyij*yin(655)
                                      zin(664) = zin(691) + dzij*zin(655)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  700

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  601

                                      ! nj =    2

                                      ! i4 = i3 =  601

                                      ! do ni = 1,    3

                                      xin(601) = xin(628) + dxij*xin(592)
                                      yin(601) = yin(628) + dyij*yin(592)
                                      zin(601) = zin(628) + dzij*zin(592)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  637

                                      ! ni =    2

                                      xin(637) = xin(664) + dxij*xin(628)
                                      yin(637) = yin(664) + dyij*yin(628)
                                      zin(637) = zin(664) + dzij*zin(628)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  673

                                      ! ni =    3

                                      xin(673) = xin(700) + dxij*xin(664)
                                      yin(673) = yin(700) + dyij*yin(664)
                                      zin(673) = zin(700) + dzij*zin(664)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  709

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  610

                                      ! nj =    3

                                      ! i4 = i3 =  610

                                      ! do ni = 1,    3

                                      xin(610) = xin(637) + dxij*xin(601)
                                      yin(610) = yin(637) + dyij*yin(601)
                                      zin(610) = zin(637) + dzij*zin(601)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  646

                                      ! ni =    2

                                      xin(646) = xin(673) + dxij*xin(637)
                                      yin(646) = yin(673) + dyij*yin(637)
                                      zin(646) = zin(673) + dzij*zin(637)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  682

                                      ! ni =    3

                                      xin(682) = xin(709) + dxij*xin(673)
                                      yin(682) = yin(709) + dyij*yin(673)
                                      zin(682) = zin(709) + dzij*zin(673)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  718

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  619

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  719

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  710

                                      xin(719) = xin(719) + dxij*xin(710)
                                      yin(719) = yin(719) + dyij*yin(710)
                                      zin(719) = zin(719) + dzij*zin(710)

                                      ! i3 = i4 =  710
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  701

                                      xin(710) = xin(710) + dxij*xin(701)
                                      yin(710) = yin(710) + dyij*yin(701)
                                      zin(710) = zin(710) + dzij*zin(701)

                                      ! i3 = i4 =  701
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  692

                                      xin(701) = xin(701) + dxij*xin(692)
                                      yin(701) = yin(701) + dyij*yin(692)
                                      zin(701) = zin(701) + dzij*zin(692)

                                      ! i3 = i4 =  692
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  719

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  710

                                      xin(719) = xin(719) + dxij*xin(710)
                                      yin(719) = yin(719) + dyij*yin(710)
                                      zin(719) = zin(719) + dzij*zin(710)

                                      ! i3 = i4 =  710
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  701

                                      xin(710) = xin(710) + dxij*xin(701)
                                      yin(710) = yin(710) + dyij*yin(701)
                                      zin(710) = zin(710) + dzij*zin(701)

                                      ! i3 = i4 =  701
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  719

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  710

                                      xin(719) = xin(719) + dxij*xin(710)
                                      yin(719) = yin(719) + dyij*yin(710)
                                      zin(719) = zin(719) + dzij*zin(710)

                                      ! i3 = i4 =  710
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  593

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  593

                                      ! do ni = 1,    3

                                      xin(593) = xin(620) + dxij*xin(584)
                                      yin(593) = yin(620) + dyij*yin(584)
                                      zin(593) = zin(620) + dzij*zin(584)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  629

                                      ! ni =    2

                                      xin(629) = xin(656) + dxij*xin(620)
                                      yin(629) = yin(656) + dyij*yin(620)
                                      zin(629) = zin(656) + dzij*zin(620)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  665

                                      ! ni =    3

                                      xin(665) = xin(692) + dxij*xin(656)
                                      yin(665) = yin(692) + dyij*yin(656)
                                      zin(665) = zin(692) + dzij*zin(656)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  701

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  602

                                      ! nj =    2

                                      ! i4 = i3 =  602

                                      ! do ni = 1,    3

                                      xin(602) = xin(629) + dxij*xin(593)
                                      yin(602) = yin(629) + dyij*yin(593)
                                      zin(602) = zin(629) + dzij*zin(593)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  638

                                      ! ni =    2

                                      xin(638) = xin(665) + dxij*xin(629)
                                      yin(638) = yin(665) + dyij*yin(629)
                                      zin(638) = zin(665) + dzij*zin(629)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  674

                                      ! ni =    3

                                      xin(674) = xin(701) + dxij*xin(665)
                                      yin(674) = yin(701) + dyij*yin(665)
                                      zin(674) = zin(701) + dzij*zin(665)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  710

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  611

                                      ! nj =    3

                                      ! i4 = i3 =  611

                                      ! do ni = 1,    3

                                      xin(611) = xin(638) + dxij*xin(602)
                                      yin(611) = yin(638) + dyij*yin(602)
                                      zin(611) = zin(638) + dzij*zin(602)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  647

                                      ! ni =    2

                                      xin(647) = xin(674) + dxij*xin(638)
                                      yin(647) = yin(674) + dyij*yin(638)
                                      zin(647) = zin(674) + dzij*zin(638)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  683

                                      ! ni =    3

                                      xin(683) = xin(710) + dxij*xin(674)
                                      yin(683) = yin(710) + dyij*yin(674)
                                      zin(683) = zin(710) + dzij*zin(674)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  719

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  620

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  720

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  711

                                      xin(720) = xin(720) + dxij*xin(711)
                                      yin(720) = yin(720) + dyij*yin(711)
                                      zin(720) = zin(720) + dzij*zin(711)

                                      ! i3 = i4 =  711
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  702

                                      xin(711) = xin(711) + dxij*xin(702)
                                      yin(711) = yin(711) + dyij*yin(702)
                                      zin(711) = zin(711) + dzij*zin(702)

                                      ! i3 = i4 =  702
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  693

                                      xin(702) = xin(702) + dxij*xin(693)
                                      yin(702) = yin(702) + dyij*yin(693)
                                      zin(702) = zin(702) + dzij*zin(693)

                                      ! i3 = i4 =  693
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  720

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  711

                                      xin(720) = xin(720) + dxij*xin(711)
                                      yin(720) = yin(720) + dyij*yin(711)
                                      zin(720) = zin(720) + dzij*zin(711)

                                      ! i3 = i4 =  711
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  702

                                      xin(711) = xin(711) + dxij*xin(702)
                                      yin(711) = yin(711) + dyij*yin(702)
                                      zin(711) = zin(711) + dzij*zin(702)

                                      ! i3 = i4 =  702
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  720

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  711

                                      xin(720) = xin(720) + dxij*xin(711)
                                      yin(720) = yin(720) + dyij*yin(711)
                                      zin(720) = zin(720) + dzij*zin(711)

                                      ! i3 = i4 =  711
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  594

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  594

                                      ! do ni = 1,    3

                                      xin(594) = xin(621) + dxij*xin(585)
                                      yin(594) = yin(621) + dyij*yin(585)
                                      zin(594) = zin(621) + dzij*zin(585)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  630

                                      ! ni =    2

                                      xin(630) = xin(657) + dxij*xin(621)
                                      yin(630) = yin(657) + dyij*yin(621)
                                      zin(630) = zin(657) + dzij*zin(621)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  666

                                      ! ni =    3

                                      xin(666) = xin(693) + dxij*xin(657)
                                      yin(666) = yin(693) + dyij*yin(657)
                                      zin(666) = zin(693) + dzij*zin(657)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  702

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  603

                                      ! nj =    2

                                      ! i4 = i3 =  603

                                      ! do ni = 1,    3

                                      xin(603) = xin(630) + dxij*xin(594)
                                      yin(603) = yin(630) + dyij*yin(594)
                                      zin(603) = zin(630) + dzij*zin(594)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  639

                                      ! ni =    2

                                      xin(639) = xin(666) + dxij*xin(630)
                                      yin(639) = yin(666) + dyij*yin(630)
                                      zin(639) = zin(666) + dzij*zin(630)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  675

                                      ! ni =    3

                                      xin(675) = xin(702) + dxij*xin(666)
                                      yin(675) = yin(702) + dyij*yin(666)
                                      zin(675) = zin(702) + dzij*zin(666)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  711

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  612

                                      ! nj =    3

                                      ! i4 = i3 =  612

                                      ! do ni = 1,    3

                                      xin(612) = xin(639) + dxij*xin(603)
                                      yin(612) = yin(639) + dyij*yin(603)
                                      zin(612) = zin(639) + dzij*zin(603)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  648

                                      ! ni =    2

                                      xin(648) = xin(675) + dxij*xin(639)
                                      yin(648) = yin(675) + dyij*yin(639)
                                      zin(648) = zin(675) + dzij*zin(639)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  684

                                      ! ni =    3

                                      xin(684) = xin(711) + dxij*xin(675)
                                      yin(684) = yin(711) + dyij*yin(675)
                                      zin(684) = zin(711) + dzij*zin(675)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  720

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  621

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  577

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  585

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  584

                                      xin(585) = xin(585) + dxkl*xin(584)
                                      yin(585) = yin(585) + dykl*yin(584)
                                      zin(585) = zin(585) + dzkl*zin(584)

                                      ! i3 = i4 =  584
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  583

                                      xin(584) = xin(584) + dxkl*xin(583)
                                      yin(584) = yin(584) + dykl*yin(583)
                                      zin(584) = zin(584) + dzkl*zin(583)

                                      ! i3 = i4 =  583
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  585

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  584

                                      xin(585) = xin(585) + dxkl*xin(584)
                                      yin(585) = yin(585) + dykl*yin(584)
                                      zin(585) = zin(585) + dzkl*zin(584)

                                      ! i3 = i4 =  584
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  578

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  578

                                      ! do nk = 1,    2

                                      xin(578) = xin(580) + dxkl*xin(577)
                                      yin(578) = yin(580) + dykl*yin(577)
                                      zin(578) = zin(580) + dzkl*zin(577)
                                      ! i4 = i4 + lang+1 =  581

                                      ! nk =    2

                                      xin(581) = xin(583) + dxkl*xin(580)
                                      yin(581) = yin(583) + dykl*yin(580)
                                      zin(581) = zin(583) + dzkl*zin(580)
                                      ! i4 = i4 + lang+1 =  584

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  579

                                      ! nl =    2

                                      ! i4 = i3 =  579

                                      ! do nk = 1,    2

                                      xin(579) = xin(581) + dxkl*xin(578)
                                      yin(579) = yin(581) + dykl*yin(578)
                                      zin(579) = zin(581) + dzkl*zin(578)
                                      ! i4 = i4 + lang+1 =  582

                                      ! nk =    2

                                      xin(582) = xin(584) + dxkl*xin(581)
                                      yin(582) = yin(584) + dykl*yin(581)
                                      zin(582) = zin(584) + dzkl*zin(581)
                                      ! i4 = i4 + lang+1 =  585

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  580

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  586

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  594

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  593

                                      xin(594) = xin(594) + dxkl*xin(593)
                                      yin(594) = yin(594) + dykl*yin(593)
                                      zin(594) = zin(594) + dzkl*zin(593)

                                      ! i3 = i4 =  593
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  592

                                      xin(593) = xin(593) + dxkl*xin(592)
                                      yin(593) = yin(593) + dykl*yin(592)
                                      zin(593) = zin(593) + dzkl*zin(592)

                                      ! i3 = i4 =  592
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  594

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  593

                                      xin(594) = xin(594) + dxkl*xin(593)
                                      yin(594) = yin(594) + dykl*yin(593)
                                      zin(594) = zin(594) + dzkl*zin(593)

                                      ! i3 = i4 =  593
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  587

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  587

                                      ! do nk = 1,    2

                                      xin(587) = xin(589) + dxkl*xin(586)
                                      yin(587) = yin(589) + dykl*yin(586)
                                      zin(587) = zin(589) + dzkl*zin(586)
                                      ! i4 = i4 + lang+1 =  590

                                      ! nk =    2

                                      xin(590) = xin(592) + dxkl*xin(589)
                                      yin(590) = yin(592) + dykl*yin(589)
                                      zin(590) = zin(592) + dzkl*zin(589)
                                      ! i4 = i4 + lang+1 =  593

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  588

                                      ! nl =    2

                                      ! i4 = i3 =  588

                                      ! do nk = 1,    2

                                      xin(588) = xin(590) + dxkl*xin(587)
                                      yin(588) = yin(590) + dykl*yin(587)
                                      zin(588) = zin(590) + dzkl*zin(587)
                                      ! i4 = i4 + lang+1 =  591

                                      ! nk =    2

                                      xin(591) = xin(593) + dxkl*xin(590)
                                      yin(591) = yin(593) + dykl*yin(590)
                                      zin(591) = zin(593) + dzkl*zin(590)
                                      ! i4 = i4 + lang+1 =  594

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  589

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  595

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  603

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  602

                                      xin(603) = xin(603) + dxkl*xin(602)
                                      yin(603) = yin(603) + dykl*yin(602)
                                      zin(603) = zin(603) + dzkl*zin(602)

                                      ! i3 = i4 =  602
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  601

                                      xin(602) = xin(602) + dxkl*xin(601)
                                      yin(602) = yin(602) + dykl*yin(601)
                                      zin(602) = zin(602) + dzkl*zin(601)

                                      ! i3 = i4 =  601
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  603

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  602

                                      xin(603) = xin(603) + dxkl*xin(602)
                                      yin(603) = yin(603) + dykl*yin(602)
                                      zin(603) = zin(603) + dzkl*zin(602)

                                      ! i3 = i4 =  602
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  596

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  596

                                      ! do nk = 1,    2

                                      xin(596) = xin(598) + dxkl*xin(595)
                                      yin(596) = yin(598) + dykl*yin(595)
                                      zin(596) = zin(598) + dzkl*zin(595)
                                      ! i4 = i4 + lang+1 =  599

                                      ! nk =    2

                                      xin(599) = xin(601) + dxkl*xin(598)
                                      yin(599) = yin(601) + dykl*yin(598)
                                      zin(599) = zin(601) + dzkl*zin(598)
                                      ! i4 = i4 + lang+1 =  602

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  597

                                      ! nl =    2

                                      ! i4 = i3 =  597

                                      ! do nk = 1,    2

                                      xin(597) = xin(599) + dxkl*xin(596)
                                      yin(597) = yin(599) + dykl*yin(596)
                                      zin(597) = zin(599) + dzkl*zin(596)
                                      ! i4 = i4 + lang+1 =  600

                                      ! nk =    2

                                      xin(600) = xin(602) + dxkl*xin(599)
                                      yin(600) = yin(602) + dykl*yin(599)
                                      zin(600) = zin(602) + dzkl*zin(599)
                                      ! i4 = i4 + lang+1 =  603

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  598

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  604

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  612

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  611

                                      xin(612) = xin(612) + dxkl*xin(611)
                                      yin(612) = yin(612) + dykl*yin(611)
                                      zin(612) = zin(612) + dzkl*zin(611)

                                      ! i3 = i4 =  611
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  610

                                      xin(611) = xin(611) + dxkl*xin(610)
                                      yin(611) = yin(611) + dykl*yin(610)
                                      zin(611) = zin(611) + dzkl*zin(610)

                                      ! i3 = i4 =  610
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  612

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  611

                                      xin(612) = xin(612) + dxkl*xin(611)
                                      yin(612) = yin(612) + dykl*yin(611)
                                      zin(612) = zin(612) + dzkl*zin(611)

                                      ! i3 = i4 =  611
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  605

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  605

                                      ! do nk = 1,    2

                                      xin(605) = xin(607) + dxkl*xin(604)
                                      yin(605) = yin(607) + dykl*yin(604)
                                      zin(605) = zin(607) + dzkl*zin(604)
                                      ! i4 = i4 + lang+1 =  608

                                      ! nk =    2

                                      xin(608) = xin(610) + dxkl*xin(607)
                                      yin(608) = yin(610) + dykl*yin(607)
                                      zin(608) = zin(610) + dzkl*zin(607)
                                      ! i4 = i4 + lang+1 =  611

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  606

                                      ! nl =    2

                                      ! i4 = i3 =  606

                                      ! do nk = 1,    2

                                      xin(606) = xin(608) + dxkl*xin(605)
                                      yin(606) = yin(608) + dykl*yin(605)
                                      zin(606) = zin(608) + dzkl*zin(605)
                                      ! i4 = i4 + lang+1 =  609

                                      ! nk =    2

                                      xin(609) = xin(611) + dxkl*xin(608)
                                      yin(609) = yin(611) + dykl*yin(608)
                                      zin(609) = zin(611) + dzkl*zin(608)
                                      ! i4 = i4 + lang+1 =  612

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  607

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  613

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  613

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  621

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  620

                                      xin(621) = xin(621) + dxkl*xin(620)
                                      yin(621) = yin(621) + dykl*yin(620)
                                      zin(621) = zin(621) + dzkl*zin(620)

                                      ! i3 = i4 =  620
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  619

                                      xin(620) = xin(620) + dxkl*xin(619)
                                      yin(620) = yin(620) + dykl*yin(619)
                                      zin(620) = zin(620) + dzkl*zin(619)

                                      ! i3 = i4 =  619
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  621

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  620

                                      xin(621) = xin(621) + dxkl*xin(620)
                                      yin(621) = yin(621) + dykl*yin(620)
                                      zin(621) = zin(621) + dzkl*zin(620)

                                      ! i3 = i4 =  620
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  614

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  614

                                      ! do nk = 1,    2

                                      xin(614) = xin(616) + dxkl*xin(613)
                                      yin(614) = yin(616) + dykl*yin(613)
                                      zin(614) = zin(616) + dzkl*zin(613)
                                      ! i4 = i4 + lang+1 =  617

                                      ! nk =    2

                                      xin(617) = xin(619) + dxkl*xin(616)
                                      yin(617) = yin(619) + dykl*yin(616)
                                      zin(617) = zin(619) + dzkl*zin(616)
                                      ! i4 = i4 + lang+1 =  620

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  615

                                      ! nl =    2

                                      ! i4 = i3 =  615

                                      ! do nk = 1,    2

                                      xin(615) = xin(617) + dxkl*xin(614)
                                      yin(615) = yin(617) + dykl*yin(614)
                                      zin(615) = zin(617) + dzkl*zin(614)
                                      ! i4 = i4 + lang+1 =  618

                                      ! nk =    2

                                      xin(618) = xin(620) + dxkl*xin(617)
                                      yin(618) = yin(620) + dykl*yin(617)
                                      zin(618) = zin(620) + dzkl*zin(617)
                                      ! i4 = i4 + lang+1 =  621

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  616

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  622

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  630

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  629

                                      xin(630) = xin(630) + dxkl*xin(629)
                                      yin(630) = yin(630) + dykl*yin(629)
                                      zin(630) = zin(630) + dzkl*zin(629)

                                      ! i3 = i4 =  629
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  628

                                      xin(629) = xin(629) + dxkl*xin(628)
                                      yin(629) = yin(629) + dykl*yin(628)
                                      zin(629) = zin(629) + dzkl*zin(628)

                                      ! i3 = i4 =  628
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  630

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  629

                                      xin(630) = xin(630) + dxkl*xin(629)
                                      yin(630) = yin(630) + dykl*yin(629)
                                      zin(630) = zin(630) + dzkl*zin(629)

                                      ! i3 = i4 =  629
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  623

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  623

                                      ! do nk = 1,    2

                                      xin(623) = xin(625) + dxkl*xin(622)
                                      yin(623) = yin(625) + dykl*yin(622)
                                      zin(623) = zin(625) + dzkl*zin(622)
                                      ! i4 = i4 + lang+1 =  626

                                      ! nk =    2

                                      xin(626) = xin(628) + dxkl*xin(625)
                                      yin(626) = yin(628) + dykl*yin(625)
                                      zin(626) = zin(628) + dzkl*zin(625)
                                      ! i4 = i4 + lang+1 =  629

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  624

                                      ! nl =    2

                                      ! i4 = i3 =  624

                                      ! do nk = 1,    2

                                      xin(624) = xin(626) + dxkl*xin(623)
                                      yin(624) = yin(626) + dykl*yin(623)
                                      zin(624) = zin(626) + dzkl*zin(623)
                                      ! i4 = i4 + lang+1 =  627

                                      ! nk =    2

                                      xin(627) = xin(629) + dxkl*xin(626)
                                      yin(627) = yin(629) + dykl*yin(626)
                                      zin(627) = zin(629) + dzkl*zin(626)
                                      ! i4 = i4 + lang+1 =  630

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  625

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  631

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  639

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  638

                                      xin(639) = xin(639) + dxkl*xin(638)
                                      yin(639) = yin(639) + dykl*yin(638)
                                      zin(639) = zin(639) + dzkl*zin(638)

                                      ! i3 = i4 =  638
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  637

                                      xin(638) = xin(638) + dxkl*xin(637)
                                      yin(638) = yin(638) + dykl*yin(637)
                                      zin(638) = zin(638) + dzkl*zin(637)

                                      ! i3 = i4 =  637
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  639

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  638

                                      xin(639) = xin(639) + dxkl*xin(638)
                                      yin(639) = yin(639) + dykl*yin(638)
                                      zin(639) = zin(639) + dzkl*zin(638)

                                      ! i3 = i4 =  638
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  632

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  632

                                      ! do nk = 1,    2

                                      xin(632) = xin(634) + dxkl*xin(631)
                                      yin(632) = yin(634) + dykl*yin(631)
                                      zin(632) = zin(634) + dzkl*zin(631)
                                      ! i4 = i4 + lang+1 =  635

                                      ! nk =    2

                                      xin(635) = xin(637) + dxkl*xin(634)
                                      yin(635) = yin(637) + dykl*yin(634)
                                      zin(635) = zin(637) + dzkl*zin(634)
                                      ! i4 = i4 + lang+1 =  638

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  633

                                      ! nl =    2

                                      ! i4 = i3 =  633

                                      ! do nk = 1,    2

                                      xin(633) = xin(635) + dxkl*xin(632)
                                      yin(633) = yin(635) + dykl*yin(632)
                                      zin(633) = zin(635) + dzkl*zin(632)
                                      ! i4 = i4 + lang+1 =  636

                                      ! nk =    2

                                      xin(636) = xin(638) + dxkl*xin(635)
                                      yin(636) = yin(638) + dykl*yin(635)
                                      zin(636) = zin(638) + dzkl*zin(635)
                                      ! i4 = i4 + lang+1 =  639

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  634

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  640

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  648

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  647

                                      xin(648) = xin(648) + dxkl*xin(647)
                                      yin(648) = yin(648) + dykl*yin(647)
                                      zin(648) = zin(648) + dzkl*zin(647)

                                      ! i3 = i4 =  647
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  646

                                      xin(647) = xin(647) + dxkl*xin(646)
                                      yin(647) = yin(647) + dykl*yin(646)
                                      zin(647) = zin(647) + dzkl*zin(646)

                                      ! i3 = i4 =  646
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  648

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  647

                                      xin(648) = xin(648) + dxkl*xin(647)
                                      yin(648) = yin(648) + dykl*yin(647)
                                      zin(648) = zin(648) + dzkl*zin(647)

                                      ! i3 = i4 =  647
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  641

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  641

                                      ! do nk = 1,    2

                                      xin(641) = xin(643) + dxkl*xin(640)
                                      yin(641) = yin(643) + dykl*yin(640)
                                      zin(641) = zin(643) + dzkl*zin(640)
                                      ! i4 = i4 + lang+1 =  644

                                      ! nk =    2

                                      xin(644) = xin(646) + dxkl*xin(643)
                                      yin(644) = yin(646) + dykl*yin(643)
                                      zin(644) = zin(646) + dzkl*zin(643)
                                      ! i4 = i4 + lang+1 =  647

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  642

                                      ! nl =    2

                                      ! i4 = i3 =  642

                                      ! do nk = 1,    2

                                      xin(642) = xin(644) + dxkl*xin(641)
                                      yin(642) = yin(644) + dykl*yin(641)
                                      zin(642) = zin(644) + dzkl*zin(641)
                                      ! i4 = i4 + lang+1 =  645

                                      ! nk =    2

                                      xin(645) = xin(647) + dxkl*xin(644)
                                      yin(645) = yin(647) + dykl*yin(644)
                                      zin(645) = zin(647) + dzkl*zin(644)
                                      ! i4 = i4 + lang+1 =  648

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  643

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  649

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  649

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  657

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  656

                                      xin(657) = xin(657) + dxkl*xin(656)
                                      yin(657) = yin(657) + dykl*yin(656)
                                      zin(657) = zin(657) + dzkl*zin(656)

                                      ! i3 = i4 =  656
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  655

                                      xin(656) = xin(656) + dxkl*xin(655)
                                      yin(656) = yin(656) + dykl*yin(655)
                                      zin(656) = zin(656) + dzkl*zin(655)

                                      ! i3 = i4 =  655
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  657

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  656

                                      xin(657) = xin(657) + dxkl*xin(656)
                                      yin(657) = yin(657) + dykl*yin(656)
                                      zin(657) = zin(657) + dzkl*zin(656)

                                      ! i3 = i4 =  656
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  650

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  650

                                      ! do nk = 1,    2

                                      xin(650) = xin(652) + dxkl*xin(649)
                                      yin(650) = yin(652) + dykl*yin(649)
                                      zin(650) = zin(652) + dzkl*zin(649)
                                      ! i4 = i4 + lang+1 =  653

                                      ! nk =    2

                                      xin(653) = xin(655) + dxkl*xin(652)
                                      yin(653) = yin(655) + dykl*yin(652)
                                      zin(653) = zin(655) + dzkl*zin(652)
                                      ! i4 = i4 + lang+1 =  656

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  651

                                      ! nl =    2

                                      ! i4 = i3 =  651

                                      ! do nk = 1,    2

                                      xin(651) = xin(653) + dxkl*xin(650)
                                      yin(651) = yin(653) + dykl*yin(650)
                                      zin(651) = zin(653) + dzkl*zin(650)
                                      ! i4 = i4 + lang+1 =  654

                                      ! nk =    2

                                      xin(654) = xin(656) + dxkl*xin(653)
                                      yin(654) = yin(656) + dykl*yin(653)
                                      zin(654) = zin(656) + dzkl*zin(653)
                                      ! i4 = i4 + lang+1 =  657

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  652

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  658

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  666

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  665

                                      xin(666) = xin(666) + dxkl*xin(665)
                                      yin(666) = yin(666) + dykl*yin(665)
                                      zin(666) = zin(666) + dzkl*zin(665)

                                      ! i3 = i4 =  665
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  664

                                      xin(665) = xin(665) + dxkl*xin(664)
                                      yin(665) = yin(665) + dykl*yin(664)
                                      zin(665) = zin(665) + dzkl*zin(664)

                                      ! i3 = i4 =  664
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  666

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  665

                                      xin(666) = xin(666) + dxkl*xin(665)
                                      yin(666) = yin(666) + dykl*yin(665)
                                      zin(666) = zin(666) + dzkl*zin(665)

                                      ! i3 = i4 =  665
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  659

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  659

                                      ! do nk = 1,    2

                                      xin(659) = xin(661) + dxkl*xin(658)
                                      yin(659) = yin(661) + dykl*yin(658)
                                      zin(659) = zin(661) + dzkl*zin(658)
                                      ! i4 = i4 + lang+1 =  662

                                      ! nk =    2

                                      xin(662) = xin(664) + dxkl*xin(661)
                                      yin(662) = yin(664) + dykl*yin(661)
                                      zin(662) = zin(664) + dzkl*zin(661)
                                      ! i4 = i4 + lang+1 =  665

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  660

                                      ! nl =    2

                                      ! i4 = i3 =  660

                                      ! do nk = 1,    2

                                      xin(660) = xin(662) + dxkl*xin(659)
                                      yin(660) = yin(662) + dykl*yin(659)
                                      zin(660) = zin(662) + dzkl*zin(659)
                                      ! i4 = i4 + lang+1 =  663

                                      ! nk =    2

                                      xin(663) = xin(665) + dxkl*xin(662)
                                      yin(663) = yin(665) + dykl*yin(662)
                                      zin(663) = zin(665) + dzkl*zin(662)
                                      ! i4 = i4 + lang+1 =  666

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  661

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  667

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  675

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  674

                                      xin(675) = xin(675) + dxkl*xin(674)
                                      yin(675) = yin(675) + dykl*yin(674)
                                      zin(675) = zin(675) + dzkl*zin(674)

                                      ! i3 = i4 =  674
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  673

                                      xin(674) = xin(674) + dxkl*xin(673)
                                      yin(674) = yin(674) + dykl*yin(673)
                                      zin(674) = zin(674) + dzkl*zin(673)

                                      ! i3 = i4 =  673
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  675

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  674

                                      xin(675) = xin(675) + dxkl*xin(674)
                                      yin(675) = yin(675) + dykl*yin(674)
                                      zin(675) = zin(675) + dzkl*zin(674)

                                      ! i3 = i4 =  674
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  668

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  668

                                      ! do nk = 1,    2

                                      xin(668) = xin(670) + dxkl*xin(667)
                                      yin(668) = yin(670) + dykl*yin(667)
                                      zin(668) = zin(670) + dzkl*zin(667)
                                      ! i4 = i4 + lang+1 =  671

                                      ! nk =    2

                                      xin(671) = xin(673) + dxkl*xin(670)
                                      yin(671) = yin(673) + dykl*yin(670)
                                      zin(671) = zin(673) + dzkl*zin(670)
                                      ! i4 = i4 + lang+1 =  674

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  669

                                      ! nl =    2

                                      ! i4 = i3 =  669

                                      ! do nk = 1,    2

                                      xin(669) = xin(671) + dxkl*xin(668)
                                      yin(669) = yin(671) + dykl*yin(668)
                                      zin(669) = zin(671) + dzkl*zin(668)
                                      ! i4 = i4 + lang+1 =  672

                                      ! nk =    2

                                      xin(672) = xin(674) + dxkl*xin(671)
                                      yin(672) = yin(674) + dykl*yin(671)
                                      zin(672) = zin(674) + dzkl*zin(671)
                                      ! i4 = i4 + lang+1 =  675

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  670

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  676

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  684

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  683

                                      xin(684) = xin(684) + dxkl*xin(683)
                                      yin(684) = yin(684) + dykl*yin(683)
                                      zin(684) = zin(684) + dzkl*zin(683)

                                      ! i3 = i4 =  683
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  682

                                      xin(683) = xin(683) + dxkl*xin(682)
                                      yin(683) = yin(683) + dykl*yin(682)
                                      zin(683) = zin(683) + dzkl*zin(682)

                                      ! i3 = i4 =  682
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  684

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  683

                                      xin(684) = xin(684) + dxkl*xin(683)
                                      yin(684) = yin(684) + dykl*yin(683)
                                      zin(684) = zin(684) + dzkl*zin(683)

                                      ! i3 = i4 =  683
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  677

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  677

                                      ! do nk = 1,    2

                                      xin(677) = xin(679) + dxkl*xin(676)
                                      yin(677) = yin(679) + dykl*yin(676)
                                      zin(677) = zin(679) + dzkl*zin(676)
                                      ! i4 = i4 + lang+1 =  680

                                      ! nk =    2

                                      xin(680) = xin(682) + dxkl*xin(679)
                                      yin(680) = yin(682) + dykl*yin(679)
                                      zin(680) = zin(682) + dzkl*zin(679)
                                      ! i4 = i4 + lang+1 =  683

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  678

                                      ! nl =    2

                                      ! i4 = i3 =  678

                                      ! do nk = 1,    2

                                      xin(678) = xin(680) + dxkl*xin(677)
                                      yin(678) = yin(680) + dykl*yin(677)
                                      zin(678) = zin(680) + dzkl*zin(677)
                                      ! i4 = i4 + lang+1 =  681

                                      ! nk =    2

                                      xin(681) = xin(683) + dxkl*xin(680)
                                      yin(681) = yin(683) + dykl*yin(680)
                                      zin(681) = zin(683) + dzkl*zin(680)
                                      ! i4 = i4 + lang+1 =  684

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  679

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  685

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  685

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  693

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  692

                                      xin(693) = xin(693) + dxkl*xin(692)
                                      yin(693) = yin(693) + dykl*yin(692)
                                      zin(693) = zin(693) + dzkl*zin(692)

                                      ! i3 = i4 =  692
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  691

                                      xin(692) = xin(692) + dxkl*xin(691)
                                      yin(692) = yin(692) + dykl*yin(691)
                                      zin(692) = zin(692) + dzkl*zin(691)

                                      ! i3 = i4 =  691
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  693

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  692

                                      xin(693) = xin(693) + dxkl*xin(692)
                                      yin(693) = yin(693) + dykl*yin(692)
                                      zin(693) = zin(693) + dzkl*zin(692)

                                      ! i3 = i4 =  692
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  686

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  686

                                      ! do nk = 1,    2

                                      xin(686) = xin(688) + dxkl*xin(685)
                                      yin(686) = yin(688) + dykl*yin(685)
                                      zin(686) = zin(688) + dzkl*zin(685)
                                      ! i4 = i4 + lang+1 =  689

                                      ! nk =    2

                                      xin(689) = xin(691) + dxkl*xin(688)
                                      yin(689) = yin(691) + dykl*yin(688)
                                      zin(689) = zin(691) + dzkl*zin(688)
                                      ! i4 = i4 + lang+1 =  692

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  687

                                      ! nl =    2

                                      ! i4 = i3 =  687

                                      ! do nk = 1,    2

                                      xin(687) = xin(689) + dxkl*xin(686)
                                      yin(687) = yin(689) + dykl*yin(686)
                                      zin(687) = zin(689) + dzkl*zin(686)
                                      ! i4 = i4 + lang+1 =  690

                                      ! nk =    2

                                      xin(690) = xin(692) + dxkl*xin(689)
                                      yin(690) = yin(692) + dykl*yin(689)
                                      zin(690) = zin(692) + dzkl*zin(689)
                                      ! i4 = i4 + lang+1 =  693

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  688

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  694

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  702

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  701

                                      xin(702) = xin(702) + dxkl*xin(701)
                                      yin(702) = yin(702) + dykl*yin(701)
                                      zin(702) = zin(702) + dzkl*zin(701)

                                      ! i3 = i4 =  701
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  700

                                      xin(701) = xin(701) + dxkl*xin(700)
                                      yin(701) = yin(701) + dykl*yin(700)
                                      zin(701) = zin(701) + dzkl*zin(700)

                                      ! i3 = i4 =  700
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  702

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  701

                                      xin(702) = xin(702) + dxkl*xin(701)
                                      yin(702) = yin(702) + dykl*yin(701)
                                      zin(702) = zin(702) + dzkl*zin(701)

                                      ! i3 = i4 =  701
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  695

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  695

                                      ! do nk = 1,    2

                                      xin(695) = xin(697) + dxkl*xin(694)
                                      yin(695) = yin(697) + dykl*yin(694)
                                      zin(695) = zin(697) + dzkl*zin(694)
                                      ! i4 = i4 + lang+1 =  698

                                      ! nk =    2

                                      xin(698) = xin(700) + dxkl*xin(697)
                                      yin(698) = yin(700) + dykl*yin(697)
                                      zin(698) = zin(700) + dzkl*zin(697)
                                      ! i4 = i4 + lang+1 =  701

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  696

                                      ! nl =    2

                                      ! i4 = i3 =  696

                                      ! do nk = 1,    2

                                      xin(696) = xin(698) + dxkl*xin(695)
                                      yin(696) = yin(698) + dykl*yin(695)
                                      zin(696) = zin(698) + dzkl*zin(695)
                                      ! i4 = i4 + lang+1 =  699

                                      ! nk =    2

                                      xin(699) = xin(701) + dxkl*xin(698)
                                      yin(699) = yin(701) + dykl*yin(698)
                                      zin(699) = zin(701) + dzkl*zin(698)
                                      ! i4 = i4 + lang+1 =  702

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  697

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  703

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  711

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  710

                                      xin(711) = xin(711) + dxkl*xin(710)
                                      yin(711) = yin(711) + dykl*yin(710)
                                      zin(711) = zin(711) + dzkl*zin(710)

                                      ! i3 = i4 =  710
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  709

                                      xin(710) = xin(710) + dxkl*xin(709)
                                      yin(710) = yin(710) + dykl*yin(709)
                                      zin(710) = zin(710) + dzkl*zin(709)

                                      ! i3 = i4 =  709
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  711

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  710

                                      xin(711) = xin(711) + dxkl*xin(710)
                                      yin(711) = yin(711) + dykl*yin(710)
                                      zin(711) = zin(711) + dzkl*zin(710)

                                      ! i3 = i4 =  710
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  704

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  704

                                      ! do nk = 1,    2

                                      xin(704) = xin(706) + dxkl*xin(703)
                                      yin(704) = yin(706) + dykl*yin(703)
                                      zin(704) = zin(706) + dzkl*zin(703)
                                      ! i4 = i4 + lang+1 =  707

                                      ! nk =    2

                                      xin(707) = xin(709) + dxkl*xin(706)
                                      yin(707) = yin(709) + dykl*yin(706)
                                      zin(707) = zin(709) + dzkl*zin(706)
                                      ! i4 = i4 + lang+1 =  710

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  705

                                      ! nl =    2

                                      ! i4 = i3 =  705

                                      ! do nk = 1,    2

                                      xin(705) = xin(707) + dxkl*xin(704)
                                      yin(705) = yin(707) + dykl*yin(704)
                                      zin(705) = zin(707) + dzkl*zin(704)
                                      ! i4 = i4 + lang+1 =  708

                                      ! nk =    2

                                      xin(708) = xin(710) + dxkl*xin(707)
                                      yin(708) = yin(710) + dykl*yin(707)
                                      zin(708) = zin(710) + dzkl*zin(707)
                                      ! i4 = i4 + lang+1 =  711

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  706

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  712

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  720

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  719

                                      xin(720) = xin(720) + dxkl*xin(719)
                                      yin(720) = yin(720) + dykl*yin(719)
                                      zin(720) = zin(720) + dzkl*zin(719)

                                      ! i3 = i4 =  719
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  718

                                      xin(719) = xin(719) + dxkl*xin(718)
                                      yin(719) = yin(719) + dykl*yin(718)
                                      zin(719) = zin(719) + dzkl*zin(718)

                                      ! i3 = i4 =  718
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  720

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  719

                                      xin(720) = xin(720) + dxkl*xin(719)
                                      yin(720) = yin(720) + dykl*yin(719)
                                      zin(720) = zin(720) + dzkl*zin(719)

                                      ! i3 = i4 =  719
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  713

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  713

                                      ! do nk = 1,    2

                                      xin(713) = xin(715) + dxkl*xin(712)
                                      yin(713) = yin(715) + dykl*yin(712)
                                      zin(713) = zin(715) + dzkl*zin(712)
                                      ! i4 = i4 + lang+1 =  716

                                      ! nk =    2

                                      xin(716) = xin(718) + dxkl*xin(715)
                                      yin(716) = yin(718) + dykl*yin(715)
                                      zin(716) = zin(718) + dzkl*zin(715)
                                      ! i4 = i4 + lang+1 =  719

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  714

                                      ! nl =    2

                                      ! i4 = i3 =  714

                                      ! do nk = 1,    2

                                      xin(714) = xin(716) + dxkl*xin(713)
                                      yin(714) = yin(716) + dykl*yin(713)
                                      zin(714) = zin(716) + dzkl*zin(713)
                                      ! i4 = i4 + lang+1 =  717

                                      ! nk =    2

                                      xin(717) = xin(719) + dxkl*xin(716)
                                      yin(717) = yin(719) + dykl*yin(716)
                                      zin(717) = zin(719) + dzkl*zin(716)
                                      ! i4 = i4 + lang+1 =  720

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  715

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  721

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  721

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  720

                                      u2 = roots(6)*rho
                                      f00 = expe*wghts(6)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  721

                                      xin(721) = 1.0_dp
                                      yin(721) = 1.0_dp
                                      zin(721) = f00

                                      ! i2 = in(2) =  757
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(757) = xc00
                                      yin(757) = yc00
                                      zin(757) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  724

                                      xin(724) = xcp00
                                      yin(724) = ycp00
                                      zin(724) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  760
                                      ! i2 =  757

                                      xin(760) = xcp00*xin(757) + cp10
                                      yin(760) = ycp00*yin(757) + cp10
                                      zin(760) = zcp00*zin(757) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  721
                                      ! i4 = i2 =  757

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  793
                                      ! i3 =  721
                                      ! i4 =  757

                                      xin(793) = c10*xin(721) + xc00*xin(757)
                                      yin(793) = c10*yin(721) + yc00*yin(757)
                                      zin(793) = c10*zin(721) + zc00*zin(757)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  796
                                      ! i5 =  793
                                      ! i4 =  757

                                      xin(796) = xcp00*xin(793) + cp10*xin(757)
                                      yin(796) = ycp00*yin(793) + cp10*yin(757)
                                      zin(796) = zcp00*zin(793) + cp10*zin(757)

                                      ! ------------------

                                      ! i3 = i4 =  757
                                      ! i4 = i5 =  793

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  829
                                      ! i3 =  757
                                      ! i4 =  793

                                      xin(829) = c10*xin(757) + xc00*xin(793)
                                      yin(829) = c10*yin(757) + yc00*yin(793)
                                      zin(829) = c10*zin(757) + zc00*zin(793)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  832
                                      ! i5 =  829
                                      ! i4 =  793

                                      xin(832) = xcp00*xin(829) + cp10*xin(793)
                                      yin(832) = ycp00*yin(829) + cp10*yin(793)
                                      zin(832) = zcp00*zin(829) + cp10*zin(793)

                                      ! ------------------

                                      ! i3 = i4 =  793
                                      ! i4 = i5 =  829

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  838
                                      ! i3 =  793
                                      ! i4 =  829

                                      xin(838) = c10*xin(793) + xc00*xin(829)
                                      yin(838) = c10*yin(793) + yc00*yin(829)
                                      zin(838) = c10*zin(793) + zc00*zin(829)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  841
                                      ! i5 =  838
                                      ! i4 =  829

                                      xin(841) = xcp00*xin(838) + cp10*xin(829)
                                      yin(841) = ycp00*yin(838) + cp10*yin(829)
                                      zin(841) = zcp00*zin(838) + cp10*zin(829)

                                      ! ------------------

                                      ! i3 = i4 =  829
                                      ! i4 = i5 =  838

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  847
                                      ! i3 =  829
                                      ! i4 =  838

                                      xin(847) = c10*xin(829) + xc00*xin(838)
                                      yin(847) = c10*yin(829) + yc00*yin(838)
                                      zin(847) = c10*zin(829) + zc00*zin(838)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  850
                                      ! i5 =  847
                                      ! i4 =  838

                                      xin(850) = xcp00*xin(847) + cp10*xin(838)
                                      yin(850) = ycp00*yin(847) + cp10*yin(838)
                                      zin(850) = zcp00*zin(847) + cp10*zin(838)

                                      ! ------------------

                                      ! i3 = i4 =  838
                                      ! i4 = i5 =  847

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  856
                                      ! i3 =  838
                                      ! i4 =  847

                                      xin(856) = c10*xin(838) + xc00*xin(847)
                                      yin(856) = c10*yin(838) + yc00*yin(847)
                                      zin(856) = c10*zin(838) + zc00*zin(847)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  859
                                      ! i5 =  856
                                      ! i4 =  847

                                      xin(859) = xcp00*xin(856) + cp10*xin(847)
                                      yin(859) = ycp00*yin(856) + cp10*yin(847)
                                      zin(859) = zcp00*zin(856) + cp10*zin(847)

                                      ! ------------------

                                      ! i3 = i4 =  847
                                      ! i4 = i5 =  856

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  721
                                      ! i4 = i1+k2 =  724

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  727
                                      ! i3 =  721
                                      ! i4 =  724

                                      xin(727) = cp01*xin(721) + xcp00*xin(724)
                                      yin(727) = cp01*yin(721) + ycp00*yin(724)
                                      zin(727) = cp01*zin(721) + zcp00*zin(724)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  763

                                      xin(763) = xc00*xin(727) + c01*xin(724)
                                      yin(763) = yc00*yin(727) + c01*yin(724)
                                      zin(763) = zc00*zin(727) + c01*zin(724)

                                      ! ------------------

                                      ! i3 = i4 =  724
                                      ! i4 = i5 =  727

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  728
                                      ! i3 =  724
                                      ! i4 =  727

                                      xin(728) = cp01*xin(724) + xcp00*xin(727)
                                      yin(728) = cp01*yin(724) + ycp00*yin(727)
                                      zin(728) = cp01*zin(724) + zcp00*zin(727)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  764

                                      xin(764) = xc00*xin(728) + c01*xin(727)
                                      yin(764) = yc00*yin(728) + c01*yin(727)
                                      zin(764) = zc00*zin(728) + c01*zin(727)

                                      ! ------------------

                                      ! i3 = i4 =  727
                                      ! i4 = i5 =  728

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  729
                                      ! i3 =  727
                                      ! i4 =  728

                                      xin(729) = cp01*xin(727) + xcp00*xin(728)
                                      yin(729) = cp01*yin(727) + ycp00*yin(728)
                                      zin(729) = cp01*zin(727) + zcp00*zin(728)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  765

                                      xin(765) = xc00*xin(729) + c01*xin(728)
                                      yin(765) = yc00*yin(729) + c01*yin(728)
                                      zin(765) = zc00*zin(729) + c01*zin(728)

                                      ! ------------------

                                      ! i3 = i4 =  728
                                      ! i4 = i5 =  729

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  721
                                      ! i4 = i2 =  757

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  793

                                      xin(799) = c10*xin(727) + xc00*xin(763) + c01*xin(760)
                                      yin(799) = c10*yin(727) + yc00*yin(763) + c01*yin(760)
                                      zin(799) = c10*zin(727) + zc00*zin(763) + c01*zin(760)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  757
                                      ! i4 = i5 =  793

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  829

                                      xin(835) = c10*xin(763) + xc00*xin(799) + c01*xin(796)
                                      yin(835) = c10*yin(763) + yc00*yin(799) + c01*yin(796)
                                      zin(835) = c10*zin(763) + zc00*zin(799) + c01*zin(796)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  793
                                      ! i4 = i5 =  829

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  838

                                      xin(844) = c10*xin(799) + xc00*xin(835) + c01*xin(832)
                                      yin(844) = c10*yin(799) + yc00*yin(835) + c01*yin(832)
                                      zin(844) = c10*zin(799) + zc00*zin(835) + c01*zin(832)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  829
                                      ! i4 = i5 =  838

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  847

                                      xin(853) = c10*xin(835) + xc00*xin(844) + c01*xin(841)
                                      yin(853) = c10*yin(835) + yc00*yin(844) + c01*yin(841)
                                      zin(853) = c10*zin(835) + zc00*zin(844) + c01*zin(841)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  838
                                      ! i4 = i5 =  847

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  856

                                      xin(862) = c10*xin(844) + xc00*xin(853) + c01*xin(850)
                                      yin(862) = c10*yin(844) + yc00*yin(853) + c01*yin(850)
                                      zin(862) = c10*zin(844) + zc00*zin(853) + c01*zin(850)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  847
                                      ! i4 = i5 =  856

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  721
                                      ! i4 = i2 =  757

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  793

                                      xin(800) = c10*xin(728) + xc00*xin(764) + c01*xin(763)
                                      yin(800) = c10*yin(728) + yc00*yin(764) + c01*yin(763)
                                      zin(800) = c10*zin(728) + zc00*zin(764) + c01*zin(763)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  757
                                      ! i4 = i5 =  793

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  829

                                      xin(836) = c10*xin(764) + xc00*xin(800) + c01*xin(799)
                                      yin(836) = c10*yin(764) + yc00*yin(800) + c01*yin(799)
                                      zin(836) = c10*zin(764) + zc00*zin(800) + c01*zin(799)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  793
                                      ! i4 = i5 =  829

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  838

                                      xin(845) = c10*xin(800) + xc00*xin(836) + c01*xin(835)
                                      yin(845) = c10*yin(800) + yc00*yin(836) + c01*yin(835)
                                      zin(845) = c10*zin(800) + zc00*zin(836) + c01*zin(835)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  829
                                      ! i4 = i5 =  838

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  847

                                      xin(854) = c10*xin(836) + xc00*xin(845) + c01*xin(844)
                                      yin(854) = c10*yin(836) + yc00*yin(845) + c01*yin(844)
                                      zin(854) = c10*zin(836) + zc00*zin(845) + c01*zin(844)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  838
                                      ! i4 = i5 =  847

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  856

                                      xin(863) = c10*xin(845) + xc00*xin(854) + c01*xin(853)
                                      yin(863) = c10*yin(845) + yc00*yin(854) + c01*yin(853)
                                      zin(863) = c10*zin(845) + zc00*zin(854) + c01*zin(853)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  847
                                      ! i4 = i5 =  856

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  721
                                      ! i4 = i2 =  757

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  793

                                      xin(801) = c10*xin(729) + xc00*xin(765) + c01*xin(764)
                                      yin(801) = c10*yin(729) + yc00*yin(765) + c01*yin(764)
                                      zin(801) = c10*zin(729) + zc00*zin(765) + c01*zin(764)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  757
                                      ! i4 = i5 =  793

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  829

                                      xin(837) = c10*xin(765) + xc00*xin(801) + c01*xin(800)
                                      yin(837) = c10*yin(765) + yc00*yin(801) + c01*yin(800)
                                      zin(837) = c10*zin(765) + zc00*zin(801) + c01*zin(800)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  793
                                      ! i4 = i5 =  829

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  838

                                      xin(846) = c10*xin(801) + xc00*xin(837) + c01*xin(836)
                                      yin(846) = c10*yin(801) + yc00*yin(837) + c01*yin(836)
                                      zin(846) = c10*zin(801) + zc00*zin(837) + c01*zin(836)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  829
                                      ! i4 = i5 =  838

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  847

                                      xin(855) = c10*xin(837) + xc00*xin(846) + c01*xin(845)
                                      yin(855) = c10*yin(837) + yc00*yin(846) + c01*yin(845)
                                      zin(855) = c10*zin(837) + zc00*zin(846) + c01*zin(845)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  838
                                      ! i4 = i5 =  847

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  856

                                      xin(864) = c10*xin(846) + xc00*xin(855) + c01*xin(854)
                                      yin(864) = c10*yin(846) + yc00*yin(855) + c01*yin(854)
                                      zin(864) = c10*zin(846) + zc00*zin(855) + c01*zin(854)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  847
                                      ! i4 = i5 =  856

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  856

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  856

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  847

                                      xin(856) = xin(856) + dxij*xin(847)
                                      yin(856) = yin(856) + dyij*yin(847)
                                      zin(856) = zin(856) + dzij*zin(847)

                                      ! i3 = i4 =  847
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  838

                                      xin(847) = xin(847) + dxij*xin(838)
                                      yin(847) = yin(847) + dyij*yin(838)
                                      zin(847) = zin(847) + dzij*zin(838)

                                      ! i3 = i4 =  838
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  829

                                      xin(838) = xin(838) + dxij*xin(829)
                                      yin(838) = yin(838) + dyij*yin(829)
                                      zin(838) = zin(838) + dzij*zin(829)

                                      ! i3 = i4 =  829
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  856

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  847

                                      xin(856) = xin(856) + dxij*xin(847)
                                      yin(856) = yin(856) + dyij*yin(847)
                                      zin(856) = zin(856) + dzij*zin(847)

                                      ! i3 = i4 =  847
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  838

                                      xin(847) = xin(847) + dxij*xin(838)
                                      yin(847) = yin(847) + dyij*yin(838)
                                      zin(847) = zin(847) + dzij*zin(838)

                                      ! i3 = i4 =  838
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  856

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  847

                                      xin(856) = xin(856) + dxij*xin(847)
                                      yin(856) = yin(856) + dyij*yin(847)
                                      zin(856) = zin(856) + dzij*zin(847)

                                      ! i3 = i4 =  847
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  730

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  730

                                      ! do ni = 1,    3

                                      xin(730) = xin(757) + dxij*xin(721)
                                      yin(730) = yin(757) + dyij*yin(721)
                                      zin(730) = zin(757) + dzij*zin(721)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  766

                                      ! ni =    2

                                      xin(766) = xin(793) + dxij*xin(757)
                                      yin(766) = yin(793) + dyij*yin(757)
                                      zin(766) = zin(793) + dzij*zin(757)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  802

                                      ! ni =    3

                                      xin(802) = xin(829) + dxij*xin(793)
                                      yin(802) = yin(829) + dyij*yin(793)
                                      zin(802) = zin(829) + dzij*zin(793)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  838

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  739

                                      ! nj =    2

                                      ! i4 = i3 =  739

                                      ! do ni = 1,    3

                                      xin(739) = xin(766) + dxij*xin(730)
                                      yin(739) = yin(766) + dyij*yin(730)
                                      zin(739) = zin(766) + dzij*zin(730)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  775

                                      ! ni =    2

                                      xin(775) = xin(802) + dxij*xin(766)
                                      yin(775) = yin(802) + dyij*yin(766)
                                      zin(775) = zin(802) + dzij*zin(766)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  811

                                      ! ni =    3

                                      xin(811) = xin(838) + dxij*xin(802)
                                      yin(811) = yin(838) + dyij*yin(802)
                                      zin(811) = zin(838) + dzij*zin(802)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  847

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  748

                                      ! nj =    3

                                      ! i4 = i3 =  748

                                      ! do ni = 1,    3

                                      xin(748) = xin(775) + dxij*xin(739)
                                      yin(748) = yin(775) + dyij*yin(739)
                                      zin(748) = zin(775) + dzij*zin(739)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  784

                                      ! ni =    2

                                      xin(784) = xin(811) + dxij*xin(775)
                                      yin(784) = yin(811) + dyij*yin(775)
                                      zin(784) = zin(811) + dzij*zin(775)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  820

                                      ! ni =    3

                                      xin(820) = xin(847) + dxij*xin(811)
                                      yin(820) = yin(847) + dyij*yin(811)
                                      zin(820) = zin(847) + dzij*zin(811)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  856

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  757

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  859

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  850

                                      xin(859) = xin(859) + dxij*xin(850)
                                      yin(859) = yin(859) + dyij*yin(850)
                                      zin(859) = zin(859) + dzij*zin(850)

                                      ! i3 = i4 =  850
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  841

                                      xin(850) = xin(850) + dxij*xin(841)
                                      yin(850) = yin(850) + dyij*yin(841)
                                      zin(850) = zin(850) + dzij*zin(841)

                                      ! i3 = i4 =  841
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  832

                                      xin(841) = xin(841) + dxij*xin(832)
                                      yin(841) = yin(841) + dyij*yin(832)
                                      zin(841) = zin(841) + dzij*zin(832)

                                      ! i3 = i4 =  832
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  859

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  850

                                      xin(859) = xin(859) + dxij*xin(850)
                                      yin(859) = yin(859) + dyij*yin(850)
                                      zin(859) = zin(859) + dzij*zin(850)

                                      ! i3 = i4 =  850
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  841

                                      xin(850) = xin(850) + dxij*xin(841)
                                      yin(850) = yin(850) + dyij*yin(841)
                                      zin(850) = zin(850) + dzij*zin(841)

                                      ! i3 = i4 =  841
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  859

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  850

                                      xin(859) = xin(859) + dxij*xin(850)
                                      yin(859) = yin(859) + dyij*yin(850)
                                      zin(859) = zin(859) + dzij*zin(850)

                                      ! i3 = i4 =  850
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  733

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  733

                                      ! do ni = 1,    3

                                      xin(733) = xin(760) + dxij*xin(724)
                                      yin(733) = yin(760) + dyij*yin(724)
                                      zin(733) = zin(760) + dzij*zin(724)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  769

                                      ! ni =    2

                                      xin(769) = xin(796) + dxij*xin(760)
                                      yin(769) = yin(796) + dyij*yin(760)
                                      zin(769) = zin(796) + dzij*zin(760)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  805

                                      ! ni =    3

                                      xin(805) = xin(832) + dxij*xin(796)
                                      yin(805) = yin(832) + dyij*yin(796)
                                      zin(805) = zin(832) + dzij*zin(796)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  841

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  742

                                      ! nj =    2

                                      ! i4 = i3 =  742

                                      ! do ni = 1,    3

                                      xin(742) = xin(769) + dxij*xin(733)
                                      yin(742) = yin(769) + dyij*yin(733)
                                      zin(742) = zin(769) + dzij*zin(733)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  778

                                      ! ni =    2

                                      xin(778) = xin(805) + dxij*xin(769)
                                      yin(778) = yin(805) + dyij*yin(769)
                                      zin(778) = zin(805) + dzij*zin(769)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  814

                                      ! ni =    3

                                      xin(814) = xin(841) + dxij*xin(805)
                                      yin(814) = yin(841) + dyij*yin(805)
                                      zin(814) = zin(841) + dzij*zin(805)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  850

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  751

                                      ! nj =    3

                                      ! i4 = i3 =  751

                                      ! do ni = 1,    3

                                      xin(751) = xin(778) + dxij*xin(742)
                                      yin(751) = yin(778) + dyij*yin(742)
                                      zin(751) = zin(778) + dzij*zin(742)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  787

                                      ! ni =    2

                                      xin(787) = xin(814) + dxij*xin(778)
                                      yin(787) = yin(814) + dyij*yin(778)
                                      zin(787) = zin(814) + dzij*zin(778)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  823

                                      ! ni =    3

                                      xin(823) = xin(850) + dxij*xin(814)
                                      yin(823) = yin(850) + dyij*yin(814)
                                      zin(823) = zin(850) + dzij*zin(814)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  859

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  760

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  862

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  853

                                      xin(862) = xin(862) + dxij*xin(853)
                                      yin(862) = yin(862) + dyij*yin(853)
                                      zin(862) = zin(862) + dzij*zin(853)

                                      ! i3 = i4 =  853
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  844

                                      xin(853) = xin(853) + dxij*xin(844)
                                      yin(853) = yin(853) + dyij*yin(844)
                                      zin(853) = zin(853) + dzij*zin(844)

                                      ! i3 = i4 =  844
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  835

                                      xin(844) = xin(844) + dxij*xin(835)
                                      yin(844) = yin(844) + dyij*yin(835)
                                      zin(844) = zin(844) + dzij*zin(835)

                                      ! i3 = i4 =  835
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  862

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  853

                                      xin(862) = xin(862) + dxij*xin(853)
                                      yin(862) = yin(862) + dyij*yin(853)
                                      zin(862) = zin(862) + dzij*zin(853)

                                      ! i3 = i4 =  853
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  844

                                      xin(853) = xin(853) + dxij*xin(844)
                                      yin(853) = yin(853) + dyij*yin(844)
                                      zin(853) = zin(853) + dzij*zin(844)

                                      ! i3 = i4 =  844
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  862

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  853

                                      xin(862) = xin(862) + dxij*xin(853)
                                      yin(862) = yin(862) + dyij*yin(853)
                                      zin(862) = zin(862) + dzij*zin(853)

                                      ! i3 = i4 =  853
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  736

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  736

                                      ! do ni = 1,    3

                                      xin(736) = xin(763) + dxij*xin(727)
                                      yin(736) = yin(763) + dyij*yin(727)
                                      zin(736) = zin(763) + dzij*zin(727)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  772

                                      ! ni =    2

                                      xin(772) = xin(799) + dxij*xin(763)
                                      yin(772) = yin(799) + dyij*yin(763)
                                      zin(772) = zin(799) + dzij*zin(763)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  808

                                      ! ni =    3

                                      xin(808) = xin(835) + dxij*xin(799)
                                      yin(808) = yin(835) + dyij*yin(799)
                                      zin(808) = zin(835) + dzij*zin(799)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  844

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  745

                                      ! nj =    2

                                      ! i4 = i3 =  745

                                      ! do ni = 1,    3

                                      xin(745) = xin(772) + dxij*xin(736)
                                      yin(745) = yin(772) + dyij*yin(736)
                                      zin(745) = zin(772) + dzij*zin(736)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  781

                                      ! ni =    2

                                      xin(781) = xin(808) + dxij*xin(772)
                                      yin(781) = yin(808) + dyij*yin(772)
                                      zin(781) = zin(808) + dzij*zin(772)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  817

                                      ! ni =    3

                                      xin(817) = xin(844) + dxij*xin(808)
                                      yin(817) = yin(844) + dyij*yin(808)
                                      zin(817) = zin(844) + dzij*zin(808)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  853

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  754

                                      ! nj =    3

                                      ! i4 = i3 =  754

                                      ! do ni = 1,    3

                                      xin(754) = xin(781) + dxij*xin(745)
                                      yin(754) = yin(781) + dyij*yin(745)
                                      zin(754) = zin(781) + dzij*zin(745)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  790

                                      ! ni =    2

                                      xin(790) = xin(817) + dxij*xin(781)
                                      yin(790) = yin(817) + dyij*yin(781)
                                      zin(790) = zin(817) + dzij*zin(781)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  826

                                      ! ni =    3

                                      xin(826) = xin(853) + dxij*xin(817)
                                      yin(826) = yin(853) + dyij*yin(817)
                                      zin(826) = zin(853) + dzij*zin(817)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  862

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  763

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  863

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  854

                                      xin(863) = xin(863) + dxij*xin(854)
                                      yin(863) = yin(863) + dyij*yin(854)
                                      zin(863) = zin(863) + dzij*zin(854)

                                      ! i3 = i4 =  854
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  845

                                      xin(854) = xin(854) + dxij*xin(845)
                                      yin(854) = yin(854) + dyij*yin(845)
                                      zin(854) = zin(854) + dzij*zin(845)

                                      ! i3 = i4 =  845
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  836

                                      xin(845) = xin(845) + dxij*xin(836)
                                      yin(845) = yin(845) + dyij*yin(836)
                                      zin(845) = zin(845) + dzij*zin(836)

                                      ! i3 = i4 =  836
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  863

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  854

                                      xin(863) = xin(863) + dxij*xin(854)
                                      yin(863) = yin(863) + dyij*yin(854)
                                      zin(863) = zin(863) + dzij*zin(854)

                                      ! i3 = i4 =  854
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  845

                                      xin(854) = xin(854) + dxij*xin(845)
                                      yin(854) = yin(854) + dyij*yin(845)
                                      zin(854) = zin(854) + dzij*zin(845)

                                      ! i3 = i4 =  845
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  863

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  854

                                      xin(863) = xin(863) + dxij*xin(854)
                                      yin(863) = yin(863) + dyij*yin(854)
                                      zin(863) = zin(863) + dzij*zin(854)

                                      ! i3 = i4 =  854
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  737

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  737

                                      ! do ni = 1,    3

                                      xin(737) = xin(764) + dxij*xin(728)
                                      yin(737) = yin(764) + dyij*yin(728)
                                      zin(737) = zin(764) + dzij*zin(728)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  773

                                      ! ni =    2

                                      xin(773) = xin(800) + dxij*xin(764)
                                      yin(773) = yin(800) + dyij*yin(764)
                                      zin(773) = zin(800) + dzij*zin(764)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  809

                                      ! ni =    3

                                      xin(809) = xin(836) + dxij*xin(800)
                                      yin(809) = yin(836) + dyij*yin(800)
                                      zin(809) = zin(836) + dzij*zin(800)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  845

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  746

                                      ! nj =    2

                                      ! i4 = i3 =  746

                                      ! do ni = 1,    3

                                      xin(746) = xin(773) + dxij*xin(737)
                                      yin(746) = yin(773) + dyij*yin(737)
                                      zin(746) = zin(773) + dzij*zin(737)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  782

                                      ! ni =    2

                                      xin(782) = xin(809) + dxij*xin(773)
                                      yin(782) = yin(809) + dyij*yin(773)
                                      zin(782) = zin(809) + dzij*zin(773)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  818

                                      ! ni =    3

                                      xin(818) = xin(845) + dxij*xin(809)
                                      yin(818) = yin(845) + dyij*yin(809)
                                      zin(818) = zin(845) + dzij*zin(809)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  854

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  755

                                      ! nj =    3

                                      ! i4 = i3 =  755

                                      ! do ni = 1,    3

                                      xin(755) = xin(782) + dxij*xin(746)
                                      yin(755) = yin(782) + dyij*yin(746)
                                      zin(755) = zin(782) + dzij*zin(746)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  791

                                      ! ni =    2

                                      xin(791) = xin(818) + dxij*xin(782)
                                      yin(791) = yin(818) + dyij*yin(782)
                                      zin(791) = zin(818) + dzij*zin(782)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  827

                                      ! ni =    3

                                      xin(827) = xin(854) + dxij*xin(818)
                                      yin(827) = yin(854) + dyij*yin(818)
                                      zin(827) = zin(854) + dzij*zin(818)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  863

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  764

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  864

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  855

                                      xin(864) = xin(864) + dxij*xin(855)
                                      yin(864) = yin(864) + dyij*yin(855)
                                      zin(864) = zin(864) + dzij*zin(855)

                                      ! i3 = i4 =  855
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  846

                                      xin(855) = xin(855) + dxij*xin(846)
                                      yin(855) = yin(855) + dyij*yin(846)
                                      zin(855) = zin(855) + dzij*zin(846)

                                      ! i3 = i4 =  846
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  837

                                      xin(846) = xin(846) + dxij*xin(837)
                                      yin(846) = yin(846) + dyij*yin(837)
                                      zin(846) = zin(846) + dzij*zin(837)

                                      ! i3 = i4 =  837
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  864

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  855

                                      xin(864) = xin(864) + dxij*xin(855)
                                      yin(864) = yin(864) + dyij*yin(855)
                                      zin(864) = zin(864) + dzij*zin(855)

                                      ! i3 = i4 =  855
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  846

                                      xin(855) = xin(855) + dxij*xin(846)
                                      yin(855) = yin(855) + dyij*yin(846)
                                      zin(855) = zin(855) + dzij*zin(846)

                                      ! i3 = i4 =  846
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  864

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  855

                                      xin(864) = xin(864) + dxij*xin(855)
                                      yin(864) = yin(864) + dyij*yin(855)
                                      zin(864) = zin(864) + dzij*zin(855)

                                      ! i3 = i4 =  855
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  738

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  738

                                      ! do ni = 1,    3

                                      xin(738) = xin(765) + dxij*xin(729)
                                      yin(738) = yin(765) + dyij*yin(729)
                                      zin(738) = zin(765) + dzij*zin(729)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  774

                                      ! ni =    2

                                      xin(774) = xin(801) + dxij*xin(765)
                                      yin(774) = yin(801) + dyij*yin(765)
                                      zin(774) = zin(801) + dzij*zin(765)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  810

                                      ! ni =    3

                                      xin(810) = xin(837) + dxij*xin(801)
                                      yin(810) = yin(837) + dyij*yin(801)
                                      zin(810) = zin(837) + dzij*zin(801)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  846

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  747

                                      ! nj =    2

                                      ! i4 = i3 =  747

                                      ! do ni = 1,    3

                                      xin(747) = xin(774) + dxij*xin(738)
                                      yin(747) = yin(774) + dyij*yin(738)
                                      zin(747) = zin(774) + dzij*zin(738)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  783

                                      ! ni =    2

                                      xin(783) = xin(810) + dxij*xin(774)
                                      yin(783) = yin(810) + dyij*yin(774)
                                      zin(783) = zin(810) + dzij*zin(774)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  819

                                      ! ni =    3

                                      xin(819) = xin(846) + dxij*xin(810)
                                      yin(819) = yin(846) + dyij*yin(810)
                                      zin(819) = zin(846) + dzij*zin(810)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  855

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  756

                                      ! nj =    3

                                      ! i4 = i3 =  756

                                      ! do ni = 1,    3

                                      xin(756) = xin(783) + dxij*xin(747)
                                      yin(756) = yin(783) + dyij*yin(747)
                                      zin(756) = zin(783) + dzij*zin(747)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  792

                                      ! ni =    2

                                      xin(792) = xin(819) + dxij*xin(783)
                                      yin(792) = yin(819) + dyij*yin(783)
                                      zin(792) = zin(819) + dzij*zin(783)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  828

                                      ! ni =    3

                                      xin(828) = xin(855) + dxij*xin(819)
                                      yin(828) = yin(855) + dyij*yin(819)
                                      zin(828) = zin(855) + dzij*zin(819)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  864

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  765

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  721

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  729

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  728

                                      xin(729) = xin(729) + dxkl*xin(728)
                                      yin(729) = yin(729) + dykl*yin(728)
                                      zin(729) = zin(729) + dzkl*zin(728)

                                      ! i3 = i4 =  728
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  727

                                      xin(728) = xin(728) + dxkl*xin(727)
                                      yin(728) = yin(728) + dykl*yin(727)
                                      zin(728) = zin(728) + dzkl*zin(727)

                                      ! i3 = i4 =  727
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  729

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  728

                                      xin(729) = xin(729) + dxkl*xin(728)
                                      yin(729) = yin(729) + dykl*yin(728)
                                      zin(729) = zin(729) + dzkl*zin(728)

                                      ! i3 = i4 =  728
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  722

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  722

                                      ! do nk = 1,    2

                                      xin(722) = xin(724) + dxkl*xin(721)
                                      yin(722) = yin(724) + dykl*yin(721)
                                      zin(722) = zin(724) + dzkl*zin(721)
                                      ! i4 = i4 + lang+1 =  725

                                      ! nk =    2

                                      xin(725) = xin(727) + dxkl*xin(724)
                                      yin(725) = yin(727) + dykl*yin(724)
                                      zin(725) = zin(727) + dzkl*zin(724)
                                      ! i4 = i4 + lang+1 =  728

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  723

                                      ! nl =    2

                                      ! i4 = i3 =  723

                                      ! do nk = 1,    2

                                      xin(723) = xin(725) + dxkl*xin(722)
                                      yin(723) = yin(725) + dykl*yin(722)
                                      zin(723) = zin(725) + dzkl*zin(722)
                                      ! i4 = i4 + lang+1 =  726

                                      ! nk =    2

                                      xin(726) = xin(728) + dxkl*xin(725)
                                      yin(726) = yin(728) + dykl*yin(725)
                                      zin(726) = zin(728) + dzkl*zin(725)
                                      ! i4 = i4 + lang+1 =  729

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  724

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  730

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  738

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  737

                                      xin(738) = xin(738) + dxkl*xin(737)
                                      yin(738) = yin(738) + dykl*yin(737)
                                      zin(738) = zin(738) + dzkl*zin(737)

                                      ! i3 = i4 =  737
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  736

                                      xin(737) = xin(737) + dxkl*xin(736)
                                      yin(737) = yin(737) + dykl*yin(736)
                                      zin(737) = zin(737) + dzkl*zin(736)

                                      ! i3 = i4 =  736
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  738

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  737

                                      xin(738) = xin(738) + dxkl*xin(737)
                                      yin(738) = yin(738) + dykl*yin(737)
                                      zin(738) = zin(738) + dzkl*zin(737)

                                      ! i3 = i4 =  737
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  731

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  731

                                      ! do nk = 1,    2

                                      xin(731) = xin(733) + dxkl*xin(730)
                                      yin(731) = yin(733) + dykl*yin(730)
                                      zin(731) = zin(733) + dzkl*zin(730)
                                      ! i4 = i4 + lang+1 =  734

                                      ! nk =    2

                                      xin(734) = xin(736) + dxkl*xin(733)
                                      yin(734) = yin(736) + dykl*yin(733)
                                      zin(734) = zin(736) + dzkl*zin(733)
                                      ! i4 = i4 + lang+1 =  737

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  732

                                      ! nl =    2

                                      ! i4 = i3 =  732

                                      ! do nk = 1,    2

                                      xin(732) = xin(734) + dxkl*xin(731)
                                      yin(732) = yin(734) + dykl*yin(731)
                                      zin(732) = zin(734) + dzkl*zin(731)
                                      ! i4 = i4 + lang+1 =  735

                                      ! nk =    2

                                      xin(735) = xin(737) + dxkl*xin(734)
                                      yin(735) = yin(737) + dykl*yin(734)
                                      zin(735) = zin(737) + dzkl*zin(734)
                                      ! i4 = i4 + lang+1 =  738

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  733

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  739

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  747

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  746

                                      xin(747) = xin(747) + dxkl*xin(746)
                                      yin(747) = yin(747) + dykl*yin(746)
                                      zin(747) = zin(747) + dzkl*zin(746)

                                      ! i3 = i4 =  746
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  745

                                      xin(746) = xin(746) + dxkl*xin(745)
                                      yin(746) = yin(746) + dykl*yin(745)
                                      zin(746) = zin(746) + dzkl*zin(745)

                                      ! i3 = i4 =  745
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  747

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  746

                                      xin(747) = xin(747) + dxkl*xin(746)
                                      yin(747) = yin(747) + dykl*yin(746)
                                      zin(747) = zin(747) + dzkl*zin(746)

                                      ! i3 = i4 =  746
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  740

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  740

                                      ! do nk = 1,    2

                                      xin(740) = xin(742) + dxkl*xin(739)
                                      yin(740) = yin(742) + dykl*yin(739)
                                      zin(740) = zin(742) + dzkl*zin(739)
                                      ! i4 = i4 + lang+1 =  743

                                      ! nk =    2

                                      xin(743) = xin(745) + dxkl*xin(742)
                                      yin(743) = yin(745) + dykl*yin(742)
                                      zin(743) = zin(745) + dzkl*zin(742)
                                      ! i4 = i4 + lang+1 =  746

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  741

                                      ! nl =    2

                                      ! i4 = i3 =  741

                                      ! do nk = 1,    2

                                      xin(741) = xin(743) + dxkl*xin(740)
                                      yin(741) = yin(743) + dykl*yin(740)
                                      zin(741) = zin(743) + dzkl*zin(740)
                                      ! i4 = i4 + lang+1 =  744

                                      ! nk =    2

                                      xin(744) = xin(746) + dxkl*xin(743)
                                      yin(744) = yin(746) + dykl*yin(743)
                                      zin(744) = zin(746) + dzkl*zin(743)
                                      ! i4 = i4 + lang+1 =  747

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  742

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  748

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  756

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  755

                                      xin(756) = xin(756) + dxkl*xin(755)
                                      yin(756) = yin(756) + dykl*yin(755)
                                      zin(756) = zin(756) + dzkl*zin(755)

                                      ! i3 = i4 =  755
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  754

                                      xin(755) = xin(755) + dxkl*xin(754)
                                      yin(755) = yin(755) + dykl*yin(754)
                                      zin(755) = zin(755) + dzkl*zin(754)

                                      ! i3 = i4 =  754
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  756

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  755

                                      xin(756) = xin(756) + dxkl*xin(755)
                                      yin(756) = yin(756) + dykl*yin(755)
                                      zin(756) = zin(756) + dzkl*zin(755)

                                      ! i3 = i4 =  755
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  749

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  749

                                      ! do nk = 1,    2

                                      xin(749) = xin(751) + dxkl*xin(748)
                                      yin(749) = yin(751) + dykl*yin(748)
                                      zin(749) = zin(751) + dzkl*zin(748)
                                      ! i4 = i4 + lang+1 =  752

                                      ! nk =    2

                                      xin(752) = xin(754) + dxkl*xin(751)
                                      yin(752) = yin(754) + dykl*yin(751)
                                      zin(752) = zin(754) + dzkl*zin(751)
                                      ! i4 = i4 + lang+1 =  755

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  750

                                      ! nl =    2

                                      ! i4 = i3 =  750

                                      ! do nk = 1,    2

                                      xin(750) = xin(752) + dxkl*xin(749)
                                      yin(750) = yin(752) + dykl*yin(749)
                                      zin(750) = zin(752) + dzkl*zin(749)
                                      ! i4 = i4 + lang+1 =  753

                                      ! nk =    2

                                      xin(753) = xin(755) + dxkl*xin(752)
                                      yin(753) = yin(755) + dykl*yin(752)
                                      zin(753) = zin(755) + dzkl*zin(752)
                                      ! i4 = i4 + lang+1 =  756

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  751

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  757

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  757

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  765

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  764

                                      xin(765) = xin(765) + dxkl*xin(764)
                                      yin(765) = yin(765) + dykl*yin(764)
                                      zin(765) = zin(765) + dzkl*zin(764)

                                      ! i3 = i4 =  764
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  763

                                      xin(764) = xin(764) + dxkl*xin(763)
                                      yin(764) = yin(764) + dykl*yin(763)
                                      zin(764) = zin(764) + dzkl*zin(763)

                                      ! i3 = i4 =  763
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  765

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  764

                                      xin(765) = xin(765) + dxkl*xin(764)
                                      yin(765) = yin(765) + dykl*yin(764)
                                      zin(765) = zin(765) + dzkl*zin(764)

                                      ! i3 = i4 =  764
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  758

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  758

                                      ! do nk = 1,    2

                                      xin(758) = xin(760) + dxkl*xin(757)
                                      yin(758) = yin(760) + dykl*yin(757)
                                      zin(758) = zin(760) + dzkl*zin(757)
                                      ! i4 = i4 + lang+1 =  761

                                      ! nk =    2

                                      xin(761) = xin(763) + dxkl*xin(760)
                                      yin(761) = yin(763) + dykl*yin(760)
                                      zin(761) = zin(763) + dzkl*zin(760)
                                      ! i4 = i4 + lang+1 =  764

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  759

                                      ! nl =    2

                                      ! i4 = i3 =  759

                                      ! do nk = 1,    2

                                      xin(759) = xin(761) + dxkl*xin(758)
                                      yin(759) = yin(761) + dykl*yin(758)
                                      zin(759) = zin(761) + dzkl*zin(758)
                                      ! i4 = i4 + lang+1 =  762

                                      ! nk =    2

                                      xin(762) = xin(764) + dxkl*xin(761)
                                      yin(762) = yin(764) + dykl*yin(761)
                                      zin(762) = zin(764) + dzkl*zin(761)
                                      ! i4 = i4 + lang+1 =  765

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  760

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  766

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  774

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  773

                                      xin(774) = xin(774) + dxkl*xin(773)
                                      yin(774) = yin(774) + dykl*yin(773)
                                      zin(774) = zin(774) + dzkl*zin(773)

                                      ! i3 = i4 =  773
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  772

                                      xin(773) = xin(773) + dxkl*xin(772)
                                      yin(773) = yin(773) + dykl*yin(772)
                                      zin(773) = zin(773) + dzkl*zin(772)

                                      ! i3 = i4 =  772
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  774

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  773

                                      xin(774) = xin(774) + dxkl*xin(773)
                                      yin(774) = yin(774) + dykl*yin(773)
                                      zin(774) = zin(774) + dzkl*zin(773)

                                      ! i3 = i4 =  773
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  767

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  767

                                      ! do nk = 1,    2

                                      xin(767) = xin(769) + dxkl*xin(766)
                                      yin(767) = yin(769) + dykl*yin(766)
                                      zin(767) = zin(769) + dzkl*zin(766)
                                      ! i4 = i4 + lang+1 =  770

                                      ! nk =    2

                                      xin(770) = xin(772) + dxkl*xin(769)
                                      yin(770) = yin(772) + dykl*yin(769)
                                      zin(770) = zin(772) + dzkl*zin(769)
                                      ! i4 = i4 + lang+1 =  773

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  768

                                      ! nl =    2

                                      ! i4 = i3 =  768

                                      ! do nk = 1,    2

                                      xin(768) = xin(770) + dxkl*xin(767)
                                      yin(768) = yin(770) + dykl*yin(767)
                                      zin(768) = zin(770) + dzkl*zin(767)
                                      ! i4 = i4 + lang+1 =  771

                                      ! nk =    2

                                      xin(771) = xin(773) + dxkl*xin(770)
                                      yin(771) = yin(773) + dykl*yin(770)
                                      zin(771) = zin(773) + dzkl*zin(770)
                                      ! i4 = i4 + lang+1 =  774

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  769

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  775

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  783

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  782

                                      xin(783) = xin(783) + dxkl*xin(782)
                                      yin(783) = yin(783) + dykl*yin(782)
                                      zin(783) = zin(783) + dzkl*zin(782)

                                      ! i3 = i4 =  782
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  781

                                      xin(782) = xin(782) + dxkl*xin(781)
                                      yin(782) = yin(782) + dykl*yin(781)
                                      zin(782) = zin(782) + dzkl*zin(781)

                                      ! i3 = i4 =  781
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  783

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  782

                                      xin(783) = xin(783) + dxkl*xin(782)
                                      yin(783) = yin(783) + dykl*yin(782)
                                      zin(783) = zin(783) + dzkl*zin(782)

                                      ! i3 = i4 =  782
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  776

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  776

                                      ! do nk = 1,    2

                                      xin(776) = xin(778) + dxkl*xin(775)
                                      yin(776) = yin(778) + dykl*yin(775)
                                      zin(776) = zin(778) + dzkl*zin(775)
                                      ! i4 = i4 + lang+1 =  779

                                      ! nk =    2

                                      xin(779) = xin(781) + dxkl*xin(778)
                                      yin(779) = yin(781) + dykl*yin(778)
                                      zin(779) = zin(781) + dzkl*zin(778)
                                      ! i4 = i4 + lang+1 =  782

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  777

                                      ! nl =    2

                                      ! i4 = i3 =  777

                                      ! do nk = 1,    2

                                      xin(777) = xin(779) + dxkl*xin(776)
                                      yin(777) = yin(779) + dykl*yin(776)
                                      zin(777) = zin(779) + dzkl*zin(776)
                                      ! i4 = i4 + lang+1 =  780

                                      ! nk =    2

                                      xin(780) = xin(782) + dxkl*xin(779)
                                      yin(780) = yin(782) + dykl*yin(779)
                                      zin(780) = zin(782) + dzkl*zin(779)
                                      ! i4 = i4 + lang+1 =  783

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  778

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  784

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  792

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  791

                                      xin(792) = xin(792) + dxkl*xin(791)
                                      yin(792) = yin(792) + dykl*yin(791)
                                      zin(792) = zin(792) + dzkl*zin(791)

                                      ! i3 = i4 =  791
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  790

                                      xin(791) = xin(791) + dxkl*xin(790)
                                      yin(791) = yin(791) + dykl*yin(790)
                                      zin(791) = zin(791) + dzkl*zin(790)

                                      ! i3 = i4 =  790
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  792

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  791

                                      xin(792) = xin(792) + dxkl*xin(791)
                                      yin(792) = yin(792) + dykl*yin(791)
                                      zin(792) = zin(792) + dzkl*zin(791)

                                      ! i3 = i4 =  791
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  785

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  785

                                      ! do nk = 1,    2

                                      xin(785) = xin(787) + dxkl*xin(784)
                                      yin(785) = yin(787) + dykl*yin(784)
                                      zin(785) = zin(787) + dzkl*zin(784)
                                      ! i4 = i4 + lang+1 =  788

                                      ! nk =    2

                                      xin(788) = xin(790) + dxkl*xin(787)
                                      yin(788) = yin(790) + dykl*yin(787)
                                      zin(788) = zin(790) + dzkl*zin(787)
                                      ! i4 = i4 + lang+1 =  791

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  786

                                      ! nl =    2

                                      ! i4 = i3 =  786

                                      ! do nk = 1,    2

                                      xin(786) = xin(788) + dxkl*xin(785)
                                      yin(786) = yin(788) + dykl*yin(785)
                                      zin(786) = zin(788) + dzkl*zin(785)
                                      ! i4 = i4 + lang+1 =  789

                                      ! nk =    2

                                      xin(789) = xin(791) + dxkl*xin(788)
                                      yin(789) = yin(791) + dykl*yin(788)
                                      zin(789) = zin(791) + dzkl*zin(788)
                                      ! i4 = i4 + lang+1 =  792

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  787

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  793

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  793

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  801

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  800

                                      xin(801) = xin(801) + dxkl*xin(800)
                                      yin(801) = yin(801) + dykl*yin(800)
                                      zin(801) = zin(801) + dzkl*zin(800)

                                      ! i3 = i4 =  800
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  799

                                      xin(800) = xin(800) + dxkl*xin(799)
                                      yin(800) = yin(800) + dykl*yin(799)
                                      zin(800) = zin(800) + dzkl*zin(799)

                                      ! i3 = i4 =  799
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  801

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  800

                                      xin(801) = xin(801) + dxkl*xin(800)
                                      yin(801) = yin(801) + dykl*yin(800)
                                      zin(801) = zin(801) + dzkl*zin(800)

                                      ! i3 = i4 =  800
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  794

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  794

                                      ! do nk = 1,    2

                                      xin(794) = xin(796) + dxkl*xin(793)
                                      yin(794) = yin(796) + dykl*yin(793)
                                      zin(794) = zin(796) + dzkl*zin(793)
                                      ! i4 = i4 + lang+1 =  797

                                      ! nk =    2

                                      xin(797) = xin(799) + dxkl*xin(796)
                                      yin(797) = yin(799) + dykl*yin(796)
                                      zin(797) = zin(799) + dzkl*zin(796)
                                      ! i4 = i4 + lang+1 =  800

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  795

                                      ! nl =    2

                                      ! i4 = i3 =  795

                                      ! do nk = 1,    2

                                      xin(795) = xin(797) + dxkl*xin(794)
                                      yin(795) = yin(797) + dykl*yin(794)
                                      zin(795) = zin(797) + dzkl*zin(794)
                                      ! i4 = i4 + lang+1 =  798

                                      ! nk =    2

                                      xin(798) = xin(800) + dxkl*xin(797)
                                      yin(798) = yin(800) + dykl*yin(797)
                                      zin(798) = zin(800) + dzkl*zin(797)
                                      ! i4 = i4 + lang+1 =  801

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  796

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  802

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  810

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  809

                                      xin(810) = xin(810) + dxkl*xin(809)
                                      yin(810) = yin(810) + dykl*yin(809)
                                      zin(810) = zin(810) + dzkl*zin(809)

                                      ! i3 = i4 =  809
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  808

                                      xin(809) = xin(809) + dxkl*xin(808)
                                      yin(809) = yin(809) + dykl*yin(808)
                                      zin(809) = zin(809) + dzkl*zin(808)

                                      ! i3 = i4 =  808
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  810

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  809

                                      xin(810) = xin(810) + dxkl*xin(809)
                                      yin(810) = yin(810) + dykl*yin(809)
                                      zin(810) = zin(810) + dzkl*zin(809)

                                      ! i3 = i4 =  809
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  803

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  803

                                      ! do nk = 1,    2

                                      xin(803) = xin(805) + dxkl*xin(802)
                                      yin(803) = yin(805) + dykl*yin(802)
                                      zin(803) = zin(805) + dzkl*zin(802)
                                      ! i4 = i4 + lang+1 =  806

                                      ! nk =    2

                                      xin(806) = xin(808) + dxkl*xin(805)
                                      yin(806) = yin(808) + dykl*yin(805)
                                      zin(806) = zin(808) + dzkl*zin(805)
                                      ! i4 = i4 + lang+1 =  809

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  804

                                      ! nl =    2

                                      ! i4 = i3 =  804

                                      ! do nk = 1,    2

                                      xin(804) = xin(806) + dxkl*xin(803)
                                      yin(804) = yin(806) + dykl*yin(803)
                                      zin(804) = zin(806) + dzkl*zin(803)
                                      ! i4 = i4 + lang+1 =  807

                                      ! nk =    2

                                      xin(807) = xin(809) + dxkl*xin(806)
                                      yin(807) = yin(809) + dykl*yin(806)
                                      zin(807) = zin(809) + dzkl*zin(806)
                                      ! i4 = i4 + lang+1 =  810

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  805

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  811

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  819

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  818

                                      xin(819) = xin(819) + dxkl*xin(818)
                                      yin(819) = yin(819) + dykl*yin(818)
                                      zin(819) = zin(819) + dzkl*zin(818)

                                      ! i3 = i4 =  818
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  817

                                      xin(818) = xin(818) + dxkl*xin(817)
                                      yin(818) = yin(818) + dykl*yin(817)
                                      zin(818) = zin(818) + dzkl*zin(817)

                                      ! i3 = i4 =  817
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  819

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  818

                                      xin(819) = xin(819) + dxkl*xin(818)
                                      yin(819) = yin(819) + dykl*yin(818)
                                      zin(819) = zin(819) + dzkl*zin(818)

                                      ! i3 = i4 =  818
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  812

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  812

                                      ! do nk = 1,    2

                                      xin(812) = xin(814) + dxkl*xin(811)
                                      yin(812) = yin(814) + dykl*yin(811)
                                      zin(812) = zin(814) + dzkl*zin(811)
                                      ! i4 = i4 + lang+1 =  815

                                      ! nk =    2

                                      xin(815) = xin(817) + dxkl*xin(814)
                                      yin(815) = yin(817) + dykl*yin(814)
                                      zin(815) = zin(817) + dzkl*zin(814)
                                      ! i4 = i4 + lang+1 =  818

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  813

                                      ! nl =    2

                                      ! i4 = i3 =  813

                                      ! do nk = 1,    2

                                      xin(813) = xin(815) + dxkl*xin(812)
                                      yin(813) = yin(815) + dykl*yin(812)
                                      zin(813) = zin(815) + dzkl*zin(812)
                                      ! i4 = i4 + lang+1 =  816

                                      ! nk =    2

                                      xin(816) = xin(818) + dxkl*xin(815)
                                      yin(816) = yin(818) + dykl*yin(815)
                                      zin(816) = zin(818) + dzkl*zin(815)
                                      ! i4 = i4 + lang+1 =  819

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  814

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  820

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  828

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  827

                                      xin(828) = xin(828) + dxkl*xin(827)
                                      yin(828) = yin(828) + dykl*yin(827)
                                      zin(828) = zin(828) + dzkl*zin(827)

                                      ! i3 = i4 =  827
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  826

                                      xin(827) = xin(827) + dxkl*xin(826)
                                      yin(827) = yin(827) + dykl*yin(826)
                                      zin(827) = zin(827) + dzkl*zin(826)

                                      ! i3 = i4 =  826
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  828

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  827

                                      xin(828) = xin(828) + dxkl*xin(827)
                                      yin(828) = yin(828) + dykl*yin(827)
                                      zin(828) = zin(828) + dzkl*zin(827)

                                      ! i3 = i4 =  827
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  821

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  821

                                      ! do nk = 1,    2

                                      xin(821) = xin(823) + dxkl*xin(820)
                                      yin(821) = yin(823) + dykl*yin(820)
                                      zin(821) = zin(823) + dzkl*zin(820)
                                      ! i4 = i4 + lang+1 =  824

                                      ! nk =    2

                                      xin(824) = xin(826) + dxkl*xin(823)
                                      yin(824) = yin(826) + dykl*yin(823)
                                      zin(824) = zin(826) + dzkl*zin(823)
                                      ! i4 = i4 + lang+1 =  827

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  822

                                      ! nl =    2

                                      ! i4 = i3 =  822

                                      ! do nk = 1,    2

                                      xin(822) = xin(824) + dxkl*xin(821)
                                      yin(822) = yin(824) + dykl*yin(821)
                                      zin(822) = zin(824) + dzkl*zin(821)
                                      ! i4 = i4 + lang+1 =  825

                                      ! nk =    2

                                      xin(825) = xin(827) + dxkl*xin(824)
                                      yin(825) = yin(827) + dykl*yin(824)
                                      zin(825) = zin(827) + dzkl*zin(824)
                                      ! i4 = i4 + lang+1 =  828

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  823

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  829

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  829

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  837

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  836

                                      xin(837) = xin(837) + dxkl*xin(836)
                                      yin(837) = yin(837) + dykl*yin(836)
                                      zin(837) = zin(837) + dzkl*zin(836)

                                      ! i3 = i4 =  836
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  835

                                      xin(836) = xin(836) + dxkl*xin(835)
                                      yin(836) = yin(836) + dykl*yin(835)
                                      zin(836) = zin(836) + dzkl*zin(835)

                                      ! i3 = i4 =  835
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  837

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  836

                                      xin(837) = xin(837) + dxkl*xin(836)
                                      yin(837) = yin(837) + dykl*yin(836)
                                      zin(837) = zin(837) + dzkl*zin(836)

                                      ! i3 = i4 =  836
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  830

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  830

                                      ! do nk = 1,    2

                                      xin(830) = xin(832) + dxkl*xin(829)
                                      yin(830) = yin(832) + dykl*yin(829)
                                      zin(830) = zin(832) + dzkl*zin(829)
                                      ! i4 = i4 + lang+1 =  833

                                      ! nk =    2

                                      xin(833) = xin(835) + dxkl*xin(832)
                                      yin(833) = yin(835) + dykl*yin(832)
                                      zin(833) = zin(835) + dzkl*zin(832)
                                      ! i4 = i4 + lang+1 =  836

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  831

                                      ! nl =    2

                                      ! i4 = i3 =  831

                                      ! do nk = 1,    2

                                      xin(831) = xin(833) + dxkl*xin(830)
                                      yin(831) = yin(833) + dykl*yin(830)
                                      zin(831) = zin(833) + dzkl*zin(830)
                                      ! i4 = i4 + lang+1 =  834

                                      ! nk =    2

                                      xin(834) = xin(836) + dxkl*xin(833)
                                      yin(834) = yin(836) + dykl*yin(833)
                                      zin(834) = zin(836) + dzkl*zin(833)
                                      ! i4 = i4 + lang+1 =  837

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  832

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  838

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  846

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  845

                                      xin(846) = xin(846) + dxkl*xin(845)
                                      yin(846) = yin(846) + dykl*yin(845)
                                      zin(846) = zin(846) + dzkl*zin(845)

                                      ! i3 = i4 =  845
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  844

                                      xin(845) = xin(845) + dxkl*xin(844)
                                      yin(845) = yin(845) + dykl*yin(844)
                                      zin(845) = zin(845) + dzkl*zin(844)

                                      ! i3 = i4 =  844
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  846

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  845

                                      xin(846) = xin(846) + dxkl*xin(845)
                                      yin(846) = yin(846) + dykl*yin(845)
                                      zin(846) = zin(846) + dzkl*zin(845)

                                      ! i3 = i4 =  845
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  839

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  839

                                      ! do nk = 1,    2

                                      xin(839) = xin(841) + dxkl*xin(838)
                                      yin(839) = yin(841) + dykl*yin(838)
                                      zin(839) = zin(841) + dzkl*zin(838)
                                      ! i4 = i4 + lang+1 =  842

                                      ! nk =    2

                                      xin(842) = xin(844) + dxkl*xin(841)
                                      yin(842) = yin(844) + dykl*yin(841)
                                      zin(842) = zin(844) + dzkl*zin(841)
                                      ! i4 = i4 + lang+1 =  845

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  840

                                      ! nl =    2

                                      ! i4 = i3 =  840

                                      ! do nk = 1,    2

                                      xin(840) = xin(842) + dxkl*xin(839)
                                      yin(840) = yin(842) + dykl*yin(839)
                                      zin(840) = zin(842) + dzkl*zin(839)
                                      ! i4 = i4 + lang+1 =  843

                                      ! nk =    2

                                      xin(843) = xin(845) + dxkl*xin(842)
                                      yin(843) = yin(845) + dykl*yin(842)
                                      zin(843) = zin(845) + dzkl*zin(842)
                                      ! i4 = i4 + lang+1 =  846

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  841

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  847

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  855

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  854

                                      xin(855) = xin(855) + dxkl*xin(854)
                                      yin(855) = yin(855) + dykl*yin(854)
                                      zin(855) = zin(855) + dzkl*zin(854)

                                      ! i3 = i4 =  854
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  853

                                      xin(854) = xin(854) + dxkl*xin(853)
                                      yin(854) = yin(854) + dykl*yin(853)
                                      zin(854) = zin(854) + dzkl*zin(853)

                                      ! i3 = i4 =  853
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  855

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  854

                                      xin(855) = xin(855) + dxkl*xin(854)
                                      yin(855) = yin(855) + dykl*yin(854)
                                      zin(855) = zin(855) + dzkl*zin(854)

                                      ! i3 = i4 =  854
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  848

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  848

                                      ! do nk = 1,    2

                                      xin(848) = xin(850) + dxkl*xin(847)
                                      yin(848) = yin(850) + dykl*yin(847)
                                      zin(848) = zin(850) + dzkl*zin(847)
                                      ! i4 = i4 + lang+1 =  851

                                      ! nk =    2

                                      xin(851) = xin(853) + dxkl*xin(850)
                                      yin(851) = yin(853) + dykl*yin(850)
                                      zin(851) = zin(853) + dzkl*zin(850)
                                      ! i4 = i4 + lang+1 =  854

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  849

                                      ! nl =    2

                                      ! i4 = i3 =  849

                                      ! do nk = 1,    2

                                      xin(849) = xin(851) + dxkl*xin(848)
                                      yin(849) = yin(851) + dykl*yin(848)
                                      zin(849) = zin(851) + dzkl*zin(848)
                                      ! i4 = i4 + lang+1 =  852

                                      ! nk =    2

                                      xin(852) = xin(854) + dxkl*xin(851)
                                      yin(852) = yin(854) + dykl*yin(851)
                                      zin(852) = zin(854) + dzkl*zin(851)
                                      ! i4 = i4 + lang+1 =  855

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  850

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  856

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  864

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  863

                                      xin(864) = xin(864) + dxkl*xin(863)
                                      yin(864) = yin(864) + dykl*yin(863)
                                      zin(864) = zin(864) + dzkl*zin(863)

                                      ! i3 = i4 =  863
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  862

                                      xin(863) = xin(863) + dxkl*xin(862)
                                      yin(863) = yin(863) + dykl*yin(862)
                                      zin(863) = zin(863) + dzkl*zin(862)

                                      ! i3 = i4 =  862
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  864

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  863

                                      xin(864) = xin(864) + dxkl*xin(863)
                                      yin(864) = yin(864) + dykl*yin(863)
                                      zin(864) = zin(864) + dzkl*zin(863)

                                      ! i3 = i4 =  863
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  857

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  857

                                      ! do nk = 1,    2

                                      xin(857) = xin(859) + dxkl*xin(856)
                                      yin(857) = yin(859) + dykl*yin(856)
                                      zin(857) = zin(859) + dzkl*zin(856)
                                      ! i4 = i4 + lang+1 =  860

                                      ! nk =    2

                                      xin(860) = xin(862) + dxkl*xin(859)
                                      yin(860) = yin(862) + dykl*yin(859)
                                      zin(860) = zin(862) + dzkl*zin(859)
                                      ! i4 = i4 + lang+1 =  863

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  858

                                      ! nl =    2

                                      ! i4 = i3 =  858

                                      ! do nk = 1,    2

                                      xin(858) = xin(860) + dxkl*xin(857)
                                      yin(858) = yin(860) + dykl*yin(857)
                                      zin(858) = zin(860) + dzkl*zin(857)
                                      ! i4 = i4 + lang+1 =  861

                                      ! nk =    2

                                      xin(861) = xin(863) + dxkl*xin(860)
                                      yin(861) = yin(863) + dykl*yin(860)
                                      zin(861) = zin(863) + dzkl*zin(860)
                                      ! i4 = i4 + lang+1 =  864

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  859

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  865

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  865

                                      ! end do

                                      ! *** Now root =    7

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  864

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 3600! loop over all integrals

                                        l = n - 36*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d22ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  2
                                                        + xin(mx + 288)*yin(my + 288)*zin(mz + 288) & ! root  3
                                                        + xin(mx + 432)*yin(my + 432)*zin(mz + 432) & ! root  4
                                                        + xin(mx + 576)*yin(my + 576)*zin(mz + 576) & ! root  5
                                                        + xin(mx + 720)*yin(my + 720)*zin(mz + 720)) ! root  6

                                        j = int(n/36) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxj2 = 10
                                  iandj = ish .eq. jsh
                                  maxl = 6
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    if (iandj) maxj2 = i

                                    ii1 = i + loci
                                    ip = (i - 1)*360 ! Stride between functions in i

                                    do j = 1, maxj2

                                      maxl2 = maxl

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*36 + ip ! Add stride between functions in j

                                      do k = 1, 6 ! # of cartesians in k

                                        if (kandl) maxl2 = k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*6 + ijp ! Add stride between functions in k

                                        do l = 1, maxl2

                                          ijklp = ijkp + l ! No stride between functions in l

                                          buff(1) = eri_value(ijklp)

                                          if (abs(buff(1)) .lt. 5.0D-11) cycle ! Skip small integrals

                                          ll1 = l + locl
                                          k2 = kk1
                                          l2 = ll1

                                          if (k2 .lt. l2) then ! sort |kl>
                                            k2 = ll1
                                            l2 = kk1
                                          end if

                                          ii = i2
                                          jj = j2
                                          kk = k2
                                          ll = l2

                                          if (ii .lt. kk) then ! sort <ij|kl>
                                            ii = k2
                                            jj = l2
                                            kk = i2
                                            ll = j2

                                          else if (ii .eq. kk .and. jj .lt. ll) then ! sort <ij|il>
                                            jj = l2
                                            ll = j2
                                          end if

                                          ii2 = res%ia(ii)
                                          jj2 = res%ia(jj)
                                          kk2 = res%ia(kk)

                                          ij = ii2 + jj
                                          ik = ii2 + kk
                                          il = ii2 + ll
                                          jk = jj2 + kk
                                          jl = jj2 + ll
                                          kl = kk2 + ll

                                          if (jj .lt. kk) jk = kk2 + jj
                                          if (jj .lt. ll) jl = res%ia(ll) + jj

                                          ! Account for identical permutations

                                          if (ii .eq. jj) buff(1) = buff(1)*0.5_dp
                                          if (kk .eq. ll) buff(1) = buff(1)*0.5_dp
                                          if (ii .eq. kk .and. jj .eq. ll) buff(1) = buff(1)*0.5_dp

                                          buff(2) = buff(1)
                                          buff(3) = buff(1)*4.0_dp

                                          buff(4) = buff(3)*density(kl)
                                          buff(5) = buff(3)*density(ij)
                                          buff(6) = -buff(2)*density(jl)
                                          buff(7) = -buff(2)*density(ik)
                                          buff(8) = -buff(2)*density(jk)
                                          buff(9) = -buff(2)*density(il)

                                          !$omp atomic update
                                          fock(ij) = fock(ij) + buff(4)
                                          !$omp atomic update
                                          fock(kl) = fock(kl) + buff(5)
                                          !$omp atomic update
                                          fock(ik) = fock(ik) + buff(6)
                                          !$omp atomic update
                                          fock(jl) = fock(jl) + buff(7)
                                          !$omp atomic update
                                          fock(il) = fock(il) + buff(8)
                                          !$omp atomic update
                                          fock(jk) = fock(jk) + buff(9)

                                        end do ! l
                                      end do ! k
                                    end do ! j
                                  end do ! i

                                  !                  --- END DIRFCK_RHF ---

                                  end if ! Screening if

                                end do ! iquart

                                !$omp end target teams distribute parallel do


                              end do ! itile

                              deallocate (n33bra)
                              deallocate (xint33bra)
                              deallocate (n22ket)
                              deallocate (xint22ket)

                              end subroutine int3322
                              end submodule
