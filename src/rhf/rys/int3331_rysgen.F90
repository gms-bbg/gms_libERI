! The total angular momentum of this class is:          10
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3331_impl
contains
  module subroutine int3331(ff_pair, pf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, pf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n13ket(:)
    real(dp), allocatable :: xint33bra(:), xint13ket(:)
    integer(kind=int64) :: nffbra, npfket
    real(dp) :: scutffbra, scutpfket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxj2
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
    real(dp) :: xin(768), yin(768), zin(768)
    real(dp) :: eri_value(3000)
    real(dp) :: d33bra(100), d13ket(30)
    integer(kind=int64) :: ix(10), jx(10), kx(10), lx(3)
    integer(kind=int64) :: iy(10), jy(10), ky(10), ly(3)
    integer(kind=int64) :: iz(10), jz(10), kz(10), lz(3)
    integer(kind=int64) :: in(7), in1(7), kn(5)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(30), kly(30), klz(30)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 33
    in1(3) = 65
    in1(4) = 97
    in1(5) = 105
    in1(6) = 113
    in1(7) = 121

    kn(1) = 0
    kn(2) = 2
    kn(3) = 4
    kn(4) = 6
    kn(5) = 7

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 1
    lx(2) = 0
    lx(3) = 0

    kx(1) = 6
    kx(2) = 0
    kx(3) = 0
    kx(4) = 4
    kx(5) = 4
    kx(6) = 2
    kx(7) = 0
    kx(8) = 2
    kx(9) = 0
    kx(10) = 2

    jx(1) = 24
    jx(2) = 0
    jx(3) = 0
    jx(4) = 16
    jx(5) = 16
    jx(6) = 8
    jx(7) = 0
    jx(8) = 8
    jx(9) = 0
    jx(10) = 8

    ix(1) = 97
    ix(2) = 1
    ix(3) = 1
    ix(4) = 65
    ix(5) = 65
    ix(6) = 33
    ix(7) = 1
    ix(8) = 33
    ix(9) = 1
    ix(10) = 33

    ! y-arrays

    ly(1) = 0
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 6
    ky(3) = 0
    ky(4) = 2
    ky(5) = 0
    ky(6) = 4
    ky(7) = 4
    ky(8) = 0
    ky(9) = 2
    ky(10) = 2

    jy(1) = 0
    jy(2) = 24
    jy(3) = 0
    jy(4) = 8
    jy(5) = 0
    jy(6) = 16
    jy(7) = 16
    jy(8) = 0
    jy(9) = 8
    jy(10) = 8

    iy(1) = 1
    iy(2) = 97
    iy(3) = 1
    iy(4) = 33
    iy(5) = 1
    iy(6) = 65
    iy(7) = 65
    iy(8) = 1
    iy(9) = 33
    iy(10) = 33

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 6
    kz(4) = 0
    kz(5) = 2
    kz(6) = 0
    kz(7) = 2
    kz(8) = 4
    kz(9) = 4
    kz(10) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 24
    jz(4) = 0
    jz(5) = 8
    jz(6) = 0
    jz(7) = 8
    jz(8) = 16
    jz(9) = 16
    jz(10) = 8

    iz(1) = 1
    iz(2) = 1
    iz(3) = 97
    iz(4) = 1
    iz(5) = 33
    iz(6) = 1
    iz(7) = 33
    iz(8) = 65
    iz(9) = 65
    iz(10) = 33

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 121
    ijx(2) = 97
    ijx(3) = 97
    ijx(4) = 113
    ijx(5) = 113
    ijx(6) = 105
    ijx(7) = 97
    ijx(8) = 105
    ijx(9) = 97
    ijx(10) = 105
    ijx(11) = 25
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 17
    ijx(15) = 17
    ijx(16) = 9
    ijx(17) = 1
    ijx(18) = 9
    ijx(19) = 1
    ijx(20) = 9
    ijx(21) = 25
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 17
    ijx(25) = 17
    ijx(26) = 9
    ijx(27) = 1
    ijx(28) = 9
    ijx(29) = 1
    ijx(30) = 9
    ijx(31) = 89
    ijx(32) = 65
    ijx(33) = 65
    ijx(34) = 81
    ijx(35) = 81
    ijx(36) = 73
    ijx(37) = 65
    ijx(38) = 73
    ijx(39) = 65
    ijx(40) = 73
    ijx(41) = 89
    ijx(42) = 65
    ijx(43) = 65
    ijx(44) = 81
    ijx(45) = 81
    ijx(46) = 73
    ijx(47) = 65
    ijx(48) = 73
    ijx(49) = 65
    ijx(50) = 73
    ijx(51) = 57
    ijx(52) = 33
    ijx(53) = 33
    ijx(54) = 49
    ijx(55) = 49
    ijx(56) = 41
    ijx(57) = 33
    ijx(58) = 41
    ijx(59) = 33
    ijx(60) = 41
    ijx(61) = 25
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 17
    ijx(65) = 17
    ijx(66) = 9
    ijx(67) = 1
    ijx(68) = 9
    ijx(69) = 1
    ijx(70) = 9
    ijx(71) = 57
    ijx(72) = 33
    ijx(73) = 33
    ijx(74) = 49
    ijx(75) = 49
    ijx(76) = 41
    ijx(77) = 33
    ijx(78) = 41
    ijx(79) = 33
    ijx(80) = 41
    ijx(81) = 25
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 17
    ijx(85) = 17
    ijx(86) = 9
    ijx(87) = 1
    ijx(88) = 9
    ijx(89) = 1
    ijx(90) = 9
    ijx(91) = 57
    ijx(92) = 33
    ijx(93) = 33
    ijx(94) = 49
    ijx(95) = 49
    ijx(96) = 41
    ijx(97) = 33
    ijx(98) = 41
    ijx(99) = 33
    ijx(100) = 41

    ijy(1) = 1
    ijy(2) = 25
    ijy(3) = 1
    ijy(4) = 9
    ijy(5) = 1
    ijy(6) = 17
    ijy(7) = 17
    ijy(8) = 1
    ijy(9) = 9
    ijy(10) = 9
    ijy(11) = 97
    ijy(12) = 121
    ijy(13) = 97
    ijy(14) = 105
    ijy(15) = 97
    ijy(16) = 113
    ijy(17) = 113
    ijy(18) = 97
    ijy(19) = 105
    ijy(20) = 105
    ijy(21) = 1
    ijy(22) = 25
    ijy(23) = 1
    ijy(24) = 9
    ijy(25) = 1
    ijy(26) = 17
    ijy(27) = 17
    ijy(28) = 1
    ijy(29) = 9
    ijy(30) = 9
    ijy(31) = 33
    ijy(32) = 57
    ijy(33) = 33
    ijy(34) = 41
    ijy(35) = 33
    ijy(36) = 49
    ijy(37) = 49
    ijy(38) = 33
    ijy(39) = 41
    ijy(40) = 41
    ijy(41) = 1
    ijy(42) = 25
    ijy(43) = 1
    ijy(44) = 9
    ijy(45) = 1
    ijy(46) = 17
    ijy(47) = 17
    ijy(48) = 1
    ijy(49) = 9
    ijy(50) = 9
    ijy(51) = 65
    ijy(52) = 89
    ijy(53) = 65
    ijy(54) = 73
    ijy(55) = 65
    ijy(56) = 81
    ijy(57) = 81
    ijy(58) = 65
    ijy(59) = 73
    ijy(60) = 73
    ijy(61) = 65
    ijy(62) = 89
    ijy(63) = 65
    ijy(64) = 73
    ijy(65) = 65
    ijy(66) = 81
    ijy(67) = 81
    ijy(68) = 65
    ijy(69) = 73
    ijy(70) = 73
    ijy(71) = 1
    ijy(72) = 25
    ijy(73) = 1
    ijy(74) = 9
    ijy(75) = 1
    ijy(76) = 17
    ijy(77) = 17
    ijy(78) = 1
    ijy(79) = 9
    ijy(80) = 9
    ijy(81) = 33
    ijy(82) = 57
    ijy(83) = 33
    ijy(84) = 41
    ijy(85) = 33
    ijy(86) = 49
    ijy(87) = 49
    ijy(88) = 33
    ijy(89) = 41
    ijy(90) = 41
    ijy(91) = 33
    ijy(92) = 57
    ijy(93) = 33
    ijy(94) = 41
    ijy(95) = 33
    ijy(96) = 49
    ijy(97) = 49
    ijy(98) = 33
    ijy(99) = 41
    ijy(100) = 41

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 25
    ijz(4) = 1
    ijz(5) = 9
    ijz(6) = 1
    ijz(7) = 9
    ijz(8) = 17
    ijz(9) = 17
    ijz(10) = 9
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 25
    ijz(14) = 1
    ijz(15) = 9
    ijz(16) = 1
    ijz(17) = 9
    ijz(18) = 17
    ijz(19) = 17
    ijz(20) = 9
    ijz(21) = 97
    ijz(22) = 97
    ijz(23) = 121
    ijz(24) = 97
    ijz(25) = 105
    ijz(26) = 97
    ijz(27) = 105
    ijz(28) = 113
    ijz(29) = 113
    ijz(30) = 105
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 25
    ijz(34) = 1
    ijz(35) = 9
    ijz(36) = 1
    ijz(37) = 9
    ijz(38) = 17
    ijz(39) = 17
    ijz(40) = 9
    ijz(41) = 33
    ijz(42) = 33
    ijz(43) = 57
    ijz(44) = 33
    ijz(45) = 41
    ijz(46) = 33
    ijz(47) = 41
    ijz(48) = 49
    ijz(49) = 49
    ijz(50) = 41
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 25
    ijz(54) = 1
    ijz(55) = 9
    ijz(56) = 1
    ijz(57) = 9
    ijz(58) = 17
    ijz(59) = 17
    ijz(60) = 9
    ijz(61) = 33
    ijz(62) = 33
    ijz(63) = 57
    ijz(64) = 33
    ijz(65) = 41
    ijz(66) = 33
    ijz(67) = 41
    ijz(68) = 49
    ijz(69) = 49
    ijz(70) = 41
    ijz(71) = 65
    ijz(72) = 65
    ijz(73) = 89
    ijz(74) = 65
    ijz(75) = 73
    ijz(76) = 65
    ijz(77) = 73
    ijz(78) = 81
    ijz(79) = 81
    ijz(80) = 73
    ijz(81) = 65
    ijz(82) = 65
    ijz(83) = 89
    ijz(84) = 65
    ijz(85) = 73
    ijz(86) = 65
    ijz(87) = 73
    ijz(88) = 81
    ijz(89) = 81
    ijz(90) = 73
    ijz(91) = 33
    ijz(92) = 33
    ijz(93) = 57
    ijz(94) = 33
    ijz(95) = 41
    ijz(96) = 33
    ijz(97) = 41
    ijz(98) = 49
    ijz(99) = 49
    ijz(100) = 41

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 7
    klx(2) = 6
    klx(3) = 6
    klx(4) = 1
    klx(5) = 0
    klx(6) = 0
    klx(7) = 1
    klx(8) = 0
    klx(9) = 0
    klx(10) = 5
    klx(11) = 4
    klx(12) = 4
    klx(13) = 5
    klx(14) = 4
    klx(15) = 4
    klx(16) = 3
    klx(17) = 2
    klx(18) = 2
    klx(19) = 1
    klx(20) = 0
    klx(21) = 0
    klx(22) = 3
    klx(23) = 2
    klx(24) = 2
    klx(25) = 1
    klx(26) = 0
    klx(27) = 0
    klx(28) = 3
    klx(29) = 2
    klx(30) = 2

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0
    kly(4) = 6
    kly(5) = 7
    kly(6) = 6
    kly(7) = 0
    kly(8) = 1
    kly(9) = 0
    kly(10) = 2
    kly(11) = 3
    kly(12) = 2
    kly(13) = 0
    kly(14) = 1
    kly(15) = 0
    kly(16) = 4
    kly(17) = 5
    kly(18) = 4
    kly(19) = 4
    kly(20) = 5
    kly(21) = 4
    kly(22) = 0
    kly(23) = 1
    kly(24) = 0
    kly(25) = 2
    kly(26) = 3
    kly(27) = 2
    kly(28) = 2
    kly(29) = 3
    kly(30) = 2

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1
    klz(4) = 0
    klz(5) = 0
    klz(6) = 1
    klz(7) = 6
    klz(8) = 6
    klz(9) = 7
    klz(10) = 0
    klz(11) = 0
    klz(12) = 1
    klz(13) = 2
    klz(14) = 2
    klz(15) = 3
    klz(16) = 0
    klz(17) = 0
    klz(18) = 1
    klz(19) = 2
    klz(20) = 2
    klz(21) = 3
    klz(22) = 4
    klz(23) = 4
    klz(24) = 5
    klz(25) = 4
    klz(26) = 4
    klz(27) = 5
    klz(28) = 2
    klz(29) = 2
    klz(30) = 3

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n13ket(res%n_p_shl*res%n_f_shl))
    allocate (xint13ket(res%n_p_shl*res%n_f_shl))

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

    scutpfket = cutoff_schwarz/maxval(pf_pair%xints)
    npfket = 0
    do ij = 1, res%n_p_shl*res%n_f_shl
      if (pf_pair%xints(ij) .ge. scutpfket) then
        npfket = npfket + 1
        xint13ket(npfket) = pf_pair%xints(ij)
        n13ket(npfket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*npfket) .le. nchunksize_int64) nchunksize_int64 = nffbra*npfket
    ntile = int(nffbra*npfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*npfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint13ket, n13ket, ff_pair, pf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d13ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d33bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,iandj)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nffbra) + 1
              kl_tmp = (iquart - 1)/nffbra + 1

              test = xint33bra(ij_tmp)*xint13ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n13ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_f_shl(ksh_tmp)
                lsh = res%i_p_shl(lsh_tmp)

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

                  t_expon_cd = pf_pair%t_expon_ab(pf_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = pf_pair%expon_b(pf_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = pf_pair%expon_a(pf_pair%pair_loc(kl) + ket_loop)
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

                  d13ket(1) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(2) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(3) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(4) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(5) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(6) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(7) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(8) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(9) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(10) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(11) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(12) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(13) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(14) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(15) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(16) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(17) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(18) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(19) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(20) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(21) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(22) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(23) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(24) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(25) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(26) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(27) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(28) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                  d13ket(29) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                  d13ket(30) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3

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

                                      ! i2 = in(2) =   33
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(33) = xc00
                                      yin(33) = yc00
                                      zin(33) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   35
                                      ! i2 =   33

                                      xin(35) = xcp00*xin(33) + cp10
                                      yin(35) = ycp00*yin(33) + cp10
                                      zin(35) = zcp00*zin(33) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   33

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   65
                                      ! i3 =    1
                                      ! i4 =   33

                                      xin(65) = c10*xin(1) + xc00*xin(33)
                                      yin(65) = c10*yin(1) + yc00*yin(33)
                                      zin(65) = c10*zin(1) + zc00*zin(33)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   67
                                      ! i5 =   65
                                      ! i4 =   33

                                      xin(67) = xcp00*xin(65) + cp10*xin(33)
                                      yin(67) = ycp00*yin(65) + cp10*yin(33)
                                      zin(67) = zcp00*zin(65) + cp10*zin(33)

                                      ! ------------------

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   65

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   97
                                      ! i3 =   33
                                      ! i4 =   65

                                      xin(97) = c10*xin(33) + xc00*xin(65)
                                      yin(97) = c10*yin(33) + yc00*yin(65)
                                      zin(97) = c10*zin(33) + zc00*zin(65)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   99
                                      ! i5 =   97
                                      ! i4 =   65

                                      xin(99) = xcp00*xin(97) + cp10*xin(65)
                                      yin(99) = ycp00*yin(97) + cp10*yin(65)
                                      zin(99) = zcp00*zin(97) + cp10*zin(65)

                                      ! ------------------

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   97

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  105
                                      ! i3 =   65
                                      ! i4 =   97

                                      xin(105) = c10*xin(65) + xc00*xin(97)
                                      yin(105) = c10*yin(65) + yc00*yin(97)
                                      zin(105) = c10*zin(65) + zc00*zin(97)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  107
                                      ! i5 =  105
                                      ! i4 =   97

                                      xin(107) = xcp00*xin(105) + cp10*xin(97)
                                      yin(107) = ycp00*yin(105) + cp10*yin(97)
                                      zin(107) = zcp00*zin(105) + cp10*zin(97)

                                      ! ------------------

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  105

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  113
                                      ! i3 =   97
                                      ! i4 =  105

                                      xin(113) = c10*xin(97) + xc00*xin(105)
                                      yin(113) = c10*yin(97) + yc00*yin(105)
                                      zin(113) = c10*zin(97) + zc00*zin(105)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  115
                                      ! i5 =  113
                                      ! i4 =  105

                                      xin(115) = xcp00*xin(113) + cp10*xin(105)
                                      yin(115) = ycp00*yin(113) + cp10*yin(105)
                                      zin(115) = zcp00*zin(113) + cp10*zin(105)

                                      ! ------------------

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =  105
                                      ! i4 =  113

                                      xin(121) = c10*xin(105) + xc00*xin(113)
                                      yin(121) = c10*yin(105) + yc00*yin(113)
                                      zin(121) = c10*zin(105) + zc00*zin(113)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  123
                                      ! i5 =  121
                                      ! i4 =  113

                                      xin(123) = xcp00*xin(121) + cp10*xin(113)
                                      yin(123) = ycp00*yin(121) + cp10*yin(113)
                                      zin(123) = zcp00*zin(121) + cp10*zin(113)

                                      ! ------------------

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    3

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    5
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(5) = cp01*xin(1) + xcp00*xin(3)
                                      yin(5) = cp01*yin(1) + ycp00*yin(3)
                                      zin(5) = cp01*zin(1) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   37

                                      xin(37) = xc00*xin(5) + c01*xin(3)
                                      yin(37) = yc00*yin(5) + c01*yin(3)
                                      zin(37) = zc00*zin(5) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    5

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    7
                                      ! i3 =    3
                                      ! i4 =    5

                                      xin(7) = cp01*xin(3) + xcp00*xin(5)
                                      yin(7) = cp01*yin(3) + ycp00*yin(5)
                                      zin(7) = cp01*zin(3) + zcp00*zin(5)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   39

                                      xin(39) = xc00*xin(7) + c01*xin(5)
                                      yin(39) = yc00*yin(7) + c01*yin(5)
                                      zin(39) = zc00*zin(7) + c01*zin(5)

                                      ! ------------------

                                      ! i3 = i4 =    5
                                      ! i4 = i5 =    7

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    8
                                      ! i3 =    5
                                      ! i4 =    7

                                      xin(8) = cp01*xin(5) + xcp00*xin(7)
                                      yin(8) = cp01*yin(5) + ycp00*yin(7)
                                      zin(8) = cp01*zin(5) + zcp00*zin(7)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   40

                                      xin(40) = xc00*xin(8) + c01*xin(7)
                                      yin(40) = yc00*yin(8) + c01*yin(7)
                                      zin(40) = zc00*zin(8) + c01*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   33

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   65

                                      xin(69) = c10*xin(5) + xc00*xin(37) + c01*xin(35)
                                      yin(69) = c10*yin(5) + yc00*yin(37) + c01*yin(35)
                                      zin(69) = c10*zin(5) + zc00*zin(37) + c01*zin(35)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   65

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   97

                                      xin(101) = c10*xin(37) + xc00*xin(69) + c01*xin(67)
                                      yin(101) = c10*yin(37) + yc00*yin(69) + c01*yin(67)
                                      zin(101) = c10*zin(37) + zc00*zin(69) + c01*zin(67)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   97

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  105

                                      xin(109) = c10*xin(69) + xc00*xin(101) + c01*xin(99)
                                      yin(109) = c10*yin(69) + yc00*yin(101) + c01*yin(99)
                                      zin(109) = c10*zin(69) + zc00*zin(101) + c01*zin(99)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  105

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  113

                                      xin(117) = c10*xin(101) + xc00*xin(109) + c01*xin(107)
                                      yin(117) = c10*yin(101) + yc00*yin(109) + c01*yin(107)
                                      zin(117) = c10*zin(101) + zc00*zin(109) + c01*zin(107)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  121

                                      xin(125) = c10*xin(109) + xc00*xin(117) + c01*xin(115)
                                      yin(125) = c10*yin(109) + yc00*yin(117) + c01*yin(115)
                                      zin(125) = c10*zin(109) + zc00*zin(117) + c01*zin(115)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   33

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   65

                                      xin(71) = c10*xin(7) + xc00*xin(39) + c01*xin(37)
                                      yin(71) = c10*yin(7) + yc00*yin(39) + c01*yin(37)
                                      zin(71) = c10*zin(7) + zc00*zin(39) + c01*zin(37)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   65

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   97

                                      xin(103) = c10*xin(39) + xc00*xin(71) + c01*xin(69)
                                      yin(103) = c10*yin(39) + yc00*yin(71) + c01*yin(69)
                                      zin(103) = c10*zin(39) + zc00*zin(71) + c01*zin(69)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   97

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  105

                                      xin(111) = c10*xin(71) + xc00*xin(103) + c01*xin(101)
                                      yin(111) = c10*yin(71) + yc00*yin(103) + c01*yin(101)
                                      zin(111) = c10*zin(71) + zc00*zin(103) + c01*zin(101)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  105

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  113

                                      xin(119) = c10*xin(103) + xc00*xin(111) + c01*xin(109)
                                      yin(119) = c10*yin(103) + yc00*yin(111) + c01*yin(109)
                                      zin(119) = c10*zin(103) + zc00*zin(111) + c01*zin(109)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  121

                                      xin(127) = c10*xin(111) + xc00*xin(119) + c01*xin(117)
                                      yin(127) = c10*yin(111) + yc00*yin(119) + c01*yin(117)
                                      zin(127) = c10*zin(111) + zc00*zin(119) + c01*zin(117)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   33

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   65

                                      xin(72) = c10*xin(8) + xc00*xin(40) + c01*xin(39)
                                      yin(72) = c10*yin(8) + yc00*yin(40) + c01*yin(39)
                                      zin(72) = c10*zin(8) + zc00*zin(40) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   65

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   97

                                      xin(104) = c10*xin(40) + xc00*xin(72) + c01*xin(71)
                                      yin(104) = c10*yin(40) + yc00*yin(72) + c01*yin(71)
                                      zin(104) = c10*zin(40) + zc00*zin(72) + c01*zin(71)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   97

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  105

                                      xin(112) = c10*xin(72) + xc00*xin(104) + c01*xin(103)
                                      yin(112) = c10*yin(72) + yc00*yin(104) + c01*yin(103)
                                      zin(112) = c10*zin(72) + zc00*zin(104) + c01*zin(103)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  105

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  113

                                      xin(120) = c10*xin(104) + xc00*xin(112) + c01*xin(111)
                                      yin(120) = c10*yin(104) + yc00*yin(112) + c01*yin(111)
                                      zin(120) = c10*zin(104) + zc00*zin(112) + c01*zin(111)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  121

                                      xin(128) = c10*xin(112) + xc00*xin(120) + c01*xin(119)
                                      yin(128) = c10*yin(112) + yc00*yin(120) + c01*yin(119)
                                      zin(128) = c10*zin(112) + zc00*zin(120) + c01*zin(119)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  121

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  121

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  113

                                      xin(121) = xin(121) + dxij*xin(113)
                                      yin(121) = yin(121) + dyij*yin(113)
                                      zin(121) = zin(121) + dzij*zin(113)

                                      ! i3 = i4 =  113
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  105

                                      xin(113) = xin(113) + dxij*xin(105)
                                      yin(113) = yin(113) + dyij*yin(105)
                                      zin(113) = zin(113) + dzij*zin(105)

                                      ! i3 = i4 =  105
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   97

                                      xin(105) = xin(105) + dxij*xin(97)
                                      yin(105) = yin(105) + dyij*yin(97)
                                      zin(105) = zin(105) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  121

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  113

                                      xin(121) = xin(121) + dxij*xin(113)
                                      yin(121) = yin(121) + dyij*yin(113)
                                      zin(121) = zin(121) + dzij*zin(113)

                                      ! i3 = i4 =  113
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  105

                                      xin(113) = xin(113) + dxij*xin(105)
                                      yin(113) = yin(113) + dyij*yin(105)
                                      zin(113) = zin(113) + dzij*zin(105)

                                      ! i3 = i4 =  105
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  121

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  113

                                      xin(121) = xin(121) + dxij*xin(113)
                                      yin(121) = yin(121) + dyij*yin(113)
                                      zin(121) = zin(121) + dzij*zin(113)

                                      ! i3 = i4 =  113
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(33) + dxij*xin(1)
                                      yin(9) = yin(33) + dyij*yin(1)
                                      zin(9) = zin(33) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    2

                                      xin(41) = xin(65) + dxij*xin(33)
                                      yin(41) = yin(65) + dyij*yin(33)
                                      zin(41) = zin(65) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! ni =    3

                                      xin(73) = xin(97) + dxij*xin(65)
                                      yin(73) = yin(97) + dyij*yin(65)
                                      zin(73) = zin(97) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    3

                                      xin(17) = xin(41) + dxij*xin(9)
                                      yin(17) = yin(41) + dyij*yin(9)
                                      zin(17) = zin(41) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! ni =    2

                                      xin(49) = xin(73) + dxij*xin(41)
                                      yin(49) = yin(73) + dyij*yin(41)
                                      zin(49) = zin(73) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(105) + dxij*xin(73)
                                      yin(81) = yin(105) + dyij*yin(73)
                                      zin(81) = zin(105) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    3

                                      ! i4 = i3 =   25

                                      ! do ni = 1,    3

                                      xin(25) = xin(49) + dxij*xin(17)
                                      yin(25) = yin(49) + dyij*yin(17)
                                      zin(25) = zin(49) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    2

                                      xin(57) = xin(81) + dxij*xin(49)
                                      yin(57) = yin(81) + dyij*yin(49)
                                      zin(57) = zin(81) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    3

                                      xin(89) = xin(113) + dxij*xin(81)
                                      yin(89) = yin(113) + dyij*yin(81)
                                      zin(89) = zin(113) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   33

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  123

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  115

                                      xin(123) = xin(123) + dxij*xin(115)
                                      yin(123) = yin(123) + dyij*yin(115)
                                      zin(123) = zin(123) + dzij*zin(115)

                                      ! i3 = i4 =  115
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  107

                                      xin(115) = xin(115) + dxij*xin(107)
                                      yin(115) = yin(115) + dyij*yin(107)
                                      zin(115) = zin(115) + dzij*zin(107)

                                      ! i3 = i4 =  107
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   99

                                      xin(107) = xin(107) + dxij*xin(99)
                                      yin(107) = yin(107) + dyij*yin(99)
                                      zin(107) = zin(107) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  123

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  115

                                      xin(123) = xin(123) + dxij*xin(115)
                                      yin(123) = yin(123) + dyij*yin(115)
                                      zin(123) = zin(123) + dzij*zin(115)

                                      ! i3 = i4 =  115
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  107

                                      xin(115) = xin(115) + dxij*xin(107)
                                      yin(115) = yin(115) + dyij*yin(107)
                                      zin(115) = zin(115) + dzij*zin(107)

                                      ! i3 = i4 =  107
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  123

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  115

                                      xin(123) = xin(123) + dxij*xin(115)
                                      yin(123) = yin(123) + dyij*yin(115)
                                      zin(123) = zin(123) + dzij*zin(115)

                                      ! i3 = i4 =  115
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(35) + dxij*xin(3)
                                      yin(11) = yin(35) + dyij*yin(3)
                                      zin(11) = zin(35) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    2

                                      xin(43) = xin(67) + dxij*xin(35)
                                      yin(43) = yin(67) + dyij*yin(35)
                                      zin(43) = zin(67) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   75

                                      ! ni =    3

                                      xin(75) = xin(99) + dxij*xin(67)
                                      yin(75) = yin(99) + dyij*yin(67)
                                      zin(75) = zin(99) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    3

                                      xin(19) = xin(43) + dxij*xin(11)
                                      yin(19) = yin(43) + dyij*yin(11)
                                      zin(19) = zin(43) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   51

                                      ! ni =    2

                                      xin(51) = xin(75) + dxij*xin(43)
                                      yin(51) = yin(75) + dyij*yin(43)
                                      zin(51) = zin(75) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    3

                                      xin(83) = xin(107) + dxij*xin(75)
                                      yin(83) = yin(107) + dyij*yin(75)
                                      zin(83) = zin(107) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    3

                                      ! i4 = i3 =   27

                                      ! do ni = 1,    3

                                      xin(27) = xin(51) + dxij*xin(19)
                                      yin(27) = yin(51) + dyij*yin(19)
                                      zin(27) = zin(51) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    2

                                      xin(59) = xin(83) + dxij*xin(51)
                                      yin(59) = yin(83) + dyij*yin(51)
                                      zin(59) = zin(83) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    3

                                      xin(91) = xin(115) + dxij*xin(83)
                                      yin(91) = yin(115) + dyij*yin(83)
                                      zin(91) = zin(115) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  123

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   35

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  117

                                      xin(125) = xin(125) + dxij*xin(117)
                                      yin(125) = yin(125) + dyij*yin(117)
                                      zin(125) = zin(125) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  109

                                      xin(117) = xin(117) + dxij*xin(109)
                                      yin(117) = yin(117) + dyij*yin(109)
                                      zin(117) = zin(117) + dzij*zin(109)

                                      ! i3 = i4 =  109
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  101

                                      xin(109) = xin(109) + dxij*xin(101)
                                      yin(109) = yin(109) + dyij*yin(101)
                                      zin(109) = zin(109) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  117

                                      xin(125) = xin(125) + dxij*xin(117)
                                      yin(125) = yin(125) + dyij*yin(117)
                                      zin(125) = zin(125) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  109

                                      xin(117) = xin(117) + dxij*xin(109)
                                      yin(117) = yin(117) + dyij*yin(109)
                                      zin(117) = zin(117) + dzij*zin(109)

                                      ! i3 = i4 =  109
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  117

                                      xin(125) = xin(125) + dxij*xin(117)
                                      yin(125) = yin(125) + dyij*yin(117)
                                      zin(125) = zin(125) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   13

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(37) + dxij*xin(5)
                                      yin(13) = yin(37) + dyij*yin(5)
                                      zin(13) = zin(37) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    2

                                      xin(45) = xin(69) + dxij*xin(37)
                                      yin(45) = yin(69) + dyij*yin(37)
                                      zin(45) = zin(69) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   77

                                      ! ni =    3

                                      xin(77) = xin(101) + dxij*xin(69)
                                      yin(77) = yin(101) + dyij*yin(69)
                                      zin(77) = zin(101) + dzij*zin(69)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    2

                                      ! i4 = i3 =   21

                                      ! do ni = 1,    3

                                      xin(21) = xin(45) + dxij*xin(13)
                                      yin(21) = yin(45) + dyij*yin(13)
                                      zin(21) = zin(45) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    2

                                      xin(53) = xin(77) + dxij*xin(45)
                                      yin(53) = yin(77) + dyij*yin(45)
                                      zin(53) = zin(77) + dzij*zin(45)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    3

                                      xin(85) = xin(109) + dxij*xin(77)
                                      yin(85) = yin(109) + dyij*yin(77)
                                      zin(85) = zin(109) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   29

                                      ! nj =    3

                                      ! i4 = i3 =   29

                                      ! do ni = 1,    3

                                      xin(29) = xin(53) + dxij*xin(21)
                                      yin(29) = yin(53) + dyij*yin(21)
                                      zin(29) = zin(53) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    2

                                      xin(61) = xin(85) + dxij*xin(53)
                                      yin(61) = yin(85) + dyij*yin(53)
                                      zin(61) = zin(85) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    3

                                      xin(93) = xin(117) + dxij*xin(85)
                                      yin(93) = yin(117) + dyij*yin(85)
                                      zin(93) = zin(117) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   37

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  119

                                      xin(127) = xin(127) + dxij*xin(119)
                                      yin(127) = yin(127) + dyij*yin(119)
                                      zin(127) = zin(127) + dzij*zin(119)

                                      ! i3 = i4 =  119
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  111

                                      xin(119) = xin(119) + dxij*xin(111)
                                      yin(119) = yin(119) + dyij*yin(111)
                                      zin(119) = zin(119) + dzij*zin(111)

                                      ! i3 = i4 =  111
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  103

                                      xin(111) = xin(111) + dxij*xin(103)
                                      yin(111) = yin(111) + dyij*yin(103)
                                      zin(111) = zin(111) + dzij*zin(103)

                                      ! i3 = i4 =  103
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  119

                                      xin(127) = xin(127) + dxij*xin(119)
                                      yin(127) = yin(127) + dyij*yin(119)
                                      zin(127) = zin(127) + dzij*zin(119)

                                      ! i3 = i4 =  119
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  111

                                      xin(119) = xin(119) + dxij*xin(111)
                                      yin(119) = yin(119) + dyij*yin(111)
                                      zin(119) = zin(119) + dzij*zin(111)

                                      ! i3 = i4 =  111
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  119

                                      xin(127) = xin(127) + dxij*xin(119)
                                      yin(127) = yin(127) + dyij*yin(119)
                                      zin(127) = zin(127) + dzij*zin(119)

                                      ! i3 = i4 =  119
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   15

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    3

                                      xin(15) = xin(39) + dxij*xin(7)
                                      yin(15) = yin(39) + dyij*yin(7)
                                      zin(15) = zin(39) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    2

                                      xin(47) = xin(71) + dxij*xin(39)
                                      yin(47) = yin(71) + dyij*yin(39)
                                      zin(47) = zin(71) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(103) + dxij*xin(71)
                                      yin(79) = yin(103) + dyij*yin(71)
                                      zin(79) = zin(103) + dzij*zin(71)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  111

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   23

                                      ! nj =    2

                                      ! i4 = i3 =   23

                                      ! do ni = 1,    3

                                      xin(23) = xin(47) + dxij*xin(15)
                                      yin(23) = yin(47) + dyij*yin(15)
                                      zin(23) = zin(47) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    2

                                      xin(55) = xin(79) + dxij*xin(47)
                                      yin(55) = yin(79) + dyij*yin(47)
                                      zin(55) = zin(79) + dzij*zin(47)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    3

                                      xin(87) = xin(111) + dxij*xin(79)
                                      yin(87) = yin(111) + dyij*yin(79)
                                      zin(87) = zin(111) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   31

                                      ! nj =    3

                                      ! i4 = i3 =   31

                                      ! do ni = 1,    3

                                      xin(31) = xin(55) + dxij*xin(23)
                                      yin(31) = yin(55) + dyij*yin(23)
                                      zin(31) = zin(55) + dzij*zin(23)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    2

                                      xin(63) = xin(87) + dxij*xin(55)
                                      yin(63) = yin(87) + dyij*yin(55)
                                      zin(63) = zin(87) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    3

                                      xin(95) = xin(119) + dxij*xin(87)
                                      yin(95) = yin(119) + dyij*yin(87)
                                      zin(95) = zin(119) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   39

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  120

                                      xin(128) = xin(128) + dxij*xin(120)
                                      yin(128) = yin(128) + dyij*yin(120)
                                      zin(128) = zin(128) + dzij*zin(120)

                                      ! i3 = i4 =  120
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  112

                                      xin(120) = xin(120) + dxij*xin(112)
                                      yin(120) = yin(120) + dyij*yin(112)
                                      zin(120) = zin(120) + dzij*zin(112)

                                      ! i3 = i4 =  112
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  104

                                      xin(112) = xin(112) + dxij*xin(104)
                                      yin(112) = yin(112) + dyij*yin(104)
                                      zin(112) = zin(112) + dzij*zin(104)

                                      ! i3 = i4 =  104
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  120

                                      xin(128) = xin(128) + dxij*xin(120)
                                      yin(128) = yin(128) + dyij*yin(120)
                                      zin(128) = zin(128) + dzij*zin(120)

                                      ! i3 = i4 =  120
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  112

                                      xin(120) = xin(120) + dxij*xin(112)
                                      yin(120) = yin(120) + dyij*yin(112)
                                      zin(120) = zin(120) + dzij*zin(112)

                                      ! i3 = i4 =  112
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  120

                                      xin(128) = xin(128) + dxij*xin(120)
                                      yin(128) = yin(128) + dyij*yin(120)
                                      zin(128) = zin(128) + dzij*zin(120)

                                      ! i3 = i4 =  120
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   16

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   16

                                      ! do ni = 1,    3

                                      xin(16) = xin(40) + dxij*xin(8)
                                      yin(16) = yin(40) + dyij*yin(8)
                                      zin(16) = zin(40) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    2

                                      xin(48) = xin(72) + dxij*xin(40)
                                      yin(48) = yin(72) + dyij*yin(40)
                                      zin(48) = zin(72) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    3

                                      xin(80) = xin(104) + dxij*xin(72)
                                      yin(80) = yin(104) + dyij*yin(72)
                                      zin(80) = zin(104) + dzij*zin(72)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  112

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   24

                                      ! nj =    2

                                      ! i4 = i3 =   24

                                      ! do ni = 1,    3

                                      xin(24) = xin(48) + dxij*xin(16)
                                      yin(24) = yin(48) + dyij*yin(16)
                                      zin(24) = zin(48) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   56

                                      ! ni =    2

                                      xin(56) = xin(80) + dxij*xin(48)
                                      yin(56) = yin(80) + dyij*yin(48)
                                      zin(56) = zin(80) + dzij*zin(48)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    3

                                      xin(88) = xin(112) + dxij*xin(80)
                                      yin(88) = yin(112) + dyij*yin(80)
                                      zin(88) = zin(112) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   32

                                      ! nj =    3

                                      ! i4 = i3 =   32

                                      ! do ni = 1,    3

                                      xin(32) = xin(56) + dxij*xin(24)
                                      yin(32) = yin(56) + dyij*yin(24)
                                      zin(32) = zin(56) + dzij*zin(24)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    2

                                      xin(64) = xin(88) + dxij*xin(56)
                                      yin(64) = yin(88) + dyij*yin(56)
                                      zin(64) = zin(88) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    3

                                      xin(96) = xin(120) + dxij*xin(88)
                                      yin(96) = yin(120) + dyij*yin(88)
                                      zin(96) = zin(120) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   40

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    8

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    3

                                      xin(2) = xin(3) + dxkl*xin(1)
                                      yin(2) = yin(3) + dykl*yin(1)
                                      zin(2) = zin(3) + dzkl*zin(1)
                                      ! i4 = i4 + lang+1 =    4

                                      ! nk =    2

                                      xin(4) = xin(5) + dxkl*xin(3)
                                      yin(4) = yin(5) + dykl*yin(3)
                                      zin(4) = zin(5) + dzkl*zin(3)
                                      ! i4 = i4 + lang+1 =    6

                                      ! nk =    3

                                      xin(6) = xin(7) + dxkl*xin(5)
                                      yin(6) = yin(7) + dykl*yin(5)
                                      zin(6) = zin(7) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    9

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   16

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   15

                                      xin(16) = xin(16) + dxkl*xin(15)
                                      yin(16) = yin(16) + dykl*yin(15)
                                      zin(16) = zin(16) + dzkl*zin(15)

                                      ! i3 = i4 =   15
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   10

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   10

                                      ! do nk = 1,    3

                                      xin(10) = xin(11) + dxkl*xin(9)
                                      yin(10) = yin(11) + dykl*yin(9)
                                      zin(10) = zin(11) + dzkl*zin(9)
                                      ! i4 = i4 + lang+1 =   12

                                      ! nk =    2

                                      xin(12) = xin(13) + dxkl*xin(11)
                                      yin(12) = yin(13) + dykl*yin(11)
                                      zin(12) = zin(13) + dzkl*zin(11)
                                      ! i4 = i4 + lang+1 =   14

                                      ! nk =    3

                                      xin(14) = xin(15) + dxkl*xin(13)
                                      yin(14) = yin(15) + dykl*yin(13)
                                      zin(14) = zin(15) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   16

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   11

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   17

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   24

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   23

                                      xin(24) = xin(24) + dxkl*xin(23)
                                      yin(24) = yin(24) + dykl*yin(23)
                                      zin(24) = zin(24) + dzkl*zin(23)

                                      ! i3 = i4 =   23
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   18

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   18

                                      ! do nk = 1,    3

                                      xin(18) = xin(19) + dxkl*xin(17)
                                      yin(18) = yin(19) + dykl*yin(17)
                                      zin(18) = zin(19) + dzkl*zin(17)
                                      ! i4 = i4 + lang+1 =   20

                                      ! nk =    2

                                      xin(20) = xin(21) + dxkl*xin(19)
                                      yin(20) = yin(21) + dykl*yin(19)
                                      zin(20) = zin(21) + dzkl*zin(19)
                                      ! i4 = i4 + lang+1 =   22

                                      ! nk =    3

                                      xin(22) = xin(23) + dxkl*xin(21)
                                      yin(22) = yin(23) + dykl*yin(21)
                                      zin(22) = zin(23) + dzkl*zin(21)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   19

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   25

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   32

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   31

                                      xin(32) = xin(32) + dxkl*xin(31)
                                      yin(32) = yin(32) + dykl*yin(31)
                                      zin(32) = zin(32) + dzkl*zin(31)

                                      ! i3 = i4 =   31
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   26

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   26

                                      ! do nk = 1,    3

                                      xin(26) = xin(27) + dxkl*xin(25)
                                      yin(26) = yin(27) + dykl*yin(25)
                                      zin(26) = zin(27) + dzkl*zin(25)
                                      ! i4 = i4 + lang+1 =   28

                                      ! nk =    2

                                      xin(28) = xin(29) + dxkl*xin(27)
                                      yin(28) = yin(29) + dykl*yin(27)
                                      zin(28) = zin(29) + dzkl*zin(27)
                                      ! i4 = i4 + lang+1 =   30

                                      ! nk =    3

                                      xin(30) = xin(31) + dxkl*xin(29)
                                      yin(30) = yin(31) + dykl*yin(29)
                                      zin(30) = zin(31) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   27

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   33

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   40

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   39

                                      xin(40) = xin(40) + dxkl*xin(39)
                                      yin(40) = yin(40) + dykl*yin(39)
                                      zin(40) = zin(40) + dzkl*zin(39)

                                      ! i3 = i4 =   39
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   34

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   34

                                      ! do nk = 1,    3

                                      xin(34) = xin(35) + dxkl*xin(33)
                                      yin(34) = yin(35) + dykl*yin(33)
                                      zin(34) = zin(35) + dzkl*zin(33)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    2

                                      xin(36) = xin(37) + dxkl*xin(35)
                                      yin(36) = yin(37) + dykl*yin(35)
                                      zin(36) = zin(37) + dzkl*zin(35)
                                      ! i4 = i4 + lang+1 =   38

                                      ! nk =    3

                                      xin(38) = xin(39) + dxkl*xin(37)
                                      yin(38) = yin(39) + dykl*yin(37)
                                      zin(38) = zin(39) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   40

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   35

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   41

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   48

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   47

                                      xin(48) = xin(48) + dxkl*xin(47)
                                      yin(48) = yin(48) + dykl*yin(47)
                                      zin(48) = zin(48) + dzkl*zin(47)

                                      ! i3 = i4 =   47
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   42

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   42

                                      ! do nk = 1,    3

                                      xin(42) = xin(43) + dxkl*xin(41)
                                      yin(42) = yin(43) + dykl*yin(41)
                                      zin(42) = zin(43) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    2

                                      xin(44) = xin(45) + dxkl*xin(43)
                                      yin(44) = yin(45) + dykl*yin(43)
                                      zin(44) = zin(45) + dzkl*zin(43)
                                      ! i4 = i4 + lang+1 =   46

                                      ! nk =    3

                                      xin(46) = xin(47) + dxkl*xin(45)
                                      yin(46) = yin(47) + dykl*yin(45)
                                      zin(46) = zin(47) + dzkl*zin(45)
                                      ! i4 = i4 + lang+1 =   48

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   43

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   49

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   56

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   55

                                      xin(56) = xin(56) + dxkl*xin(55)
                                      yin(56) = yin(56) + dykl*yin(55)
                                      zin(56) = zin(56) + dzkl*zin(55)

                                      ! i3 = i4 =   55
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   50

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   50

                                      ! do nk = 1,    3

                                      xin(50) = xin(51) + dxkl*xin(49)
                                      yin(50) = yin(51) + dykl*yin(49)
                                      zin(50) = zin(51) + dzkl*zin(49)
                                      ! i4 = i4 + lang+1 =   52

                                      ! nk =    2

                                      xin(52) = xin(53) + dxkl*xin(51)
                                      yin(52) = yin(53) + dykl*yin(51)
                                      zin(52) = zin(53) + dzkl*zin(51)
                                      ! i4 = i4 + lang+1 =   54

                                      ! nk =    3

                                      xin(54) = xin(55) + dxkl*xin(53)
                                      yin(54) = yin(55) + dykl*yin(53)
                                      zin(54) = zin(55) + dzkl*zin(53)
                                      ! i4 = i4 + lang+1 =   56

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   51

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   57

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   64

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   63

                                      xin(64) = xin(64) + dxkl*xin(63)
                                      yin(64) = yin(64) + dykl*yin(63)
                                      zin(64) = zin(64) + dzkl*zin(63)

                                      ! i3 = i4 =   63
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   58

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   58

                                      ! do nk = 1,    3

                                      xin(58) = xin(59) + dxkl*xin(57)
                                      yin(58) = yin(59) + dykl*yin(57)
                                      zin(58) = zin(59) + dzkl*zin(57)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      xin(60) = xin(61) + dxkl*xin(59)
                                      yin(60) = yin(61) + dykl*yin(59)
                                      zin(60) = zin(61) + dzkl*zin(59)
                                      ! i4 = i4 + lang+1 =   62

                                      ! nk =    3

                                      xin(62) = xin(63) + dxkl*xin(61)
                                      yin(62) = yin(63) + dykl*yin(61)
                                      zin(62) = zin(63) + dzkl*zin(61)
                                      ! i4 = i4 + lang+1 =   64

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   59

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   65

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   66

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    3

                                      xin(66) = xin(67) + dxkl*xin(65)
                                      yin(66) = yin(67) + dykl*yin(65)
                                      zin(66) = zin(67) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      xin(68) = xin(69) + dxkl*xin(67)
                                      yin(68) = yin(69) + dykl*yin(67)
                                      zin(68) = zin(69) + dzkl*zin(67)
                                      ! i4 = i4 + lang+1 =   70

                                      ! nk =    3

                                      xin(70) = xin(71) + dxkl*xin(69)
                                      yin(70) = yin(71) + dykl*yin(69)
                                      zin(70) = zin(71) + dzkl*zin(69)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   80

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    3

                                      xin(74) = xin(75) + dxkl*xin(73)
                                      yin(74) = yin(75) + dykl*yin(73)
                                      zin(74) = zin(75) + dzkl*zin(73)
                                      ! i4 = i4 + lang+1 =   76

                                      ! nk =    2

                                      xin(76) = xin(77) + dxkl*xin(75)
                                      yin(76) = yin(77) + dykl*yin(75)
                                      zin(76) = zin(77) + dzkl*zin(75)
                                      ! i4 = i4 + lang+1 =   78

                                      ! nk =    3

                                      xin(78) = xin(79) + dxkl*xin(77)
                                      yin(78) = yin(79) + dykl*yin(77)
                                      zin(78) = zin(79) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   81

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   88

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   87

                                      xin(88) = xin(88) + dxkl*xin(87)
                                      yin(88) = yin(88) + dykl*yin(87)
                                      zin(88) = zin(88) + dzkl*zin(87)

                                      ! i3 = i4 =   87
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   82

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   82

                                      ! do nk = 1,    3

                                      xin(82) = xin(83) + dxkl*xin(81)
                                      yin(82) = yin(83) + dykl*yin(81)
                                      zin(82) = zin(83) + dzkl*zin(81)
                                      ! i4 = i4 + lang+1 =   84

                                      ! nk =    2

                                      xin(84) = xin(85) + dxkl*xin(83)
                                      yin(84) = yin(85) + dykl*yin(83)
                                      zin(84) = zin(85) + dzkl*zin(83)
                                      ! i4 = i4 + lang+1 =   86

                                      ! nk =    3

                                      xin(86) = xin(87) + dxkl*xin(85)
                                      yin(86) = yin(87) + dykl*yin(85)
                                      zin(86) = zin(87) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   88

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   83

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   89

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   96

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   95

                                      xin(96) = xin(96) + dxkl*xin(95)
                                      yin(96) = yin(96) + dykl*yin(95)
                                      zin(96) = zin(96) + dzkl*zin(95)

                                      ! i3 = i4 =   95
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   90

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   90

                                      ! do nk = 1,    3

                                      xin(90) = xin(91) + dxkl*xin(89)
                                      yin(90) = yin(91) + dykl*yin(89)
                                      zin(90) = zin(91) + dzkl*zin(89)
                                      ! i4 = i4 + lang+1 =   92

                                      ! nk =    2

                                      xin(92) = xin(93) + dxkl*xin(91)
                                      yin(92) = yin(93) + dykl*yin(91)
                                      zin(92) = zin(93) + dzkl*zin(91)
                                      ! i4 = i4 + lang+1 =   94

                                      ! nk =    3

                                      xin(94) = xin(95) + dxkl*xin(93)
                                      yin(94) = yin(95) + dykl*yin(93)
                                      zin(94) = zin(95) + dzkl*zin(93)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   91

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   97

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  104

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  103

                                      xin(104) = xin(104) + dxkl*xin(103)
                                      yin(104) = yin(104) + dykl*yin(103)
                                      zin(104) = zin(104) + dzkl*zin(103)

                                      ! i3 = i4 =  103
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   98

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   98

                                      ! do nk = 1,    3

                                      xin(98) = xin(99) + dxkl*xin(97)
                                      yin(98) = yin(99) + dykl*yin(97)
                                      zin(98) = zin(99) + dzkl*zin(97)
                                      ! i4 = i4 + lang+1 =  100

                                      ! nk =    2

                                      xin(100) = xin(101) + dxkl*xin(99)
                                      yin(100) = yin(101) + dykl*yin(99)
                                      zin(100) = zin(101) + dzkl*zin(99)
                                      ! i4 = i4 + lang+1 =  102

                                      ! nk =    3

                                      xin(102) = xin(103) + dxkl*xin(101)
                                      yin(102) = yin(103) + dykl*yin(101)
                                      zin(102) = zin(103) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   99

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  105

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  112

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  111

                                      xin(112) = xin(112) + dxkl*xin(111)
                                      yin(112) = yin(112) + dykl*yin(111)
                                      zin(112) = zin(112) + dzkl*zin(111)

                                      ! i3 = i4 =  111
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  106

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  106

                                      ! do nk = 1,    3

                                      xin(106) = xin(107) + dxkl*xin(105)
                                      yin(106) = yin(107) + dykl*yin(105)
                                      zin(106) = zin(107) + dzkl*zin(105)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    2

                                      xin(108) = xin(109) + dxkl*xin(107)
                                      yin(108) = yin(109) + dykl*yin(107)
                                      zin(108) = zin(109) + dzkl*zin(107)
                                      ! i4 = i4 + lang+1 =  110

                                      ! nk =    3

                                      xin(110) = xin(111) + dxkl*xin(109)
                                      yin(110) = yin(111) + dykl*yin(109)
                                      zin(110) = zin(111) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  112

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  107

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  113

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  120

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  119

                                      xin(120) = xin(120) + dxkl*xin(119)
                                      yin(120) = yin(120) + dykl*yin(119)
                                      zin(120) = zin(120) + dzkl*zin(119)

                                      ! i3 = i4 =  119
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  114

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  114

                                      ! do nk = 1,    3

                                      xin(114) = xin(115) + dxkl*xin(113)
                                      yin(114) = yin(115) + dykl*yin(113)
                                      zin(114) = zin(115) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    2

                                      xin(116) = xin(117) + dxkl*xin(115)
                                      yin(116) = yin(117) + dykl*yin(115)
                                      zin(116) = zin(117) + dzkl*zin(115)
                                      ! i4 = i4 + lang+1 =  118

                                      ! nk =    3

                                      xin(118) = xin(119) + dxkl*xin(117)
                                      yin(118) = yin(119) + dykl*yin(117)
                                      zin(118) = zin(119) + dzkl*zin(117)
                                      ! i4 = i4 + lang+1 =  120

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  115

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  121

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  128

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  127

                                      xin(128) = xin(128) + dxkl*xin(127)
                                      yin(128) = yin(128) + dykl*yin(127)
                                      zin(128) = zin(128) + dzkl*zin(127)

                                      ! i3 = i4 =  127
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  122

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  122

                                      ! do nk = 1,    3

                                      xin(122) = xin(123) + dxkl*xin(121)
                                      yin(122) = yin(123) + dykl*yin(121)
                                      zin(122) = zin(123) + dzkl*zin(121)
                                      ! i4 = i4 + lang+1 =  124

                                      ! nk =    2

                                      xin(124) = xin(125) + dxkl*xin(123)
                                      yin(124) = yin(125) + dykl*yin(123)
                                      zin(124) = zin(125) + dzkl*zin(123)
                                      ! i4 = i4 + lang+1 =  126

                                      ! nk =    3

                                      xin(126) = xin(127) + dxkl*xin(125)
                                      yin(126) = yin(127) + dykl*yin(125)
                                      zin(126) = zin(127) + dzkl*zin(125)
                                      ! i4 = i4 + lang+1 =  128

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  123

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  129

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  128

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

                                      ! i1 = in(1) =  129

                                      xin(129) = 1.0_dp
                                      yin(129) = 1.0_dp
                                      zin(129) = f00

                                      ! i2 = in(2) =  161
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(161) = xc00
                                      yin(161) = yc00
                                      zin(161) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  131

                                      xin(131) = xcp00
                                      yin(131) = ycp00
                                      zin(131) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  163
                                      ! i2 =  161

                                      xin(163) = xcp00*xin(161) + cp10
                                      yin(163) = ycp00*yin(161) + cp10
                                      zin(163) = zcp00*zin(161) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  161

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  193
                                      ! i3 =  129
                                      ! i4 =  161

                                      xin(193) = c10*xin(129) + xc00*xin(161)
                                      yin(193) = c10*yin(129) + yc00*yin(161)
                                      zin(193) = c10*zin(129) + zc00*zin(161)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  195
                                      ! i5 =  193
                                      ! i4 =  161

                                      xin(195) = xcp00*xin(193) + cp10*xin(161)
                                      yin(195) = ycp00*yin(193) + cp10*yin(161)
                                      zin(195) = zcp00*zin(193) + cp10*zin(161)

                                      ! ------------------

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  193

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  225
                                      ! i3 =  161
                                      ! i4 =  193

                                      xin(225) = c10*xin(161) + xc00*xin(193)
                                      yin(225) = c10*yin(161) + yc00*yin(193)
                                      zin(225) = c10*zin(161) + zc00*zin(193)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  227
                                      ! i5 =  225
                                      ! i4 =  193

                                      xin(227) = xcp00*xin(225) + cp10*xin(193)
                                      yin(227) = ycp00*yin(225) + cp10*yin(193)
                                      zin(227) = zcp00*zin(225) + cp10*zin(193)

                                      ! ------------------

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  225

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  233
                                      ! i3 =  193
                                      ! i4 =  225

                                      xin(233) = c10*xin(193) + xc00*xin(225)
                                      yin(233) = c10*yin(193) + yc00*yin(225)
                                      zin(233) = c10*zin(193) + zc00*zin(225)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  235
                                      ! i5 =  233
                                      ! i4 =  225

                                      xin(235) = xcp00*xin(233) + cp10*xin(225)
                                      yin(235) = ycp00*yin(233) + cp10*yin(225)
                                      zin(235) = zcp00*zin(233) + cp10*zin(225)

                                      ! ------------------

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  233

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  241
                                      ! i3 =  225
                                      ! i4 =  233

                                      xin(241) = c10*xin(225) + xc00*xin(233)
                                      yin(241) = c10*yin(225) + yc00*yin(233)
                                      zin(241) = c10*zin(225) + zc00*zin(233)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  243
                                      ! i5 =  241
                                      ! i4 =  233

                                      xin(243) = xcp00*xin(241) + cp10*xin(233)
                                      yin(243) = ycp00*yin(241) + cp10*yin(233)
                                      zin(243) = zcp00*zin(241) + cp10*zin(233)

                                      ! ------------------

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  241

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  249
                                      ! i3 =  233
                                      ! i4 =  241

                                      xin(249) = c10*xin(233) + xc00*xin(241)
                                      yin(249) = c10*yin(233) + yc00*yin(241)
                                      zin(249) = c10*zin(233) + zc00*zin(241)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  251
                                      ! i5 =  249
                                      ! i4 =  241

                                      xin(251) = xcp00*xin(249) + cp10*xin(241)
                                      yin(251) = ycp00*yin(249) + cp10*yin(241)
                                      zin(251) = zcp00*zin(249) + cp10*zin(241)

                                      ! ------------------

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  249

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  129
                                      ! i4 = i1+k2 =  131

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  133
                                      ! i3 =  129
                                      ! i4 =  131

                                      xin(133) = cp01*xin(129) + xcp00*xin(131)
                                      yin(133) = cp01*yin(129) + ycp00*yin(131)
                                      zin(133) = cp01*zin(129) + zcp00*zin(131)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  165

                                      xin(165) = xc00*xin(133) + c01*xin(131)
                                      yin(165) = yc00*yin(133) + c01*yin(131)
                                      zin(165) = zc00*zin(133) + c01*zin(131)

                                      ! ------------------

                                      ! i3 = i4 =  131
                                      ! i4 = i5 =  133

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  135
                                      ! i3 =  131
                                      ! i4 =  133

                                      xin(135) = cp01*xin(131) + xcp00*xin(133)
                                      yin(135) = cp01*yin(131) + ycp00*yin(133)
                                      zin(135) = cp01*zin(131) + zcp00*zin(133)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  167

                                      xin(167) = xc00*xin(135) + c01*xin(133)
                                      yin(167) = yc00*yin(135) + c01*yin(133)
                                      zin(167) = zc00*zin(135) + c01*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  135

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  136
                                      ! i3 =  133
                                      ! i4 =  135

                                      xin(136) = cp01*xin(133) + xcp00*xin(135)
                                      yin(136) = cp01*yin(133) + ycp00*yin(135)
                                      zin(136) = cp01*zin(133) + zcp00*zin(135)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  168

                                      xin(168) = xc00*xin(136) + c01*xin(135)
                                      yin(168) = yc00*yin(136) + c01*yin(135)
                                      zin(168) = zc00*zin(136) + c01*zin(135)

                                      ! ------------------

                                      ! i3 = i4 =  135
                                      ! i4 = i5 =  136

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  161

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  193

                                      xin(197) = c10*xin(133) + xc00*xin(165) + c01*xin(163)
                                      yin(197) = c10*yin(133) + yc00*yin(165) + c01*yin(163)
                                      zin(197) = c10*zin(133) + zc00*zin(165) + c01*zin(163)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  225

                                      xin(229) = c10*xin(165) + xc00*xin(197) + c01*xin(195)
                                      yin(229) = c10*yin(165) + yc00*yin(197) + c01*yin(195)
                                      zin(229) = c10*zin(165) + zc00*zin(197) + c01*zin(195)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  225

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  233

                                      xin(237) = c10*xin(197) + xc00*xin(229) + c01*xin(227)
                                      yin(237) = c10*yin(197) + yc00*yin(229) + c01*yin(227)
                                      zin(237) = c10*zin(197) + zc00*zin(229) + c01*zin(227)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  233

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  241

                                      xin(245) = c10*xin(229) + xc00*xin(237) + c01*xin(235)
                                      yin(245) = c10*yin(229) + yc00*yin(237) + c01*yin(235)
                                      zin(245) = c10*zin(229) + zc00*zin(237) + c01*zin(235)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  241

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  249

                                      xin(253) = c10*xin(237) + xc00*xin(245) + c01*xin(243)
                                      yin(253) = c10*yin(237) + yc00*yin(245) + c01*yin(243)
                                      zin(253) = c10*zin(237) + zc00*zin(245) + c01*zin(243)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  249

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  161

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  193

                                      xin(199) = c10*xin(135) + xc00*xin(167) + c01*xin(165)
                                      yin(199) = c10*yin(135) + yc00*yin(167) + c01*yin(165)
                                      zin(199) = c10*zin(135) + zc00*zin(167) + c01*zin(165)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  225

                                      xin(231) = c10*xin(167) + xc00*xin(199) + c01*xin(197)
                                      yin(231) = c10*yin(167) + yc00*yin(199) + c01*yin(197)
                                      zin(231) = c10*zin(167) + zc00*zin(199) + c01*zin(197)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  225

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  233

                                      xin(239) = c10*xin(199) + xc00*xin(231) + c01*xin(229)
                                      yin(239) = c10*yin(199) + yc00*yin(231) + c01*yin(229)
                                      zin(239) = c10*zin(199) + zc00*zin(231) + c01*zin(229)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  233

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  241

                                      xin(247) = c10*xin(231) + xc00*xin(239) + c01*xin(237)
                                      yin(247) = c10*yin(231) + yc00*yin(239) + c01*yin(237)
                                      zin(247) = c10*zin(231) + zc00*zin(239) + c01*zin(237)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  241

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  249

                                      xin(255) = c10*xin(239) + xc00*xin(247) + c01*xin(245)
                                      yin(255) = c10*yin(239) + yc00*yin(247) + c01*yin(245)
                                      zin(255) = c10*zin(239) + zc00*zin(247) + c01*zin(245)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  249

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  161

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  193

                                      xin(200) = c10*xin(136) + xc00*xin(168) + c01*xin(167)
                                      yin(200) = c10*yin(136) + yc00*yin(168) + c01*yin(167)
                                      zin(200) = c10*zin(136) + zc00*zin(168) + c01*zin(167)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  225

                                      xin(232) = c10*xin(168) + xc00*xin(200) + c01*xin(199)
                                      yin(232) = c10*yin(168) + yc00*yin(200) + c01*yin(199)
                                      zin(232) = c10*zin(168) + zc00*zin(200) + c01*zin(199)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  225

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  233

                                      xin(240) = c10*xin(200) + xc00*xin(232) + c01*xin(231)
                                      yin(240) = c10*yin(200) + yc00*yin(232) + c01*yin(231)
                                      zin(240) = c10*zin(200) + zc00*zin(232) + c01*zin(231)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  233

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  241

                                      xin(248) = c10*xin(232) + xc00*xin(240) + c01*xin(239)
                                      yin(248) = c10*yin(232) + yc00*yin(240) + c01*yin(239)
                                      zin(248) = c10*zin(232) + zc00*zin(240) + c01*zin(239)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  241

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  249

                                      xin(256) = c10*xin(240) + xc00*xin(248) + c01*xin(247)
                                      yin(256) = c10*yin(240) + yc00*yin(248) + c01*yin(247)
                                      zin(256) = c10*zin(240) + zc00*zin(248) + c01*zin(247)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  249

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  249

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  249

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  241

                                      xin(249) = xin(249) + dxij*xin(241)
                                      yin(249) = yin(249) + dyij*yin(241)
                                      zin(249) = zin(249) + dzij*zin(241)

                                      ! i3 = i4 =  241
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  233

                                      xin(241) = xin(241) + dxij*xin(233)
                                      yin(241) = yin(241) + dyij*yin(233)
                                      zin(241) = zin(241) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  225

                                      xin(233) = xin(233) + dxij*xin(225)
                                      yin(233) = yin(233) + dyij*yin(225)
                                      zin(233) = zin(233) + dzij*zin(225)

                                      ! i3 = i4 =  225
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  249

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  241

                                      xin(249) = xin(249) + dxij*xin(241)
                                      yin(249) = yin(249) + dyij*yin(241)
                                      zin(249) = zin(249) + dzij*zin(241)

                                      ! i3 = i4 =  241
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  233

                                      xin(241) = xin(241) + dxij*xin(233)
                                      yin(241) = yin(241) + dyij*yin(233)
                                      zin(241) = zin(241) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  249

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  241

                                      xin(249) = xin(249) + dxij*xin(241)
                                      yin(249) = yin(249) + dyij*yin(241)
                                      zin(249) = zin(249) + dzij*zin(241)

                                      ! i3 = i4 =  241
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  137

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  137

                                      ! do ni = 1,    3

                                      xin(137) = xin(161) + dxij*xin(129)
                                      yin(137) = yin(161) + dyij*yin(129)
                                      zin(137) = zin(161) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! ni =    2

                                      xin(169) = xin(193) + dxij*xin(161)
                                      yin(169) = yin(193) + dyij*yin(161)
                                      zin(169) = zin(193) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  201

                                      ! ni =    3

                                      xin(201) = xin(225) + dxij*xin(193)
                                      yin(201) = yin(225) + dyij*yin(193)
                                      zin(201) = zin(225) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  145

                                      ! nj =    2

                                      ! i4 = i3 =  145

                                      ! do ni = 1,    3

                                      xin(145) = xin(169) + dxij*xin(137)
                                      yin(145) = yin(169) + dyij*yin(137)
                                      zin(145) = zin(169) + dzij*zin(137)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    2

                                      xin(177) = xin(201) + dxij*xin(169)
                                      yin(177) = yin(201) + dyij*yin(169)
                                      zin(177) = zin(201) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    3

                                      xin(209) = xin(233) + dxij*xin(201)
                                      yin(209) = yin(233) + dyij*yin(201)
                                      zin(209) = zin(233) + dzij*zin(201)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  153

                                      ! nj =    3

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(177) + dxij*xin(145)
                                      yin(153) = yin(177) + dyij*yin(145)
                                      zin(153) = zin(177) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    2

                                      xin(185) = xin(209) + dxij*xin(177)
                                      yin(185) = yin(209) + dyij*yin(177)
                                      zin(185) = zin(209) + dzij*zin(177)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! ni =    3

                                      xin(217) = xin(241) + dxij*xin(209)
                                      yin(217) = yin(241) + dyij*yin(209)
                                      zin(217) = zin(241) + dzij*zin(209)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  249

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  161

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  251

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  243

                                      xin(251) = xin(251) + dxij*xin(243)
                                      yin(251) = yin(251) + dyij*yin(243)
                                      zin(251) = zin(251) + dzij*zin(243)

                                      ! i3 = i4 =  243
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  235

                                      xin(243) = xin(243) + dxij*xin(235)
                                      yin(243) = yin(243) + dyij*yin(235)
                                      zin(243) = zin(243) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  227

                                      xin(235) = xin(235) + dxij*xin(227)
                                      yin(235) = yin(235) + dyij*yin(227)
                                      zin(235) = zin(235) + dzij*zin(227)

                                      ! i3 = i4 =  227
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  251

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  243

                                      xin(251) = xin(251) + dxij*xin(243)
                                      yin(251) = yin(251) + dyij*yin(243)
                                      zin(251) = zin(251) + dzij*zin(243)

                                      ! i3 = i4 =  243
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  235

                                      xin(243) = xin(243) + dxij*xin(235)
                                      yin(243) = yin(243) + dyij*yin(235)
                                      zin(243) = zin(243) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  251

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  243

                                      xin(251) = xin(251) + dxij*xin(243)
                                      yin(251) = yin(251) + dyij*yin(243)
                                      zin(251) = zin(251) + dzij*zin(243)

                                      ! i3 = i4 =  243
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  139

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  139

                                      ! do ni = 1,    3

                                      xin(139) = xin(163) + dxij*xin(131)
                                      yin(139) = yin(163) + dyij*yin(131)
                                      zin(139) = zin(163) + dzij*zin(131)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  171

                                      ! ni =    2

                                      xin(171) = xin(195) + dxij*xin(163)
                                      yin(171) = yin(195) + dyij*yin(163)
                                      zin(171) = zin(195) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  203

                                      ! ni =    3

                                      xin(203) = xin(227) + dxij*xin(195)
                                      yin(203) = yin(227) + dyij*yin(195)
                                      zin(203) = zin(227) + dzij*zin(195)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  147

                                      ! nj =    2

                                      ! i4 = i3 =  147

                                      ! do ni = 1,    3

                                      xin(147) = xin(171) + dxij*xin(139)
                                      yin(147) = yin(171) + dyij*yin(139)
                                      zin(147) = zin(171) + dzij*zin(139)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    2

                                      xin(179) = xin(203) + dxij*xin(171)
                                      yin(179) = yin(203) + dyij*yin(171)
                                      zin(179) = zin(203) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    3

                                      xin(211) = xin(235) + dxij*xin(203)
                                      yin(211) = yin(235) + dyij*yin(203)
                                      zin(211) = zin(235) + dzij*zin(203)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  243

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  155

                                      ! nj =    3

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(179) + dxij*xin(147)
                                      yin(155) = yin(179) + dyij*yin(147)
                                      zin(155) = zin(179) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    2

                                      xin(187) = xin(211) + dxij*xin(179)
                                      yin(187) = yin(211) + dyij*yin(179)
                                      zin(187) = zin(211) + dzij*zin(179)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  219

                                      ! ni =    3

                                      xin(219) = xin(243) + dxij*xin(211)
                                      yin(219) = yin(243) + dyij*yin(211)
                                      zin(219) = zin(243) + dzij*zin(211)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  163

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  245

                                      xin(253) = xin(253) + dxij*xin(245)
                                      yin(253) = yin(253) + dyij*yin(245)
                                      zin(253) = zin(253) + dzij*zin(245)

                                      ! i3 = i4 =  245
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  237

                                      xin(245) = xin(245) + dxij*xin(237)
                                      yin(245) = yin(245) + dyij*yin(237)
                                      zin(245) = zin(245) + dzij*zin(237)

                                      ! i3 = i4 =  237
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  229

                                      xin(237) = xin(237) + dxij*xin(229)
                                      yin(237) = yin(237) + dyij*yin(229)
                                      zin(237) = zin(237) + dzij*zin(229)

                                      ! i3 = i4 =  229
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  245

                                      xin(253) = xin(253) + dxij*xin(245)
                                      yin(253) = yin(253) + dyij*yin(245)
                                      zin(253) = zin(253) + dzij*zin(245)

                                      ! i3 = i4 =  245
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  237

                                      xin(245) = xin(245) + dxij*xin(237)
                                      yin(245) = yin(245) + dyij*yin(237)
                                      zin(245) = zin(245) + dzij*zin(237)

                                      ! i3 = i4 =  237
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  245

                                      xin(253) = xin(253) + dxij*xin(245)
                                      yin(253) = yin(253) + dyij*yin(245)
                                      zin(253) = zin(253) + dzij*zin(245)

                                      ! i3 = i4 =  245
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  141

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  141

                                      ! do ni = 1,    3

                                      xin(141) = xin(165) + dxij*xin(133)
                                      yin(141) = yin(165) + dyij*yin(133)
                                      zin(141) = zin(165) + dzij*zin(133)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    2

                                      xin(173) = xin(197) + dxij*xin(165)
                                      yin(173) = yin(197) + dyij*yin(165)
                                      zin(173) = zin(197) + dzij*zin(165)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    3

                                      xin(205) = xin(229) + dxij*xin(197)
                                      yin(205) = yin(229) + dyij*yin(197)
                                      zin(205) = zin(229) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  237

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  149

                                      ! nj =    2

                                      ! i4 = i3 =  149

                                      ! do ni = 1,    3

                                      xin(149) = xin(173) + dxij*xin(141)
                                      yin(149) = yin(173) + dyij*yin(141)
                                      zin(149) = zin(173) + dzij*zin(141)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! ni =    2

                                      xin(181) = xin(205) + dxij*xin(173)
                                      yin(181) = yin(205) + dyij*yin(173)
                                      zin(181) = zin(205) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    3

                                      xin(213) = xin(237) + dxij*xin(205)
                                      yin(213) = yin(237) + dyij*yin(205)
                                      zin(213) = zin(237) + dzij*zin(205)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  245

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    3

                                      ! i4 = i3 =  157

                                      ! do ni = 1,    3

                                      xin(157) = xin(181) + dxij*xin(149)
                                      yin(157) = yin(181) + dyij*yin(149)
                                      zin(157) = zin(181) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    2

                                      xin(189) = xin(213) + dxij*xin(181)
                                      yin(189) = yin(213) + dyij*yin(181)
                                      zin(189) = zin(213) + dzij*zin(181)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  221

                                      ! ni =    3

                                      xin(221) = xin(245) + dxij*xin(213)
                                      yin(221) = yin(245) + dyij*yin(213)
                                      zin(221) = zin(245) + dzij*zin(213)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  165

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  247

                                      xin(255) = xin(255) + dxij*xin(247)
                                      yin(255) = yin(255) + dyij*yin(247)
                                      zin(255) = zin(255) + dzij*zin(247)

                                      ! i3 = i4 =  247
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  239

                                      xin(247) = xin(247) + dxij*xin(239)
                                      yin(247) = yin(247) + dyij*yin(239)
                                      zin(247) = zin(247) + dzij*zin(239)

                                      ! i3 = i4 =  239
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  231

                                      xin(239) = xin(239) + dxij*xin(231)
                                      yin(239) = yin(239) + dyij*yin(231)
                                      zin(239) = zin(239) + dzij*zin(231)

                                      ! i3 = i4 =  231
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  247

                                      xin(255) = xin(255) + dxij*xin(247)
                                      yin(255) = yin(255) + dyij*yin(247)
                                      zin(255) = zin(255) + dzij*zin(247)

                                      ! i3 = i4 =  247
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  239

                                      xin(247) = xin(247) + dxij*xin(239)
                                      yin(247) = yin(247) + dyij*yin(239)
                                      zin(247) = zin(247) + dzij*zin(239)

                                      ! i3 = i4 =  239
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  247

                                      xin(255) = xin(255) + dxij*xin(247)
                                      yin(255) = yin(255) + dyij*yin(247)
                                      zin(255) = zin(255) + dzij*zin(247)

                                      ! i3 = i4 =  247
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  143

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  143

                                      ! do ni = 1,    3

                                      xin(143) = xin(167) + dxij*xin(135)
                                      yin(143) = yin(167) + dyij*yin(135)
                                      zin(143) = zin(167) + dzij*zin(135)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    2

                                      xin(175) = xin(199) + dxij*xin(167)
                                      yin(175) = yin(199) + dyij*yin(167)
                                      zin(175) = zin(199) + dzij*zin(167)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    3

                                      xin(207) = xin(231) + dxij*xin(199)
                                      yin(207) = yin(231) + dyij*yin(199)
                                      zin(207) = zin(231) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  151

                                      ! nj =    2

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(175) + dxij*xin(143)
                                      yin(151) = yin(175) + dyij*yin(143)
                                      zin(151) = zin(175) + dzij*zin(143)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  183

                                      ! ni =    2

                                      xin(183) = xin(207) + dxij*xin(175)
                                      yin(183) = yin(207) + dyij*yin(175)
                                      zin(183) = zin(207) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    3

                                      xin(215) = xin(239) + dxij*xin(207)
                                      yin(215) = yin(239) + dyij*yin(207)
                                      zin(215) = zin(239) + dzij*zin(207)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    3

                                      ! i4 = i3 =  159

                                      ! do ni = 1,    3

                                      xin(159) = xin(183) + dxij*xin(151)
                                      yin(159) = yin(183) + dyij*yin(151)
                                      zin(159) = zin(183) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    2

                                      xin(191) = xin(215) + dxij*xin(183)
                                      yin(191) = yin(215) + dyij*yin(183)
                                      zin(191) = zin(215) + dzij*zin(183)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  223

                                      ! ni =    3

                                      xin(223) = xin(247) + dxij*xin(215)
                                      yin(223) = yin(247) + dyij*yin(215)
                                      zin(223) = zin(247) + dzij*zin(215)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  255

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  167

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  248

                                      xin(256) = xin(256) + dxij*xin(248)
                                      yin(256) = yin(256) + dyij*yin(248)
                                      zin(256) = zin(256) + dzij*zin(248)

                                      ! i3 = i4 =  248
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  240

                                      xin(248) = xin(248) + dxij*xin(240)
                                      yin(248) = yin(248) + dyij*yin(240)
                                      zin(248) = zin(248) + dzij*zin(240)

                                      ! i3 = i4 =  240
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  232

                                      xin(240) = xin(240) + dxij*xin(232)
                                      yin(240) = yin(240) + dyij*yin(232)
                                      zin(240) = zin(240) + dzij*zin(232)

                                      ! i3 = i4 =  232
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  248

                                      xin(256) = xin(256) + dxij*xin(248)
                                      yin(256) = yin(256) + dyij*yin(248)
                                      zin(256) = zin(256) + dzij*zin(248)

                                      ! i3 = i4 =  248
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  240

                                      xin(248) = xin(248) + dxij*xin(240)
                                      yin(248) = yin(248) + dyij*yin(240)
                                      zin(248) = zin(248) + dzij*zin(240)

                                      ! i3 = i4 =  240
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  248

                                      xin(256) = xin(256) + dxij*xin(248)
                                      yin(256) = yin(256) + dyij*yin(248)
                                      zin(256) = zin(256) + dzij*zin(248)

                                      ! i3 = i4 =  248
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  144

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  144

                                      ! do ni = 1,    3

                                      xin(144) = xin(168) + dxij*xin(136)
                                      yin(144) = yin(168) + dyij*yin(136)
                                      zin(144) = zin(168) + dzij*zin(136)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  176

                                      ! ni =    2

                                      xin(176) = xin(200) + dxij*xin(168)
                                      yin(176) = yin(200) + dyij*yin(168)
                                      zin(176) = zin(200) + dzij*zin(168)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                      ! ni =    3

                                      xin(208) = xin(232) + dxij*xin(200)
                                      yin(208) = yin(232) + dyij*yin(200)
                                      zin(208) = zin(232) + dzij*zin(200)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  152

                                      ! nj =    2

                                      ! i4 = i3 =  152

                                      ! do ni = 1,    3

                                      xin(152) = xin(176) + dxij*xin(144)
                                      yin(152) = yin(176) + dyij*yin(144)
                                      zin(152) = zin(176) + dzij*zin(144)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  184

                                      ! ni =    2

                                      xin(184) = xin(208) + dxij*xin(176)
                                      yin(184) = yin(208) + dyij*yin(176)
                                      zin(184) = zin(208) + dzij*zin(176)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    3

                                      xin(216) = xin(240) + dxij*xin(208)
                                      yin(216) = yin(240) + dyij*yin(208)
                                      zin(216) = zin(240) + dzij*zin(208)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  248

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  160

                                      ! nj =    3

                                      ! i4 = i3 =  160

                                      ! do ni = 1,    3

                                      xin(160) = xin(184) + dxij*xin(152)
                                      yin(160) = yin(184) + dyij*yin(152)
                                      zin(160) = zin(184) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    2

                                      xin(192) = xin(216) + dxij*xin(184)
                                      yin(192) = yin(216) + dyij*yin(184)
                                      zin(192) = zin(216) + dzij*zin(184)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  224

                                      ! ni =    3

                                      xin(224) = xin(248) + dxij*xin(216)
                                      yin(224) = yin(248) + dyij*yin(216)
                                      zin(224) = zin(248) + dzij*zin(216)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  256

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  168

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  129

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  136

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  135

                                      xin(136) = xin(136) + dxkl*xin(135)
                                      yin(136) = yin(136) + dykl*yin(135)
                                      zin(136) = zin(136) + dzkl*zin(135)

                                      ! i3 = i4 =  135
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  130

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  130

                                      ! do nk = 1,    3

                                      xin(130) = xin(131) + dxkl*xin(129)
                                      yin(130) = yin(131) + dykl*yin(129)
                                      zin(130) = zin(131) + dzkl*zin(129)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    2

                                      xin(132) = xin(133) + dxkl*xin(131)
                                      yin(132) = yin(133) + dykl*yin(131)
                                      zin(132) = zin(133) + dzkl*zin(131)
                                      ! i4 = i4 + lang+1 =  134

                                      ! nk =    3

                                      xin(134) = xin(135) + dxkl*xin(133)
                                      yin(134) = yin(135) + dykl*yin(133)
                                      zin(134) = zin(135) + dzkl*zin(133)
                                      ! i4 = i4 + lang+1 =  136

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  131

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  137

                                      ! nj = nj + 1 =    1

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  138

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  138

                                      ! do nk = 1,    3

                                      xin(138) = xin(139) + dxkl*xin(137)
                                      yin(138) = yin(139) + dykl*yin(137)
                                      zin(138) = zin(139) + dzkl*zin(137)
                                      ! i4 = i4 + lang+1 =  140

                                      ! nk =    2

                                      xin(140) = xin(141) + dxkl*xin(139)
                                      yin(140) = yin(141) + dykl*yin(139)
                                      zin(140) = zin(141) + dzkl*zin(139)
                                      ! i4 = i4 + lang+1 =  142

                                      ! nk =    3

                                      xin(142) = xin(143) + dxkl*xin(141)
                                      yin(142) = yin(143) + dykl*yin(141)
                                      zin(142) = zin(143) + dzkl*zin(141)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  139

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  152

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  151

                                      xin(152) = xin(152) + dxkl*xin(151)
                                      yin(152) = yin(152) + dykl*yin(151)
                                      zin(152) = zin(152) + dzkl*zin(151)

                                      ! i3 = i4 =  151
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    3

                                      xin(146) = xin(147) + dxkl*xin(145)
                                      yin(146) = yin(147) + dykl*yin(145)
                                      zin(146) = zin(147) + dzkl*zin(145)
                                      ! i4 = i4 + lang+1 =  148

                                      ! nk =    2

                                      xin(148) = xin(149) + dxkl*xin(147)
                                      yin(148) = yin(149) + dykl*yin(147)
                                      zin(148) = zin(149) + dzkl*zin(147)
                                      ! i4 = i4 + lang+1 =  150

                                      ! nk =    3

                                      xin(150) = xin(151) + dxkl*xin(149)
                                      yin(150) = yin(151) + dykl*yin(149)
                                      zin(150) = zin(151) + dzkl*zin(149)
                                      ! i4 = i4 + lang+1 =  152

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  153

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  160

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  159

                                      xin(160) = xin(160) + dxkl*xin(159)
                                      yin(160) = yin(160) + dykl*yin(159)
                                      zin(160) = zin(160) + dzkl*zin(159)

                                      ! i3 = i4 =  159
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  154

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  154

                                      ! do nk = 1,    3

                                      xin(154) = xin(155) + dxkl*xin(153)
                                      yin(154) = yin(155) + dykl*yin(153)
                                      zin(154) = zin(155) + dzkl*zin(153)
                                      ! i4 = i4 + lang+1 =  156

                                      ! nk =    2

                                      xin(156) = xin(157) + dxkl*xin(155)
                                      yin(156) = yin(157) + dykl*yin(155)
                                      zin(156) = zin(157) + dzkl*zin(155)
                                      ! i4 = i4 + lang+1 =  158

                                      ! nk =    3

                                      xin(158) = xin(159) + dxkl*xin(157)
                                      yin(158) = yin(159) + dykl*yin(157)
                                      zin(158) = zin(159) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  160

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  155

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  161

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  168

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  167

                                      xin(168) = xin(168) + dxkl*xin(167)
                                      yin(168) = yin(168) + dykl*yin(167)
                                      zin(168) = zin(168) + dzkl*zin(167)

                                      ! i3 = i4 =  167
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  162

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  162

                                      ! do nk = 1,    3

                                      xin(162) = xin(163) + dxkl*xin(161)
                                      yin(162) = yin(163) + dykl*yin(161)
                                      zin(162) = zin(163) + dzkl*zin(161)
                                      ! i4 = i4 + lang+1 =  164

                                      ! nk =    2

                                      xin(164) = xin(165) + dxkl*xin(163)
                                      yin(164) = yin(165) + dykl*yin(163)
                                      zin(164) = zin(165) + dzkl*zin(163)
                                      ! i4 = i4 + lang+1 =  166

                                      ! nk =    3

                                      xin(166) = xin(167) + dxkl*xin(165)
                                      yin(166) = yin(167) + dykl*yin(165)
                                      zin(166) = zin(167) + dzkl*zin(165)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  163

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  169

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  176

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  175

                                      xin(176) = xin(176) + dxkl*xin(175)
                                      yin(176) = yin(176) + dykl*yin(175)
                                      zin(176) = zin(176) + dzkl*zin(175)

                                      ! i3 = i4 =  175
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  170

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  170

                                      ! do nk = 1,    3

                                      xin(170) = xin(171) + dxkl*xin(169)
                                      yin(170) = yin(171) + dykl*yin(169)
                                      zin(170) = zin(171) + dzkl*zin(169)
                                      ! i4 = i4 + lang+1 =  172

                                      ! nk =    2

                                      xin(172) = xin(173) + dxkl*xin(171)
                                      yin(172) = yin(173) + dykl*yin(171)
                                      zin(172) = zin(173) + dzkl*zin(171)
                                      ! i4 = i4 + lang+1 =  174

                                      ! nk =    3

                                      xin(174) = xin(175) + dxkl*xin(173)
                                      yin(174) = yin(175) + dykl*yin(173)
                                      zin(174) = zin(175) + dzkl*zin(173)
                                      ! i4 = i4 + lang+1 =  176

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  171

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  177

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  184

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  183

                                      xin(184) = xin(184) + dxkl*xin(183)
                                      yin(184) = yin(184) + dykl*yin(183)
                                      zin(184) = zin(184) + dzkl*zin(183)

                                      ! i3 = i4 =  183
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  178

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  178

                                      ! do nk = 1,    3

                                      xin(178) = xin(179) + dxkl*xin(177)
                                      yin(178) = yin(179) + dykl*yin(177)
                                      zin(178) = zin(179) + dzkl*zin(177)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    2

                                      xin(180) = xin(181) + dxkl*xin(179)
                                      yin(180) = yin(181) + dykl*yin(179)
                                      zin(180) = zin(181) + dzkl*zin(179)
                                      ! i4 = i4 + lang+1 =  182

                                      ! nk =    3

                                      xin(182) = xin(183) + dxkl*xin(181)
                                      yin(182) = yin(183) + dykl*yin(181)
                                      zin(182) = zin(183) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  184

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  179

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  185

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  192

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  191

                                      xin(192) = xin(192) + dxkl*xin(191)
                                      yin(192) = yin(192) + dykl*yin(191)
                                      zin(192) = zin(192) + dzkl*zin(191)

                                      ! i3 = i4 =  191
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  186

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  186

                                      ! do nk = 1,    3

                                      xin(186) = xin(187) + dxkl*xin(185)
                                      yin(186) = yin(187) + dykl*yin(185)
                                      zin(186) = zin(187) + dzkl*zin(185)
                                      ! i4 = i4 + lang+1 =  188

                                      ! nk =    2

                                      xin(188) = xin(189) + dxkl*xin(187)
                                      yin(188) = yin(189) + dykl*yin(187)
                                      zin(188) = zin(189) + dzkl*zin(187)
                                      ! i4 = i4 + lang+1 =  190

                                      ! nk =    3

                                      xin(190) = xin(191) + dxkl*xin(189)
                                      yin(190) = yin(191) + dykl*yin(189)
                                      zin(190) = zin(191) + dzkl*zin(189)
                                      ! i4 = i4 + lang+1 =  192

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  187

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  193

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  200

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  199

                                      xin(200) = xin(200) + dxkl*xin(199)
                                      yin(200) = yin(200) + dykl*yin(199)
                                      zin(200) = zin(200) + dzkl*zin(199)

                                      ! i3 = i4 =  199
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  194

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  194

                                      ! do nk = 1,    3

                                      xin(194) = xin(195) + dxkl*xin(193)
                                      yin(194) = yin(195) + dykl*yin(193)
                                      zin(194) = zin(195) + dzkl*zin(193)
                                      ! i4 = i4 + lang+1 =  196

                                      ! nk =    2

                                      xin(196) = xin(197) + dxkl*xin(195)
                                      yin(196) = yin(197) + dykl*yin(195)
                                      zin(196) = zin(197) + dzkl*zin(195)
                                      ! i4 = i4 + lang+1 =  198

                                      ! nk =    3

                                      xin(198) = xin(199) + dxkl*xin(197)
                                      yin(198) = yin(199) + dykl*yin(197)
                                      zin(198) = zin(199) + dzkl*zin(197)
                                      ! i4 = i4 + lang+1 =  200

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  195

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  201

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  208

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  207

                                      xin(208) = xin(208) + dxkl*xin(207)
                                      yin(208) = yin(208) + dykl*yin(207)
                                      zin(208) = zin(208) + dzkl*zin(207)

                                      ! i3 = i4 =  207
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  202

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  202

                                      ! do nk = 1,    3

                                      xin(202) = xin(203) + dxkl*xin(201)
                                      yin(202) = yin(203) + dykl*yin(201)
                                      zin(202) = zin(203) + dzkl*zin(201)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    2

                                      xin(204) = xin(205) + dxkl*xin(203)
                                      yin(204) = yin(205) + dykl*yin(203)
                                      zin(204) = zin(205) + dzkl*zin(203)
                                      ! i4 = i4 + lang+1 =  206

                                      ! nk =    3

                                      xin(206) = xin(207) + dxkl*xin(205)
                                      yin(206) = yin(207) + dykl*yin(205)
                                      zin(206) = zin(207) + dzkl*zin(205)
                                      ! i4 = i4 + lang+1 =  208

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  203

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  209

                                      ! nj = nj + 1 =    2

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  210

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  210

                                      ! do nk = 1,    3

                                      xin(210) = xin(211) + dxkl*xin(209)
                                      yin(210) = yin(211) + dykl*yin(209)
                                      zin(210) = zin(211) + dzkl*zin(209)
                                      ! i4 = i4 + lang+1 =  212

                                      ! nk =    2

                                      xin(212) = xin(213) + dxkl*xin(211)
                                      yin(212) = yin(213) + dykl*yin(211)
                                      zin(212) = zin(213) + dzkl*zin(211)
                                      ! i4 = i4 + lang+1 =  214

                                      ! nk =    3

                                      xin(214) = xin(215) + dxkl*xin(213)
                                      yin(214) = yin(215) + dykl*yin(213)
                                      zin(214) = zin(215) + dzkl*zin(213)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  211

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  224

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  223

                                      xin(224) = xin(224) + dxkl*xin(223)
                                      yin(224) = yin(224) + dykl*yin(223)
                                      zin(224) = zin(224) + dzkl*zin(223)

                                      ! i3 = i4 =  223
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    3

                                      xin(218) = xin(219) + dxkl*xin(217)
                                      yin(218) = yin(219) + dykl*yin(217)
                                      zin(218) = zin(219) + dzkl*zin(217)
                                      ! i4 = i4 + lang+1 =  220

                                      ! nk =    2

                                      xin(220) = xin(221) + dxkl*xin(219)
                                      yin(220) = yin(221) + dykl*yin(219)
                                      zin(220) = zin(221) + dzkl*zin(219)
                                      ! i4 = i4 + lang+1 =  222

                                      ! nk =    3

                                      xin(222) = xin(223) + dxkl*xin(221)
                                      yin(222) = yin(223) + dykl*yin(221)
                                      zin(222) = zin(223) + dzkl*zin(221)
                                      ! i4 = i4 + lang+1 =  224

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  225

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  225

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  232

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  231

                                      xin(232) = xin(232) + dxkl*xin(231)
                                      yin(232) = yin(232) + dykl*yin(231)
                                      zin(232) = zin(232) + dzkl*zin(231)

                                      ! i3 = i4 =  231
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  226

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  226

                                      ! do nk = 1,    3

                                      xin(226) = xin(227) + dxkl*xin(225)
                                      yin(226) = yin(227) + dykl*yin(225)
                                      zin(226) = zin(227) + dzkl*zin(225)
                                      ! i4 = i4 + lang+1 =  228

                                      ! nk =    2

                                      xin(228) = xin(229) + dxkl*xin(227)
                                      yin(228) = yin(229) + dykl*yin(227)
                                      zin(228) = zin(229) + dzkl*zin(227)
                                      ! i4 = i4 + lang+1 =  230

                                      ! nk =    3

                                      xin(230) = xin(231) + dxkl*xin(229)
                                      yin(230) = yin(231) + dykl*yin(229)
                                      zin(230) = zin(231) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  232

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  227

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  233

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  240

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  239

                                      xin(240) = xin(240) + dxkl*xin(239)
                                      yin(240) = yin(240) + dykl*yin(239)
                                      zin(240) = zin(240) + dzkl*zin(239)

                                      ! i3 = i4 =  239
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  234

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  234

                                      ! do nk = 1,    3

                                      xin(234) = xin(235) + dxkl*xin(233)
                                      yin(234) = yin(235) + dykl*yin(233)
                                      zin(234) = zin(235) + dzkl*zin(233)
                                      ! i4 = i4 + lang+1 =  236

                                      ! nk =    2

                                      xin(236) = xin(237) + dxkl*xin(235)
                                      yin(236) = yin(237) + dykl*yin(235)
                                      zin(236) = zin(237) + dzkl*zin(235)
                                      ! i4 = i4 + lang+1 =  238

                                      ! nk =    3

                                      xin(238) = xin(239) + dxkl*xin(237)
                                      yin(238) = yin(239) + dykl*yin(237)
                                      zin(238) = zin(239) + dzkl*zin(237)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  235

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  241

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  248

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  247

                                      xin(248) = xin(248) + dxkl*xin(247)
                                      yin(248) = yin(248) + dykl*yin(247)
                                      zin(248) = zin(248) + dzkl*zin(247)

                                      ! i3 = i4 =  247
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  242

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  242

                                      ! do nk = 1,    3

                                      xin(242) = xin(243) + dxkl*xin(241)
                                      yin(242) = yin(243) + dykl*yin(241)
                                      zin(242) = zin(243) + dzkl*zin(241)
                                      ! i4 = i4 + lang+1 =  244

                                      ! nk =    2

                                      xin(244) = xin(245) + dxkl*xin(243)
                                      yin(244) = yin(245) + dykl*yin(243)
                                      zin(244) = zin(245) + dzkl*zin(243)
                                      ! i4 = i4 + lang+1 =  246

                                      ! nk =    3

                                      xin(246) = xin(247) + dxkl*xin(245)
                                      yin(246) = yin(247) + dykl*yin(245)
                                      zin(246) = zin(247) + dzkl*zin(245)
                                      ! i4 = i4 + lang+1 =  248

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  243

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  249

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  256

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  255

                                      xin(256) = xin(256) + dxkl*xin(255)
                                      yin(256) = yin(256) + dykl*yin(255)
                                      zin(256) = zin(256) + dzkl*zin(255)

                                      ! i3 = i4 =  255
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  250

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  250

                                      ! do nk = 1,    3

                                      xin(250) = xin(251) + dxkl*xin(249)
                                      yin(250) = yin(251) + dykl*yin(249)
                                      zin(250) = zin(251) + dzkl*zin(249)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    2

                                      xin(252) = xin(253) + dxkl*xin(251)
                                      yin(252) = yin(253) + dykl*yin(251)
                                      zin(252) = zin(253) + dzkl*zin(251)
                                      ! i4 = i4 + lang+1 =  254

                                      ! nk =    3

                                      xin(254) = xin(255) + dxkl*xin(253)
                                      yin(254) = yin(255) + dykl*yin(253)
                                      zin(254) = zin(255) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  256

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  251

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  257

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  257

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  256

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

                                      ! i1 = in(1) =  257

                                      xin(257) = 1.0_dp
                                      yin(257) = 1.0_dp
                                      zin(257) = f00

                                      ! i2 = in(2) =  289
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(289) = xc00
                                      yin(289) = yc00
                                      zin(289) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  259

                                      xin(259) = xcp00
                                      yin(259) = ycp00
                                      zin(259) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  291
                                      ! i2 =  289

                                      xin(291) = xcp00*xin(289) + cp10
                                      yin(291) = ycp00*yin(289) + cp10
                                      zin(291) = zcp00*zin(289) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  289

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  321
                                      ! i3 =  257
                                      ! i4 =  289

                                      xin(321) = c10*xin(257) + xc00*xin(289)
                                      yin(321) = c10*yin(257) + yc00*yin(289)
                                      zin(321) = c10*zin(257) + zc00*zin(289)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  323
                                      ! i5 =  321
                                      ! i4 =  289

                                      xin(323) = xcp00*xin(321) + cp10*xin(289)
                                      yin(323) = ycp00*yin(321) + cp10*yin(289)
                                      zin(323) = zcp00*zin(321) + cp10*zin(289)

                                      ! ------------------

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  321

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  353
                                      ! i3 =  289
                                      ! i4 =  321

                                      xin(353) = c10*xin(289) + xc00*xin(321)
                                      yin(353) = c10*yin(289) + yc00*yin(321)
                                      zin(353) = c10*zin(289) + zc00*zin(321)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  355
                                      ! i5 =  353
                                      ! i4 =  321

                                      xin(355) = xcp00*xin(353) + cp10*xin(321)
                                      yin(355) = ycp00*yin(353) + cp10*yin(321)
                                      zin(355) = zcp00*zin(353) + cp10*zin(321)

                                      ! ------------------

                                      ! i3 = i4 =  321
                                      ! i4 = i5 =  353

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  361
                                      ! i3 =  321
                                      ! i4 =  353

                                      xin(361) = c10*xin(321) + xc00*xin(353)
                                      yin(361) = c10*yin(321) + yc00*yin(353)
                                      zin(361) = c10*zin(321) + zc00*zin(353)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  363
                                      ! i5 =  361
                                      ! i4 =  353

                                      xin(363) = xcp00*xin(361) + cp10*xin(353)
                                      yin(363) = ycp00*yin(361) + cp10*yin(353)
                                      zin(363) = zcp00*zin(361) + cp10*zin(353)

                                      ! ------------------

                                      ! i3 = i4 =  353
                                      ! i4 = i5 =  361

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  369
                                      ! i3 =  353
                                      ! i4 =  361

                                      xin(369) = c10*xin(353) + xc00*xin(361)
                                      yin(369) = c10*yin(353) + yc00*yin(361)
                                      zin(369) = c10*zin(353) + zc00*zin(361)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  371
                                      ! i5 =  369
                                      ! i4 =  361

                                      xin(371) = xcp00*xin(369) + cp10*xin(361)
                                      yin(371) = ycp00*yin(369) + cp10*yin(361)
                                      zin(371) = zcp00*zin(369) + cp10*zin(361)

                                      ! ------------------

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  369

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  377
                                      ! i3 =  361
                                      ! i4 =  369

                                      xin(377) = c10*xin(361) + xc00*xin(369)
                                      yin(377) = c10*yin(361) + yc00*yin(369)
                                      zin(377) = c10*zin(361) + zc00*zin(369)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  379
                                      ! i5 =  377
                                      ! i4 =  369

                                      xin(379) = xcp00*xin(377) + cp10*xin(369)
                                      yin(379) = ycp00*yin(377) + cp10*yin(369)
                                      zin(379) = zcp00*zin(377) + cp10*zin(369)

                                      ! ------------------

                                      ! i3 = i4 =  369
                                      ! i4 = i5 =  377

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  257
                                      ! i4 = i1+k2 =  259

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  261
                                      ! i3 =  257
                                      ! i4 =  259

                                      xin(261) = cp01*xin(257) + xcp00*xin(259)
                                      yin(261) = cp01*yin(257) + ycp00*yin(259)
                                      zin(261) = cp01*zin(257) + zcp00*zin(259)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  293

                                      xin(293) = xc00*xin(261) + c01*xin(259)
                                      yin(293) = yc00*yin(261) + c01*yin(259)
                                      zin(293) = zc00*zin(261) + c01*zin(259)

                                      ! ------------------

                                      ! i3 = i4 =  259
                                      ! i4 = i5 =  261

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  263
                                      ! i3 =  259
                                      ! i4 =  261

                                      xin(263) = cp01*xin(259) + xcp00*xin(261)
                                      yin(263) = cp01*yin(259) + ycp00*yin(261)
                                      zin(263) = cp01*zin(259) + zcp00*zin(261)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  295

                                      xin(295) = xc00*xin(263) + c01*xin(261)
                                      yin(295) = yc00*yin(263) + c01*yin(261)
                                      zin(295) = zc00*zin(263) + c01*zin(261)

                                      ! ------------------

                                      ! i3 = i4 =  261
                                      ! i4 = i5 =  263

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  264
                                      ! i3 =  261
                                      ! i4 =  263

                                      xin(264) = cp01*xin(261) + xcp00*xin(263)
                                      yin(264) = cp01*yin(261) + ycp00*yin(263)
                                      zin(264) = cp01*zin(261) + zcp00*zin(263)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  296

                                      xin(296) = xc00*xin(264) + c01*xin(263)
                                      yin(296) = yc00*yin(264) + c01*yin(263)
                                      zin(296) = zc00*zin(264) + c01*zin(263)

                                      ! ------------------

                                      ! i3 = i4 =  263
                                      ! i4 = i5 =  264

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  289

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  321

                                      xin(325) = c10*xin(261) + xc00*xin(293) + c01*xin(291)
                                      yin(325) = c10*yin(261) + yc00*yin(293) + c01*yin(291)
                                      zin(325) = c10*zin(261) + zc00*zin(293) + c01*zin(291)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  321

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  353

                                      xin(357) = c10*xin(293) + xc00*xin(325) + c01*xin(323)
                                      yin(357) = c10*yin(293) + yc00*yin(325) + c01*yin(323)
                                      zin(357) = c10*zin(293) + zc00*zin(325) + c01*zin(323)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  321
                                      ! i4 = i5 =  353

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  361

                                      xin(365) = c10*xin(325) + xc00*xin(357) + c01*xin(355)
                                      yin(365) = c10*yin(325) + yc00*yin(357) + c01*yin(355)
                                      zin(365) = c10*zin(325) + zc00*zin(357) + c01*zin(355)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  353
                                      ! i4 = i5 =  361

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  369

                                      xin(373) = c10*xin(357) + xc00*xin(365) + c01*xin(363)
                                      yin(373) = c10*yin(357) + yc00*yin(365) + c01*yin(363)
                                      zin(373) = c10*zin(357) + zc00*zin(365) + c01*zin(363)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  369

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  377

                                      xin(381) = c10*xin(365) + xc00*xin(373) + c01*xin(371)
                                      yin(381) = c10*yin(365) + yc00*yin(373) + c01*yin(371)
                                      zin(381) = c10*zin(365) + zc00*zin(373) + c01*zin(371)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  369
                                      ! i4 = i5 =  377

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  289

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  321

                                      xin(327) = c10*xin(263) + xc00*xin(295) + c01*xin(293)
                                      yin(327) = c10*yin(263) + yc00*yin(295) + c01*yin(293)
                                      zin(327) = c10*zin(263) + zc00*zin(295) + c01*zin(293)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  321

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  353

                                      xin(359) = c10*xin(295) + xc00*xin(327) + c01*xin(325)
                                      yin(359) = c10*yin(295) + yc00*yin(327) + c01*yin(325)
                                      zin(359) = c10*zin(295) + zc00*zin(327) + c01*zin(325)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  321
                                      ! i4 = i5 =  353

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  361

                                      xin(367) = c10*xin(327) + xc00*xin(359) + c01*xin(357)
                                      yin(367) = c10*yin(327) + yc00*yin(359) + c01*yin(357)
                                      zin(367) = c10*zin(327) + zc00*zin(359) + c01*zin(357)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  353
                                      ! i4 = i5 =  361

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  369

                                      xin(375) = c10*xin(359) + xc00*xin(367) + c01*xin(365)
                                      yin(375) = c10*yin(359) + yc00*yin(367) + c01*yin(365)
                                      zin(375) = c10*zin(359) + zc00*zin(367) + c01*zin(365)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  369

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  377

                                      xin(383) = c10*xin(367) + xc00*xin(375) + c01*xin(373)
                                      yin(383) = c10*yin(367) + yc00*yin(375) + c01*yin(373)
                                      zin(383) = c10*zin(367) + zc00*zin(375) + c01*zin(373)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  369
                                      ! i4 = i5 =  377

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  289

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  321

                                      xin(328) = c10*xin(264) + xc00*xin(296) + c01*xin(295)
                                      yin(328) = c10*yin(264) + yc00*yin(296) + c01*yin(295)
                                      zin(328) = c10*zin(264) + zc00*zin(296) + c01*zin(295)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  321

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  353

                                      xin(360) = c10*xin(296) + xc00*xin(328) + c01*xin(327)
                                      yin(360) = c10*yin(296) + yc00*yin(328) + c01*yin(327)
                                      zin(360) = c10*zin(296) + zc00*zin(328) + c01*zin(327)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  321
                                      ! i4 = i5 =  353

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  361

                                      xin(368) = c10*xin(328) + xc00*xin(360) + c01*xin(359)
                                      yin(368) = c10*yin(328) + yc00*yin(360) + c01*yin(359)
                                      zin(368) = c10*zin(328) + zc00*zin(360) + c01*zin(359)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  353
                                      ! i4 = i5 =  361

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  369

                                      xin(376) = c10*xin(360) + xc00*xin(368) + c01*xin(367)
                                      yin(376) = c10*yin(360) + yc00*yin(368) + c01*yin(367)
                                      zin(376) = c10*zin(360) + zc00*zin(368) + c01*zin(367)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  369

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  377

                                      xin(384) = c10*xin(368) + xc00*xin(376) + c01*xin(375)
                                      yin(384) = c10*yin(368) + yc00*yin(376) + c01*yin(375)
                                      zin(384) = c10*zin(368) + zc00*zin(376) + c01*zin(375)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  369
                                      ! i4 = i5 =  377

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  377

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  377

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  369

                                      xin(377) = xin(377) + dxij*xin(369)
                                      yin(377) = yin(377) + dyij*yin(369)
                                      zin(377) = zin(377) + dzij*zin(369)

                                      ! i3 = i4 =  369
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  361

                                      xin(369) = xin(369) + dxij*xin(361)
                                      yin(369) = yin(369) + dyij*yin(361)
                                      zin(369) = zin(369) + dzij*zin(361)

                                      ! i3 = i4 =  361
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  353

                                      xin(361) = xin(361) + dxij*xin(353)
                                      yin(361) = yin(361) + dyij*yin(353)
                                      zin(361) = zin(361) + dzij*zin(353)

                                      ! i3 = i4 =  353
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  377

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  369

                                      xin(377) = xin(377) + dxij*xin(369)
                                      yin(377) = yin(377) + dyij*yin(369)
                                      zin(377) = zin(377) + dzij*zin(369)

                                      ! i3 = i4 =  369
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  361

                                      xin(369) = xin(369) + dxij*xin(361)
                                      yin(369) = yin(369) + dyij*yin(361)
                                      zin(369) = zin(369) + dzij*zin(361)

                                      ! i3 = i4 =  361
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  377

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  369

                                      xin(377) = xin(377) + dxij*xin(369)
                                      yin(377) = yin(377) + dyij*yin(369)
                                      zin(377) = zin(377) + dzij*zin(369)

                                      ! i3 = i4 =  369
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  265

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  265

                                      ! do ni = 1,    3

                                      xin(265) = xin(289) + dxij*xin(257)
                                      yin(265) = yin(289) + dyij*yin(257)
                                      zin(265) = zin(289) + dzij*zin(257)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  297

                                      ! ni =    2

                                      xin(297) = xin(321) + dxij*xin(289)
                                      yin(297) = yin(321) + dyij*yin(289)
                                      zin(297) = zin(321) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  329

                                      ! ni =    3

                                      xin(329) = xin(353) + dxij*xin(321)
                                      yin(329) = yin(353) + dyij*yin(321)
                                      zin(329) = zin(353) + dzij*zin(321)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  273

                                      ! nj =    2

                                      ! i4 = i3 =  273

                                      ! do ni = 1,    3

                                      xin(273) = xin(297) + dxij*xin(265)
                                      yin(273) = yin(297) + dyij*yin(265)
                                      zin(273) = zin(297) + dzij*zin(265)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  305

                                      ! ni =    2

                                      xin(305) = xin(329) + dxij*xin(297)
                                      yin(305) = yin(329) + dyij*yin(297)
                                      zin(305) = zin(329) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  337

                                      ! ni =    3

                                      xin(337) = xin(361) + dxij*xin(329)
                                      yin(337) = yin(361) + dyij*yin(329)
                                      zin(337) = zin(361) + dzij*zin(329)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  369

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  281

                                      ! nj =    3

                                      ! i4 = i3 =  281

                                      ! do ni = 1,    3

                                      xin(281) = xin(305) + dxij*xin(273)
                                      yin(281) = yin(305) + dyij*yin(273)
                                      zin(281) = zin(305) + dzij*zin(273)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! ni =    2

                                      xin(313) = xin(337) + dxij*xin(305)
                                      yin(313) = yin(337) + dyij*yin(305)
                                      zin(313) = zin(337) + dzij*zin(305)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  345

                                      ! ni =    3

                                      xin(345) = xin(369) + dxij*xin(337)
                                      yin(345) = yin(369) + dyij*yin(337)
                                      zin(345) = zin(369) + dzij*zin(337)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  377

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  289

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  371

                                      xin(379) = xin(379) + dxij*xin(371)
                                      yin(379) = yin(379) + dyij*yin(371)
                                      zin(379) = zin(379) + dzij*zin(371)

                                      ! i3 = i4 =  371
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  363

                                      xin(371) = xin(371) + dxij*xin(363)
                                      yin(371) = yin(371) + dyij*yin(363)
                                      zin(371) = zin(371) + dzij*zin(363)

                                      ! i3 = i4 =  363
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  355

                                      xin(363) = xin(363) + dxij*xin(355)
                                      yin(363) = yin(363) + dyij*yin(355)
                                      zin(363) = zin(363) + dzij*zin(355)

                                      ! i3 = i4 =  355
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  371

                                      xin(379) = xin(379) + dxij*xin(371)
                                      yin(379) = yin(379) + dyij*yin(371)
                                      zin(379) = zin(379) + dzij*zin(371)

                                      ! i3 = i4 =  371
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  363

                                      xin(371) = xin(371) + dxij*xin(363)
                                      yin(371) = yin(371) + dyij*yin(363)
                                      zin(371) = zin(371) + dzij*zin(363)

                                      ! i3 = i4 =  363
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  371

                                      xin(379) = xin(379) + dxij*xin(371)
                                      yin(379) = yin(379) + dyij*yin(371)
                                      zin(379) = zin(379) + dzij*zin(371)

                                      ! i3 = i4 =  371
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  267

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  267

                                      ! do ni = 1,    3

                                      xin(267) = xin(291) + dxij*xin(259)
                                      yin(267) = yin(291) + dyij*yin(259)
                                      zin(267) = zin(291) + dzij*zin(259)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  299

                                      ! ni =    2

                                      xin(299) = xin(323) + dxij*xin(291)
                                      yin(299) = yin(323) + dyij*yin(291)
                                      zin(299) = zin(323) + dzij*zin(291)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  331

                                      ! ni =    3

                                      xin(331) = xin(355) + dxij*xin(323)
                                      yin(331) = yin(355) + dyij*yin(323)
                                      zin(331) = zin(355) + dzij*zin(323)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  363

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  275

                                      ! nj =    2

                                      ! i4 = i3 =  275

                                      ! do ni = 1,    3

                                      xin(275) = xin(299) + dxij*xin(267)
                                      yin(275) = yin(299) + dyij*yin(267)
                                      zin(275) = zin(299) + dzij*zin(267)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  307

                                      ! ni =    2

                                      xin(307) = xin(331) + dxij*xin(299)
                                      yin(307) = yin(331) + dyij*yin(299)
                                      zin(307) = zin(331) + dzij*zin(299)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  339

                                      ! ni =    3

                                      xin(339) = xin(363) + dxij*xin(331)
                                      yin(339) = yin(363) + dyij*yin(331)
                                      zin(339) = zin(363) + dzij*zin(331)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  371

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  283

                                      ! nj =    3

                                      ! i4 = i3 =  283

                                      ! do ni = 1,    3

                                      xin(283) = xin(307) + dxij*xin(275)
                                      yin(283) = yin(307) + dyij*yin(275)
                                      zin(283) = zin(307) + dzij*zin(275)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  315

                                      ! ni =    2

                                      xin(315) = xin(339) + dxij*xin(307)
                                      yin(315) = yin(339) + dyij*yin(307)
                                      zin(315) = zin(339) + dzij*zin(307)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  347

                                      ! ni =    3

                                      xin(347) = xin(371) + dxij*xin(339)
                                      yin(347) = yin(371) + dyij*yin(339)
                                      zin(347) = zin(371) + dzij*zin(339)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  379

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  291

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(381) = xin(381) + dxij*xin(373)
                                      yin(381) = yin(381) + dyij*yin(373)
                                      zin(381) = zin(381) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  365

                                      xin(373) = xin(373) + dxij*xin(365)
                                      yin(373) = yin(373) + dyij*yin(365)
                                      zin(373) = zin(373) + dzij*zin(365)

                                      ! i3 = i4 =  365
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  357

                                      xin(365) = xin(365) + dxij*xin(357)
                                      yin(365) = yin(365) + dyij*yin(357)
                                      zin(365) = zin(365) + dzij*zin(357)

                                      ! i3 = i4 =  357
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(381) = xin(381) + dxij*xin(373)
                                      yin(381) = yin(381) + dyij*yin(373)
                                      zin(381) = zin(381) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  365

                                      xin(373) = xin(373) + dxij*xin(365)
                                      yin(373) = yin(373) + dyij*yin(365)
                                      zin(373) = zin(373) + dzij*zin(365)

                                      ! i3 = i4 =  365
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(381) = xin(381) + dxij*xin(373)
                                      yin(381) = yin(381) + dyij*yin(373)
                                      zin(381) = zin(381) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  269

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  269

                                      ! do ni = 1,    3

                                      xin(269) = xin(293) + dxij*xin(261)
                                      yin(269) = yin(293) + dyij*yin(261)
                                      zin(269) = zin(293) + dzij*zin(261)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  301

                                      ! ni =    2

                                      xin(301) = xin(325) + dxij*xin(293)
                                      yin(301) = yin(325) + dyij*yin(293)
                                      zin(301) = zin(325) + dzij*zin(293)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  333

                                      ! ni =    3

                                      xin(333) = xin(357) + dxij*xin(325)
                                      yin(333) = yin(357) + dyij*yin(325)
                                      zin(333) = zin(357) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  365

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  277

                                      ! nj =    2

                                      ! i4 = i3 =  277

                                      ! do ni = 1,    3

                                      xin(277) = xin(301) + dxij*xin(269)
                                      yin(277) = yin(301) + dyij*yin(269)
                                      zin(277) = zin(301) + dzij*zin(269)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  309

                                      ! ni =    2

                                      xin(309) = xin(333) + dxij*xin(301)
                                      yin(309) = yin(333) + dyij*yin(301)
                                      zin(309) = zin(333) + dzij*zin(301)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  341

                                      ! ni =    3

                                      xin(341) = xin(365) + dxij*xin(333)
                                      yin(341) = yin(365) + dyij*yin(333)
                                      zin(341) = zin(365) + dzij*zin(333)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  373

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  285

                                      ! nj =    3

                                      ! i4 = i3 =  285

                                      ! do ni = 1,    3

                                      xin(285) = xin(309) + dxij*xin(277)
                                      yin(285) = yin(309) + dyij*yin(277)
                                      zin(285) = zin(309) + dzij*zin(277)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  317

                                      ! ni =    2

                                      xin(317) = xin(341) + dxij*xin(309)
                                      yin(317) = yin(341) + dyij*yin(309)
                                      zin(317) = zin(341) + dzij*zin(309)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                      ! ni =    3

                                      xin(349) = xin(373) + dxij*xin(341)
                                      yin(349) = yin(373) + dyij*yin(341)
                                      zin(349) = zin(373) + dzij*zin(341)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  381

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  293

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(383) = xin(383) + dxij*xin(375)
                                      yin(383) = yin(383) + dyij*yin(375)
                                      zin(383) = zin(383) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  367

                                      xin(375) = xin(375) + dxij*xin(367)
                                      yin(375) = yin(375) + dyij*yin(367)
                                      zin(375) = zin(375) + dzij*zin(367)

                                      ! i3 = i4 =  367
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  359

                                      xin(367) = xin(367) + dxij*xin(359)
                                      yin(367) = yin(367) + dyij*yin(359)
                                      zin(367) = zin(367) + dzij*zin(359)

                                      ! i3 = i4 =  359
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(383) = xin(383) + dxij*xin(375)
                                      yin(383) = yin(383) + dyij*yin(375)
                                      zin(383) = zin(383) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  367

                                      xin(375) = xin(375) + dxij*xin(367)
                                      yin(375) = yin(375) + dyij*yin(367)
                                      zin(375) = zin(375) + dzij*zin(367)

                                      ! i3 = i4 =  367
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(383) = xin(383) + dxij*xin(375)
                                      yin(383) = yin(383) + dyij*yin(375)
                                      zin(383) = zin(383) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  271

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  271

                                      ! do ni = 1,    3

                                      xin(271) = xin(295) + dxij*xin(263)
                                      yin(271) = yin(295) + dyij*yin(263)
                                      zin(271) = zin(295) + dzij*zin(263)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  303

                                      ! ni =    2

                                      xin(303) = xin(327) + dxij*xin(295)
                                      yin(303) = yin(327) + dyij*yin(295)
                                      zin(303) = zin(327) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  335

                                      ! ni =    3

                                      xin(335) = xin(359) + dxij*xin(327)
                                      yin(335) = yin(359) + dyij*yin(327)
                                      zin(335) = zin(359) + dzij*zin(327)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  367

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  279

                                      ! nj =    2

                                      ! i4 = i3 =  279

                                      ! do ni = 1,    3

                                      xin(279) = xin(303) + dxij*xin(271)
                                      yin(279) = yin(303) + dyij*yin(271)
                                      zin(279) = zin(303) + dzij*zin(271)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  311

                                      ! ni =    2

                                      xin(311) = xin(335) + dxij*xin(303)
                                      yin(311) = yin(335) + dyij*yin(303)
                                      zin(311) = zin(335) + dzij*zin(303)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  343

                                      ! ni =    3

                                      xin(343) = xin(367) + dxij*xin(335)
                                      yin(343) = yin(367) + dyij*yin(335)
                                      zin(343) = zin(367) + dzij*zin(335)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  375

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  287

                                      ! nj =    3

                                      ! i4 = i3 =  287

                                      ! do ni = 1,    3

                                      xin(287) = xin(311) + dxij*xin(279)
                                      yin(287) = yin(311) + dyij*yin(279)
                                      zin(287) = zin(311) + dzij*zin(279)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  319

                                      ! ni =    2

                                      xin(319) = xin(343) + dxij*xin(311)
                                      yin(319) = yin(343) + dyij*yin(311)
                                      zin(319) = zin(343) + dzij*zin(311)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  351

                                      ! ni =    3

                                      xin(351) = xin(375) + dxij*xin(343)
                                      yin(351) = yin(375) + dyij*yin(343)
                                      zin(351) = zin(375) + dzij*zin(343)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  383

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  295

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  376

                                      xin(384) = xin(384) + dxij*xin(376)
                                      yin(384) = yin(384) + dyij*yin(376)
                                      zin(384) = zin(384) + dzij*zin(376)

                                      ! i3 = i4 =  376
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  368

                                      xin(376) = xin(376) + dxij*xin(368)
                                      yin(376) = yin(376) + dyij*yin(368)
                                      zin(376) = zin(376) + dzij*zin(368)

                                      ! i3 = i4 =  368
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  360

                                      xin(368) = xin(368) + dxij*xin(360)
                                      yin(368) = yin(368) + dyij*yin(360)
                                      zin(368) = zin(368) + dzij*zin(360)

                                      ! i3 = i4 =  360
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  376

                                      xin(384) = xin(384) + dxij*xin(376)
                                      yin(384) = yin(384) + dyij*yin(376)
                                      zin(384) = zin(384) + dzij*zin(376)

                                      ! i3 = i4 =  376
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  368

                                      xin(376) = xin(376) + dxij*xin(368)
                                      yin(376) = yin(376) + dyij*yin(368)
                                      zin(376) = zin(376) + dzij*zin(368)

                                      ! i3 = i4 =  368
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  376

                                      xin(384) = xin(384) + dxij*xin(376)
                                      yin(384) = yin(384) + dyij*yin(376)
                                      zin(384) = zin(384) + dzij*zin(376)

                                      ! i3 = i4 =  376
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  272

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  272

                                      ! do ni = 1,    3

                                      xin(272) = xin(296) + dxij*xin(264)
                                      yin(272) = yin(296) + dyij*yin(264)
                                      zin(272) = zin(296) + dzij*zin(264)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  304

                                      ! ni =    2

                                      xin(304) = xin(328) + dxij*xin(296)
                                      yin(304) = yin(328) + dyij*yin(296)
                                      zin(304) = zin(328) + dzij*zin(296)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  336

                                      ! ni =    3

                                      xin(336) = xin(360) + dxij*xin(328)
                                      yin(336) = yin(360) + dyij*yin(328)
                                      zin(336) = zin(360) + dzij*zin(328)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  368

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  280

                                      ! nj =    2

                                      ! i4 = i3 =  280

                                      ! do ni = 1,    3

                                      xin(280) = xin(304) + dxij*xin(272)
                                      yin(280) = yin(304) + dyij*yin(272)
                                      zin(280) = zin(304) + dzij*zin(272)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  312

                                      ! ni =    2

                                      xin(312) = xin(336) + dxij*xin(304)
                                      yin(312) = yin(336) + dyij*yin(304)
                                      zin(312) = zin(336) + dzij*zin(304)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  344

                                      ! ni =    3

                                      xin(344) = xin(368) + dxij*xin(336)
                                      yin(344) = yin(368) + dyij*yin(336)
                                      zin(344) = zin(368) + dzij*zin(336)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  376

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  288

                                      ! nj =    3

                                      ! i4 = i3 =  288

                                      ! do ni = 1,    3

                                      xin(288) = xin(312) + dxij*xin(280)
                                      yin(288) = yin(312) + dyij*yin(280)
                                      zin(288) = zin(312) + dzij*zin(280)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  320

                                      ! ni =    2

                                      xin(320) = xin(344) + dxij*xin(312)
                                      yin(320) = yin(344) + dyij*yin(312)
                                      zin(320) = zin(344) + dzij*zin(312)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  352

                                      ! ni =    3

                                      xin(352) = xin(376) + dxij*xin(344)
                                      yin(352) = yin(376) + dyij*yin(344)
                                      zin(352) = zin(376) + dzij*zin(344)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  384

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  296

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  257

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  264

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  263

                                      xin(264) = xin(264) + dxkl*xin(263)
                                      yin(264) = yin(264) + dykl*yin(263)
                                      zin(264) = zin(264) + dzkl*zin(263)

                                      ! i3 = i4 =  263
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  258

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  258

                                      ! do nk = 1,    3

                                      xin(258) = xin(259) + dxkl*xin(257)
                                      yin(258) = yin(259) + dykl*yin(257)
                                      zin(258) = zin(259) + dzkl*zin(257)
                                      ! i4 = i4 + lang+1 =  260

                                      ! nk =    2

                                      xin(260) = xin(261) + dxkl*xin(259)
                                      yin(260) = yin(261) + dykl*yin(259)
                                      zin(260) = zin(261) + dzkl*zin(259)
                                      ! i4 = i4 + lang+1 =  262

                                      ! nk =    3

                                      xin(262) = xin(263) + dxkl*xin(261)
                                      yin(262) = yin(263) + dykl*yin(261)
                                      zin(262) = zin(263) + dzkl*zin(261)
                                      ! i4 = i4 + lang+1 =  264

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  259

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  265

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  272

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  271

                                      xin(272) = xin(272) + dxkl*xin(271)
                                      yin(272) = yin(272) + dykl*yin(271)
                                      zin(272) = zin(272) + dzkl*zin(271)

                                      ! i3 = i4 =  271
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  266

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  266

                                      ! do nk = 1,    3

                                      xin(266) = xin(267) + dxkl*xin(265)
                                      yin(266) = yin(267) + dykl*yin(265)
                                      zin(266) = zin(267) + dzkl*zin(265)
                                      ! i4 = i4 + lang+1 =  268

                                      ! nk =    2

                                      xin(268) = xin(269) + dxkl*xin(267)
                                      yin(268) = yin(269) + dykl*yin(267)
                                      zin(268) = zin(269) + dzkl*zin(267)
                                      ! i4 = i4 + lang+1 =  270

                                      ! nk =    3

                                      xin(270) = xin(271) + dxkl*xin(269)
                                      yin(270) = yin(271) + dykl*yin(269)
                                      zin(270) = zin(271) + dzkl*zin(269)
                                      ! i4 = i4 + lang+1 =  272

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  267

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  273

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  280

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  279

                                      xin(280) = xin(280) + dxkl*xin(279)
                                      yin(280) = yin(280) + dykl*yin(279)
                                      zin(280) = zin(280) + dzkl*zin(279)

                                      ! i3 = i4 =  279
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  274

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  274

                                      ! do nk = 1,    3

                                      xin(274) = xin(275) + dxkl*xin(273)
                                      yin(274) = yin(275) + dykl*yin(273)
                                      zin(274) = zin(275) + dzkl*zin(273)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    2

                                      xin(276) = xin(277) + dxkl*xin(275)
                                      yin(276) = yin(277) + dykl*yin(275)
                                      zin(276) = zin(277) + dzkl*zin(275)
                                      ! i4 = i4 + lang+1 =  278

                                      ! nk =    3

                                      xin(278) = xin(279) + dxkl*xin(277)
                                      yin(278) = yin(279) + dykl*yin(277)
                                      zin(278) = zin(279) + dzkl*zin(277)
                                      ! i4 = i4 + lang+1 =  280

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  275

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  281

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  282

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  282

                                      ! do nk = 1,    3

                                      xin(282) = xin(283) + dxkl*xin(281)
                                      yin(282) = yin(283) + dykl*yin(281)
                                      zin(282) = zin(283) + dzkl*zin(281)
                                      ! i4 = i4 + lang+1 =  284

                                      ! nk =    2

                                      xin(284) = xin(285) + dxkl*xin(283)
                                      yin(284) = yin(285) + dykl*yin(283)
                                      zin(284) = zin(285) + dzkl*zin(283)
                                      ! i4 = i4 + lang+1 =  286

                                      ! nk =    3

                                      xin(286) = xin(287) + dxkl*xin(285)
                                      yin(286) = yin(287) + dykl*yin(285)
                                      zin(286) = zin(287) + dzkl*zin(285)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  283

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  296

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  295

                                      xin(296) = xin(296) + dxkl*xin(295)
                                      yin(296) = yin(296) + dykl*yin(295)
                                      zin(296) = zin(296) + dzkl*zin(295)

                                      ! i3 = i4 =  295
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    3

                                      xin(290) = xin(291) + dxkl*xin(289)
                                      yin(290) = yin(291) + dykl*yin(289)
                                      zin(290) = zin(291) + dzkl*zin(289)
                                      ! i4 = i4 + lang+1 =  292

                                      ! nk =    2

                                      xin(292) = xin(293) + dxkl*xin(291)
                                      yin(292) = yin(293) + dykl*yin(291)
                                      zin(292) = zin(293) + dzkl*zin(291)
                                      ! i4 = i4 + lang+1 =  294

                                      ! nk =    3

                                      xin(294) = xin(295) + dxkl*xin(293)
                                      yin(294) = yin(295) + dykl*yin(293)
                                      zin(294) = zin(295) + dzkl*zin(293)
                                      ! i4 = i4 + lang+1 =  296

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  297

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  304

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  303

                                      xin(304) = xin(304) + dxkl*xin(303)
                                      yin(304) = yin(304) + dykl*yin(303)
                                      zin(304) = zin(304) + dzkl*zin(303)

                                      ! i3 = i4 =  303
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  298

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  298

                                      ! do nk = 1,    3

                                      xin(298) = xin(299) + dxkl*xin(297)
                                      yin(298) = yin(299) + dykl*yin(297)
                                      zin(298) = zin(299) + dzkl*zin(297)
                                      ! i4 = i4 + lang+1 =  300

                                      ! nk =    2

                                      xin(300) = xin(301) + dxkl*xin(299)
                                      yin(300) = yin(301) + dykl*yin(299)
                                      zin(300) = zin(301) + dzkl*zin(299)
                                      ! i4 = i4 + lang+1 =  302

                                      ! nk =    3

                                      xin(302) = xin(303) + dxkl*xin(301)
                                      yin(302) = yin(303) + dykl*yin(301)
                                      zin(302) = zin(303) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  304

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  299

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  305

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  312

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  311

                                      xin(312) = xin(312) + dxkl*xin(311)
                                      yin(312) = yin(312) + dykl*yin(311)
                                      zin(312) = zin(312) + dzkl*zin(311)

                                      ! i3 = i4 =  311
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  306

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  306

                                      ! do nk = 1,    3

                                      xin(306) = xin(307) + dxkl*xin(305)
                                      yin(306) = yin(307) + dykl*yin(305)
                                      zin(306) = zin(307) + dzkl*zin(305)
                                      ! i4 = i4 + lang+1 =  308

                                      ! nk =    2

                                      xin(308) = xin(309) + dxkl*xin(307)
                                      yin(308) = yin(309) + dykl*yin(307)
                                      zin(308) = zin(309) + dzkl*zin(307)
                                      ! i4 = i4 + lang+1 =  310

                                      ! nk =    3

                                      xin(310) = xin(311) + dxkl*xin(309)
                                      yin(310) = yin(311) + dykl*yin(309)
                                      zin(310) = zin(311) + dzkl*zin(309)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  307

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  313

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  320

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  319

                                      xin(320) = xin(320) + dxkl*xin(319)
                                      yin(320) = yin(320) + dykl*yin(319)
                                      zin(320) = zin(320) + dzkl*zin(319)

                                      ! i3 = i4 =  319
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  314

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  314

                                      ! do nk = 1,    3

                                      xin(314) = xin(315) + dxkl*xin(313)
                                      yin(314) = yin(315) + dykl*yin(313)
                                      zin(314) = zin(315) + dzkl*zin(313)
                                      ! i4 = i4 + lang+1 =  316

                                      ! nk =    2

                                      xin(316) = xin(317) + dxkl*xin(315)
                                      yin(316) = yin(317) + dykl*yin(315)
                                      zin(316) = zin(317) + dzkl*zin(315)
                                      ! i4 = i4 + lang+1 =  318

                                      ! nk =    3

                                      xin(318) = xin(319) + dxkl*xin(317)
                                      yin(318) = yin(319) + dykl*yin(317)
                                      zin(318) = zin(319) + dzkl*zin(317)
                                      ! i4 = i4 + lang+1 =  320

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  315

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  321

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  321

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  328

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  327

                                      xin(328) = xin(328) + dxkl*xin(327)
                                      yin(328) = yin(328) + dykl*yin(327)
                                      zin(328) = zin(328) + dzkl*zin(327)

                                      ! i3 = i4 =  327
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  322

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  322

                                      ! do nk = 1,    3

                                      xin(322) = xin(323) + dxkl*xin(321)
                                      yin(322) = yin(323) + dykl*yin(321)
                                      zin(322) = zin(323) + dzkl*zin(321)
                                      ! i4 = i4 + lang+1 =  324

                                      ! nk =    2

                                      xin(324) = xin(325) + dxkl*xin(323)
                                      yin(324) = yin(325) + dykl*yin(323)
                                      zin(324) = zin(325) + dzkl*zin(323)
                                      ! i4 = i4 + lang+1 =  326

                                      ! nk =    3

                                      xin(326) = xin(327) + dxkl*xin(325)
                                      yin(326) = yin(327) + dykl*yin(325)
                                      zin(326) = zin(327) + dzkl*zin(325)
                                      ! i4 = i4 + lang+1 =  328

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  323

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  329

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  336

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  335

                                      xin(336) = xin(336) + dxkl*xin(335)
                                      yin(336) = yin(336) + dykl*yin(335)
                                      zin(336) = zin(336) + dzkl*zin(335)

                                      ! i3 = i4 =  335
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  330

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  330

                                      ! do nk = 1,    3

                                      xin(330) = xin(331) + dxkl*xin(329)
                                      yin(330) = yin(331) + dykl*yin(329)
                                      zin(330) = zin(331) + dzkl*zin(329)
                                      ! i4 = i4 + lang+1 =  332

                                      ! nk =    2

                                      xin(332) = xin(333) + dxkl*xin(331)
                                      yin(332) = yin(333) + dykl*yin(331)
                                      zin(332) = zin(333) + dzkl*zin(331)
                                      ! i4 = i4 + lang+1 =  334

                                      ! nk =    3

                                      xin(334) = xin(335) + dxkl*xin(333)
                                      yin(334) = yin(335) + dykl*yin(333)
                                      zin(334) = zin(335) + dzkl*zin(333)
                                      ! i4 = i4 + lang+1 =  336

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  331

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  337

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  344

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  343

                                      xin(344) = xin(344) + dxkl*xin(343)
                                      yin(344) = yin(344) + dykl*yin(343)
                                      zin(344) = zin(344) + dzkl*zin(343)

                                      ! i3 = i4 =  343
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  338

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  338

                                      ! do nk = 1,    3

                                      xin(338) = xin(339) + dxkl*xin(337)
                                      yin(338) = yin(339) + dykl*yin(337)
                                      zin(338) = zin(339) + dzkl*zin(337)
                                      ! i4 = i4 + lang+1 =  340

                                      ! nk =    2

                                      xin(340) = xin(341) + dxkl*xin(339)
                                      yin(340) = yin(341) + dykl*yin(339)
                                      zin(340) = zin(341) + dzkl*zin(339)
                                      ! i4 = i4 + lang+1 =  342

                                      ! nk =    3

                                      xin(342) = xin(343) + dxkl*xin(341)
                                      yin(342) = yin(343) + dykl*yin(341)
                                      zin(342) = zin(343) + dzkl*zin(341)
                                      ! i4 = i4 + lang+1 =  344

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  339

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  345

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  352

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  351

                                      xin(352) = xin(352) + dxkl*xin(351)
                                      yin(352) = yin(352) + dykl*yin(351)
                                      zin(352) = zin(352) + dzkl*zin(351)

                                      ! i3 = i4 =  351
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  346

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  346

                                      ! do nk = 1,    3

                                      xin(346) = xin(347) + dxkl*xin(345)
                                      yin(346) = yin(347) + dykl*yin(345)
                                      zin(346) = zin(347) + dzkl*zin(345)
                                      ! i4 = i4 + lang+1 =  348

                                      ! nk =    2

                                      xin(348) = xin(349) + dxkl*xin(347)
                                      yin(348) = yin(349) + dykl*yin(347)
                                      zin(348) = zin(349) + dzkl*zin(347)
                                      ! i4 = i4 + lang+1 =  350

                                      ! nk =    3

                                      xin(350) = xin(351) + dxkl*xin(349)
                                      yin(350) = yin(351) + dykl*yin(349)
                                      zin(350) = zin(351) + dzkl*zin(349)
                                      ! i4 = i4 + lang+1 =  352

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  347

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  353

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  353

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  354

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  354

                                      ! do nk = 1,    3

                                      xin(354) = xin(355) + dxkl*xin(353)
                                      yin(354) = yin(355) + dykl*yin(353)
                                      zin(354) = zin(355) + dzkl*zin(353)
                                      ! i4 = i4 + lang+1 =  356

                                      ! nk =    2

                                      xin(356) = xin(357) + dxkl*xin(355)
                                      yin(356) = yin(357) + dykl*yin(355)
                                      zin(356) = zin(357) + dzkl*zin(355)
                                      ! i4 = i4 + lang+1 =  358

                                      ! nk =    3

                                      xin(358) = xin(359) + dxkl*xin(357)
                                      yin(358) = yin(359) + dykl*yin(357)
                                      zin(358) = zin(359) + dzkl*zin(357)
                                      ! i4 = i4 + lang+1 =  360

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  355

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  361

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  368

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  367

                                      xin(368) = xin(368) + dxkl*xin(367)
                                      yin(368) = yin(368) + dykl*yin(367)
                                      zin(368) = zin(368) + dzkl*zin(367)

                                      ! i3 = i4 =  367
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  362

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  362

                                      ! do nk = 1,    3

                                      xin(362) = xin(363) + dxkl*xin(361)
                                      yin(362) = yin(363) + dykl*yin(361)
                                      zin(362) = zin(363) + dzkl*zin(361)
                                      ! i4 = i4 + lang+1 =  364

                                      ! nk =    2

                                      xin(364) = xin(365) + dxkl*xin(363)
                                      yin(364) = yin(365) + dykl*yin(363)
                                      zin(364) = zin(365) + dzkl*zin(363)
                                      ! i4 = i4 + lang+1 =  366

                                      ! nk =    3

                                      xin(366) = xin(367) + dxkl*xin(365)
                                      yin(366) = yin(367) + dykl*yin(365)
                                      zin(366) = zin(367) + dzkl*zin(365)
                                      ! i4 = i4 + lang+1 =  368

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  363

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  369

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  376

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  375

                                      xin(376) = xin(376) + dxkl*xin(375)
                                      yin(376) = yin(376) + dykl*yin(375)
                                      zin(376) = zin(376) + dzkl*zin(375)

                                      ! i3 = i4 =  375
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  370

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  370

                                      ! do nk = 1,    3

                                      xin(370) = xin(371) + dxkl*xin(369)
                                      yin(370) = yin(371) + dykl*yin(369)
                                      zin(370) = zin(371) + dzkl*zin(369)
                                      ! i4 = i4 + lang+1 =  372

                                      ! nk =    2

                                      xin(372) = xin(373) + dxkl*xin(371)
                                      yin(372) = yin(373) + dykl*yin(371)
                                      zin(372) = zin(373) + dzkl*zin(371)
                                      ! i4 = i4 + lang+1 =  374

                                      ! nk =    3

                                      xin(374) = xin(375) + dxkl*xin(373)
                                      yin(374) = yin(375) + dykl*yin(373)
                                      zin(374) = zin(375) + dzkl*zin(373)
                                      ! i4 = i4 + lang+1 =  376

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  371

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  377

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  384

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  383

                                      xin(384) = xin(384) + dxkl*xin(383)
                                      yin(384) = yin(384) + dykl*yin(383)
                                      zin(384) = zin(384) + dzkl*zin(383)

                                      ! i3 = i4 =  383
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  378

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  378

                                      ! do nk = 1,    3

                                      xin(378) = xin(379) + dxkl*xin(377)
                                      yin(378) = yin(379) + dykl*yin(377)
                                      zin(378) = zin(379) + dzkl*zin(377)
                                      ! i4 = i4 + lang+1 =  380

                                      ! nk =    2

                                      xin(380) = xin(381) + dxkl*xin(379)
                                      yin(380) = yin(381) + dykl*yin(379)
                                      zin(380) = zin(381) + dzkl*zin(379)
                                      ! i4 = i4 + lang+1 =  382

                                      ! nk =    3

                                      xin(382) = xin(383) + dxkl*xin(381)
                                      yin(382) = yin(383) + dykl*yin(381)
                                      zin(382) = zin(383) + dzkl*zin(381)
                                      ! i4 = i4 + lang+1 =  384

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  379

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  385

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  385

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  384

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

                                      ! i1 = in(1) =  385

                                      xin(385) = 1.0_dp
                                      yin(385) = 1.0_dp
                                      zin(385) = f00

                                      ! i2 = in(2) =  417
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(417) = xc00
                                      yin(417) = yc00
                                      zin(417) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  387

                                      xin(387) = xcp00
                                      yin(387) = ycp00
                                      zin(387) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  419
                                      ! i2 =  417

                                      xin(419) = xcp00*xin(417) + cp10
                                      yin(419) = ycp00*yin(417) + cp10
                                      zin(419) = zcp00*zin(417) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  417

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  449
                                      ! i3 =  385
                                      ! i4 =  417

                                      xin(449) = c10*xin(385) + xc00*xin(417)
                                      yin(449) = c10*yin(385) + yc00*yin(417)
                                      zin(449) = c10*zin(385) + zc00*zin(417)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  451
                                      ! i5 =  449
                                      ! i4 =  417

                                      xin(451) = xcp00*xin(449) + cp10*xin(417)
                                      yin(451) = ycp00*yin(449) + cp10*yin(417)
                                      zin(451) = zcp00*zin(449) + cp10*zin(417)

                                      ! ------------------

                                      ! i3 = i4 =  417
                                      ! i4 = i5 =  449

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  481
                                      ! i3 =  417
                                      ! i4 =  449

                                      xin(481) = c10*xin(417) + xc00*xin(449)
                                      yin(481) = c10*yin(417) + yc00*yin(449)
                                      zin(481) = c10*zin(417) + zc00*zin(449)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  483
                                      ! i5 =  481
                                      ! i4 =  449

                                      xin(483) = xcp00*xin(481) + cp10*xin(449)
                                      yin(483) = ycp00*yin(481) + cp10*yin(449)
                                      zin(483) = zcp00*zin(481) + cp10*zin(449)

                                      ! ------------------

                                      ! i3 = i4 =  449
                                      ! i4 = i5 =  481

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  489
                                      ! i3 =  449
                                      ! i4 =  481

                                      xin(489) = c10*xin(449) + xc00*xin(481)
                                      yin(489) = c10*yin(449) + yc00*yin(481)
                                      zin(489) = c10*zin(449) + zc00*zin(481)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  491
                                      ! i5 =  489
                                      ! i4 =  481

                                      xin(491) = xcp00*xin(489) + cp10*xin(481)
                                      yin(491) = ycp00*yin(489) + cp10*yin(481)
                                      zin(491) = zcp00*zin(489) + cp10*zin(481)

                                      ! ------------------

                                      ! i3 = i4 =  481
                                      ! i4 = i5 =  489

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  497
                                      ! i3 =  481
                                      ! i4 =  489

                                      xin(497) = c10*xin(481) + xc00*xin(489)
                                      yin(497) = c10*yin(481) + yc00*yin(489)
                                      zin(497) = c10*zin(481) + zc00*zin(489)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  499
                                      ! i5 =  497
                                      ! i4 =  489

                                      xin(499) = xcp00*xin(497) + cp10*xin(489)
                                      yin(499) = ycp00*yin(497) + cp10*yin(489)
                                      zin(499) = zcp00*zin(497) + cp10*zin(489)

                                      ! ------------------

                                      ! i3 = i4 =  489
                                      ! i4 = i5 =  497

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  505
                                      ! i3 =  489
                                      ! i4 =  497

                                      xin(505) = c10*xin(489) + xc00*xin(497)
                                      yin(505) = c10*yin(489) + yc00*yin(497)
                                      zin(505) = c10*zin(489) + zc00*zin(497)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  507
                                      ! i5 =  505
                                      ! i4 =  497

                                      xin(507) = xcp00*xin(505) + cp10*xin(497)
                                      yin(507) = ycp00*yin(505) + cp10*yin(497)
                                      zin(507) = zcp00*zin(505) + cp10*zin(497)

                                      ! ------------------

                                      ! i3 = i4 =  497
                                      ! i4 = i5 =  505

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  385
                                      ! i4 = i1+k2 =  387

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  389
                                      ! i3 =  385
                                      ! i4 =  387

                                      xin(389) = cp01*xin(385) + xcp00*xin(387)
                                      yin(389) = cp01*yin(385) + ycp00*yin(387)
                                      zin(389) = cp01*zin(385) + zcp00*zin(387)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  421

                                      xin(421) = xc00*xin(389) + c01*xin(387)
                                      yin(421) = yc00*yin(389) + c01*yin(387)
                                      zin(421) = zc00*zin(389) + c01*zin(387)

                                      ! ------------------

                                      ! i3 = i4 =  387
                                      ! i4 = i5 =  389

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  391
                                      ! i3 =  387
                                      ! i4 =  389

                                      xin(391) = cp01*xin(387) + xcp00*xin(389)
                                      yin(391) = cp01*yin(387) + ycp00*yin(389)
                                      zin(391) = cp01*zin(387) + zcp00*zin(389)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  423

                                      xin(423) = xc00*xin(391) + c01*xin(389)
                                      yin(423) = yc00*yin(391) + c01*yin(389)
                                      zin(423) = zc00*zin(391) + c01*zin(389)

                                      ! ------------------

                                      ! i3 = i4 =  389
                                      ! i4 = i5 =  391

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  392
                                      ! i3 =  389
                                      ! i4 =  391

                                      xin(392) = cp01*xin(389) + xcp00*xin(391)
                                      yin(392) = cp01*yin(389) + ycp00*yin(391)
                                      zin(392) = cp01*zin(389) + zcp00*zin(391)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  424

                                      xin(424) = xc00*xin(392) + c01*xin(391)
                                      yin(424) = yc00*yin(392) + c01*yin(391)
                                      zin(424) = zc00*zin(392) + c01*zin(391)

                                      ! ------------------

                                      ! i3 = i4 =  391
                                      ! i4 = i5 =  392

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  417

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  449

                                      xin(453) = c10*xin(389) + xc00*xin(421) + c01*xin(419)
                                      yin(453) = c10*yin(389) + yc00*yin(421) + c01*yin(419)
                                      zin(453) = c10*zin(389) + zc00*zin(421) + c01*zin(419)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  417
                                      ! i4 = i5 =  449

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  481

                                      xin(485) = c10*xin(421) + xc00*xin(453) + c01*xin(451)
                                      yin(485) = c10*yin(421) + yc00*yin(453) + c01*yin(451)
                                      zin(485) = c10*zin(421) + zc00*zin(453) + c01*zin(451)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  449
                                      ! i4 = i5 =  481

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  489

                                      xin(493) = c10*xin(453) + xc00*xin(485) + c01*xin(483)
                                      yin(493) = c10*yin(453) + yc00*yin(485) + c01*yin(483)
                                      zin(493) = c10*zin(453) + zc00*zin(485) + c01*zin(483)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  481
                                      ! i4 = i5 =  489

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  497

                                      xin(501) = c10*xin(485) + xc00*xin(493) + c01*xin(491)
                                      yin(501) = c10*yin(485) + yc00*yin(493) + c01*yin(491)
                                      zin(501) = c10*zin(485) + zc00*zin(493) + c01*zin(491)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  489
                                      ! i4 = i5 =  497

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  505

                                      xin(509) = c10*xin(493) + xc00*xin(501) + c01*xin(499)
                                      yin(509) = c10*yin(493) + yc00*yin(501) + c01*yin(499)
                                      zin(509) = c10*zin(493) + zc00*zin(501) + c01*zin(499)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  497
                                      ! i4 = i5 =  505

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  417

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  449

                                      xin(455) = c10*xin(391) + xc00*xin(423) + c01*xin(421)
                                      yin(455) = c10*yin(391) + yc00*yin(423) + c01*yin(421)
                                      zin(455) = c10*zin(391) + zc00*zin(423) + c01*zin(421)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  417
                                      ! i4 = i5 =  449

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  481

                                      xin(487) = c10*xin(423) + xc00*xin(455) + c01*xin(453)
                                      yin(487) = c10*yin(423) + yc00*yin(455) + c01*yin(453)
                                      zin(487) = c10*zin(423) + zc00*zin(455) + c01*zin(453)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  449
                                      ! i4 = i5 =  481

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  489

                                      xin(495) = c10*xin(455) + xc00*xin(487) + c01*xin(485)
                                      yin(495) = c10*yin(455) + yc00*yin(487) + c01*yin(485)
                                      zin(495) = c10*zin(455) + zc00*zin(487) + c01*zin(485)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  481
                                      ! i4 = i5 =  489

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  497

                                      xin(503) = c10*xin(487) + xc00*xin(495) + c01*xin(493)
                                      yin(503) = c10*yin(487) + yc00*yin(495) + c01*yin(493)
                                      zin(503) = c10*zin(487) + zc00*zin(495) + c01*zin(493)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  489
                                      ! i4 = i5 =  497

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  505

                                      xin(511) = c10*xin(495) + xc00*xin(503) + c01*xin(501)
                                      yin(511) = c10*yin(495) + yc00*yin(503) + c01*yin(501)
                                      zin(511) = c10*zin(495) + zc00*zin(503) + c01*zin(501)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  497
                                      ! i4 = i5 =  505

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  417

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  449

                                      xin(456) = c10*xin(392) + xc00*xin(424) + c01*xin(423)
                                      yin(456) = c10*yin(392) + yc00*yin(424) + c01*yin(423)
                                      zin(456) = c10*zin(392) + zc00*zin(424) + c01*zin(423)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  417
                                      ! i4 = i5 =  449

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  481

                                      xin(488) = c10*xin(424) + xc00*xin(456) + c01*xin(455)
                                      yin(488) = c10*yin(424) + yc00*yin(456) + c01*yin(455)
                                      zin(488) = c10*zin(424) + zc00*zin(456) + c01*zin(455)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  449
                                      ! i4 = i5 =  481

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  489

                                      xin(496) = c10*xin(456) + xc00*xin(488) + c01*xin(487)
                                      yin(496) = c10*yin(456) + yc00*yin(488) + c01*yin(487)
                                      zin(496) = c10*zin(456) + zc00*zin(488) + c01*zin(487)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  481
                                      ! i4 = i5 =  489

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  497

                                      xin(504) = c10*xin(488) + xc00*xin(496) + c01*xin(495)
                                      yin(504) = c10*yin(488) + yc00*yin(496) + c01*yin(495)
                                      zin(504) = c10*zin(488) + zc00*zin(496) + c01*zin(495)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  489
                                      ! i4 = i5 =  497

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  505

                                      xin(512) = c10*xin(496) + xc00*xin(504) + c01*xin(503)
                                      yin(512) = c10*yin(496) + yc00*yin(504) + c01*yin(503)
                                      zin(512) = c10*zin(496) + zc00*zin(504) + c01*zin(503)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  497
                                      ! i4 = i5 =  505

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  505

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  505

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  497

                                      xin(505) = xin(505) + dxij*xin(497)
                                      yin(505) = yin(505) + dyij*yin(497)
                                      zin(505) = zin(505) + dzij*zin(497)

                                      ! i3 = i4 =  497
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  489

                                      xin(497) = xin(497) + dxij*xin(489)
                                      yin(497) = yin(497) + dyij*yin(489)
                                      zin(497) = zin(497) + dzij*zin(489)

                                      ! i3 = i4 =  489
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  481

                                      xin(489) = xin(489) + dxij*xin(481)
                                      yin(489) = yin(489) + dyij*yin(481)
                                      zin(489) = zin(489) + dzij*zin(481)

                                      ! i3 = i4 =  481
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  505

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  497

                                      xin(505) = xin(505) + dxij*xin(497)
                                      yin(505) = yin(505) + dyij*yin(497)
                                      zin(505) = zin(505) + dzij*zin(497)

                                      ! i3 = i4 =  497
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  489

                                      xin(497) = xin(497) + dxij*xin(489)
                                      yin(497) = yin(497) + dyij*yin(489)
                                      zin(497) = zin(497) + dzij*zin(489)

                                      ! i3 = i4 =  489
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  505

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  497

                                      xin(505) = xin(505) + dxij*xin(497)
                                      yin(505) = yin(505) + dyij*yin(497)
                                      zin(505) = zin(505) + dzij*zin(497)

                                      ! i3 = i4 =  497
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  393

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  393

                                      ! do ni = 1,    3

                                      xin(393) = xin(417) + dxij*xin(385)
                                      yin(393) = yin(417) + dyij*yin(385)
                                      zin(393) = zin(417) + dzij*zin(385)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  425

                                      ! ni =    2

                                      xin(425) = xin(449) + dxij*xin(417)
                                      yin(425) = yin(449) + dyij*yin(417)
                                      zin(425) = zin(449) + dzij*zin(417)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  457

                                      ! ni =    3

                                      xin(457) = xin(481) + dxij*xin(449)
                                      yin(457) = yin(481) + dyij*yin(449)
                                      zin(457) = zin(481) + dzij*zin(449)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  489

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  401

                                      ! nj =    2

                                      ! i4 = i3 =  401

                                      ! do ni = 1,    3

                                      xin(401) = xin(425) + dxij*xin(393)
                                      yin(401) = yin(425) + dyij*yin(393)
                                      zin(401) = zin(425) + dzij*zin(393)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  433

                                      ! ni =    2

                                      xin(433) = xin(457) + dxij*xin(425)
                                      yin(433) = yin(457) + dyij*yin(425)
                                      zin(433) = zin(457) + dzij*zin(425)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  465

                                      ! ni =    3

                                      xin(465) = xin(489) + dxij*xin(457)
                                      yin(465) = yin(489) + dyij*yin(457)
                                      zin(465) = zin(489) + dzij*zin(457)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  497

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  409

                                      ! nj =    3

                                      ! i4 = i3 =  409

                                      ! do ni = 1,    3

                                      xin(409) = xin(433) + dxij*xin(401)
                                      yin(409) = yin(433) + dyij*yin(401)
                                      zin(409) = zin(433) + dzij*zin(401)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  441

                                      ! ni =    2

                                      xin(441) = xin(465) + dxij*xin(433)
                                      yin(441) = yin(465) + dyij*yin(433)
                                      zin(441) = zin(465) + dzij*zin(433)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  473

                                      ! ni =    3

                                      xin(473) = xin(497) + dxij*xin(465)
                                      yin(473) = yin(497) + dyij*yin(465)
                                      zin(473) = zin(497) + dzij*zin(465)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  505

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  417

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  507

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  499

                                      xin(507) = xin(507) + dxij*xin(499)
                                      yin(507) = yin(507) + dyij*yin(499)
                                      zin(507) = zin(507) + dzij*zin(499)

                                      ! i3 = i4 =  499
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  491

                                      xin(499) = xin(499) + dxij*xin(491)
                                      yin(499) = yin(499) + dyij*yin(491)
                                      zin(499) = zin(499) + dzij*zin(491)

                                      ! i3 = i4 =  491
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  483

                                      xin(491) = xin(491) + dxij*xin(483)
                                      yin(491) = yin(491) + dyij*yin(483)
                                      zin(491) = zin(491) + dzij*zin(483)

                                      ! i3 = i4 =  483
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  507

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  499

                                      xin(507) = xin(507) + dxij*xin(499)
                                      yin(507) = yin(507) + dyij*yin(499)
                                      zin(507) = zin(507) + dzij*zin(499)

                                      ! i3 = i4 =  499
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  491

                                      xin(499) = xin(499) + dxij*xin(491)
                                      yin(499) = yin(499) + dyij*yin(491)
                                      zin(499) = zin(499) + dzij*zin(491)

                                      ! i3 = i4 =  491
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  507

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  499

                                      xin(507) = xin(507) + dxij*xin(499)
                                      yin(507) = yin(507) + dyij*yin(499)
                                      zin(507) = zin(507) + dzij*zin(499)

                                      ! i3 = i4 =  499
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  395

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  395

                                      ! do ni = 1,    3

                                      xin(395) = xin(419) + dxij*xin(387)
                                      yin(395) = yin(419) + dyij*yin(387)
                                      zin(395) = zin(419) + dzij*zin(387)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  427

                                      ! ni =    2

                                      xin(427) = xin(451) + dxij*xin(419)
                                      yin(427) = yin(451) + dyij*yin(419)
                                      zin(427) = zin(451) + dzij*zin(419)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  459

                                      ! ni =    3

                                      xin(459) = xin(483) + dxij*xin(451)
                                      yin(459) = yin(483) + dyij*yin(451)
                                      zin(459) = zin(483) + dzij*zin(451)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  491

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  403

                                      ! nj =    2

                                      ! i4 = i3 =  403

                                      ! do ni = 1,    3

                                      xin(403) = xin(427) + dxij*xin(395)
                                      yin(403) = yin(427) + dyij*yin(395)
                                      zin(403) = zin(427) + dzij*zin(395)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  435

                                      ! ni =    2

                                      xin(435) = xin(459) + dxij*xin(427)
                                      yin(435) = yin(459) + dyij*yin(427)
                                      zin(435) = zin(459) + dzij*zin(427)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  467

                                      ! ni =    3

                                      xin(467) = xin(491) + dxij*xin(459)
                                      yin(467) = yin(491) + dyij*yin(459)
                                      zin(467) = zin(491) + dzij*zin(459)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  499

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  411

                                      ! nj =    3

                                      ! i4 = i3 =  411

                                      ! do ni = 1,    3

                                      xin(411) = xin(435) + dxij*xin(403)
                                      yin(411) = yin(435) + dyij*yin(403)
                                      zin(411) = zin(435) + dzij*zin(403)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  443

                                      ! ni =    2

                                      xin(443) = xin(467) + dxij*xin(435)
                                      yin(443) = yin(467) + dyij*yin(435)
                                      zin(443) = zin(467) + dzij*zin(435)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  475

                                      ! ni =    3

                                      xin(475) = xin(499) + dxij*xin(467)
                                      yin(475) = yin(499) + dyij*yin(467)
                                      zin(475) = zin(499) + dzij*zin(467)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  507

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  419

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  509

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  501

                                      xin(509) = xin(509) + dxij*xin(501)
                                      yin(509) = yin(509) + dyij*yin(501)
                                      zin(509) = zin(509) + dzij*zin(501)

                                      ! i3 = i4 =  501
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  493

                                      xin(501) = xin(501) + dxij*xin(493)
                                      yin(501) = yin(501) + dyij*yin(493)
                                      zin(501) = zin(501) + dzij*zin(493)

                                      ! i3 = i4 =  493
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  485

                                      xin(493) = xin(493) + dxij*xin(485)
                                      yin(493) = yin(493) + dyij*yin(485)
                                      zin(493) = zin(493) + dzij*zin(485)

                                      ! i3 = i4 =  485
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  509

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  501

                                      xin(509) = xin(509) + dxij*xin(501)
                                      yin(509) = yin(509) + dyij*yin(501)
                                      zin(509) = zin(509) + dzij*zin(501)

                                      ! i3 = i4 =  501
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  493

                                      xin(501) = xin(501) + dxij*xin(493)
                                      yin(501) = yin(501) + dyij*yin(493)
                                      zin(501) = zin(501) + dzij*zin(493)

                                      ! i3 = i4 =  493
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  509

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  501

                                      xin(509) = xin(509) + dxij*xin(501)
                                      yin(509) = yin(509) + dyij*yin(501)
                                      zin(509) = zin(509) + dzij*zin(501)

                                      ! i3 = i4 =  501
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  397

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  397

                                      ! do ni = 1,    3

                                      xin(397) = xin(421) + dxij*xin(389)
                                      yin(397) = yin(421) + dyij*yin(389)
                                      zin(397) = zin(421) + dzij*zin(389)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  429

                                      ! ni =    2

                                      xin(429) = xin(453) + dxij*xin(421)
                                      yin(429) = yin(453) + dyij*yin(421)
                                      zin(429) = zin(453) + dzij*zin(421)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  461

                                      ! ni =    3

                                      xin(461) = xin(485) + dxij*xin(453)
                                      yin(461) = yin(485) + dyij*yin(453)
                                      zin(461) = zin(485) + dzij*zin(453)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  493

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  405

                                      ! nj =    2

                                      ! i4 = i3 =  405

                                      ! do ni = 1,    3

                                      xin(405) = xin(429) + dxij*xin(397)
                                      yin(405) = yin(429) + dyij*yin(397)
                                      zin(405) = zin(429) + dzij*zin(397)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  437

                                      ! ni =    2

                                      xin(437) = xin(461) + dxij*xin(429)
                                      yin(437) = yin(461) + dyij*yin(429)
                                      zin(437) = zin(461) + dzij*zin(429)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  469

                                      ! ni =    3

                                      xin(469) = xin(493) + dxij*xin(461)
                                      yin(469) = yin(493) + dyij*yin(461)
                                      zin(469) = zin(493) + dzij*zin(461)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  501

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  413

                                      ! nj =    3

                                      ! i4 = i3 =  413

                                      ! do ni = 1,    3

                                      xin(413) = xin(437) + dxij*xin(405)
                                      yin(413) = yin(437) + dyij*yin(405)
                                      zin(413) = zin(437) + dzij*zin(405)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  445

                                      ! ni =    2

                                      xin(445) = xin(469) + dxij*xin(437)
                                      yin(445) = yin(469) + dyij*yin(437)
                                      zin(445) = zin(469) + dzij*zin(437)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  477

                                      ! ni =    3

                                      xin(477) = xin(501) + dxij*xin(469)
                                      yin(477) = yin(501) + dyij*yin(469)
                                      zin(477) = zin(501) + dzij*zin(469)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  509

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  421

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  511

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  503

                                      xin(511) = xin(511) + dxij*xin(503)
                                      yin(511) = yin(511) + dyij*yin(503)
                                      zin(511) = zin(511) + dzij*zin(503)

                                      ! i3 = i4 =  503
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  495

                                      xin(503) = xin(503) + dxij*xin(495)
                                      yin(503) = yin(503) + dyij*yin(495)
                                      zin(503) = zin(503) + dzij*zin(495)

                                      ! i3 = i4 =  495
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  487

                                      xin(495) = xin(495) + dxij*xin(487)
                                      yin(495) = yin(495) + dyij*yin(487)
                                      zin(495) = zin(495) + dzij*zin(487)

                                      ! i3 = i4 =  487
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  511

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  503

                                      xin(511) = xin(511) + dxij*xin(503)
                                      yin(511) = yin(511) + dyij*yin(503)
                                      zin(511) = zin(511) + dzij*zin(503)

                                      ! i3 = i4 =  503
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  495

                                      xin(503) = xin(503) + dxij*xin(495)
                                      yin(503) = yin(503) + dyij*yin(495)
                                      zin(503) = zin(503) + dzij*zin(495)

                                      ! i3 = i4 =  495
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  511

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  503

                                      xin(511) = xin(511) + dxij*xin(503)
                                      yin(511) = yin(511) + dyij*yin(503)
                                      zin(511) = zin(511) + dzij*zin(503)

                                      ! i3 = i4 =  503
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  399

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  399

                                      ! do ni = 1,    3

                                      xin(399) = xin(423) + dxij*xin(391)
                                      yin(399) = yin(423) + dyij*yin(391)
                                      zin(399) = zin(423) + dzij*zin(391)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  431

                                      ! ni =    2

                                      xin(431) = xin(455) + dxij*xin(423)
                                      yin(431) = yin(455) + dyij*yin(423)
                                      zin(431) = zin(455) + dzij*zin(423)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  463

                                      ! ni =    3

                                      xin(463) = xin(487) + dxij*xin(455)
                                      yin(463) = yin(487) + dyij*yin(455)
                                      zin(463) = zin(487) + dzij*zin(455)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  495

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  407

                                      ! nj =    2

                                      ! i4 = i3 =  407

                                      ! do ni = 1,    3

                                      xin(407) = xin(431) + dxij*xin(399)
                                      yin(407) = yin(431) + dyij*yin(399)
                                      zin(407) = zin(431) + dzij*zin(399)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  439

                                      ! ni =    2

                                      xin(439) = xin(463) + dxij*xin(431)
                                      yin(439) = yin(463) + dyij*yin(431)
                                      zin(439) = zin(463) + dzij*zin(431)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  471

                                      ! ni =    3

                                      xin(471) = xin(495) + dxij*xin(463)
                                      yin(471) = yin(495) + dyij*yin(463)
                                      zin(471) = zin(495) + dzij*zin(463)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  503

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  415

                                      ! nj =    3

                                      ! i4 = i3 =  415

                                      ! do ni = 1,    3

                                      xin(415) = xin(439) + dxij*xin(407)
                                      yin(415) = yin(439) + dyij*yin(407)
                                      zin(415) = zin(439) + dzij*zin(407)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  447

                                      ! ni =    2

                                      xin(447) = xin(471) + dxij*xin(439)
                                      yin(447) = yin(471) + dyij*yin(439)
                                      zin(447) = zin(471) + dzij*zin(439)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  479

                                      ! ni =    3

                                      xin(479) = xin(503) + dxij*xin(471)
                                      yin(479) = yin(503) + dyij*yin(471)
                                      zin(479) = zin(503) + dzij*zin(471)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  511

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  423

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  512

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  504

                                      xin(512) = xin(512) + dxij*xin(504)
                                      yin(512) = yin(512) + dyij*yin(504)
                                      zin(512) = zin(512) + dzij*zin(504)

                                      ! i3 = i4 =  504
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  496

                                      xin(504) = xin(504) + dxij*xin(496)
                                      yin(504) = yin(504) + dyij*yin(496)
                                      zin(504) = zin(504) + dzij*zin(496)

                                      ! i3 = i4 =  496
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  488

                                      xin(496) = xin(496) + dxij*xin(488)
                                      yin(496) = yin(496) + dyij*yin(488)
                                      zin(496) = zin(496) + dzij*zin(488)

                                      ! i3 = i4 =  488
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  512

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  504

                                      xin(512) = xin(512) + dxij*xin(504)
                                      yin(512) = yin(512) + dyij*yin(504)
                                      zin(512) = zin(512) + dzij*zin(504)

                                      ! i3 = i4 =  504
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  496

                                      xin(504) = xin(504) + dxij*xin(496)
                                      yin(504) = yin(504) + dyij*yin(496)
                                      zin(504) = zin(504) + dzij*zin(496)

                                      ! i3 = i4 =  496
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  512

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  504

                                      xin(512) = xin(512) + dxij*xin(504)
                                      yin(512) = yin(512) + dyij*yin(504)
                                      zin(512) = zin(512) + dzij*zin(504)

                                      ! i3 = i4 =  504
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  400

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  400

                                      ! do ni = 1,    3

                                      xin(400) = xin(424) + dxij*xin(392)
                                      yin(400) = yin(424) + dyij*yin(392)
                                      zin(400) = zin(424) + dzij*zin(392)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  432

                                      ! ni =    2

                                      xin(432) = xin(456) + dxij*xin(424)
                                      yin(432) = yin(456) + dyij*yin(424)
                                      zin(432) = zin(456) + dzij*zin(424)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  464

                                      ! ni =    3

                                      xin(464) = xin(488) + dxij*xin(456)
                                      yin(464) = yin(488) + dyij*yin(456)
                                      zin(464) = zin(488) + dzij*zin(456)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  496

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  408

                                      ! nj =    2

                                      ! i4 = i3 =  408

                                      ! do ni = 1,    3

                                      xin(408) = xin(432) + dxij*xin(400)
                                      yin(408) = yin(432) + dyij*yin(400)
                                      zin(408) = zin(432) + dzij*zin(400)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  440

                                      ! ni =    2

                                      xin(440) = xin(464) + dxij*xin(432)
                                      yin(440) = yin(464) + dyij*yin(432)
                                      zin(440) = zin(464) + dzij*zin(432)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  472

                                      ! ni =    3

                                      xin(472) = xin(496) + dxij*xin(464)
                                      yin(472) = yin(496) + dyij*yin(464)
                                      zin(472) = zin(496) + dzij*zin(464)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  504

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  416

                                      ! nj =    3

                                      ! i4 = i3 =  416

                                      ! do ni = 1,    3

                                      xin(416) = xin(440) + dxij*xin(408)
                                      yin(416) = yin(440) + dyij*yin(408)
                                      zin(416) = zin(440) + dzij*zin(408)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  448

                                      ! ni =    2

                                      xin(448) = xin(472) + dxij*xin(440)
                                      yin(448) = yin(472) + dyij*yin(440)
                                      zin(448) = zin(472) + dzij*zin(440)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  480

                                      ! ni =    3

                                      xin(480) = xin(504) + dxij*xin(472)
                                      yin(480) = yin(504) + dyij*yin(472)
                                      zin(480) = zin(504) + dzij*zin(472)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  512

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  424

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  385

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  392

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  391

                                      xin(392) = xin(392) + dxkl*xin(391)
                                      yin(392) = yin(392) + dykl*yin(391)
                                      zin(392) = zin(392) + dzkl*zin(391)

                                      ! i3 = i4 =  391
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  386

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  386

                                      ! do nk = 1,    3

                                      xin(386) = xin(387) + dxkl*xin(385)
                                      yin(386) = yin(387) + dykl*yin(385)
                                      zin(386) = zin(387) + dzkl*zin(385)
                                      ! i4 = i4 + lang+1 =  388

                                      ! nk =    2

                                      xin(388) = xin(389) + dxkl*xin(387)
                                      yin(388) = yin(389) + dykl*yin(387)
                                      zin(388) = zin(389) + dzkl*zin(387)
                                      ! i4 = i4 + lang+1 =  390

                                      ! nk =    3

                                      xin(390) = xin(391) + dxkl*xin(389)
                                      yin(390) = yin(391) + dykl*yin(389)
                                      zin(390) = zin(391) + dzkl*zin(389)
                                      ! i4 = i4 + lang+1 =  392

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  387

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  393

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  400

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  399

                                      xin(400) = xin(400) + dxkl*xin(399)
                                      yin(400) = yin(400) + dykl*yin(399)
                                      zin(400) = zin(400) + dzkl*zin(399)

                                      ! i3 = i4 =  399
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  394

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  394

                                      ! do nk = 1,    3

                                      xin(394) = xin(395) + dxkl*xin(393)
                                      yin(394) = yin(395) + dykl*yin(393)
                                      zin(394) = zin(395) + dzkl*zin(393)
                                      ! i4 = i4 + lang+1 =  396

                                      ! nk =    2

                                      xin(396) = xin(397) + dxkl*xin(395)
                                      yin(396) = yin(397) + dykl*yin(395)
                                      zin(396) = zin(397) + dzkl*zin(395)
                                      ! i4 = i4 + lang+1 =  398

                                      ! nk =    3

                                      xin(398) = xin(399) + dxkl*xin(397)
                                      yin(398) = yin(399) + dykl*yin(397)
                                      zin(398) = zin(399) + dzkl*zin(397)
                                      ! i4 = i4 + lang+1 =  400

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  395

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  401

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  408

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  407

                                      xin(408) = xin(408) + dxkl*xin(407)
                                      yin(408) = yin(408) + dykl*yin(407)
                                      zin(408) = zin(408) + dzkl*zin(407)

                                      ! i3 = i4 =  407
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  402

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  402

                                      ! do nk = 1,    3

                                      xin(402) = xin(403) + dxkl*xin(401)
                                      yin(402) = yin(403) + dykl*yin(401)
                                      zin(402) = zin(403) + dzkl*zin(401)
                                      ! i4 = i4 + lang+1 =  404

                                      ! nk =    2

                                      xin(404) = xin(405) + dxkl*xin(403)
                                      yin(404) = yin(405) + dykl*yin(403)
                                      zin(404) = zin(405) + dzkl*zin(403)
                                      ! i4 = i4 + lang+1 =  406

                                      ! nk =    3

                                      xin(406) = xin(407) + dxkl*xin(405)
                                      yin(406) = yin(407) + dykl*yin(405)
                                      zin(406) = zin(407) + dzkl*zin(405)
                                      ! i4 = i4 + lang+1 =  408

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  403

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  409

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  416

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  415

                                      xin(416) = xin(416) + dxkl*xin(415)
                                      yin(416) = yin(416) + dykl*yin(415)
                                      zin(416) = zin(416) + dzkl*zin(415)

                                      ! i3 = i4 =  415
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  410

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  410

                                      ! do nk = 1,    3

                                      xin(410) = xin(411) + dxkl*xin(409)
                                      yin(410) = yin(411) + dykl*yin(409)
                                      zin(410) = zin(411) + dzkl*zin(409)
                                      ! i4 = i4 + lang+1 =  412

                                      ! nk =    2

                                      xin(412) = xin(413) + dxkl*xin(411)
                                      yin(412) = yin(413) + dykl*yin(411)
                                      zin(412) = zin(413) + dzkl*zin(411)
                                      ! i4 = i4 + lang+1 =  414

                                      ! nk =    3

                                      xin(414) = xin(415) + dxkl*xin(413)
                                      yin(414) = yin(415) + dykl*yin(413)
                                      zin(414) = zin(415) + dzkl*zin(413)
                                      ! i4 = i4 + lang+1 =  416

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  411

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  417

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  417

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  424

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  423

                                      xin(424) = xin(424) + dxkl*xin(423)
                                      yin(424) = yin(424) + dykl*yin(423)
                                      zin(424) = zin(424) + dzkl*zin(423)

                                      ! i3 = i4 =  423
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  418

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  418

                                      ! do nk = 1,    3

                                      xin(418) = xin(419) + dxkl*xin(417)
                                      yin(418) = yin(419) + dykl*yin(417)
                                      zin(418) = zin(419) + dzkl*zin(417)
                                      ! i4 = i4 + lang+1 =  420

                                      ! nk =    2

                                      xin(420) = xin(421) + dxkl*xin(419)
                                      yin(420) = yin(421) + dykl*yin(419)
                                      zin(420) = zin(421) + dzkl*zin(419)
                                      ! i4 = i4 + lang+1 =  422

                                      ! nk =    3

                                      xin(422) = xin(423) + dxkl*xin(421)
                                      yin(422) = yin(423) + dykl*yin(421)
                                      zin(422) = zin(423) + dzkl*zin(421)
                                      ! i4 = i4 + lang+1 =  424

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  419

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  425

                                      ! nj = nj + 1 =    1

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  426

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  426

                                      ! do nk = 1,    3

                                      xin(426) = xin(427) + dxkl*xin(425)
                                      yin(426) = yin(427) + dykl*yin(425)
                                      zin(426) = zin(427) + dzkl*zin(425)
                                      ! i4 = i4 + lang+1 =  428

                                      ! nk =    2

                                      xin(428) = xin(429) + dxkl*xin(427)
                                      yin(428) = yin(429) + dykl*yin(427)
                                      zin(428) = zin(429) + dzkl*zin(427)
                                      ! i4 = i4 + lang+1 =  430

                                      ! nk =    3

                                      xin(430) = xin(431) + dxkl*xin(429)
                                      yin(430) = yin(431) + dykl*yin(429)
                                      zin(430) = zin(431) + dzkl*zin(429)
                                      ! i4 = i4 + lang+1 =  432

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  427

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  433

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  440

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  439

                                      xin(440) = xin(440) + dxkl*xin(439)
                                      yin(440) = yin(440) + dykl*yin(439)
                                      zin(440) = zin(440) + dzkl*zin(439)

                                      ! i3 = i4 =  439
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  434

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  434

                                      ! do nk = 1,    3

                                      xin(434) = xin(435) + dxkl*xin(433)
                                      yin(434) = yin(435) + dykl*yin(433)
                                      zin(434) = zin(435) + dzkl*zin(433)
                                      ! i4 = i4 + lang+1 =  436

                                      ! nk =    2

                                      xin(436) = xin(437) + dxkl*xin(435)
                                      yin(436) = yin(437) + dykl*yin(435)
                                      zin(436) = zin(437) + dzkl*zin(435)
                                      ! i4 = i4 + lang+1 =  438

                                      ! nk =    3

                                      xin(438) = xin(439) + dxkl*xin(437)
                                      yin(438) = yin(439) + dykl*yin(437)
                                      zin(438) = zin(439) + dzkl*zin(437)
                                      ! i4 = i4 + lang+1 =  440

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  435

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  441

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  448

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  447

                                      xin(448) = xin(448) + dxkl*xin(447)
                                      yin(448) = yin(448) + dykl*yin(447)
                                      zin(448) = zin(448) + dzkl*zin(447)

                                      ! i3 = i4 =  447
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  442

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  442

                                      ! do nk = 1,    3

                                      xin(442) = xin(443) + dxkl*xin(441)
                                      yin(442) = yin(443) + dykl*yin(441)
                                      zin(442) = zin(443) + dzkl*zin(441)
                                      ! i4 = i4 + lang+1 =  444

                                      ! nk =    2

                                      xin(444) = xin(445) + dxkl*xin(443)
                                      yin(444) = yin(445) + dykl*yin(443)
                                      zin(444) = zin(445) + dzkl*zin(443)
                                      ! i4 = i4 + lang+1 =  446

                                      ! nk =    3

                                      xin(446) = xin(447) + dxkl*xin(445)
                                      yin(446) = yin(447) + dykl*yin(445)
                                      zin(446) = zin(447) + dzkl*zin(445)
                                      ! i4 = i4 + lang+1 =  448

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  443

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  449

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  449

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  456

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  455

                                      xin(456) = xin(456) + dxkl*xin(455)
                                      yin(456) = yin(456) + dykl*yin(455)
                                      zin(456) = zin(456) + dzkl*zin(455)

                                      ! i3 = i4 =  455
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  450

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  450

                                      ! do nk = 1,    3

                                      xin(450) = xin(451) + dxkl*xin(449)
                                      yin(450) = yin(451) + dykl*yin(449)
                                      zin(450) = zin(451) + dzkl*zin(449)
                                      ! i4 = i4 + lang+1 =  452

                                      ! nk =    2

                                      xin(452) = xin(453) + dxkl*xin(451)
                                      yin(452) = yin(453) + dykl*yin(451)
                                      zin(452) = zin(453) + dzkl*zin(451)
                                      ! i4 = i4 + lang+1 =  454

                                      ! nk =    3

                                      xin(454) = xin(455) + dxkl*xin(453)
                                      yin(454) = yin(455) + dykl*yin(453)
                                      zin(454) = zin(455) + dzkl*zin(453)
                                      ! i4 = i4 + lang+1 =  456

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  451

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  457

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  464

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  463

                                      xin(464) = xin(464) + dxkl*xin(463)
                                      yin(464) = yin(464) + dykl*yin(463)
                                      zin(464) = zin(464) + dzkl*zin(463)

                                      ! i3 = i4 =  463
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  458

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  458

                                      ! do nk = 1,    3

                                      xin(458) = xin(459) + dxkl*xin(457)
                                      yin(458) = yin(459) + dykl*yin(457)
                                      zin(458) = zin(459) + dzkl*zin(457)
                                      ! i4 = i4 + lang+1 =  460

                                      ! nk =    2

                                      xin(460) = xin(461) + dxkl*xin(459)
                                      yin(460) = yin(461) + dykl*yin(459)
                                      zin(460) = zin(461) + dzkl*zin(459)
                                      ! i4 = i4 + lang+1 =  462

                                      ! nk =    3

                                      xin(462) = xin(463) + dxkl*xin(461)
                                      yin(462) = yin(463) + dykl*yin(461)
                                      zin(462) = zin(463) + dzkl*zin(461)
                                      ! i4 = i4 + lang+1 =  464

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  459

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  465

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  472

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  471

                                      xin(472) = xin(472) + dxkl*xin(471)
                                      yin(472) = yin(472) + dykl*yin(471)
                                      zin(472) = zin(472) + dzkl*zin(471)

                                      ! i3 = i4 =  471
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  466

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  466

                                      ! do nk = 1,    3

                                      xin(466) = xin(467) + dxkl*xin(465)
                                      yin(466) = yin(467) + dykl*yin(465)
                                      zin(466) = zin(467) + dzkl*zin(465)
                                      ! i4 = i4 + lang+1 =  468

                                      ! nk =    2

                                      xin(468) = xin(469) + dxkl*xin(467)
                                      yin(468) = yin(469) + dykl*yin(467)
                                      zin(468) = zin(469) + dzkl*zin(467)
                                      ! i4 = i4 + lang+1 =  470

                                      ! nk =    3

                                      xin(470) = xin(471) + dxkl*xin(469)
                                      yin(470) = yin(471) + dykl*yin(469)
                                      zin(470) = zin(471) + dzkl*zin(469)
                                      ! i4 = i4 + lang+1 =  472

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  467

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  473

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  480

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  479

                                      xin(480) = xin(480) + dxkl*xin(479)
                                      yin(480) = yin(480) + dykl*yin(479)
                                      zin(480) = zin(480) + dzkl*zin(479)

                                      ! i3 = i4 =  479
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  474

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  474

                                      ! do nk = 1,    3

                                      xin(474) = xin(475) + dxkl*xin(473)
                                      yin(474) = yin(475) + dykl*yin(473)
                                      zin(474) = zin(475) + dzkl*zin(473)
                                      ! i4 = i4 + lang+1 =  476

                                      ! nk =    2

                                      xin(476) = xin(477) + dxkl*xin(475)
                                      yin(476) = yin(477) + dykl*yin(475)
                                      zin(476) = zin(477) + dzkl*zin(475)
                                      ! i4 = i4 + lang+1 =  478

                                      ! nk =    3

                                      xin(478) = xin(479) + dxkl*xin(477)
                                      yin(478) = yin(479) + dykl*yin(477)
                                      zin(478) = zin(479) + dzkl*zin(477)
                                      ! i4 = i4 + lang+1 =  480

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  475

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  481

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  481

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  488

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  487

                                      xin(488) = xin(488) + dxkl*xin(487)
                                      yin(488) = yin(488) + dykl*yin(487)
                                      zin(488) = zin(488) + dzkl*zin(487)

                                      ! i3 = i4 =  487
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  482

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  482

                                      ! do nk = 1,    3

                                      xin(482) = xin(483) + dxkl*xin(481)
                                      yin(482) = yin(483) + dykl*yin(481)
                                      zin(482) = zin(483) + dzkl*zin(481)
                                      ! i4 = i4 + lang+1 =  484

                                      ! nk =    2

                                      xin(484) = xin(485) + dxkl*xin(483)
                                      yin(484) = yin(485) + dykl*yin(483)
                                      zin(484) = zin(485) + dzkl*zin(483)
                                      ! i4 = i4 + lang+1 =  486

                                      ! nk =    3

                                      xin(486) = xin(487) + dxkl*xin(485)
                                      yin(486) = yin(487) + dykl*yin(485)
                                      zin(486) = zin(487) + dzkl*zin(485)
                                      ! i4 = i4 + lang+1 =  488

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  483

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  489

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  496

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  495

                                      xin(496) = xin(496) + dxkl*xin(495)
                                      yin(496) = yin(496) + dykl*yin(495)
                                      zin(496) = zin(496) + dzkl*zin(495)

                                      ! i3 = i4 =  495
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  490

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  490

                                      ! do nk = 1,    3

                                      xin(490) = xin(491) + dxkl*xin(489)
                                      yin(490) = yin(491) + dykl*yin(489)
                                      zin(490) = zin(491) + dzkl*zin(489)
                                      ! i4 = i4 + lang+1 =  492

                                      ! nk =    2

                                      xin(492) = xin(493) + dxkl*xin(491)
                                      yin(492) = yin(493) + dykl*yin(491)
                                      zin(492) = zin(493) + dzkl*zin(491)
                                      ! i4 = i4 + lang+1 =  494

                                      ! nk =    3

                                      xin(494) = xin(495) + dxkl*xin(493)
                                      yin(494) = yin(495) + dykl*yin(493)
                                      zin(494) = zin(495) + dzkl*zin(493)
                                      ! i4 = i4 + lang+1 =  496

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  491

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  497

                                      ! nj = nj + 1 =    2

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  498

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  498

                                      ! do nk = 1,    3

                                      xin(498) = xin(499) + dxkl*xin(497)
                                      yin(498) = yin(499) + dykl*yin(497)
                                      zin(498) = zin(499) + dzkl*zin(497)
                                      ! i4 = i4 + lang+1 =  500

                                      ! nk =    2

                                      xin(500) = xin(501) + dxkl*xin(499)
                                      yin(500) = yin(501) + dykl*yin(499)
                                      zin(500) = zin(501) + dzkl*zin(499)
                                      ! i4 = i4 + lang+1 =  502

                                      ! nk =    3

                                      xin(502) = xin(503) + dxkl*xin(501)
                                      yin(502) = yin(503) + dykl*yin(501)
                                      zin(502) = zin(503) + dzkl*zin(501)
                                      ! i4 = i4 + lang+1 =  504

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  499

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  505

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  512

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  511

                                      xin(512) = xin(512) + dxkl*xin(511)
                                      yin(512) = yin(512) + dykl*yin(511)
                                      zin(512) = zin(512) + dzkl*zin(511)

                                      ! i3 = i4 =  511
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  506

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  506

                                      ! do nk = 1,    3

                                      xin(506) = xin(507) + dxkl*xin(505)
                                      yin(506) = yin(507) + dykl*yin(505)
                                      zin(506) = zin(507) + dzkl*zin(505)
                                      ! i4 = i4 + lang+1 =  508

                                      ! nk =    2

                                      xin(508) = xin(509) + dxkl*xin(507)
                                      yin(508) = yin(509) + dykl*yin(507)
                                      zin(508) = zin(509) + dzkl*zin(507)
                                      ! i4 = i4 + lang+1 =  510

                                      ! nk =    3

                                      xin(510) = xin(511) + dxkl*xin(509)
                                      yin(510) = yin(511) + dykl*yin(509)
                                      zin(510) = zin(511) + dzkl*zin(509)
                                      ! i4 = i4 + lang+1 =  512

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  507

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  513

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  513

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  512

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

                                      ! i1 = in(1) =  513

                                      xin(513) = 1.0_dp
                                      yin(513) = 1.0_dp
                                      zin(513) = f00

                                      ! i2 = in(2) =  545
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(545) = xc00
                                      yin(545) = yc00
                                      zin(545) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  515

                                      xin(515) = xcp00
                                      yin(515) = ycp00
                                      zin(515) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  547
                                      ! i2 =  545

                                      xin(547) = xcp00*xin(545) + cp10
                                      yin(547) = ycp00*yin(545) + cp10
                                      zin(547) = zcp00*zin(545) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  513
                                      ! i4 = i2 =  545

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  577
                                      ! i3 =  513
                                      ! i4 =  545

                                      xin(577) = c10*xin(513) + xc00*xin(545)
                                      yin(577) = c10*yin(513) + yc00*yin(545)
                                      zin(577) = c10*zin(513) + zc00*zin(545)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  579
                                      ! i5 =  577
                                      ! i4 =  545

                                      xin(579) = xcp00*xin(577) + cp10*xin(545)
                                      yin(579) = ycp00*yin(577) + cp10*yin(545)
                                      zin(579) = zcp00*zin(577) + cp10*zin(545)

                                      ! ------------------

                                      ! i3 = i4 =  545
                                      ! i4 = i5 =  577

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  609
                                      ! i3 =  545
                                      ! i4 =  577

                                      xin(609) = c10*xin(545) + xc00*xin(577)
                                      yin(609) = c10*yin(545) + yc00*yin(577)
                                      zin(609) = c10*zin(545) + zc00*zin(577)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  611
                                      ! i5 =  609
                                      ! i4 =  577

                                      xin(611) = xcp00*xin(609) + cp10*xin(577)
                                      yin(611) = ycp00*yin(609) + cp10*yin(577)
                                      zin(611) = zcp00*zin(609) + cp10*zin(577)

                                      ! ------------------

                                      ! i3 = i4 =  577
                                      ! i4 = i5 =  609

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  617
                                      ! i3 =  577
                                      ! i4 =  609

                                      xin(617) = c10*xin(577) + xc00*xin(609)
                                      yin(617) = c10*yin(577) + yc00*yin(609)
                                      zin(617) = c10*zin(577) + zc00*zin(609)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  619
                                      ! i5 =  617
                                      ! i4 =  609

                                      xin(619) = xcp00*xin(617) + cp10*xin(609)
                                      yin(619) = ycp00*yin(617) + cp10*yin(609)
                                      zin(619) = zcp00*zin(617) + cp10*zin(609)

                                      ! ------------------

                                      ! i3 = i4 =  609
                                      ! i4 = i5 =  617

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  625
                                      ! i3 =  609
                                      ! i4 =  617

                                      xin(625) = c10*xin(609) + xc00*xin(617)
                                      yin(625) = c10*yin(609) + yc00*yin(617)
                                      zin(625) = c10*zin(609) + zc00*zin(617)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  627
                                      ! i5 =  625
                                      ! i4 =  617

                                      xin(627) = xcp00*xin(625) + cp10*xin(617)
                                      yin(627) = ycp00*yin(625) + cp10*yin(617)
                                      zin(627) = zcp00*zin(625) + cp10*zin(617)

                                      ! ------------------

                                      ! i3 = i4 =  617
                                      ! i4 = i5 =  625

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  633
                                      ! i3 =  617
                                      ! i4 =  625

                                      xin(633) = c10*xin(617) + xc00*xin(625)
                                      yin(633) = c10*yin(617) + yc00*yin(625)
                                      zin(633) = c10*zin(617) + zc00*zin(625)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  635
                                      ! i5 =  633
                                      ! i4 =  625

                                      xin(635) = xcp00*xin(633) + cp10*xin(625)
                                      yin(635) = ycp00*yin(633) + cp10*yin(625)
                                      zin(635) = zcp00*zin(633) + cp10*zin(625)

                                      ! ------------------

                                      ! i3 = i4 =  625
                                      ! i4 = i5 =  633

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  513
                                      ! i4 = i1+k2 =  515

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  517
                                      ! i3 =  513
                                      ! i4 =  515

                                      xin(517) = cp01*xin(513) + xcp00*xin(515)
                                      yin(517) = cp01*yin(513) + ycp00*yin(515)
                                      zin(517) = cp01*zin(513) + zcp00*zin(515)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  549

                                      xin(549) = xc00*xin(517) + c01*xin(515)
                                      yin(549) = yc00*yin(517) + c01*yin(515)
                                      zin(549) = zc00*zin(517) + c01*zin(515)

                                      ! ------------------

                                      ! i3 = i4 =  515
                                      ! i4 = i5 =  517

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  519
                                      ! i3 =  515
                                      ! i4 =  517

                                      xin(519) = cp01*xin(515) + xcp00*xin(517)
                                      yin(519) = cp01*yin(515) + ycp00*yin(517)
                                      zin(519) = cp01*zin(515) + zcp00*zin(517)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  551

                                      xin(551) = xc00*xin(519) + c01*xin(517)
                                      yin(551) = yc00*yin(519) + c01*yin(517)
                                      zin(551) = zc00*zin(519) + c01*zin(517)

                                      ! ------------------

                                      ! i3 = i4 =  517
                                      ! i4 = i5 =  519

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  520
                                      ! i3 =  517
                                      ! i4 =  519

                                      xin(520) = cp01*xin(517) + xcp00*xin(519)
                                      yin(520) = cp01*yin(517) + ycp00*yin(519)
                                      zin(520) = cp01*zin(517) + zcp00*zin(519)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  552

                                      xin(552) = xc00*xin(520) + c01*xin(519)
                                      yin(552) = yc00*yin(520) + c01*yin(519)
                                      zin(552) = zc00*zin(520) + c01*zin(519)

                                      ! ------------------

                                      ! i3 = i4 =  519
                                      ! i4 = i5 =  520

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  513
                                      ! i4 = i2 =  545

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  577

                                      xin(581) = c10*xin(517) + xc00*xin(549) + c01*xin(547)
                                      yin(581) = c10*yin(517) + yc00*yin(549) + c01*yin(547)
                                      zin(581) = c10*zin(517) + zc00*zin(549) + c01*zin(547)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  545
                                      ! i4 = i5 =  577

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  609

                                      xin(613) = c10*xin(549) + xc00*xin(581) + c01*xin(579)
                                      yin(613) = c10*yin(549) + yc00*yin(581) + c01*yin(579)
                                      zin(613) = c10*zin(549) + zc00*zin(581) + c01*zin(579)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  577
                                      ! i4 = i5 =  609

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  617

                                      xin(621) = c10*xin(581) + xc00*xin(613) + c01*xin(611)
                                      yin(621) = c10*yin(581) + yc00*yin(613) + c01*yin(611)
                                      zin(621) = c10*zin(581) + zc00*zin(613) + c01*zin(611)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  609
                                      ! i4 = i5 =  617

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  625

                                      xin(629) = c10*xin(613) + xc00*xin(621) + c01*xin(619)
                                      yin(629) = c10*yin(613) + yc00*yin(621) + c01*yin(619)
                                      zin(629) = c10*zin(613) + zc00*zin(621) + c01*zin(619)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  617
                                      ! i4 = i5 =  625

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  633

                                      xin(637) = c10*xin(621) + xc00*xin(629) + c01*xin(627)
                                      yin(637) = c10*yin(621) + yc00*yin(629) + c01*yin(627)
                                      zin(637) = c10*zin(621) + zc00*zin(629) + c01*zin(627)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  625
                                      ! i4 = i5 =  633

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  513
                                      ! i4 = i2 =  545

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  577

                                      xin(583) = c10*xin(519) + xc00*xin(551) + c01*xin(549)
                                      yin(583) = c10*yin(519) + yc00*yin(551) + c01*yin(549)
                                      zin(583) = c10*zin(519) + zc00*zin(551) + c01*zin(549)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  545
                                      ! i4 = i5 =  577

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  609

                                      xin(615) = c10*xin(551) + xc00*xin(583) + c01*xin(581)
                                      yin(615) = c10*yin(551) + yc00*yin(583) + c01*yin(581)
                                      zin(615) = c10*zin(551) + zc00*zin(583) + c01*zin(581)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  577
                                      ! i4 = i5 =  609

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  617

                                      xin(623) = c10*xin(583) + xc00*xin(615) + c01*xin(613)
                                      yin(623) = c10*yin(583) + yc00*yin(615) + c01*yin(613)
                                      zin(623) = c10*zin(583) + zc00*zin(615) + c01*zin(613)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  609
                                      ! i4 = i5 =  617

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  625

                                      xin(631) = c10*xin(615) + xc00*xin(623) + c01*xin(621)
                                      yin(631) = c10*yin(615) + yc00*yin(623) + c01*yin(621)
                                      zin(631) = c10*zin(615) + zc00*zin(623) + c01*zin(621)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  617
                                      ! i4 = i5 =  625

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  633

                                      xin(639) = c10*xin(623) + xc00*xin(631) + c01*xin(629)
                                      yin(639) = c10*yin(623) + yc00*yin(631) + c01*yin(629)
                                      zin(639) = c10*zin(623) + zc00*zin(631) + c01*zin(629)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  625
                                      ! i4 = i5 =  633

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  513
                                      ! i4 = i2 =  545

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  577

                                      xin(584) = c10*xin(520) + xc00*xin(552) + c01*xin(551)
                                      yin(584) = c10*yin(520) + yc00*yin(552) + c01*yin(551)
                                      zin(584) = c10*zin(520) + zc00*zin(552) + c01*zin(551)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  545
                                      ! i4 = i5 =  577

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  609

                                      xin(616) = c10*xin(552) + xc00*xin(584) + c01*xin(583)
                                      yin(616) = c10*yin(552) + yc00*yin(584) + c01*yin(583)
                                      zin(616) = c10*zin(552) + zc00*zin(584) + c01*zin(583)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  577
                                      ! i4 = i5 =  609

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  617

                                      xin(624) = c10*xin(584) + xc00*xin(616) + c01*xin(615)
                                      yin(624) = c10*yin(584) + yc00*yin(616) + c01*yin(615)
                                      zin(624) = c10*zin(584) + zc00*zin(616) + c01*zin(615)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  609
                                      ! i4 = i5 =  617

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  625

                                      xin(632) = c10*xin(616) + xc00*xin(624) + c01*xin(623)
                                      yin(632) = c10*yin(616) + yc00*yin(624) + c01*yin(623)
                                      zin(632) = c10*zin(616) + zc00*zin(624) + c01*zin(623)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  617
                                      ! i4 = i5 =  625

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  633

                                      xin(640) = c10*xin(624) + xc00*xin(632) + c01*xin(631)
                                      yin(640) = c10*yin(624) + yc00*yin(632) + c01*yin(631)
                                      zin(640) = c10*zin(624) + zc00*zin(632) + c01*zin(631)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  625
                                      ! i4 = i5 =  633

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  633

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  633

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  625

                                      xin(633) = xin(633) + dxij*xin(625)
                                      yin(633) = yin(633) + dyij*yin(625)
                                      zin(633) = zin(633) + dzij*zin(625)

                                      ! i3 = i4 =  625
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  617

                                      xin(625) = xin(625) + dxij*xin(617)
                                      yin(625) = yin(625) + dyij*yin(617)
                                      zin(625) = zin(625) + dzij*zin(617)

                                      ! i3 = i4 =  617
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  609

                                      xin(617) = xin(617) + dxij*xin(609)
                                      yin(617) = yin(617) + dyij*yin(609)
                                      zin(617) = zin(617) + dzij*zin(609)

                                      ! i3 = i4 =  609
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  633

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  625

                                      xin(633) = xin(633) + dxij*xin(625)
                                      yin(633) = yin(633) + dyij*yin(625)
                                      zin(633) = zin(633) + dzij*zin(625)

                                      ! i3 = i4 =  625
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  617

                                      xin(625) = xin(625) + dxij*xin(617)
                                      yin(625) = yin(625) + dyij*yin(617)
                                      zin(625) = zin(625) + dzij*zin(617)

                                      ! i3 = i4 =  617
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  633

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  625

                                      xin(633) = xin(633) + dxij*xin(625)
                                      yin(633) = yin(633) + dyij*yin(625)
                                      zin(633) = zin(633) + dzij*zin(625)

                                      ! i3 = i4 =  625
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  521

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  521

                                      ! do ni = 1,    3

                                      xin(521) = xin(545) + dxij*xin(513)
                                      yin(521) = yin(545) + dyij*yin(513)
                                      zin(521) = zin(545) + dzij*zin(513)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  553

                                      ! ni =    2

                                      xin(553) = xin(577) + dxij*xin(545)
                                      yin(553) = yin(577) + dyij*yin(545)
                                      zin(553) = zin(577) + dzij*zin(545)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  585

                                      ! ni =    3

                                      xin(585) = xin(609) + dxij*xin(577)
                                      yin(585) = yin(609) + dyij*yin(577)
                                      zin(585) = zin(609) + dzij*zin(577)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  617

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  529

                                      ! nj =    2

                                      ! i4 = i3 =  529

                                      ! do ni = 1,    3

                                      xin(529) = xin(553) + dxij*xin(521)
                                      yin(529) = yin(553) + dyij*yin(521)
                                      zin(529) = zin(553) + dzij*zin(521)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  561

                                      ! ni =    2

                                      xin(561) = xin(585) + dxij*xin(553)
                                      yin(561) = yin(585) + dyij*yin(553)
                                      zin(561) = zin(585) + dzij*zin(553)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  593

                                      ! ni =    3

                                      xin(593) = xin(617) + dxij*xin(585)
                                      yin(593) = yin(617) + dyij*yin(585)
                                      zin(593) = zin(617) + dzij*zin(585)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  625

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  537

                                      ! nj =    3

                                      ! i4 = i3 =  537

                                      ! do ni = 1,    3

                                      xin(537) = xin(561) + dxij*xin(529)
                                      yin(537) = yin(561) + dyij*yin(529)
                                      zin(537) = zin(561) + dzij*zin(529)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  569

                                      ! ni =    2

                                      xin(569) = xin(593) + dxij*xin(561)
                                      yin(569) = yin(593) + dyij*yin(561)
                                      zin(569) = zin(593) + dzij*zin(561)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  601

                                      ! ni =    3

                                      xin(601) = xin(625) + dxij*xin(593)
                                      yin(601) = yin(625) + dyij*yin(593)
                                      zin(601) = zin(625) + dzij*zin(593)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  633

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  545

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  635

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  627

                                      xin(635) = xin(635) + dxij*xin(627)
                                      yin(635) = yin(635) + dyij*yin(627)
                                      zin(635) = zin(635) + dzij*zin(627)

                                      ! i3 = i4 =  627
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  619

                                      xin(627) = xin(627) + dxij*xin(619)
                                      yin(627) = yin(627) + dyij*yin(619)
                                      zin(627) = zin(627) + dzij*zin(619)

                                      ! i3 = i4 =  619
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  611

                                      xin(619) = xin(619) + dxij*xin(611)
                                      yin(619) = yin(619) + dyij*yin(611)
                                      zin(619) = zin(619) + dzij*zin(611)

                                      ! i3 = i4 =  611
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  635

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  627

                                      xin(635) = xin(635) + dxij*xin(627)
                                      yin(635) = yin(635) + dyij*yin(627)
                                      zin(635) = zin(635) + dzij*zin(627)

                                      ! i3 = i4 =  627
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  619

                                      xin(627) = xin(627) + dxij*xin(619)
                                      yin(627) = yin(627) + dyij*yin(619)
                                      zin(627) = zin(627) + dzij*zin(619)

                                      ! i3 = i4 =  619
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  635

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  627

                                      xin(635) = xin(635) + dxij*xin(627)
                                      yin(635) = yin(635) + dyij*yin(627)
                                      zin(635) = zin(635) + dzij*zin(627)

                                      ! i3 = i4 =  627
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  523

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  523

                                      ! do ni = 1,    3

                                      xin(523) = xin(547) + dxij*xin(515)
                                      yin(523) = yin(547) + dyij*yin(515)
                                      zin(523) = zin(547) + dzij*zin(515)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  555

                                      ! ni =    2

                                      xin(555) = xin(579) + dxij*xin(547)
                                      yin(555) = yin(579) + dyij*yin(547)
                                      zin(555) = zin(579) + dzij*zin(547)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  587

                                      ! ni =    3

                                      xin(587) = xin(611) + dxij*xin(579)
                                      yin(587) = yin(611) + dyij*yin(579)
                                      zin(587) = zin(611) + dzij*zin(579)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  619

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  531

                                      ! nj =    2

                                      ! i4 = i3 =  531

                                      ! do ni = 1,    3

                                      xin(531) = xin(555) + dxij*xin(523)
                                      yin(531) = yin(555) + dyij*yin(523)
                                      zin(531) = zin(555) + dzij*zin(523)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  563

                                      ! ni =    2

                                      xin(563) = xin(587) + dxij*xin(555)
                                      yin(563) = yin(587) + dyij*yin(555)
                                      zin(563) = zin(587) + dzij*zin(555)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  595

                                      ! ni =    3

                                      xin(595) = xin(619) + dxij*xin(587)
                                      yin(595) = yin(619) + dyij*yin(587)
                                      zin(595) = zin(619) + dzij*zin(587)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  627

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  539

                                      ! nj =    3

                                      ! i4 = i3 =  539

                                      ! do ni = 1,    3

                                      xin(539) = xin(563) + dxij*xin(531)
                                      yin(539) = yin(563) + dyij*yin(531)
                                      zin(539) = zin(563) + dzij*zin(531)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  571

                                      ! ni =    2

                                      xin(571) = xin(595) + dxij*xin(563)
                                      yin(571) = yin(595) + dyij*yin(563)
                                      zin(571) = zin(595) + dzij*zin(563)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  603

                                      ! ni =    3

                                      xin(603) = xin(627) + dxij*xin(595)
                                      yin(603) = yin(627) + dyij*yin(595)
                                      zin(603) = zin(627) + dzij*zin(595)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  635

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  547

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  637

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  629

                                      xin(637) = xin(637) + dxij*xin(629)
                                      yin(637) = yin(637) + dyij*yin(629)
                                      zin(637) = zin(637) + dzij*zin(629)

                                      ! i3 = i4 =  629
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  621

                                      xin(629) = xin(629) + dxij*xin(621)
                                      yin(629) = yin(629) + dyij*yin(621)
                                      zin(629) = zin(629) + dzij*zin(621)

                                      ! i3 = i4 =  621
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  613

                                      xin(621) = xin(621) + dxij*xin(613)
                                      yin(621) = yin(621) + dyij*yin(613)
                                      zin(621) = zin(621) + dzij*zin(613)

                                      ! i3 = i4 =  613
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  637

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  629

                                      xin(637) = xin(637) + dxij*xin(629)
                                      yin(637) = yin(637) + dyij*yin(629)
                                      zin(637) = zin(637) + dzij*zin(629)

                                      ! i3 = i4 =  629
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  621

                                      xin(629) = xin(629) + dxij*xin(621)
                                      yin(629) = yin(629) + dyij*yin(621)
                                      zin(629) = zin(629) + dzij*zin(621)

                                      ! i3 = i4 =  621
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  637

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  629

                                      xin(637) = xin(637) + dxij*xin(629)
                                      yin(637) = yin(637) + dyij*yin(629)
                                      zin(637) = zin(637) + dzij*zin(629)

                                      ! i3 = i4 =  629
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  525

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  525

                                      ! do ni = 1,    3

                                      xin(525) = xin(549) + dxij*xin(517)
                                      yin(525) = yin(549) + dyij*yin(517)
                                      zin(525) = zin(549) + dzij*zin(517)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  557

                                      ! ni =    2

                                      xin(557) = xin(581) + dxij*xin(549)
                                      yin(557) = yin(581) + dyij*yin(549)
                                      zin(557) = zin(581) + dzij*zin(549)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  589

                                      ! ni =    3

                                      xin(589) = xin(613) + dxij*xin(581)
                                      yin(589) = yin(613) + dyij*yin(581)
                                      zin(589) = zin(613) + dzij*zin(581)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  621

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  533

                                      ! nj =    2

                                      ! i4 = i3 =  533

                                      ! do ni = 1,    3

                                      xin(533) = xin(557) + dxij*xin(525)
                                      yin(533) = yin(557) + dyij*yin(525)
                                      zin(533) = zin(557) + dzij*zin(525)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  565

                                      ! ni =    2

                                      xin(565) = xin(589) + dxij*xin(557)
                                      yin(565) = yin(589) + dyij*yin(557)
                                      zin(565) = zin(589) + dzij*zin(557)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  597

                                      ! ni =    3

                                      xin(597) = xin(621) + dxij*xin(589)
                                      yin(597) = yin(621) + dyij*yin(589)
                                      zin(597) = zin(621) + dzij*zin(589)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  629

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  541

                                      ! nj =    3

                                      ! i4 = i3 =  541

                                      ! do ni = 1,    3

                                      xin(541) = xin(565) + dxij*xin(533)
                                      yin(541) = yin(565) + dyij*yin(533)
                                      zin(541) = zin(565) + dzij*zin(533)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  573

                                      ! ni =    2

                                      xin(573) = xin(597) + dxij*xin(565)
                                      yin(573) = yin(597) + dyij*yin(565)
                                      zin(573) = zin(597) + dzij*zin(565)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  605

                                      ! ni =    3

                                      xin(605) = xin(629) + dxij*xin(597)
                                      yin(605) = yin(629) + dyij*yin(597)
                                      zin(605) = zin(629) + dzij*zin(597)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  637

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  549

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  639

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  631

                                      xin(639) = xin(639) + dxij*xin(631)
                                      yin(639) = yin(639) + dyij*yin(631)
                                      zin(639) = zin(639) + dzij*zin(631)

                                      ! i3 = i4 =  631
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  623

                                      xin(631) = xin(631) + dxij*xin(623)
                                      yin(631) = yin(631) + dyij*yin(623)
                                      zin(631) = zin(631) + dzij*zin(623)

                                      ! i3 = i4 =  623
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  615

                                      xin(623) = xin(623) + dxij*xin(615)
                                      yin(623) = yin(623) + dyij*yin(615)
                                      zin(623) = zin(623) + dzij*zin(615)

                                      ! i3 = i4 =  615
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  639

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  631

                                      xin(639) = xin(639) + dxij*xin(631)
                                      yin(639) = yin(639) + dyij*yin(631)
                                      zin(639) = zin(639) + dzij*zin(631)

                                      ! i3 = i4 =  631
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  623

                                      xin(631) = xin(631) + dxij*xin(623)
                                      yin(631) = yin(631) + dyij*yin(623)
                                      zin(631) = zin(631) + dzij*zin(623)

                                      ! i3 = i4 =  623
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  639

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  631

                                      xin(639) = xin(639) + dxij*xin(631)
                                      yin(639) = yin(639) + dyij*yin(631)
                                      zin(639) = zin(639) + dzij*zin(631)

                                      ! i3 = i4 =  631
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  527

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  527

                                      ! do ni = 1,    3

                                      xin(527) = xin(551) + dxij*xin(519)
                                      yin(527) = yin(551) + dyij*yin(519)
                                      zin(527) = zin(551) + dzij*zin(519)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  559

                                      ! ni =    2

                                      xin(559) = xin(583) + dxij*xin(551)
                                      yin(559) = yin(583) + dyij*yin(551)
                                      zin(559) = zin(583) + dzij*zin(551)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  591

                                      ! ni =    3

                                      xin(591) = xin(615) + dxij*xin(583)
                                      yin(591) = yin(615) + dyij*yin(583)
                                      zin(591) = zin(615) + dzij*zin(583)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  623

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  535

                                      ! nj =    2

                                      ! i4 = i3 =  535

                                      ! do ni = 1,    3

                                      xin(535) = xin(559) + dxij*xin(527)
                                      yin(535) = yin(559) + dyij*yin(527)
                                      zin(535) = zin(559) + dzij*zin(527)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  567

                                      ! ni =    2

                                      xin(567) = xin(591) + dxij*xin(559)
                                      yin(567) = yin(591) + dyij*yin(559)
                                      zin(567) = zin(591) + dzij*zin(559)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  599

                                      ! ni =    3

                                      xin(599) = xin(623) + dxij*xin(591)
                                      yin(599) = yin(623) + dyij*yin(591)
                                      zin(599) = zin(623) + dzij*zin(591)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  631

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  543

                                      ! nj =    3

                                      ! i4 = i3 =  543

                                      ! do ni = 1,    3

                                      xin(543) = xin(567) + dxij*xin(535)
                                      yin(543) = yin(567) + dyij*yin(535)
                                      zin(543) = zin(567) + dzij*zin(535)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  575

                                      ! ni =    2

                                      xin(575) = xin(599) + dxij*xin(567)
                                      yin(575) = yin(599) + dyij*yin(567)
                                      zin(575) = zin(599) + dzij*zin(567)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  607

                                      ! ni =    3

                                      xin(607) = xin(631) + dxij*xin(599)
                                      yin(607) = yin(631) + dyij*yin(599)
                                      zin(607) = zin(631) + dzij*zin(599)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  639

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  551

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  640

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  632

                                      xin(640) = xin(640) + dxij*xin(632)
                                      yin(640) = yin(640) + dyij*yin(632)
                                      zin(640) = zin(640) + dzij*zin(632)

                                      ! i3 = i4 =  632
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  624

                                      xin(632) = xin(632) + dxij*xin(624)
                                      yin(632) = yin(632) + dyij*yin(624)
                                      zin(632) = zin(632) + dzij*zin(624)

                                      ! i3 = i4 =  624
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  616

                                      xin(624) = xin(624) + dxij*xin(616)
                                      yin(624) = yin(624) + dyij*yin(616)
                                      zin(624) = zin(624) + dzij*zin(616)

                                      ! i3 = i4 =  616
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  640

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  632

                                      xin(640) = xin(640) + dxij*xin(632)
                                      yin(640) = yin(640) + dyij*yin(632)
                                      zin(640) = zin(640) + dzij*zin(632)

                                      ! i3 = i4 =  632
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  624

                                      xin(632) = xin(632) + dxij*xin(624)
                                      yin(632) = yin(632) + dyij*yin(624)
                                      zin(632) = zin(632) + dzij*zin(624)

                                      ! i3 = i4 =  624
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  640

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  632

                                      xin(640) = xin(640) + dxij*xin(632)
                                      yin(640) = yin(640) + dyij*yin(632)
                                      zin(640) = zin(640) + dzij*zin(632)

                                      ! i3 = i4 =  632
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  528

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  528

                                      ! do ni = 1,    3

                                      xin(528) = xin(552) + dxij*xin(520)
                                      yin(528) = yin(552) + dyij*yin(520)
                                      zin(528) = zin(552) + dzij*zin(520)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  560

                                      ! ni =    2

                                      xin(560) = xin(584) + dxij*xin(552)
                                      yin(560) = yin(584) + dyij*yin(552)
                                      zin(560) = zin(584) + dzij*zin(552)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  592

                                      ! ni =    3

                                      xin(592) = xin(616) + dxij*xin(584)
                                      yin(592) = yin(616) + dyij*yin(584)
                                      zin(592) = zin(616) + dzij*zin(584)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  624

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  536

                                      ! nj =    2

                                      ! i4 = i3 =  536

                                      ! do ni = 1,    3

                                      xin(536) = xin(560) + dxij*xin(528)
                                      yin(536) = yin(560) + dyij*yin(528)
                                      zin(536) = zin(560) + dzij*zin(528)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  568

                                      ! ni =    2

                                      xin(568) = xin(592) + dxij*xin(560)
                                      yin(568) = yin(592) + dyij*yin(560)
                                      zin(568) = zin(592) + dzij*zin(560)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  600

                                      ! ni =    3

                                      xin(600) = xin(624) + dxij*xin(592)
                                      yin(600) = yin(624) + dyij*yin(592)
                                      zin(600) = zin(624) + dzij*zin(592)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  632

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  544

                                      ! nj =    3

                                      ! i4 = i3 =  544

                                      ! do ni = 1,    3

                                      xin(544) = xin(568) + dxij*xin(536)
                                      yin(544) = yin(568) + dyij*yin(536)
                                      zin(544) = zin(568) + dzij*zin(536)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  576

                                      ! ni =    2

                                      xin(576) = xin(600) + dxij*xin(568)
                                      yin(576) = yin(600) + dyij*yin(568)
                                      zin(576) = zin(600) + dzij*zin(568)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  608

                                      ! ni =    3

                                      xin(608) = xin(632) + dxij*xin(600)
                                      yin(608) = yin(632) + dyij*yin(600)
                                      zin(608) = zin(632) + dzij*zin(600)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  640

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  552

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  513

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  520

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  519

                                      xin(520) = xin(520) + dxkl*xin(519)
                                      yin(520) = yin(520) + dykl*yin(519)
                                      zin(520) = zin(520) + dzkl*zin(519)

                                      ! i3 = i4 =  519
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  514

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  514

                                      ! do nk = 1,    3

                                      xin(514) = xin(515) + dxkl*xin(513)
                                      yin(514) = yin(515) + dykl*yin(513)
                                      zin(514) = zin(515) + dzkl*zin(513)
                                      ! i4 = i4 + lang+1 =  516

                                      ! nk =    2

                                      xin(516) = xin(517) + dxkl*xin(515)
                                      yin(516) = yin(517) + dykl*yin(515)
                                      zin(516) = zin(517) + dzkl*zin(515)
                                      ! i4 = i4 + lang+1 =  518

                                      ! nk =    3

                                      xin(518) = xin(519) + dxkl*xin(517)
                                      yin(518) = yin(519) + dykl*yin(517)
                                      zin(518) = zin(519) + dzkl*zin(517)
                                      ! i4 = i4 + lang+1 =  520

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  515

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  521

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  528

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  527

                                      xin(528) = xin(528) + dxkl*xin(527)
                                      yin(528) = yin(528) + dykl*yin(527)
                                      zin(528) = zin(528) + dzkl*zin(527)

                                      ! i3 = i4 =  527
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  522

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  522

                                      ! do nk = 1,    3

                                      xin(522) = xin(523) + dxkl*xin(521)
                                      yin(522) = yin(523) + dykl*yin(521)
                                      zin(522) = zin(523) + dzkl*zin(521)
                                      ! i4 = i4 + lang+1 =  524

                                      ! nk =    2

                                      xin(524) = xin(525) + dxkl*xin(523)
                                      yin(524) = yin(525) + dykl*yin(523)
                                      zin(524) = zin(525) + dzkl*zin(523)
                                      ! i4 = i4 + lang+1 =  526

                                      ! nk =    3

                                      xin(526) = xin(527) + dxkl*xin(525)
                                      yin(526) = yin(527) + dykl*yin(525)
                                      zin(526) = zin(527) + dzkl*zin(525)
                                      ! i4 = i4 + lang+1 =  528

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  523

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  529

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  536

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  535

                                      xin(536) = xin(536) + dxkl*xin(535)
                                      yin(536) = yin(536) + dykl*yin(535)
                                      zin(536) = zin(536) + dzkl*zin(535)

                                      ! i3 = i4 =  535
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  530

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  530

                                      ! do nk = 1,    3

                                      xin(530) = xin(531) + dxkl*xin(529)
                                      yin(530) = yin(531) + dykl*yin(529)
                                      zin(530) = zin(531) + dzkl*zin(529)
                                      ! i4 = i4 + lang+1 =  532

                                      ! nk =    2

                                      xin(532) = xin(533) + dxkl*xin(531)
                                      yin(532) = yin(533) + dykl*yin(531)
                                      zin(532) = zin(533) + dzkl*zin(531)
                                      ! i4 = i4 + lang+1 =  534

                                      ! nk =    3

                                      xin(534) = xin(535) + dxkl*xin(533)
                                      yin(534) = yin(535) + dykl*yin(533)
                                      zin(534) = zin(535) + dzkl*zin(533)
                                      ! i4 = i4 + lang+1 =  536

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  531

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  537

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  544

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  543

                                      xin(544) = xin(544) + dxkl*xin(543)
                                      yin(544) = yin(544) + dykl*yin(543)
                                      zin(544) = zin(544) + dzkl*zin(543)

                                      ! i3 = i4 =  543
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  538

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  538

                                      ! do nk = 1,    3

                                      xin(538) = xin(539) + dxkl*xin(537)
                                      yin(538) = yin(539) + dykl*yin(537)
                                      zin(538) = zin(539) + dzkl*zin(537)
                                      ! i4 = i4 + lang+1 =  540

                                      ! nk =    2

                                      xin(540) = xin(541) + dxkl*xin(539)
                                      yin(540) = yin(541) + dykl*yin(539)
                                      zin(540) = zin(541) + dzkl*zin(539)
                                      ! i4 = i4 + lang+1 =  542

                                      ! nk =    3

                                      xin(542) = xin(543) + dxkl*xin(541)
                                      yin(542) = yin(543) + dykl*yin(541)
                                      zin(542) = zin(543) + dzkl*zin(541)
                                      ! i4 = i4 + lang+1 =  544

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  539

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  545

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  545

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  552

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  551

                                      xin(552) = xin(552) + dxkl*xin(551)
                                      yin(552) = yin(552) + dykl*yin(551)
                                      zin(552) = zin(552) + dzkl*zin(551)

                                      ! i3 = i4 =  551
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  546

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  546

                                      ! do nk = 1,    3

                                      xin(546) = xin(547) + dxkl*xin(545)
                                      yin(546) = yin(547) + dykl*yin(545)
                                      zin(546) = zin(547) + dzkl*zin(545)
                                      ! i4 = i4 + lang+1 =  548

                                      ! nk =    2

                                      xin(548) = xin(549) + dxkl*xin(547)
                                      yin(548) = yin(549) + dykl*yin(547)
                                      zin(548) = zin(549) + dzkl*zin(547)
                                      ! i4 = i4 + lang+1 =  550

                                      ! nk =    3

                                      xin(550) = xin(551) + dxkl*xin(549)
                                      yin(550) = yin(551) + dykl*yin(549)
                                      zin(550) = zin(551) + dzkl*zin(549)
                                      ! i4 = i4 + lang+1 =  552

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  547

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  553

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  560

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  559

                                      xin(560) = xin(560) + dxkl*xin(559)
                                      yin(560) = yin(560) + dykl*yin(559)
                                      zin(560) = zin(560) + dzkl*zin(559)

                                      ! i3 = i4 =  559
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  554

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  554

                                      ! do nk = 1,    3

                                      xin(554) = xin(555) + dxkl*xin(553)
                                      yin(554) = yin(555) + dykl*yin(553)
                                      zin(554) = zin(555) + dzkl*zin(553)
                                      ! i4 = i4 + lang+1 =  556

                                      ! nk =    2

                                      xin(556) = xin(557) + dxkl*xin(555)
                                      yin(556) = yin(557) + dykl*yin(555)
                                      zin(556) = zin(557) + dzkl*zin(555)
                                      ! i4 = i4 + lang+1 =  558

                                      ! nk =    3

                                      xin(558) = xin(559) + dxkl*xin(557)
                                      yin(558) = yin(559) + dykl*yin(557)
                                      zin(558) = zin(559) + dzkl*zin(557)
                                      ! i4 = i4 + lang+1 =  560

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  555

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  561

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  568

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  567

                                      xin(568) = xin(568) + dxkl*xin(567)
                                      yin(568) = yin(568) + dykl*yin(567)
                                      zin(568) = zin(568) + dzkl*zin(567)

                                      ! i3 = i4 =  567
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  562

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  562

                                      ! do nk = 1,    3

                                      xin(562) = xin(563) + dxkl*xin(561)
                                      yin(562) = yin(563) + dykl*yin(561)
                                      zin(562) = zin(563) + dzkl*zin(561)
                                      ! i4 = i4 + lang+1 =  564

                                      ! nk =    2

                                      xin(564) = xin(565) + dxkl*xin(563)
                                      yin(564) = yin(565) + dykl*yin(563)
                                      zin(564) = zin(565) + dzkl*zin(563)
                                      ! i4 = i4 + lang+1 =  566

                                      ! nk =    3

                                      xin(566) = xin(567) + dxkl*xin(565)
                                      yin(566) = yin(567) + dykl*yin(565)
                                      zin(566) = zin(567) + dzkl*zin(565)
                                      ! i4 = i4 + lang+1 =  568

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  563

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  569

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  570

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  570

                                      ! do nk = 1,    3

                                      xin(570) = xin(571) + dxkl*xin(569)
                                      yin(570) = yin(571) + dykl*yin(569)
                                      zin(570) = zin(571) + dzkl*zin(569)
                                      ! i4 = i4 + lang+1 =  572

                                      ! nk =    2

                                      xin(572) = xin(573) + dxkl*xin(571)
                                      yin(572) = yin(573) + dykl*yin(571)
                                      zin(572) = zin(573) + dzkl*zin(571)
                                      ! i4 = i4 + lang+1 =  574

                                      ! nk =    3

                                      xin(574) = xin(575) + dxkl*xin(573)
                                      yin(574) = yin(575) + dykl*yin(573)
                                      zin(574) = zin(575) + dzkl*zin(573)
                                      ! i4 = i4 + lang+1 =  576

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  571

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  577

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  577

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  584

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  583

                                      xin(584) = xin(584) + dxkl*xin(583)
                                      yin(584) = yin(584) + dykl*yin(583)
                                      zin(584) = zin(584) + dzkl*zin(583)

                                      ! i3 = i4 =  583
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  578

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  578

                                      ! do nk = 1,    3

                                      xin(578) = xin(579) + dxkl*xin(577)
                                      yin(578) = yin(579) + dykl*yin(577)
                                      zin(578) = zin(579) + dzkl*zin(577)
                                      ! i4 = i4 + lang+1 =  580

                                      ! nk =    2

                                      xin(580) = xin(581) + dxkl*xin(579)
                                      yin(580) = yin(581) + dykl*yin(579)
                                      zin(580) = zin(581) + dzkl*zin(579)
                                      ! i4 = i4 + lang+1 =  582

                                      ! nk =    3

                                      xin(582) = xin(583) + dxkl*xin(581)
                                      yin(582) = yin(583) + dykl*yin(581)
                                      zin(582) = zin(583) + dzkl*zin(581)
                                      ! i4 = i4 + lang+1 =  584

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  579

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  585

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  592

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  591

                                      xin(592) = xin(592) + dxkl*xin(591)
                                      yin(592) = yin(592) + dykl*yin(591)
                                      zin(592) = zin(592) + dzkl*zin(591)

                                      ! i3 = i4 =  591
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  586

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  586

                                      ! do nk = 1,    3

                                      xin(586) = xin(587) + dxkl*xin(585)
                                      yin(586) = yin(587) + dykl*yin(585)
                                      zin(586) = zin(587) + dzkl*zin(585)
                                      ! i4 = i4 + lang+1 =  588

                                      ! nk =    2

                                      xin(588) = xin(589) + dxkl*xin(587)
                                      yin(588) = yin(589) + dykl*yin(587)
                                      zin(588) = zin(589) + dzkl*zin(587)
                                      ! i4 = i4 + lang+1 =  590

                                      ! nk =    3

                                      xin(590) = xin(591) + dxkl*xin(589)
                                      yin(590) = yin(591) + dykl*yin(589)
                                      zin(590) = zin(591) + dzkl*zin(589)
                                      ! i4 = i4 + lang+1 =  592

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  587

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  593

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  600

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  599

                                      xin(600) = xin(600) + dxkl*xin(599)
                                      yin(600) = yin(600) + dykl*yin(599)
                                      zin(600) = zin(600) + dzkl*zin(599)

                                      ! i3 = i4 =  599
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  594

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  594

                                      ! do nk = 1,    3

                                      xin(594) = xin(595) + dxkl*xin(593)
                                      yin(594) = yin(595) + dykl*yin(593)
                                      zin(594) = zin(595) + dzkl*zin(593)
                                      ! i4 = i4 + lang+1 =  596

                                      ! nk =    2

                                      xin(596) = xin(597) + dxkl*xin(595)
                                      yin(596) = yin(597) + dykl*yin(595)
                                      zin(596) = zin(597) + dzkl*zin(595)
                                      ! i4 = i4 + lang+1 =  598

                                      ! nk =    3

                                      xin(598) = xin(599) + dxkl*xin(597)
                                      yin(598) = yin(599) + dykl*yin(597)
                                      zin(598) = zin(599) + dzkl*zin(597)
                                      ! i4 = i4 + lang+1 =  600

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  595

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  601

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  608

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  607

                                      xin(608) = xin(608) + dxkl*xin(607)
                                      yin(608) = yin(608) + dykl*yin(607)
                                      zin(608) = zin(608) + dzkl*zin(607)

                                      ! i3 = i4 =  607
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  602

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  602

                                      ! do nk = 1,    3

                                      xin(602) = xin(603) + dxkl*xin(601)
                                      yin(602) = yin(603) + dykl*yin(601)
                                      zin(602) = zin(603) + dzkl*zin(601)
                                      ! i4 = i4 + lang+1 =  604

                                      ! nk =    2

                                      xin(604) = xin(605) + dxkl*xin(603)
                                      yin(604) = yin(605) + dykl*yin(603)
                                      zin(604) = zin(605) + dzkl*zin(603)
                                      ! i4 = i4 + lang+1 =  606

                                      ! nk =    3

                                      xin(606) = xin(607) + dxkl*xin(605)
                                      yin(606) = yin(607) + dykl*yin(605)
                                      zin(606) = zin(607) + dzkl*zin(605)
                                      ! i4 = i4 + lang+1 =  608

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  603

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  609

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  609

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  616

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  615

                                      xin(616) = xin(616) + dxkl*xin(615)
                                      yin(616) = yin(616) + dykl*yin(615)
                                      zin(616) = zin(616) + dzkl*zin(615)

                                      ! i3 = i4 =  615
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  610

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  610

                                      ! do nk = 1,    3

                                      xin(610) = xin(611) + dxkl*xin(609)
                                      yin(610) = yin(611) + dykl*yin(609)
                                      zin(610) = zin(611) + dzkl*zin(609)
                                      ! i4 = i4 + lang+1 =  612

                                      ! nk =    2

                                      xin(612) = xin(613) + dxkl*xin(611)
                                      yin(612) = yin(613) + dykl*yin(611)
                                      zin(612) = zin(613) + dzkl*zin(611)
                                      ! i4 = i4 + lang+1 =  614

                                      ! nk =    3

                                      xin(614) = xin(615) + dxkl*xin(613)
                                      yin(614) = yin(615) + dykl*yin(613)
                                      zin(614) = zin(615) + dzkl*zin(613)
                                      ! i4 = i4 + lang+1 =  616

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  611

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  617

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  624

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  623

                                      xin(624) = xin(624) + dxkl*xin(623)
                                      yin(624) = yin(624) + dykl*yin(623)
                                      zin(624) = zin(624) + dzkl*zin(623)

                                      ! i3 = i4 =  623
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  618

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  618

                                      ! do nk = 1,    3

                                      xin(618) = xin(619) + dxkl*xin(617)
                                      yin(618) = yin(619) + dykl*yin(617)
                                      zin(618) = zin(619) + dzkl*zin(617)
                                      ! i4 = i4 + lang+1 =  620

                                      ! nk =    2

                                      xin(620) = xin(621) + dxkl*xin(619)
                                      yin(620) = yin(621) + dykl*yin(619)
                                      zin(620) = zin(621) + dzkl*zin(619)
                                      ! i4 = i4 + lang+1 =  622

                                      ! nk =    3

                                      xin(622) = xin(623) + dxkl*xin(621)
                                      yin(622) = yin(623) + dykl*yin(621)
                                      zin(622) = zin(623) + dzkl*zin(621)
                                      ! i4 = i4 + lang+1 =  624

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  619

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  625

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  632

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  631

                                      xin(632) = xin(632) + dxkl*xin(631)
                                      yin(632) = yin(632) + dykl*yin(631)
                                      zin(632) = zin(632) + dzkl*zin(631)

                                      ! i3 = i4 =  631
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  626

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  626

                                      ! do nk = 1,    3

                                      xin(626) = xin(627) + dxkl*xin(625)
                                      yin(626) = yin(627) + dykl*yin(625)
                                      zin(626) = zin(627) + dzkl*zin(625)
                                      ! i4 = i4 + lang+1 =  628

                                      ! nk =    2

                                      xin(628) = xin(629) + dxkl*xin(627)
                                      yin(628) = yin(629) + dykl*yin(627)
                                      zin(628) = zin(629) + dzkl*zin(627)
                                      ! i4 = i4 + lang+1 =  630

                                      ! nk =    3

                                      xin(630) = xin(631) + dxkl*xin(629)
                                      yin(630) = yin(631) + dykl*yin(629)
                                      zin(630) = zin(631) + dzkl*zin(629)
                                      ! i4 = i4 + lang+1 =  632

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  627

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  633

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  640

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  639

                                      xin(640) = xin(640) + dxkl*xin(639)
                                      yin(640) = yin(640) + dykl*yin(639)
                                      zin(640) = zin(640) + dzkl*zin(639)

                                      ! i3 = i4 =  639
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  634

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  634

                                      ! do nk = 1,    3

                                      xin(634) = xin(635) + dxkl*xin(633)
                                      yin(634) = yin(635) + dykl*yin(633)
                                      zin(634) = zin(635) + dzkl*zin(633)
                                      ! i4 = i4 + lang+1 =  636

                                      ! nk =    2

                                      xin(636) = xin(637) + dxkl*xin(635)
                                      yin(636) = yin(637) + dykl*yin(635)
                                      zin(636) = zin(637) + dzkl*zin(635)
                                      ! i4 = i4 + lang+1 =  638

                                      ! nk =    3

                                      xin(638) = xin(639) + dxkl*xin(637)
                                      yin(638) = yin(639) + dykl*yin(637)
                                      zin(638) = zin(639) + dzkl*zin(637)
                                      ! i4 = i4 + lang+1 =  640

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  635

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  641

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  641

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  640

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

                                      ! i1 = in(1) =  641

                                      xin(641) = 1.0_dp
                                      yin(641) = 1.0_dp
                                      zin(641) = f00

                                      ! i2 = in(2) =  673
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(673) = xc00
                                      yin(673) = yc00
                                      zin(673) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  643

                                      xin(643) = xcp00
                                      yin(643) = ycp00
                                      zin(643) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  675
                                      ! i2 =  673

                                      xin(675) = xcp00*xin(673) + cp10
                                      yin(675) = ycp00*yin(673) + cp10
                                      zin(675) = zcp00*zin(673) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  641
                                      ! i4 = i2 =  673

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  705
                                      ! i3 =  641
                                      ! i4 =  673

                                      xin(705) = c10*xin(641) + xc00*xin(673)
                                      yin(705) = c10*yin(641) + yc00*yin(673)
                                      zin(705) = c10*zin(641) + zc00*zin(673)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  707
                                      ! i5 =  705
                                      ! i4 =  673

                                      xin(707) = xcp00*xin(705) + cp10*xin(673)
                                      yin(707) = ycp00*yin(705) + cp10*yin(673)
                                      zin(707) = zcp00*zin(705) + cp10*zin(673)

                                      ! ------------------

                                      ! i3 = i4 =  673
                                      ! i4 = i5 =  705

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  737
                                      ! i3 =  673
                                      ! i4 =  705

                                      xin(737) = c10*xin(673) + xc00*xin(705)
                                      yin(737) = c10*yin(673) + yc00*yin(705)
                                      zin(737) = c10*zin(673) + zc00*zin(705)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  739
                                      ! i5 =  737
                                      ! i4 =  705

                                      xin(739) = xcp00*xin(737) + cp10*xin(705)
                                      yin(739) = ycp00*yin(737) + cp10*yin(705)
                                      zin(739) = zcp00*zin(737) + cp10*zin(705)

                                      ! ------------------

                                      ! i3 = i4 =  705
                                      ! i4 = i5 =  737

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  745
                                      ! i3 =  705
                                      ! i4 =  737

                                      xin(745) = c10*xin(705) + xc00*xin(737)
                                      yin(745) = c10*yin(705) + yc00*yin(737)
                                      zin(745) = c10*zin(705) + zc00*zin(737)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  747
                                      ! i5 =  745
                                      ! i4 =  737

                                      xin(747) = xcp00*xin(745) + cp10*xin(737)
                                      yin(747) = ycp00*yin(745) + cp10*yin(737)
                                      zin(747) = zcp00*zin(745) + cp10*zin(737)

                                      ! ------------------

                                      ! i3 = i4 =  737
                                      ! i4 = i5 =  745

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  753
                                      ! i3 =  737
                                      ! i4 =  745

                                      xin(753) = c10*xin(737) + xc00*xin(745)
                                      yin(753) = c10*yin(737) + yc00*yin(745)
                                      zin(753) = c10*zin(737) + zc00*zin(745)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  755
                                      ! i5 =  753
                                      ! i4 =  745

                                      xin(755) = xcp00*xin(753) + cp10*xin(745)
                                      yin(755) = ycp00*yin(753) + cp10*yin(745)
                                      zin(755) = zcp00*zin(753) + cp10*zin(745)

                                      ! ------------------

                                      ! i3 = i4 =  745
                                      ! i4 = i5 =  753

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  761
                                      ! i3 =  745
                                      ! i4 =  753

                                      xin(761) = c10*xin(745) + xc00*xin(753)
                                      yin(761) = c10*yin(745) + yc00*yin(753)
                                      zin(761) = c10*zin(745) + zc00*zin(753)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  763
                                      ! i5 =  761
                                      ! i4 =  753

                                      xin(763) = xcp00*xin(761) + cp10*xin(753)
                                      yin(763) = ycp00*yin(761) + cp10*yin(753)
                                      zin(763) = zcp00*zin(761) + cp10*zin(753)

                                      ! ------------------

                                      ! i3 = i4 =  753
                                      ! i4 = i5 =  761

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  641
                                      ! i4 = i1+k2 =  643

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  645
                                      ! i3 =  641
                                      ! i4 =  643

                                      xin(645) = cp01*xin(641) + xcp00*xin(643)
                                      yin(645) = cp01*yin(641) + ycp00*yin(643)
                                      zin(645) = cp01*zin(641) + zcp00*zin(643)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  677

                                      xin(677) = xc00*xin(645) + c01*xin(643)
                                      yin(677) = yc00*yin(645) + c01*yin(643)
                                      zin(677) = zc00*zin(645) + c01*zin(643)

                                      ! ------------------

                                      ! i3 = i4 =  643
                                      ! i4 = i5 =  645

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  647
                                      ! i3 =  643
                                      ! i4 =  645

                                      xin(647) = cp01*xin(643) + xcp00*xin(645)
                                      yin(647) = cp01*yin(643) + ycp00*yin(645)
                                      zin(647) = cp01*zin(643) + zcp00*zin(645)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  679

                                      xin(679) = xc00*xin(647) + c01*xin(645)
                                      yin(679) = yc00*yin(647) + c01*yin(645)
                                      zin(679) = zc00*zin(647) + c01*zin(645)

                                      ! ------------------

                                      ! i3 = i4 =  645
                                      ! i4 = i5 =  647

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  648
                                      ! i3 =  645
                                      ! i4 =  647

                                      xin(648) = cp01*xin(645) + xcp00*xin(647)
                                      yin(648) = cp01*yin(645) + ycp00*yin(647)
                                      zin(648) = cp01*zin(645) + zcp00*zin(647)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  680

                                      xin(680) = xc00*xin(648) + c01*xin(647)
                                      yin(680) = yc00*yin(648) + c01*yin(647)
                                      zin(680) = zc00*zin(648) + c01*zin(647)

                                      ! ------------------

                                      ! i3 = i4 =  647
                                      ! i4 = i5 =  648

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  641
                                      ! i4 = i2 =  673

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  705

                                      xin(709) = c10*xin(645) + xc00*xin(677) + c01*xin(675)
                                      yin(709) = c10*yin(645) + yc00*yin(677) + c01*yin(675)
                                      zin(709) = c10*zin(645) + zc00*zin(677) + c01*zin(675)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  673
                                      ! i4 = i5 =  705

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  737

                                      xin(741) = c10*xin(677) + xc00*xin(709) + c01*xin(707)
                                      yin(741) = c10*yin(677) + yc00*yin(709) + c01*yin(707)
                                      zin(741) = c10*zin(677) + zc00*zin(709) + c01*zin(707)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  705
                                      ! i4 = i5 =  737

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  745

                                      xin(749) = c10*xin(709) + xc00*xin(741) + c01*xin(739)
                                      yin(749) = c10*yin(709) + yc00*yin(741) + c01*yin(739)
                                      zin(749) = c10*zin(709) + zc00*zin(741) + c01*zin(739)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  737
                                      ! i4 = i5 =  745

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  753

                                      xin(757) = c10*xin(741) + xc00*xin(749) + c01*xin(747)
                                      yin(757) = c10*yin(741) + yc00*yin(749) + c01*yin(747)
                                      zin(757) = c10*zin(741) + zc00*zin(749) + c01*zin(747)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  745
                                      ! i4 = i5 =  753

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  761

                                      xin(765) = c10*xin(749) + xc00*xin(757) + c01*xin(755)
                                      yin(765) = c10*yin(749) + yc00*yin(757) + c01*yin(755)
                                      zin(765) = c10*zin(749) + zc00*zin(757) + c01*zin(755)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  753
                                      ! i4 = i5 =  761

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  641
                                      ! i4 = i2 =  673

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  705

                                      xin(711) = c10*xin(647) + xc00*xin(679) + c01*xin(677)
                                      yin(711) = c10*yin(647) + yc00*yin(679) + c01*yin(677)
                                      zin(711) = c10*zin(647) + zc00*zin(679) + c01*zin(677)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  673
                                      ! i4 = i5 =  705

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  737

                                      xin(743) = c10*xin(679) + xc00*xin(711) + c01*xin(709)
                                      yin(743) = c10*yin(679) + yc00*yin(711) + c01*yin(709)
                                      zin(743) = c10*zin(679) + zc00*zin(711) + c01*zin(709)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  705
                                      ! i4 = i5 =  737

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  745

                                      xin(751) = c10*xin(711) + xc00*xin(743) + c01*xin(741)
                                      yin(751) = c10*yin(711) + yc00*yin(743) + c01*yin(741)
                                      zin(751) = c10*zin(711) + zc00*zin(743) + c01*zin(741)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  737
                                      ! i4 = i5 =  745

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  753

                                      xin(759) = c10*xin(743) + xc00*xin(751) + c01*xin(749)
                                      yin(759) = c10*yin(743) + yc00*yin(751) + c01*yin(749)
                                      zin(759) = c10*zin(743) + zc00*zin(751) + c01*zin(749)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  745
                                      ! i4 = i5 =  753

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  761

                                      xin(767) = c10*xin(751) + xc00*xin(759) + c01*xin(757)
                                      yin(767) = c10*yin(751) + yc00*yin(759) + c01*yin(757)
                                      zin(767) = c10*zin(751) + zc00*zin(759) + c01*zin(757)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  753
                                      ! i4 = i5 =  761

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  641
                                      ! i4 = i2 =  673

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  705

                                      xin(712) = c10*xin(648) + xc00*xin(680) + c01*xin(679)
                                      yin(712) = c10*yin(648) + yc00*yin(680) + c01*yin(679)
                                      zin(712) = c10*zin(648) + zc00*zin(680) + c01*zin(679)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  673
                                      ! i4 = i5 =  705

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  737

                                      xin(744) = c10*xin(680) + xc00*xin(712) + c01*xin(711)
                                      yin(744) = c10*yin(680) + yc00*yin(712) + c01*yin(711)
                                      zin(744) = c10*zin(680) + zc00*zin(712) + c01*zin(711)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  705
                                      ! i4 = i5 =  737

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  745

                                      xin(752) = c10*xin(712) + xc00*xin(744) + c01*xin(743)
                                      yin(752) = c10*yin(712) + yc00*yin(744) + c01*yin(743)
                                      zin(752) = c10*zin(712) + zc00*zin(744) + c01*zin(743)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  737
                                      ! i4 = i5 =  745

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  753

                                      xin(760) = c10*xin(744) + xc00*xin(752) + c01*xin(751)
                                      yin(760) = c10*yin(744) + yc00*yin(752) + c01*yin(751)
                                      zin(760) = c10*zin(744) + zc00*zin(752) + c01*zin(751)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  745
                                      ! i4 = i5 =  753

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  761

                                      xin(768) = c10*xin(752) + xc00*xin(760) + c01*xin(759)
                                      yin(768) = c10*yin(752) + yc00*yin(760) + c01*yin(759)
                                      zin(768) = c10*zin(752) + zc00*zin(760) + c01*zin(759)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  753
                                      ! i4 = i5 =  761

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  761

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  761

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  753

                                      xin(761) = xin(761) + dxij*xin(753)
                                      yin(761) = yin(761) + dyij*yin(753)
                                      zin(761) = zin(761) + dzij*zin(753)

                                      ! i3 = i4 =  753
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  745

                                      xin(753) = xin(753) + dxij*xin(745)
                                      yin(753) = yin(753) + dyij*yin(745)
                                      zin(753) = zin(753) + dzij*zin(745)

                                      ! i3 = i4 =  745
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  737

                                      xin(745) = xin(745) + dxij*xin(737)
                                      yin(745) = yin(745) + dyij*yin(737)
                                      zin(745) = zin(745) + dzij*zin(737)

                                      ! i3 = i4 =  737
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  761

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  753

                                      xin(761) = xin(761) + dxij*xin(753)
                                      yin(761) = yin(761) + dyij*yin(753)
                                      zin(761) = zin(761) + dzij*zin(753)

                                      ! i3 = i4 =  753
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  745

                                      xin(753) = xin(753) + dxij*xin(745)
                                      yin(753) = yin(753) + dyij*yin(745)
                                      zin(753) = zin(753) + dzij*zin(745)

                                      ! i3 = i4 =  745
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  761

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  753

                                      xin(761) = xin(761) + dxij*xin(753)
                                      yin(761) = yin(761) + dyij*yin(753)
                                      zin(761) = zin(761) + dzij*zin(753)

                                      ! i3 = i4 =  753
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  649

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  649

                                      ! do ni = 1,    3

                                      xin(649) = xin(673) + dxij*xin(641)
                                      yin(649) = yin(673) + dyij*yin(641)
                                      zin(649) = zin(673) + dzij*zin(641)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  681

                                      ! ni =    2

                                      xin(681) = xin(705) + dxij*xin(673)
                                      yin(681) = yin(705) + dyij*yin(673)
                                      zin(681) = zin(705) + dzij*zin(673)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  713

                                      ! ni =    3

                                      xin(713) = xin(737) + dxij*xin(705)
                                      yin(713) = yin(737) + dyij*yin(705)
                                      zin(713) = zin(737) + dzij*zin(705)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  745

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  657

                                      ! nj =    2

                                      ! i4 = i3 =  657

                                      ! do ni = 1,    3

                                      xin(657) = xin(681) + dxij*xin(649)
                                      yin(657) = yin(681) + dyij*yin(649)
                                      zin(657) = zin(681) + dzij*zin(649)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  689

                                      ! ni =    2

                                      xin(689) = xin(713) + dxij*xin(681)
                                      yin(689) = yin(713) + dyij*yin(681)
                                      zin(689) = zin(713) + dzij*zin(681)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  721

                                      ! ni =    3

                                      xin(721) = xin(745) + dxij*xin(713)
                                      yin(721) = yin(745) + dyij*yin(713)
                                      zin(721) = zin(745) + dzij*zin(713)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  753

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  665

                                      ! nj =    3

                                      ! i4 = i3 =  665

                                      ! do ni = 1,    3

                                      xin(665) = xin(689) + dxij*xin(657)
                                      yin(665) = yin(689) + dyij*yin(657)
                                      zin(665) = zin(689) + dzij*zin(657)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  697

                                      ! ni =    2

                                      xin(697) = xin(721) + dxij*xin(689)
                                      yin(697) = yin(721) + dyij*yin(689)
                                      zin(697) = zin(721) + dzij*zin(689)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  729

                                      ! ni =    3

                                      xin(729) = xin(753) + dxij*xin(721)
                                      yin(729) = yin(753) + dyij*yin(721)
                                      zin(729) = zin(753) + dzij*zin(721)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  761

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  673

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  763

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  755

                                      xin(763) = xin(763) + dxij*xin(755)
                                      yin(763) = yin(763) + dyij*yin(755)
                                      zin(763) = zin(763) + dzij*zin(755)

                                      ! i3 = i4 =  755
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  747

                                      xin(755) = xin(755) + dxij*xin(747)
                                      yin(755) = yin(755) + dyij*yin(747)
                                      zin(755) = zin(755) + dzij*zin(747)

                                      ! i3 = i4 =  747
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  739

                                      xin(747) = xin(747) + dxij*xin(739)
                                      yin(747) = yin(747) + dyij*yin(739)
                                      zin(747) = zin(747) + dzij*zin(739)

                                      ! i3 = i4 =  739
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  763

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  755

                                      xin(763) = xin(763) + dxij*xin(755)
                                      yin(763) = yin(763) + dyij*yin(755)
                                      zin(763) = zin(763) + dzij*zin(755)

                                      ! i3 = i4 =  755
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  747

                                      xin(755) = xin(755) + dxij*xin(747)
                                      yin(755) = yin(755) + dyij*yin(747)
                                      zin(755) = zin(755) + dzij*zin(747)

                                      ! i3 = i4 =  747
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  763

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  755

                                      xin(763) = xin(763) + dxij*xin(755)
                                      yin(763) = yin(763) + dyij*yin(755)
                                      zin(763) = zin(763) + dzij*zin(755)

                                      ! i3 = i4 =  755
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  651

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  651

                                      ! do ni = 1,    3

                                      xin(651) = xin(675) + dxij*xin(643)
                                      yin(651) = yin(675) + dyij*yin(643)
                                      zin(651) = zin(675) + dzij*zin(643)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  683

                                      ! ni =    2

                                      xin(683) = xin(707) + dxij*xin(675)
                                      yin(683) = yin(707) + dyij*yin(675)
                                      zin(683) = zin(707) + dzij*zin(675)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  715

                                      ! ni =    3

                                      xin(715) = xin(739) + dxij*xin(707)
                                      yin(715) = yin(739) + dyij*yin(707)
                                      zin(715) = zin(739) + dzij*zin(707)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  747

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  659

                                      ! nj =    2

                                      ! i4 = i3 =  659

                                      ! do ni = 1,    3

                                      xin(659) = xin(683) + dxij*xin(651)
                                      yin(659) = yin(683) + dyij*yin(651)
                                      zin(659) = zin(683) + dzij*zin(651)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  691

                                      ! ni =    2

                                      xin(691) = xin(715) + dxij*xin(683)
                                      yin(691) = yin(715) + dyij*yin(683)
                                      zin(691) = zin(715) + dzij*zin(683)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  723

                                      ! ni =    3

                                      xin(723) = xin(747) + dxij*xin(715)
                                      yin(723) = yin(747) + dyij*yin(715)
                                      zin(723) = zin(747) + dzij*zin(715)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  755

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  667

                                      ! nj =    3

                                      ! i4 = i3 =  667

                                      ! do ni = 1,    3

                                      xin(667) = xin(691) + dxij*xin(659)
                                      yin(667) = yin(691) + dyij*yin(659)
                                      zin(667) = zin(691) + dzij*zin(659)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  699

                                      ! ni =    2

                                      xin(699) = xin(723) + dxij*xin(691)
                                      yin(699) = yin(723) + dyij*yin(691)
                                      zin(699) = zin(723) + dzij*zin(691)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  731

                                      ! ni =    3

                                      xin(731) = xin(755) + dxij*xin(723)
                                      yin(731) = yin(755) + dyij*yin(723)
                                      zin(731) = zin(755) + dzij*zin(723)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  763

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  675

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  765

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  757

                                      xin(765) = xin(765) + dxij*xin(757)
                                      yin(765) = yin(765) + dyij*yin(757)
                                      zin(765) = zin(765) + dzij*zin(757)

                                      ! i3 = i4 =  757
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  749

                                      xin(757) = xin(757) + dxij*xin(749)
                                      yin(757) = yin(757) + dyij*yin(749)
                                      zin(757) = zin(757) + dzij*zin(749)

                                      ! i3 = i4 =  749
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  741

                                      xin(749) = xin(749) + dxij*xin(741)
                                      yin(749) = yin(749) + dyij*yin(741)
                                      zin(749) = zin(749) + dzij*zin(741)

                                      ! i3 = i4 =  741
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  765

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  757

                                      xin(765) = xin(765) + dxij*xin(757)
                                      yin(765) = yin(765) + dyij*yin(757)
                                      zin(765) = zin(765) + dzij*zin(757)

                                      ! i3 = i4 =  757
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  749

                                      xin(757) = xin(757) + dxij*xin(749)
                                      yin(757) = yin(757) + dyij*yin(749)
                                      zin(757) = zin(757) + dzij*zin(749)

                                      ! i3 = i4 =  749
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  765

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  757

                                      xin(765) = xin(765) + dxij*xin(757)
                                      yin(765) = yin(765) + dyij*yin(757)
                                      zin(765) = zin(765) + dzij*zin(757)

                                      ! i3 = i4 =  757
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  653

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  653

                                      ! do ni = 1,    3

                                      xin(653) = xin(677) + dxij*xin(645)
                                      yin(653) = yin(677) + dyij*yin(645)
                                      zin(653) = zin(677) + dzij*zin(645)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  685

                                      ! ni =    2

                                      xin(685) = xin(709) + dxij*xin(677)
                                      yin(685) = yin(709) + dyij*yin(677)
                                      zin(685) = zin(709) + dzij*zin(677)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  717

                                      ! ni =    3

                                      xin(717) = xin(741) + dxij*xin(709)
                                      yin(717) = yin(741) + dyij*yin(709)
                                      zin(717) = zin(741) + dzij*zin(709)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  749

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  661

                                      ! nj =    2

                                      ! i4 = i3 =  661

                                      ! do ni = 1,    3

                                      xin(661) = xin(685) + dxij*xin(653)
                                      yin(661) = yin(685) + dyij*yin(653)
                                      zin(661) = zin(685) + dzij*zin(653)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  693

                                      ! ni =    2

                                      xin(693) = xin(717) + dxij*xin(685)
                                      yin(693) = yin(717) + dyij*yin(685)
                                      zin(693) = zin(717) + dzij*zin(685)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  725

                                      ! ni =    3

                                      xin(725) = xin(749) + dxij*xin(717)
                                      yin(725) = yin(749) + dyij*yin(717)
                                      zin(725) = zin(749) + dzij*zin(717)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  757

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  669

                                      ! nj =    3

                                      ! i4 = i3 =  669

                                      ! do ni = 1,    3

                                      xin(669) = xin(693) + dxij*xin(661)
                                      yin(669) = yin(693) + dyij*yin(661)
                                      zin(669) = zin(693) + dzij*zin(661)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  701

                                      ! ni =    2

                                      xin(701) = xin(725) + dxij*xin(693)
                                      yin(701) = yin(725) + dyij*yin(693)
                                      zin(701) = zin(725) + dzij*zin(693)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  733

                                      ! ni =    3

                                      xin(733) = xin(757) + dxij*xin(725)
                                      yin(733) = yin(757) + dyij*yin(725)
                                      zin(733) = zin(757) + dzij*zin(725)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  765

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  677

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  767

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  759

                                      xin(767) = xin(767) + dxij*xin(759)
                                      yin(767) = yin(767) + dyij*yin(759)
                                      zin(767) = zin(767) + dzij*zin(759)

                                      ! i3 = i4 =  759
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  751

                                      xin(759) = xin(759) + dxij*xin(751)
                                      yin(759) = yin(759) + dyij*yin(751)
                                      zin(759) = zin(759) + dzij*zin(751)

                                      ! i3 = i4 =  751
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  743

                                      xin(751) = xin(751) + dxij*xin(743)
                                      yin(751) = yin(751) + dyij*yin(743)
                                      zin(751) = zin(751) + dzij*zin(743)

                                      ! i3 = i4 =  743
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  767

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  759

                                      xin(767) = xin(767) + dxij*xin(759)
                                      yin(767) = yin(767) + dyij*yin(759)
                                      zin(767) = zin(767) + dzij*zin(759)

                                      ! i3 = i4 =  759
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  751

                                      xin(759) = xin(759) + dxij*xin(751)
                                      yin(759) = yin(759) + dyij*yin(751)
                                      zin(759) = zin(759) + dzij*zin(751)

                                      ! i3 = i4 =  751
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  767

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  759

                                      xin(767) = xin(767) + dxij*xin(759)
                                      yin(767) = yin(767) + dyij*yin(759)
                                      zin(767) = zin(767) + dzij*zin(759)

                                      ! i3 = i4 =  759
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  655

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  655

                                      ! do ni = 1,    3

                                      xin(655) = xin(679) + dxij*xin(647)
                                      yin(655) = yin(679) + dyij*yin(647)
                                      zin(655) = zin(679) + dzij*zin(647)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  687

                                      ! ni =    2

                                      xin(687) = xin(711) + dxij*xin(679)
                                      yin(687) = yin(711) + dyij*yin(679)
                                      zin(687) = zin(711) + dzij*zin(679)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  719

                                      ! ni =    3

                                      xin(719) = xin(743) + dxij*xin(711)
                                      yin(719) = yin(743) + dyij*yin(711)
                                      zin(719) = zin(743) + dzij*zin(711)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  751

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  663

                                      ! nj =    2

                                      ! i4 = i3 =  663

                                      ! do ni = 1,    3

                                      xin(663) = xin(687) + dxij*xin(655)
                                      yin(663) = yin(687) + dyij*yin(655)
                                      zin(663) = zin(687) + dzij*zin(655)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  695

                                      ! ni =    2

                                      xin(695) = xin(719) + dxij*xin(687)
                                      yin(695) = yin(719) + dyij*yin(687)
                                      zin(695) = zin(719) + dzij*zin(687)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  727

                                      ! ni =    3

                                      xin(727) = xin(751) + dxij*xin(719)
                                      yin(727) = yin(751) + dyij*yin(719)
                                      zin(727) = zin(751) + dzij*zin(719)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  759

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  671

                                      ! nj =    3

                                      ! i4 = i3 =  671

                                      ! do ni = 1,    3

                                      xin(671) = xin(695) + dxij*xin(663)
                                      yin(671) = yin(695) + dyij*yin(663)
                                      zin(671) = zin(695) + dzij*zin(663)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  703

                                      ! ni =    2

                                      xin(703) = xin(727) + dxij*xin(695)
                                      yin(703) = yin(727) + dyij*yin(695)
                                      zin(703) = zin(727) + dzij*zin(695)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  735

                                      ! ni =    3

                                      xin(735) = xin(759) + dxij*xin(727)
                                      yin(735) = yin(759) + dyij*yin(727)
                                      zin(735) = zin(759) + dzij*zin(727)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  767

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  679

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  768

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  760

                                      xin(768) = xin(768) + dxij*xin(760)
                                      yin(768) = yin(768) + dyij*yin(760)
                                      zin(768) = zin(768) + dzij*zin(760)

                                      ! i3 = i4 =  760
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  752

                                      xin(760) = xin(760) + dxij*xin(752)
                                      yin(760) = yin(760) + dyij*yin(752)
                                      zin(760) = zin(760) + dzij*zin(752)

                                      ! i3 = i4 =  752
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  744

                                      xin(752) = xin(752) + dxij*xin(744)
                                      yin(752) = yin(752) + dyij*yin(744)
                                      zin(752) = zin(752) + dzij*zin(744)

                                      ! i3 = i4 =  744
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  768

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  760

                                      xin(768) = xin(768) + dxij*xin(760)
                                      yin(768) = yin(768) + dyij*yin(760)
                                      zin(768) = zin(768) + dzij*zin(760)

                                      ! i3 = i4 =  760
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  752

                                      xin(760) = xin(760) + dxij*xin(752)
                                      yin(760) = yin(760) + dyij*yin(752)
                                      zin(760) = zin(760) + dzij*zin(752)

                                      ! i3 = i4 =  752
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  768

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  760

                                      xin(768) = xin(768) + dxij*xin(760)
                                      yin(768) = yin(768) + dyij*yin(760)
                                      zin(768) = zin(768) + dzij*zin(760)

                                      ! i3 = i4 =  760
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  656

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  656

                                      ! do ni = 1,    3

                                      xin(656) = xin(680) + dxij*xin(648)
                                      yin(656) = yin(680) + dyij*yin(648)
                                      zin(656) = zin(680) + dzij*zin(648)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  688

                                      ! ni =    2

                                      xin(688) = xin(712) + dxij*xin(680)
                                      yin(688) = yin(712) + dyij*yin(680)
                                      zin(688) = zin(712) + dzij*zin(680)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  720

                                      ! ni =    3

                                      xin(720) = xin(744) + dxij*xin(712)
                                      yin(720) = yin(744) + dyij*yin(712)
                                      zin(720) = zin(744) + dzij*zin(712)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  752

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  664

                                      ! nj =    2

                                      ! i4 = i3 =  664

                                      ! do ni = 1,    3

                                      xin(664) = xin(688) + dxij*xin(656)
                                      yin(664) = yin(688) + dyij*yin(656)
                                      zin(664) = zin(688) + dzij*zin(656)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  696

                                      ! ni =    2

                                      xin(696) = xin(720) + dxij*xin(688)
                                      yin(696) = yin(720) + dyij*yin(688)
                                      zin(696) = zin(720) + dzij*zin(688)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  728

                                      ! ni =    3

                                      xin(728) = xin(752) + dxij*xin(720)
                                      yin(728) = yin(752) + dyij*yin(720)
                                      zin(728) = zin(752) + dzij*zin(720)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  760

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  672

                                      ! nj =    3

                                      ! i4 = i3 =  672

                                      ! do ni = 1,    3

                                      xin(672) = xin(696) + dxij*xin(664)
                                      yin(672) = yin(696) + dyij*yin(664)
                                      zin(672) = zin(696) + dzij*zin(664)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  704

                                      ! ni =    2

                                      xin(704) = xin(728) + dxij*xin(696)
                                      yin(704) = yin(728) + dyij*yin(696)
                                      zin(704) = zin(728) + dzij*zin(696)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  736

                                      ! ni =    3

                                      xin(736) = xin(760) + dxij*xin(728)
                                      yin(736) = yin(760) + dyij*yin(728)
                                      zin(736) = zin(760) + dzij*zin(728)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  768

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  680

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  641

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  642

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  642

                                      ! do nk = 1,    3

                                      xin(642) = xin(643) + dxkl*xin(641)
                                      yin(642) = yin(643) + dykl*yin(641)
                                      zin(642) = zin(643) + dzkl*zin(641)
                                      ! i4 = i4 + lang+1 =  644

                                      ! nk =    2

                                      xin(644) = xin(645) + dxkl*xin(643)
                                      yin(644) = yin(645) + dykl*yin(643)
                                      zin(644) = zin(645) + dzkl*zin(643)
                                      ! i4 = i4 + lang+1 =  646

                                      ! nk =    3

                                      xin(646) = xin(647) + dxkl*xin(645)
                                      yin(646) = yin(647) + dykl*yin(645)
                                      zin(646) = zin(647) + dzkl*zin(645)
                                      ! i4 = i4 + lang+1 =  648

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  643

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  649

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  656

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  655

                                      xin(656) = xin(656) + dxkl*xin(655)
                                      yin(656) = yin(656) + dykl*yin(655)
                                      zin(656) = zin(656) + dzkl*zin(655)

                                      ! i3 = i4 =  655
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  650

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  650

                                      ! do nk = 1,    3

                                      xin(650) = xin(651) + dxkl*xin(649)
                                      yin(650) = yin(651) + dykl*yin(649)
                                      zin(650) = zin(651) + dzkl*zin(649)
                                      ! i4 = i4 + lang+1 =  652

                                      ! nk =    2

                                      xin(652) = xin(653) + dxkl*xin(651)
                                      yin(652) = yin(653) + dykl*yin(651)
                                      zin(652) = zin(653) + dzkl*zin(651)
                                      ! i4 = i4 + lang+1 =  654

                                      ! nk =    3

                                      xin(654) = xin(655) + dxkl*xin(653)
                                      yin(654) = yin(655) + dykl*yin(653)
                                      zin(654) = zin(655) + dzkl*zin(653)
                                      ! i4 = i4 + lang+1 =  656

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  651

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  657

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  664

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  663

                                      xin(664) = xin(664) + dxkl*xin(663)
                                      yin(664) = yin(664) + dykl*yin(663)
                                      zin(664) = zin(664) + dzkl*zin(663)

                                      ! i3 = i4 =  663
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  658

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  658

                                      ! do nk = 1,    3

                                      xin(658) = xin(659) + dxkl*xin(657)
                                      yin(658) = yin(659) + dykl*yin(657)
                                      zin(658) = zin(659) + dzkl*zin(657)
                                      ! i4 = i4 + lang+1 =  660

                                      ! nk =    2

                                      xin(660) = xin(661) + dxkl*xin(659)
                                      yin(660) = yin(661) + dykl*yin(659)
                                      zin(660) = zin(661) + dzkl*zin(659)
                                      ! i4 = i4 + lang+1 =  662

                                      ! nk =    3

                                      xin(662) = xin(663) + dxkl*xin(661)
                                      yin(662) = yin(663) + dykl*yin(661)
                                      zin(662) = zin(663) + dzkl*zin(661)
                                      ! i4 = i4 + lang+1 =  664

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  659

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  665

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  672

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  671

                                      xin(672) = xin(672) + dxkl*xin(671)
                                      yin(672) = yin(672) + dykl*yin(671)
                                      zin(672) = zin(672) + dzkl*zin(671)

                                      ! i3 = i4 =  671
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  666

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  666

                                      ! do nk = 1,    3

                                      xin(666) = xin(667) + dxkl*xin(665)
                                      yin(666) = yin(667) + dykl*yin(665)
                                      zin(666) = zin(667) + dzkl*zin(665)
                                      ! i4 = i4 + lang+1 =  668

                                      ! nk =    2

                                      xin(668) = xin(669) + dxkl*xin(667)
                                      yin(668) = yin(669) + dykl*yin(667)
                                      zin(668) = zin(669) + dzkl*zin(667)
                                      ! i4 = i4 + lang+1 =  670

                                      ! nk =    3

                                      xin(670) = xin(671) + dxkl*xin(669)
                                      yin(670) = yin(671) + dykl*yin(669)
                                      zin(670) = zin(671) + dzkl*zin(669)
                                      ! i4 = i4 + lang+1 =  672

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  667

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  673

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  673

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  680

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  679

                                      xin(680) = xin(680) + dxkl*xin(679)
                                      yin(680) = yin(680) + dykl*yin(679)
                                      zin(680) = zin(680) + dzkl*zin(679)

                                      ! i3 = i4 =  679
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  674

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  674

                                      ! do nk = 1,    3

                                      xin(674) = xin(675) + dxkl*xin(673)
                                      yin(674) = yin(675) + dykl*yin(673)
                                      zin(674) = zin(675) + dzkl*zin(673)
                                      ! i4 = i4 + lang+1 =  676

                                      ! nk =    2

                                      xin(676) = xin(677) + dxkl*xin(675)
                                      yin(676) = yin(677) + dykl*yin(675)
                                      zin(676) = zin(677) + dzkl*zin(675)
                                      ! i4 = i4 + lang+1 =  678

                                      ! nk =    3

                                      xin(678) = xin(679) + dxkl*xin(677)
                                      yin(678) = yin(679) + dykl*yin(677)
                                      zin(678) = zin(679) + dzkl*zin(677)
                                      ! i4 = i4 + lang+1 =  680

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  675

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  681

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  688

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  687

                                      xin(688) = xin(688) + dxkl*xin(687)
                                      yin(688) = yin(688) + dykl*yin(687)
                                      zin(688) = zin(688) + dzkl*zin(687)

                                      ! i3 = i4 =  687
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  682

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  682

                                      ! do nk = 1,    3

                                      xin(682) = xin(683) + dxkl*xin(681)
                                      yin(682) = yin(683) + dykl*yin(681)
                                      zin(682) = zin(683) + dzkl*zin(681)
                                      ! i4 = i4 + lang+1 =  684

                                      ! nk =    2

                                      xin(684) = xin(685) + dxkl*xin(683)
                                      yin(684) = yin(685) + dykl*yin(683)
                                      zin(684) = zin(685) + dzkl*zin(683)
                                      ! i4 = i4 + lang+1 =  686

                                      ! nk =    3

                                      xin(686) = xin(687) + dxkl*xin(685)
                                      yin(686) = yin(687) + dykl*yin(685)
                                      zin(686) = zin(687) + dzkl*zin(685)
                                      ! i4 = i4 + lang+1 =  688

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  683

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  689

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  696

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  695

                                      xin(696) = xin(696) + dxkl*xin(695)
                                      yin(696) = yin(696) + dykl*yin(695)
                                      zin(696) = zin(696) + dzkl*zin(695)

                                      ! i3 = i4 =  695
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  690

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  690

                                      ! do nk = 1,    3

                                      xin(690) = xin(691) + dxkl*xin(689)
                                      yin(690) = yin(691) + dykl*yin(689)
                                      zin(690) = zin(691) + dzkl*zin(689)
                                      ! i4 = i4 + lang+1 =  692

                                      ! nk =    2

                                      xin(692) = xin(693) + dxkl*xin(691)
                                      yin(692) = yin(693) + dykl*yin(691)
                                      zin(692) = zin(693) + dzkl*zin(691)
                                      ! i4 = i4 + lang+1 =  694

                                      ! nk =    3

                                      xin(694) = xin(695) + dxkl*xin(693)
                                      yin(694) = yin(695) + dykl*yin(693)
                                      zin(694) = zin(695) + dzkl*zin(693)
                                      ! i4 = i4 + lang+1 =  696

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  691

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  697

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  704

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  703

                                      xin(704) = xin(704) + dxkl*xin(703)
                                      yin(704) = yin(704) + dykl*yin(703)
                                      zin(704) = zin(704) + dzkl*zin(703)

                                      ! i3 = i4 =  703
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  698

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  698

                                      ! do nk = 1,    3

                                      xin(698) = xin(699) + dxkl*xin(697)
                                      yin(698) = yin(699) + dykl*yin(697)
                                      zin(698) = zin(699) + dzkl*zin(697)
                                      ! i4 = i4 + lang+1 =  700

                                      ! nk =    2

                                      xin(700) = xin(701) + dxkl*xin(699)
                                      yin(700) = yin(701) + dykl*yin(699)
                                      zin(700) = zin(701) + dzkl*zin(699)
                                      ! i4 = i4 + lang+1 =  702

                                      ! nk =    3

                                      xin(702) = xin(703) + dxkl*xin(701)
                                      yin(702) = yin(703) + dykl*yin(701)
                                      zin(702) = zin(703) + dzkl*zin(701)
                                      ! i4 = i4 + lang+1 =  704

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  699

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  705

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  705

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  712

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  711

                                      xin(712) = xin(712) + dxkl*xin(711)
                                      yin(712) = yin(712) + dykl*yin(711)
                                      zin(712) = zin(712) + dzkl*zin(711)

                                      ! i3 = i4 =  711
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  706

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  706

                                      ! do nk = 1,    3

                                      xin(706) = xin(707) + dxkl*xin(705)
                                      yin(706) = yin(707) + dykl*yin(705)
                                      zin(706) = zin(707) + dzkl*zin(705)
                                      ! i4 = i4 + lang+1 =  708

                                      ! nk =    2

                                      xin(708) = xin(709) + dxkl*xin(707)
                                      yin(708) = yin(709) + dykl*yin(707)
                                      zin(708) = zin(709) + dzkl*zin(707)
                                      ! i4 = i4 + lang+1 =  710

                                      ! nk =    3

                                      xin(710) = xin(711) + dxkl*xin(709)
                                      yin(710) = yin(711) + dykl*yin(709)
                                      zin(710) = zin(711) + dzkl*zin(709)
                                      ! i4 = i4 + lang+1 =  712

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  707

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  713

                                      ! nj = nj + 1 =    1

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

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  714

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  714

                                      ! do nk = 1,    3

                                      xin(714) = xin(715) + dxkl*xin(713)
                                      yin(714) = yin(715) + dykl*yin(713)
                                      zin(714) = zin(715) + dzkl*zin(713)
                                      ! i4 = i4 + lang+1 =  716

                                      ! nk =    2

                                      xin(716) = xin(717) + dxkl*xin(715)
                                      yin(716) = yin(717) + dykl*yin(715)
                                      zin(716) = zin(717) + dzkl*zin(715)
                                      ! i4 = i4 + lang+1 =  718

                                      ! nk =    3

                                      xin(718) = xin(719) + dxkl*xin(717)
                                      yin(718) = yin(719) + dykl*yin(717)
                                      zin(718) = zin(719) + dzkl*zin(717)
                                      ! i4 = i4 + lang+1 =  720

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  715

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  721

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  728

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  727

                                      xin(728) = xin(728) + dxkl*xin(727)
                                      yin(728) = yin(728) + dykl*yin(727)
                                      zin(728) = zin(728) + dzkl*zin(727)

                                      ! i3 = i4 =  727
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  722

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  722

                                      ! do nk = 1,    3

                                      xin(722) = xin(723) + dxkl*xin(721)
                                      yin(722) = yin(723) + dykl*yin(721)
                                      zin(722) = zin(723) + dzkl*zin(721)
                                      ! i4 = i4 + lang+1 =  724

                                      ! nk =    2

                                      xin(724) = xin(725) + dxkl*xin(723)
                                      yin(724) = yin(725) + dykl*yin(723)
                                      zin(724) = zin(725) + dzkl*zin(723)
                                      ! i4 = i4 + lang+1 =  726

                                      ! nk =    3

                                      xin(726) = xin(727) + dxkl*xin(725)
                                      yin(726) = yin(727) + dykl*yin(725)
                                      zin(726) = zin(727) + dzkl*zin(725)
                                      ! i4 = i4 + lang+1 =  728

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  723

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  729

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  736

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  735

                                      xin(736) = xin(736) + dxkl*xin(735)
                                      yin(736) = yin(736) + dykl*yin(735)
                                      zin(736) = zin(736) + dzkl*zin(735)

                                      ! i3 = i4 =  735
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  730

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  730

                                      ! do nk = 1,    3

                                      xin(730) = xin(731) + dxkl*xin(729)
                                      yin(730) = yin(731) + dykl*yin(729)
                                      zin(730) = zin(731) + dzkl*zin(729)
                                      ! i4 = i4 + lang+1 =  732

                                      ! nk =    2

                                      xin(732) = xin(733) + dxkl*xin(731)
                                      yin(732) = yin(733) + dykl*yin(731)
                                      zin(732) = zin(733) + dzkl*zin(731)
                                      ! i4 = i4 + lang+1 =  734

                                      ! nk =    3

                                      xin(734) = xin(735) + dxkl*xin(733)
                                      yin(734) = yin(735) + dykl*yin(733)
                                      zin(734) = zin(735) + dzkl*zin(733)
                                      ! i4 = i4 + lang+1 =  736

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  731

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  737

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  737

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  744

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  743

                                      xin(744) = xin(744) + dxkl*xin(743)
                                      yin(744) = yin(744) + dykl*yin(743)
                                      zin(744) = zin(744) + dzkl*zin(743)

                                      ! i3 = i4 =  743
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  738

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  738

                                      ! do nk = 1,    3

                                      xin(738) = xin(739) + dxkl*xin(737)
                                      yin(738) = yin(739) + dykl*yin(737)
                                      zin(738) = zin(739) + dzkl*zin(737)
                                      ! i4 = i4 + lang+1 =  740

                                      ! nk =    2

                                      xin(740) = xin(741) + dxkl*xin(739)
                                      yin(740) = yin(741) + dykl*yin(739)
                                      zin(740) = zin(741) + dzkl*zin(739)
                                      ! i4 = i4 + lang+1 =  742

                                      ! nk =    3

                                      xin(742) = xin(743) + dxkl*xin(741)
                                      yin(742) = yin(743) + dykl*yin(741)
                                      zin(742) = zin(743) + dzkl*zin(741)
                                      ! i4 = i4 + lang+1 =  744

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  739

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  745

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  752

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  751

                                      xin(752) = xin(752) + dxkl*xin(751)
                                      yin(752) = yin(752) + dykl*yin(751)
                                      zin(752) = zin(752) + dzkl*zin(751)

                                      ! i3 = i4 =  751
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  746

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  746

                                      ! do nk = 1,    3

                                      xin(746) = xin(747) + dxkl*xin(745)
                                      yin(746) = yin(747) + dykl*yin(745)
                                      zin(746) = zin(747) + dzkl*zin(745)
                                      ! i4 = i4 + lang+1 =  748

                                      ! nk =    2

                                      xin(748) = xin(749) + dxkl*xin(747)
                                      yin(748) = yin(749) + dykl*yin(747)
                                      zin(748) = zin(749) + dzkl*zin(747)
                                      ! i4 = i4 + lang+1 =  750

                                      ! nk =    3

                                      xin(750) = xin(751) + dxkl*xin(749)
                                      yin(750) = yin(751) + dykl*yin(749)
                                      zin(750) = zin(751) + dzkl*zin(749)
                                      ! i4 = i4 + lang+1 =  752

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  747

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  753

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  760

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  759

                                      xin(760) = xin(760) + dxkl*xin(759)
                                      yin(760) = yin(760) + dykl*yin(759)
                                      zin(760) = zin(760) + dzkl*zin(759)

                                      ! i3 = i4 =  759
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  754

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  754

                                      ! do nk = 1,    3

                                      xin(754) = xin(755) + dxkl*xin(753)
                                      yin(754) = yin(755) + dykl*yin(753)
                                      zin(754) = zin(755) + dzkl*zin(753)
                                      ! i4 = i4 + lang+1 =  756

                                      ! nk =    2

                                      xin(756) = xin(757) + dxkl*xin(755)
                                      yin(756) = yin(757) + dykl*yin(755)
                                      zin(756) = zin(757) + dzkl*zin(755)
                                      ! i4 = i4 + lang+1 =  758

                                      ! nk =    3

                                      xin(758) = xin(759) + dxkl*xin(757)
                                      yin(758) = yin(759) + dykl*yin(757)
                                      zin(758) = zin(759) + dzkl*zin(757)
                                      ! i4 = i4 + lang+1 =  760

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  755

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  761

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  768

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  767

                                      xin(768) = xin(768) + dxkl*xin(767)
                                      yin(768) = yin(768) + dykl*yin(767)
                                      zin(768) = zin(768) + dzkl*zin(767)

                                      ! i3 = i4 =  767
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  762

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  762

                                      ! do nk = 1,    3

                                      xin(762) = xin(763) + dxkl*xin(761)
                                      yin(762) = yin(763) + dykl*yin(761)
                                      zin(762) = zin(763) + dzkl*zin(761)
                                      ! i4 = i4 + lang+1 =  764

                                      ! nk =    2

                                      xin(764) = xin(765) + dxkl*xin(763)
                                      yin(764) = yin(765) + dykl*yin(763)
                                      zin(764) = zin(765) + dzkl*zin(763)
                                      ! i4 = i4 + lang+1 =  766

                                      ! nk =    3

                                      xin(766) = xin(767) + dxkl*xin(765)
                                      yin(766) = yin(767) + dykl*yin(765)
                                      zin(766) = zin(767) + dzkl*zin(765)
                                      ! i4 = i4 + lang+1 =  768

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  763

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  769

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  769

                                      ! end do

                                      ! *** Now root =    7

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  768

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 3000! loop over all integrals

                                        l = n - 30*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d13ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 128)*yin(my + 128)*zin(mz + 128) & ! root  2
                                                        + xin(mx + 256)*yin(my + 256)*zin(mz + 256) & ! root  3
                                                        + xin(mx + 384)*yin(my + 384)*zin(mz + 384) & ! root  4
                                                        + xin(mx + 512)*yin(my + 512)*zin(mz + 512) & ! root  5
                                                        + xin(mx + 640)*yin(my + 640)*zin(mz + 640)) ! root  6

                                        j = int(n/30) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxj2 = 10
                                  iandj = ish .eq. jsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    if (iandj) maxj2 = i

                                    ii1 = i + loci
                                    ip = (i - 1)*300 ! Stride between functions in i

                                    do j = 1, maxj2

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*30 + ip ! Add stride between functions in j

                                      do k = 1, 10 ! # of cartesians in k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*3 + ijp ! Add stride between functions in k

                                        do l = 1, 3! # of cartesians in l

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
                              deallocate (n13ket)
                              deallocate (xint13ket)

                              end subroutine int3331
                              end submodule
