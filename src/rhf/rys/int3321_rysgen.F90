! The total angular momentum of this class is:           9
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3321_impl
contains
  module subroutine int3321(ff_pair, pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n12ket(:)
    real(dp), allocatable :: xint33bra(:), xint12ket(:)
    integer(kind=int64) :: nffbra, npdket
    real(dp) :: scutffbra, scutpdket, test
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
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(480), yin(480), zin(480)
    real(dp) :: eri_value(1800)
    real(dp) :: d33bra(100), d12ket(18)
    integer(kind=int64) :: ix(10), jx(10), kx(6), lx(3)
    integer(kind=int64) :: iy(10), jy(10), ky(6), ly(3)
    integer(kind=int64) :: iz(10), jz(10), kz(6), lz(3)
    integer(kind=int64) :: in(7), in1(7), kn(4)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(18), kly(18), klz(18)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 25
    in1(3) = 49
    in1(4) = 73
    in1(5) = 79
    in1(6) = 85
    in1(7) = 91

    kn(1) = 0
    kn(2) = 2
    kn(3) = 4
    kn(4) = 5

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 1
    lx(2) = 0
    lx(3) = 0

    kx(1) = 4
    kx(2) = 0
    kx(3) = 0
    kx(4) = 2
    kx(5) = 2
    kx(6) = 0

    jx(1) = 18
    jx(2) = 0
    jx(3) = 0
    jx(4) = 12
    jx(5) = 12
    jx(6) = 6
    jx(7) = 0
    jx(8) = 6
    jx(9) = 0
    jx(10) = 6

    ix(1) = 73
    ix(2) = 1
    ix(3) = 1
    ix(4) = 49
    ix(5) = 49
    ix(6) = 25
    ix(7) = 1
    ix(8) = 25
    ix(9) = 1
    ix(10) = 25

    ! y-arrays

    ly(1) = 0
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 4
    ky(3) = 0
    ky(4) = 2
    ky(5) = 0
    ky(6) = 2

    jy(1) = 0
    jy(2) = 18
    jy(3) = 0
    jy(4) = 6
    jy(5) = 0
    jy(6) = 12
    jy(7) = 12
    jy(8) = 0
    jy(9) = 6
    jy(10) = 6

    iy(1) = 1
    iy(2) = 73
    iy(3) = 1
    iy(4) = 25
    iy(5) = 1
    iy(6) = 49
    iy(7) = 49
    iy(8) = 1
    iy(9) = 25
    iy(10) = 25

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 4
    kz(4) = 0
    kz(5) = 2
    kz(6) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 18
    jz(4) = 0
    jz(5) = 6
    jz(6) = 0
    jz(7) = 6
    jz(8) = 12
    jz(9) = 12
    jz(10) = 6

    iz(1) = 1
    iz(2) = 1
    iz(3) = 73
    iz(4) = 1
    iz(5) = 25
    iz(6) = 1
    iz(7) = 25
    iz(8) = 49
    iz(9) = 49
    iz(10) = 25

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 91
    ijx(2) = 73
    ijx(3) = 73
    ijx(4) = 85
    ijx(5) = 85
    ijx(6) = 79
    ijx(7) = 73
    ijx(8) = 79
    ijx(9) = 73
    ijx(10) = 79
    ijx(11) = 19
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 13
    ijx(15) = 13
    ijx(16) = 7
    ijx(17) = 1
    ijx(18) = 7
    ijx(19) = 1
    ijx(20) = 7
    ijx(21) = 19
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 13
    ijx(25) = 13
    ijx(26) = 7
    ijx(27) = 1
    ijx(28) = 7
    ijx(29) = 1
    ijx(30) = 7
    ijx(31) = 67
    ijx(32) = 49
    ijx(33) = 49
    ijx(34) = 61
    ijx(35) = 61
    ijx(36) = 55
    ijx(37) = 49
    ijx(38) = 55
    ijx(39) = 49
    ijx(40) = 55
    ijx(41) = 67
    ijx(42) = 49
    ijx(43) = 49
    ijx(44) = 61
    ijx(45) = 61
    ijx(46) = 55
    ijx(47) = 49
    ijx(48) = 55
    ijx(49) = 49
    ijx(50) = 55
    ijx(51) = 43
    ijx(52) = 25
    ijx(53) = 25
    ijx(54) = 37
    ijx(55) = 37
    ijx(56) = 31
    ijx(57) = 25
    ijx(58) = 31
    ijx(59) = 25
    ijx(60) = 31
    ijx(61) = 19
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 13
    ijx(65) = 13
    ijx(66) = 7
    ijx(67) = 1
    ijx(68) = 7
    ijx(69) = 1
    ijx(70) = 7
    ijx(71) = 43
    ijx(72) = 25
    ijx(73) = 25
    ijx(74) = 37
    ijx(75) = 37
    ijx(76) = 31
    ijx(77) = 25
    ijx(78) = 31
    ijx(79) = 25
    ijx(80) = 31
    ijx(81) = 19
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 13
    ijx(85) = 13
    ijx(86) = 7
    ijx(87) = 1
    ijx(88) = 7
    ijx(89) = 1
    ijx(90) = 7
    ijx(91) = 43
    ijx(92) = 25
    ijx(93) = 25
    ijx(94) = 37
    ijx(95) = 37
    ijx(96) = 31
    ijx(97) = 25
    ijx(98) = 31
    ijx(99) = 25
    ijx(100) = 31

    ijy(1) = 1
    ijy(2) = 19
    ijy(3) = 1
    ijy(4) = 7
    ijy(5) = 1
    ijy(6) = 13
    ijy(7) = 13
    ijy(8) = 1
    ijy(9) = 7
    ijy(10) = 7
    ijy(11) = 73
    ijy(12) = 91
    ijy(13) = 73
    ijy(14) = 79
    ijy(15) = 73
    ijy(16) = 85
    ijy(17) = 85
    ijy(18) = 73
    ijy(19) = 79
    ijy(20) = 79
    ijy(21) = 1
    ijy(22) = 19
    ijy(23) = 1
    ijy(24) = 7
    ijy(25) = 1
    ijy(26) = 13
    ijy(27) = 13
    ijy(28) = 1
    ijy(29) = 7
    ijy(30) = 7
    ijy(31) = 25
    ijy(32) = 43
    ijy(33) = 25
    ijy(34) = 31
    ijy(35) = 25
    ijy(36) = 37
    ijy(37) = 37
    ijy(38) = 25
    ijy(39) = 31
    ijy(40) = 31
    ijy(41) = 1
    ijy(42) = 19
    ijy(43) = 1
    ijy(44) = 7
    ijy(45) = 1
    ijy(46) = 13
    ijy(47) = 13
    ijy(48) = 1
    ijy(49) = 7
    ijy(50) = 7
    ijy(51) = 49
    ijy(52) = 67
    ijy(53) = 49
    ijy(54) = 55
    ijy(55) = 49
    ijy(56) = 61
    ijy(57) = 61
    ijy(58) = 49
    ijy(59) = 55
    ijy(60) = 55
    ijy(61) = 49
    ijy(62) = 67
    ijy(63) = 49
    ijy(64) = 55
    ijy(65) = 49
    ijy(66) = 61
    ijy(67) = 61
    ijy(68) = 49
    ijy(69) = 55
    ijy(70) = 55
    ijy(71) = 1
    ijy(72) = 19
    ijy(73) = 1
    ijy(74) = 7
    ijy(75) = 1
    ijy(76) = 13
    ijy(77) = 13
    ijy(78) = 1
    ijy(79) = 7
    ijy(80) = 7
    ijy(81) = 25
    ijy(82) = 43
    ijy(83) = 25
    ijy(84) = 31
    ijy(85) = 25
    ijy(86) = 37
    ijy(87) = 37
    ijy(88) = 25
    ijy(89) = 31
    ijy(90) = 31
    ijy(91) = 25
    ijy(92) = 43
    ijy(93) = 25
    ijy(94) = 31
    ijy(95) = 25
    ijy(96) = 37
    ijy(97) = 37
    ijy(98) = 25
    ijy(99) = 31
    ijy(100) = 31

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 19
    ijz(4) = 1
    ijz(5) = 7
    ijz(6) = 1
    ijz(7) = 7
    ijz(8) = 13
    ijz(9) = 13
    ijz(10) = 7
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 19
    ijz(14) = 1
    ijz(15) = 7
    ijz(16) = 1
    ijz(17) = 7
    ijz(18) = 13
    ijz(19) = 13
    ijz(20) = 7
    ijz(21) = 73
    ijz(22) = 73
    ijz(23) = 91
    ijz(24) = 73
    ijz(25) = 79
    ijz(26) = 73
    ijz(27) = 79
    ijz(28) = 85
    ijz(29) = 85
    ijz(30) = 79
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 19
    ijz(34) = 1
    ijz(35) = 7
    ijz(36) = 1
    ijz(37) = 7
    ijz(38) = 13
    ijz(39) = 13
    ijz(40) = 7
    ijz(41) = 25
    ijz(42) = 25
    ijz(43) = 43
    ijz(44) = 25
    ijz(45) = 31
    ijz(46) = 25
    ijz(47) = 31
    ijz(48) = 37
    ijz(49) = 37
    ijz(50) = 31
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 19
    ijz(54) = 1
    ijz(55) = 7
    ijz(56) = 1
    ijz(57) = 7
    ijz(58) = 13
    ijz(59) = 13
    ijz(60) = 7
    ijz(61) = 25
    ijz(62) = 25
    ijz(63) = 43
    ijz(64) = 25
    ijz(65) = 31
    ijz(66) = 25
    ijz(67) = 31
    ijz(68) = 37
    ijz(69) = 37
    ijz(70) = 31
    ijz(71) = 49
    ijz(72) = 49
    ijz(73) = 67
    ijz(74) = 49
    ijz(75) = 55
    ijz(76) = 49
    ijz(77) = 55
    ijz(78) = 61
    ijz(79) = 61
    ijz(80) = 55
    ijz(81) = 49
    ijz(82) = 49
    ijz(83) = 67
    ijz(84) = 49
    ijz(85) = 55
    ijz(86) = 49
    ijz(87) = 55
    ijz(88) = 61
    ijz(89) = 61
    ijz(90) = 55
    ijz(91) = 25
    ijz(92) = 25
    ijz(93) = 43
    ijz(94) = 25
    ijz(95) = 31
    ijz(96) = 25
    ijz(97) = 31
    ijz(98) = 37
    ijz(99) = 37
    ijz(100) = 31

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 5
    klx(2) = 4
    klx(3) = 4
    klx(4) = 1
    klx(5) = 0
    klx(6) = 0
    klx(7) = 1
    klx(8) = 0
    klx(9) = 0
    klx(10) = 3
    klx(11) = 2
    klx(12) = 2
    klx(13) = 3
    klx(14) = 2
    klx(15) = 2
    klx(16) = 1
    klx(17) = 0
    klx(18) = 0

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0
    kly(4) = 4
    kly(5) = 5
    kly(6) = 4
    kly(7) = 0
    kly(8) = 1
    kly(9) = 0
    kly(10) = 2
    kly(11) = 3
    kly(12) = 2
    kly(13) = 0
    kly(14) = 1
    kly(15) = 0
    kly(16) = 2
    kly(17) = 3
    kly(18) = 2

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1
    klz(4) = 0
    klz(5) = 0
    klz(6) = 1
    klz(7) = 4
    klz(8) = 4
    klz(9) = 5
    klz(10) = 0
    klz(11) = 0
    klz(12) = 1
    klz(13) = 2
    klz(14) = 2
    klz(15) = 3
    klz(16) = 2
    klz(17) = 2
    klz(18) = 3

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n12ket(res%n_p_shl*res%n_d_shl))
    allocate (xint12ket(res%n_p_shl*res%n_d_shl))

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

    scutpdket = cutoff_schwarz/maxval(pd_pair%xints)
    npdket = 0
    do ij = 1, res%n_p_shl*res%n_d_shl
      if (pd_pair%xints(ij) .ge. scutpdket) then
        npdket = npdket + 1
        xint12ket(npdket) = pd_pair%xints(ij)
        n12ket(npdket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*npdket) .le. nchunksize_int64) nchunksize_int64 = nffbra*npdket
    ntile = int(nffbra*npdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*npdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint12ket, n12ket, ff_pair, pd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d12ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
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

              test = xint33bra(ij_tmp)*xint12ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n12ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_d_shl(ksh_tmp)
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

                  t_expon_cd = pd_pair%t_expon_ab(pd_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = pd_pair%expon_b(pd_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = pd_pair%expon_a(pd_pair%pair_loc(kl) + ket_loop)
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

                  d12ket(1) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(2) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(3) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(4) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(5) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(6) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(7) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(8) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(9) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d12ket(10) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(11) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(12) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(13) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(14) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(15) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(16) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(17) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d12ket(18) = pd_pair%d_coeff_alt(pd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3

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

                    if (xx .ge. 55.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.1175813202117781D+00
                      rts(2) = factr*0.1074562012436902D+01
                      rts(3) = factr*0.3085937443717543D+01
                      rts(4) = factr*0.6414729733662035D+01
                      rts(5) = factr*0.1180718948997174D+02

                      wts(1) = factw*0.6108626337353259D+00
                      wts(2) = factw*0.2401386110823149D+00
                      wts(3) = factw*0.3387439445548109D-01
                      wts(4) = factw*0.1343645746781227D-02
                      wts(5) = factw*0.7640432855232614D-05

                    else ! "regular" evaluation

                      rgrid(1) = 0.1314956323727580D-05
                      rgrid(2) = 0.3638644488856184D-04
                      rgrid(3) = 0.2184836363610052D-03
                      rgrid(4) = 0.7467882061060832D-03
                      rgrid(5) = 0.1898594940807536D-02
                      rgrid(6) = 0.4018347564842567D-02
                      rgrid(7) = 0.7503899366138534D-02
                      rgrid(8) = 0.1279045049262753D-01
                      rgrid(9) = 0.2033269213059855D-01
                      rgrid(10) = 0.3058574876834499D-01
                      rgrid(11) = 0.4398555195277359D-01
                      rgrid(12) = 0.6092930105499368D-01
                      rgrid(13) = 0.8175666785532240D-01
                      rgrid(14) = 0.1067323821680571D+00
                      rgrid(15) = 0.1360307958356441D+00
                      rgrid(16) = 0.1697229634514230D+00
                      rgrid(17) = 0.2077667019402563D+00
                      rgrid(18) = 0.2499999999999998D+00
                      rgrid(19) = 0.2961380452159158D+00
                      rgrid(20) = 0.3457740246174122D+00
                      rgrid(21) = 0.3983837370449402D+00
                      rgrid(22) = 0.4533339365988709D+00
                      rgrid(23) = 0.5098942093731368D+00
                      rgrid(24) = 0.5672520742964827D+00
                      rgrid(25) = 0.6245308967025374D+00
                      rgrid(26) = 0.6808101134342349D+00
                      rgrid(27) = 0.7351471936872270D+00
                      rgrid(28) = 0.7866007027795398D+00
                      rgrid(29) = 0.8342537984583634D+00
                      rgrid(30) = 0.8772374725900651D+00
                      rgrid(31) = 0.9147528563001249D+00
                      rgrid(32) = 0.9460919364139329D+00
                      rgrid(33) = 0.9706560996755904D+00
                      rgrid(34) = 0.9879721508887392D+00
                      rgrid(35) = 0.9977078840559234D+00

                      wgrid(1) = 0.2941716710221880D-02*exp(-xx*0.1314956323727580D-05)
                      wgrid(2) = 0.6825414174180503D-02*exp(-xx*0.3638644488856184D-04)
                      wgrid(3) = 0.1066148995574242D-01*exp(-xx*0.2184836363610052D-03)
                      wgrid(4) = 0.1441463005444643D-01*exp(-xx*0.7467882061060832D-03)
                      wgrid(5) = 0.1805505793173160D-01*exp(-xx*0.1898594940807536D-02)
                      wgrid(6) = 0.2155421116308512D-01*exp(-xx*0.4018347564842567D-02)
                      wgrid(7) = 0.2488468520067677D-01*exp(-xx*0.7503899366138534D-02)
                      wgrid(8) = 0.2802040810618503D-01*exp(-xx*0.1279045049262753D-01)
                      wgrid(9) = 0.3093683598304008D-01*exp(-xx*0.2033269213059855D-01)
                      wgrid(10) = 0.3361114263454359D-01*exp(-xx*0.3058574876834499D-01)
                      wgrid(11) = 0.3602239738628011D-01*exp(-xx*0.4398555195277359D-01)
                      wgrid(12) = 0.3815172857772094D-01*exp(-xx*0.6092930105499368D-01)
                      wgrid(13) = 0.3998247112116213D-01*exp(-xx*0.8175666785532240D-01)
                      wgrid(14) = 0.4150029686442836D-01*exp(-xx*0.1067323821680571D+00)
                      wgrid(15) = 0.4269332669604951D-01*exp(-xx*0.1360307958356441D+00)
                      wgrid(16) = 0.4355222349859141D-01*exp(-xx*0.1697229634514230D+00)
                      wgrid(17) = 0.4407026521513794D-01*exp(-xx*0.2077667019402563D+00)
                      wgrid(18) = 0.4424339745355207D-01*exp(-xx*0.2499999999999998D+00)
                      wgrid(19) = 0.4407026521513754D-01*exp(-xx*0.2961380452159158D+00)
                      wgrid(20) = 0.4355222349859152D-01*exp(-xx*0.3457740246174122D+00)
                      wgrid(21) = 0.4269332669604975D-01*exp(-xx*0.3983837370449402D+00)
                      wgrid(22) = 0.4150029686442849D-01*exp(-xx*0.4533339365988709D+00)
                      wgrid(23) = 0.3998247112116278D-01*exp(-xx*0.5098942093731368D+00)
                      wgrid(24) = 0.3815172857772100D-01*exp(-xx*0.5672520742964827D+00)
                      wgrid(25) = 0.3602239738628051D-01*exp(-xx*0.6245308967025374D+00)
                      wgrid(26) = 0.3361114263454346D-01*exp(-xx*0.6808101134342349D+00)
                      wgrid(27) = 0.3093683598304035D-01*exp(-xx*0.7351471936872270D+00)
                      wgrid(28) = 0.2802040810618515D-01*exp(-xx*0.7866007027795398D+00)
                      wgrid(29) = 0.2488468520067681D-01*exp(-xx*0.8342537984583634D+00)
                      wgrid(30) = 0.2155421116308484D-01*exp(-xx*0.8772374725900651D+00)
                      wgrid(31) = 0.1805505793173228D-01*exp(-xx*0.9147528563001249D+00)
                      wgrid(32) = 0.1441463005444692D-01*exp(-xx*0.9460919364139329D+00)
                      wgrid(33) = 0.1066148995574185D-01*exp(-xx*0.9706560996755904D+00)
                      wgrid(34) = 0.6825414174181070D-02*exp(-xx*0.9879721508887392D+00)
                      wgrid(35) = 0.2941716710221409D-02*exp(-xx*0.9977078840559234D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 35
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 35
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 4

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 35

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
                        wrk(5) = 0.0D+00
                        do 100 kk = 2, 5

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 5

                            jj = 0

105                         do 110 m = l, 5
                              if (m .eq. 5) go to 120
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

                                do 300 ii = 2, 5

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 5
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

                                    do 310 kk = 1, 5
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 5
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

                                      ! i2 = in(2) =   25
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(25) = xc00
                                      yin(25) = yc00
                                      zin(25) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   27
                                      ! i2 =   25

                                      xin(27) = xcp00*xin(25) + cp10
                                      yin(27) = ycp00*yin(25) + cp10
                                      zin(27) = zcp00*zin(25) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =    1
                                      ! i4 =   25

                                      xin(49) = c10*xin(1) + xc00*xin(25)
                                      yin(49) = c10*yin(1) + yc00*yin(25)
                                      zin(49) = c10*zin(1) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   51
                                      ! i5 =   49
                                      ! i4 =   25

                                      xin(51) = xcp00*xin(49) + cp10*xin(25)
                                      yin(51) = ycp00*yin(49) + cp10*yin(25)
                                      zin(51) = zcp00*zin(49) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   73
                                      ! i3 =   25
                                      ! i4 =   49

                                      xin(73) = c10*xin(25) + xc00*xin(49)
                                      yin(73) = c10*yin(25) + yc00*yin(49)
                                      zin(73) = c10*zin(25) + zc00*zin(49)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   75
                                      ! i5 =   73
                                      ! i4 =   49

                                      xin(75) = xcp00*xin(73) + cp10*xin(49)
                                      yin(75) = ycp00*yin(73) + cp10*yin(49)
                                      zin(75) = zcp00*zin(73) + cp10*zin(49)

                                      ! ------------------

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   73

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   79
                                      ! i3 =   49
                                      ! i4 =   73

                                      xin(79) = c10*xin(49) + xc00*xin(73)
                                      yin(79) = c10*yin(49) + yc00*yin(73)
                                      zin(79) = c10*zin(49) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   81
                                      ! i5 =   79
                                      ! i4 =   73

                                      xin(81) = xcp00*xin(79) + cp10*xin(73)
                                      yin(81) = ycp00*yin(79) + cp10*yin(73)
                                      zin(81) = zcp00*zin(79) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   79

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   85
                                      ! i3 =   73
                                      ! i4 =   79

                                      xin(85) = c10*xin(73) + xc00*xin(79)
                                      yin(85) = c10*yin(73) + yc00*yin(79)
                                      zin(85) = c10*zin(73) + zc00*zin(79)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   87
                                      ! i5 =   85
                                      ! i4 =   79

                                      xin(87) = xcp00*xin(85) + cp10*xin(79)
                                      yin(87) = ycp00*yin(85) + cp10*yin(79)
                                      zin(87) = zcp00*zin(85) + cp10*zin(79)

                                      ! ------------------

                                      ! i3 = i4 =   79
                                      ! i4 = i5 =   85

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   79
                                      ! i4 =   85

                                      xin(91) = c10*xin(79) + xc00*xin(85)
                                      yin(91) = c10*yin(79) + yc00*yin(85)
                                      zin(91) = c10*zin(79) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   93
                                      ! i5 =   91
                                      ! i4 =   85

                                      xin(93) = xcp00*xin(91) + cp10*xin(85)
                                      yin(93) = ycp00*yin(91) + cp10*yin(85)
                                      zin(93) = zcp00*zin(91) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    3

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    5
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(5) = cp01*xin(1) + xcp00*xin(3)
                                      yin(5) = cp01*yin(1) + ycp00*yin(3)
                                      zin(5) = cp01*zin(1) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   29

                                      xin(29) = xc00*xin(5) + c01*xin(3)
                                      yin(29) = yc00*yin(5) + c01*yin(3)
                                      zin(29) = zc00*zin(5) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    5

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    6
                                      ! i3 =    3
                                      ! i4 =    5

                                      xin(6) = cp01*xin(3) + xcp00*xin(5)
                                      yin(6) = cp01*yin(3) + ycp00*yin(5)
                                      zin(6) = cp01*zin(3) + zcp00*zin(5)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   30

                                      xin(30) = xc00*xin(6) + c01*xin(5)
                                      yin(30) = yc00*yin(6) + c01*yin(5)
                                      zin(30) = zc00*zin(6) + c01*zin(5)

                                      ! ------------------

                                      ! i3 = i4 =    5
                                      ! i4 = i5 =    6

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   49

                                      xin(53) = c10*xin(5) + xc00*xin(29) + c01*xin(27)
                                      yin(53) = c10*yin(5) + yc00*yin(29) + c01*yin(27)
                                      zin(53) = c10*zin(5) + zc00*zin(29) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   73

                                      xin(77) = c10*xin(29) + xc00*xin(53) + c01*xin(51)
                                      yin(77) = c10*yin(29) + yc00*yin(53) + c01*yin(51)
                                      zin(77) = c10*zin(29) + zc00*zin(53) + c01*zin(51)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   73

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   79

                                      xin(83) = c10*xin(53) + xc00*xin(77) + c01*xin(75)
                                      yin(83) = c10*yin(53) + yc00*yin(77) + c01*yin(75)
                                      zin(83) = c10*zin(53) + zc00*zin(77) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   79

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   85

                                      xin(89) = c10*xin(77) + xc00*xin(83) + c01*xin(81)
                                      yin(89) = c10*yin(77) + yc00*yin(83) + c01*yin(81)
                                      zin(89) = c10*zin(77) + zc00*zin(83) + c01*zin(81)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   79
                                      ! i4 = i5 =   85

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   91

                                      xin(95) = c10*xin(83) + xc00*xin(89) + c01*xin(87)
                                      yin(95) = c10*yin(83) + yc00*yin(89) + c01*yin(87)
                                      zin(95) = c10*zin(83) + zc00*zin(89) + c01*zin(87)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   49

                                      xin(54) = c10*xin(6) + xc00*xin(30) + c01*xin(29)
                                      yin(54) = c10*yin(6) + yc00*yin(30) + c01*yin(29)
                                      zin(54) = c10*zin(6) + zc00*zin(30) + c01*zin(29)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   73

                                      xin(78) = c10*xin(30) + xc00*xin(54) + c01*xin(53)
                                      yin(78) = c10*yin(30) + yc00*yin(54) + c01*yin(53)
                                      zin(78) = c10*zin(30) + zc00*zin(54) + c01*zin(53)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   73

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   79

                                      xin(84) = c10*xin(54) + xc00*xin(78) + c01*xin(77)
                                      yin(84) = c10*yin(54) + yc00*yin(78) + c01*yin(77)
                                      zin(84) = c10*zin(54) + zc00*zin(78) + c01*zin(77)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   79

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   85

                                      xin(90) = c10*xin(78) + xc00*xin(84) + c01*xin(83)
                                      yin(90) = c10*yin(78) + yc00*yin(84) + c01*yin(83)
                                      zin(90) = c10*zin(78) + zc00*zin(84) + c01*zin(83)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   79
                                      ! i4 = i5 =   85

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   91

                                      xin(96) = c10*xin(84) + xc00*xin(90) + c01*xin(89)
                                      yin(96) = c10*yin(84) + yc00*yin(90) + c01*yin(89)
                                      zin(96) = c10*zin(84) + zc00*zin(90) + c01*zin(89)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   91

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   91

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   85

                                      xin(91) = xin(91) + dxij*xin(85)
                                      yin(91) = yin(91) + dyij*yin(85)
                                      zin(91) = zin(91) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   79

                                      xin(85) = xin(85) + dxij*xin(79)
                                      yin(85) = yin(85) + dyij*yin(79)
                                      zin(85) = zin(85) + dzij*zin(79)

                                      ! i3 = i4 =   79
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   73

                                      xin(79) = xin(79) + dxij*xin(73)
                                      yin(79) = yin(79) + dyij*yin(73)
                                      zin(79) = zin(79) + dzij*zin(73)

                                      ! i3 = i4 =   73
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   91

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   85

                                      xin(91) = xin(91) + dxij*xin(85)
                                      yin(91) = yin(91) + dyij*yin(85)
                                      zin(91) = zin(91) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   79

                                      xin(85) = xin(85) + dxij*xin(79)
                                      yin(85) = yin(85) + dyij*yin(79)
                                      zin(85) = zin(85) + dzij*zin(79)

                                      ! i3 = i4 =   79
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   91

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   85

                                      xin(91) = xin(91) + dxij*xin(85)
                                      yin(91) = yin(91) + dyij*yin(85)
                                      zin(91) = zin(91) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(25) + dxij*xin(1)
                                      yin(7) = yin(25) + dyij*yin(1)
                                      zin(7) = zin(25) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    2

                                      xin(31) = xin(49) + dxij*xin(25)
                                      yin(31) = yin(49) + dyij*yin(25)
                                      zin(31) = zin(49) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    3

                                      xin(55) = xin(73) + dxij*xin(49)
                                      yin(55) = yin(73) + dyij*yin(49)
                                      zin(55) = zin(73) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    2

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(31) + dxij*xin(7)
                                      yin(13) = yin(31) + dyij*yin(7)
                                      zin(13) = zin(31) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! ni =    2

                                      xin(37) = xin(55) + dxij*xin(31)
                                      yin(37) = yin(55) + dyij*yin(31)
                                      zin(37) = zin(55) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    3

                                      xin(61) = xin(79) + dxij*xin(55)
                                      yin(61) = yin(79) + dyij*yin(55)
                                      zin(61) = zin(79) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    3

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    3

                                      xin(19) = xin(37) + dxij*xin(13)
                                      yin(19) = yin(37) + dyij*yin(13)
                                      zin(19) = zin(37) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    2

                                      xin(43) = xin(61) + dxij*xin(37)
                                      yin(43) = yin(61) + dyij*yin(37)
                                      zin(43) = zin(61) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    3

                                      xin(67) = xin(85) + dxij*xin(61)
                                      yin(67) = yin(85) + dyij*yin(61)
                                      zin(67) = zin(85) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   87

                                      xin(93) = xin(93) + dxij*xin(87)
                                      yin(93) = yin(93) + dyij*yin(87)
                                      zin(93) = zin(93) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   81

                                      xin(87) = xin(87) + dxij*xin(81)
                                      yin(87) = yin(87) + dyij*yin(81)
                                      zin(87) = zin(87) + dzij*zin(81)

                                      ! i3 = i4 =   81
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   75

                                      xin(81) = xin(81) + dxij*xin(75)
                                      yin(81) = yin(81) + dyij*yin(75)
                                      zin(81) = zin(81) + dzij*zin(75)

                                      ! i3 = i4 =   75
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   87

                                      xin(93) = xin(93) + dxij*xin(87)
                                      yin(93) = yin(93) + dyij*yin(87)
                                      zin(93) = zin(93) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   81

                                      xin(87) = xin(87) + dxij*xin(81)
                                      yin(87) = yin(87) + dyij*yin(81)
                                      zin(87) = zin(87) + dzij*zin(81)

                                      ! i3 = i4 =   81
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   87

                                      xin(93) = xin(93) + dxij*xin(87)
                                      yin(93) = yin(93) + dyij*yin(87)
                                      zin(93) = zin(93) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(27) + dxij*xin(3)
                                      yin(9) = yin(27) + dyij*yin(3)
                                      zin(9) = zin(27) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    2

                                      xin(33) = xin(51) + dxij*xin(27)
                                      yin(33) = yin(51) + dyij*yin(27)
                                      zin(33) = zin(51) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    3

                                      xin(57) = xin(75) + dxij*xin(51)
                                      yin(57) = yin(75) + dyij*yin(51)
                                      zin(57) = zin(75) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    2

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    3

                                      xin(15) = xin(33) + dxij*xin(9)
                                      yin(15) = yin(33) + dyij*yin(9)
                                      zin(15) = zin(33) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    2

                                      xin(39) = xin(57) + dxij*xin(33)
                                      yin(39) = yin(57) + dyij*yin(33)
                                      zin(39) = zin(57) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    3

                                      xin(63) = xin(81) + dxij*xin(57)
                                      yin(63) = yin(81) + dyij*yin(57)
                                      zin(63) = zin(81) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    3

                                      ! i4 = i3 =   21

                                      ! do ni = 1,    3

                                      xin(21) = xin(39) + dxij*xin(15)
                                      yin(21) = yin(39) + dyij*yin(15)
                                      zin(21) = zin(39) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    2

                                      xin(45) = xin(63) + dxij*xin(39)
                                      yin(45) = yin(63) + dyij*yin(39)
                                      zin(45) = zin(63) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    3

                                      xin(69) = xin(87) + dxij*xin(63)
                                      yin(69) = yin(87) + dyij*yin(63)
                                      zin(69) = zin(87) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(95) = xin(95) + dxij*xin(89)
                                      yin(95) = yin(95) + dyij*yin(89)
                                      zin(95) = zin(95) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   83

                                      xin(89) = xin(89) + dxij*xin(83)
                                      yin(89) = yin(89) + dyij*yin(83)
                                      zin(89) = zin(89) + dzij*zin(83)

                                      ! i3 = i4 =   83
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   77

                                      xin(83) = xin(83) + dxij*xin(77)
                                      yin(83) = yin(83) + dyij*yin(77)
                                      zin(83) = zin(83) + dzij*zin(77)

                                      ! i3 = i4 =   77
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(95) = xin(95) + dxij*xin(89)
                                      yin(95) = yin(95) + dyij*yin(89)
                                      zin(95) = zin(95) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   83

                                      xin(89) = xin(89) + dxij*xin(83)
                                      yin(89) = yin(89) + dyij*yin(83)
                                      zin(89) = zin(89) + dzij*zin(83)

                                      ! i3 = i4 =   83
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(95) = xin(95) + dxij*xin(89)
                                      yin(95) = yin(95) + dyij*yin(89)
                                      zin(95) = zin(95) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(29) + dxij*xin(5)
                                      yin(11) = yin(29) + dyij*yin(5)
                                      zin(11) = zin(29) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    2

                                      xin(35) = xin(53) + dxij*xin(29)
                                      yin(35) = yin(53) + dyij*yin(29)
                                      zin(35) = zin(53) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    3

                                      xin(59) = xin(77) + dxij*xin(53)
                                      yin(59) = yin(77) + dyij*yin(53)
                                      zin(59) = zin(77) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    3

                                      xin(17) = xin(35) + dxij*xin(11)
                                      yin(17) = yin(35) + dyij*yin(11)
                                      zin(17) = zin(35) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    2

                                      xin(41) = xin(59) + dxij*xin(35)
                                      yin(41) = yin(59) + dyij*yin(35)
                                      zin(41) = zin(59) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    3

                                      xin(65) = xin(83) + dxij*xin(59)
                                      yin(65) = yin(83) + dyij*yin(59)
                                      zin(65) = zin(83) + dzij*zin(59)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   23

                                      ! nj =    3

                                      ! i4 = i3 =   23

                                      ! do ni = 1,    3

                                      xin(23) = xin(41) + dxij*xin(17)
                                      yin(23) = yin(41) + dyij*yin(17)
                                      zin(23) = zin(41) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    2

                                      xin(47) = xin(65) + dxij*xin(41)
                                      yin(47) = yin(65) + dyij*yin(41)
                                      zin(47) = zin(65) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    3

                                      xin(71) = xin(89) + dxij*xin(65)
                                      yin(71) = yin(89) + dyij*yin(65)
                                      zin(71) = zin(89) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   29

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(96) = xin(96) + dxij*xin(90)
                                      yin(96) = yin(96) + dyij*yin(90)
                                      zin(96) = zin(96) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   84

                                      xin(90) = xin(90) + dxij*xin(84)
                                      yin(90) = yin(90) + dyij*yin(84)
                                      zin(90) = zin(90) + dzij*zin(84)

                                      ! i3 = i4 =   84
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   78

                                      xin(84) = xin(84) + dxij*xin(78)
                                      yin(84) = yin(84) + dyij*yin(78)
                                      zin(84) = zin(84) + dzij*zin(78)

                                      ! i3 = i4 =   78
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(96) = xin(96) + dxij*xin(90)
                                      yin(96) = yin(96) + dyij*yin(90)
                                      zin(96) = zin(96) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   84

                                      xin(90) = xin(90) + dxij*xin(84)
                                      yin(90) = yin(90) + dyij*yin(84)
                                      zin(90) = zin(90) + dzij*zin(84)

                                      ! i3 = i4 =   84
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(96) = xin(96) + dxij*xin(90)
                                      yin(96) = yin(96) + dyij*yin(90)
                                      zin(96) = zin(96) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   12

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(30) + dxij*xin(6)
                                      yin(12) = yin(30) + dyij*yin(6)
                                      zin(12) = zin(30) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    2

                                      xin(36) = xin(54) + dxij*xin(30)
                                      yin(36) = yin(54) + dyij*yin(30)
                                      zin(36) = zin(54) + dzij*zin(30)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    3

                                      xin(60) = xin(78) + dxij*xin(54)
                                      yin(60) = yin(78) + dyij*yin(54)
                                      zin(60) = zin(78) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   18

                                      ! nj =    2

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(36) + dxij*xin(12)
                                      yin(18) = yin(36) + dyij*yin(12)
                                      zin(18) = zin(36) + dzij*zin(12)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    2

                                      xin(42) = xin(60) + dxij*xin(36)
                                      yin(42) = yin(60) + dyij*yin(36)
                                      zin(42) = zin(60) + dzij*zin(36)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    3

                                      xin(66) = xin(84) + dxij*xin(60)
                                      yin(66) = yin(84) + dyij*yin(60)
                                      zin(66) = zin(84) + dzij*zin(60)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   24

                                      ! nj =    3

                                      ! i4 = i3 =   24

                                      ! do ni = 1,    3

                                      xin(24) = xin(42) + dxij*xin(18)
                                      yin(24) = yin(42) + dyij*yin(18)
                                      zin(24) = zin(42) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    2

                                      xin(48) = xin(66) + dxij*xin(42)
                                      yin(48) = yin(66) + dyij*yin(42)
                                      zin(48) = zin(66) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    3

                                      xin(72) = xin(90) + dxij*xin(66)
                                      yin(72) = yin(90) + dyij*yin(66)
                                      zin(72) = zin(90) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   30

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =    6

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    5

                                      xin(6) = xin(6) + dxkl*xin(5)
                                      yin(6) = yin(6) + dykl*yin(5)
                                      zin(6) = zin(6) + dzkl*zin(5)

                                      ! i3 = i4 =    5
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    7

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   12

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   11

                                      xin(12) = xin(12) + dxkl*xin(11)
                                      yin(12) = yin(12) + dykl*yin(11)
                                      zin(12) = zin(12) + dzkl*zin(11)

                                      ! i3 = i4 =   11
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =    8

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    8

                                      ! do nk = 1,    2

                                      xin(8) = xin(9) + dxkl*xin(7)
                                      yin(8) = yin(9) + dykl*yin(7)
                                      zin(8) = zin(9) + dzkl*zin(7)
                                      ! i4 = i4 + lang+1 =   10

                                      ! nk =    2

                                      xin(10) = xin(11) + dxkl*xin(9)
                                      yin(10) = yin(11) + dykl*yin(9)
                                      zin(10) = zin(11) + dzkl*zin(9)
                                      ! i4 = i4 + lang+1 =   12

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =    9

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   13

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   18

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   17

                                      xin(18) = xin(18) + dxkl*xin(17)
                                      yin(18) = yin(18) + dykl*yin(17)
                                      zin(18) = zin(18) + dzkl*zin(17)

                                      ! i3 = i4 =   17
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   14

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   14

                                      ! do nk = 1,    2

                                      xin(14) = xin(15) + dxkl*xin(13)
                                      yin(14) = yin(15) + dykl*yin(13)
                                      zin(14) = zin(15) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   16

                                      ! nk =    2

                                      xin(16) = xin(17) + dxkl*xin(15)
                                      yin(16) = yin(17) + dykl*yin(15)
                                      zin(16) = zin(17) + dzkl*zin(15)
                                      ! i4 = i4 + lang+1 =   18

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   15

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   19

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   24

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   23

                                      xin(24) = xin(24) + dxkl*xin(23)
                                      yin(24) = yin(24) + dykl*yin(23)
                                      zin(24) = zin(24) + dzkl*zin(23)

                                      ! i3 = i4 =   23
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   20

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   20

                                      ! do nk = 1,    2

                                      xin(20) = xin(21) + dxkl*xin(19)
                                      yin(20) = yin(21) + dykl*yin(19)
                                      zin(20) = zin(21) + dzkl*zin(19)
                                      ! i4 = i4 + lang+1 =   22

                                      ! nk =    2

                                      xin(22) = xin(23) + dxkl*xin(21)
                                      yin(22) = yin(23) + dykl*yin(21)
                                      zin(22) = zin(23) + dzkl*zin(21)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   21

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   25

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   30

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   29

                                      xin(30) = xin(30) + dxkl*xin(29)
                                      yin(30) = yin(30) + dykl*yin(29)
                                      zin(30) = zin(30) + dzkl*zin(29)

                                      ! i3 = i4 =   29
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   26

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   26

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =   27

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   31

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   32

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   32

                                      ! do nk = 1,    2

                                      xin(32) = xin(33) + dxkl*xin(31)
                                      yin(32) = yin(33) + dykl*yin(31)
                                      zin(32) = zin(33) + dzkl*zin(31)
                                      ! i4 = i4 + lang+1 =   34

                                      ! nk =    2

                                      xin(34) = xin(35) + dxkl*xin(33)
                                      yin(34) = yin(35) + dykl*yin(33)
                                      zin(34) = zin(35) + dzkl*zin(33)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   33

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   37

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   42

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   41

                                      xin(42) = xin(42) + dxkl*xin(41)
                                      yin(42) = yin(42) + dykl*yin(41)
                                      zin(42) = zin(42) + dzkl*zin(41)

                                      ! i3 = i4 =   41
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   38

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   38

                                      ! do nk = 1,    2

                                      xin(38) = xin(39) + dxkl*xin(37)
                                      yin(38) = yin(39) + dykl*yin(37)
                                      zin(38) = zin(39) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   40

                                      ! nk =    2

                                      xin(40) = xin(41) + dxkl*xin(39)
                                      yin(40) = yin(41) + dykl*yin(39)
                                      zin(40) = zin(41) + dzkl*zin(39)
                                      ! i4 = i4 + lang+1 =   42

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   39

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   43

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   48

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   47

                                      xin(48) = xin(48) + dxkl*xin(47)
                                      yin(48) = yin(48) + dykl*yin(47)
                                      zin(48) = zin(48) + dzkl*zin(47)

                                      ! i3 = i4 =   47
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   44

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   44

                                      ! do nk = 1,    2

                                      xin(44) = xin(45) + dxkl*xin(43)
                                      yin(44) = yin(45) + dykl*yin(43)
                                      zin(44) = zin(45) + dzkl*zin(43)
                                      ! i4 = i4 + lang+1 =   46

                                      ! nk =    2

                                      xin(46) = xin(47) + dxkl*xin(45)
                                      yin(46) = yin(47) + dykl*yin(45)
                                      zin(46) = zin(47) + dzkl*zin(45)
                                      ! i4 = i4 + lang+1 =   48

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   45

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   49

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   54

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   53

                                      xin(54) = xin(54) + dxkl*xin(53)
                                      yin(54) = yin(54) + dykl*yin(53)
                                      zin(54) = zin(54) + dzkl*zin(53)

                                      ! i3 = i4 =   53
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   50

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   50

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =   51

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   55

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   60

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   59

                                      xin(60) = xin(60) + dxkl*xin(59)
                                      yin(60) = yin(60) + dykl*yin(59)
                                      zin(60) = zin(60) + dzkl*zin(59)

                                      ! i3 = i4 =   59
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   56

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   56

                                      ! do nk = 1,    2

                                      xin(56) = xin(57) + dxkl*xin(55)
                                      yin(56) = yin(57) + dykl*yin(55)
                                      zin(56) = zin(57) + dzkl*zin(55)
                                      ! i4 = i4 + lang+1 =   58

                                      ! nk =    2

                                      xin(58) = xin(59) + dxkl*xin(57)
                                      yin(58) = yin(59) + dykl*yin(57)
                                      zin(58) = zin(59) + dzkl*zin(57)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   57

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   61

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   66

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   65

                                      xin(66) = xin(66) + dxkl*xin(65)
                                      yin(66) = yin(66) + dykl*yin(65)
                                      zin(66) = zin(66) + dzkl*zin(65)

                                      ! i3 = i4 =   65
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   62

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   62

                                      ! do nk = 1,    2

                                      xin(62) = xin(63) + dxkl*xin(61)
                                      yin(62) = yin(63) + dykl*yin(61)
                                      zin(62) = zin(63) + dzkl*zin(61)
                                      ! i4 = i4 + lang+1 =   64

                                      ! nk =    2

                                      xin(64) = xin(65) + dxkl*xin(63)
                                      yin(64) = yin(65) + dykl*yin(63)
                                      zin(64) = zin(65) + dzkl*zin(63)
                                      ! i4 = i4 + lang+1 =   66

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   63

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   67

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   68

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   68

                                      ! do nk = 1,    2

                                      xin(68) = xin(69) + dxkl*xin(67)
                                      yin(68) = yin(69) + dykl*yin(67)
                                      zin(68) = zin(69) + dzkl*zin(67)
                                      ! i4 = i4 + lang+1 =   70

                                      ! nk =    2

                                      xin(70) = xin(71) + dxkl*xin(69)
                                      yin(70) = yin(71) + dykl*yin(69)
                                      zin(70) = zin(71) + dzkl*zin(69)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   69

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   78

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   77

                                      xin(78) = xin(78) + dxkl*xin(77)
                                      yin(78) = yin(78) + dykl*yin(77)
                                      zin(78) = zin(78) + dzkl*zin(77)

                                      ! i3 = i4 =   77
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   79

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   84

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   83

                                      xin(84) = xin(84) + dxkl*xin(83)
                                      yin(84) = yin(84) + dykl*yin(83)
                                      zin(84) = zin(84) + dzkl*zin(83)

                                      ! i3 = i4 =   83
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   80

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   80

                                      ! do nk = 1,    2

                                      xin(80) = xin(81) + dxkl*xin(79)
                                      yin(80) = yin(81) + dykl*yin(79)
                                      zin(80) = zin(81) + dzkl*zin(79)
                                      ! i4 = i4 + lang+1 =   82

                                      ! nk =    2

                                      xin(82) = xin(83) + dxkl*xin(81)
                                      yin(82) = yin(83) + dykl*yin(81)
                                      zin(82) = zin(83) + dzkl*zin(81)
                                      ! i4 = i4 + lang+1 =   84

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   81

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   85

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   90

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   89

                                      xin(90) = xin(90) + dxkl*xin(89)
                                      yin(90) = yin(90) + dykl*yin(89)
                                      zin(90) = zin(90) + dzkl*zin(89)

                                      ! i3 = i4 =   89
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   86

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   86

                                      ! do nk = 1,    2

                                      xin(86) = xin(87) + dxkl*xin(85)
                                      yin(86) = yin(87) + dykl*yin(85)
                                      zin(86) = zin(87) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   88

                                      ! nk =    2

                                      xin(88) = xin(89) + dxkl*xin(87)
                                      yin(88) = yin(89) + dykl*yin(87)
                                      zin(88) = zin(89) + dzkl*zin(87)
                                      ! i4 = i4 + lang+1 =   90

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   87

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   91

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =   96

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   95

                                      xin(96) = xin(96) + dxkl*xin(95)
                                      yin(96) = yin(96) + dykl*yin(95)
                                      zin(96) = zin(96) + dzkl*zin(95)

                                      ! i3 = i4 =   95
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   92

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   92

                                      ! do nk = 1,    2

                                      xin(92) = xin(93) + dxkl*xin(91)
                                      yin(92) = yin(93) + dykl*yin(91)
                                      zin(92) = zin(93) + dzkl*zin(91)
                                      ! i4 = i4 + lang+1 =   94

                                      ! nk =    2

                                      xin(94) = xin(95) + dxkl*xin(93)
                                      yin(94) = yin(95) + dykl*yin(93)
                                      zin(94) = zin(95) + dzkl*zin(93)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   93

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   97

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

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

                                      ! i1 = in(1) =   97

                                      xin(97) = 1.0_dp
                                      yin(97) = 1.0_dp
                                      zin(97) = f00

                                      ! i2 = in(2) =  121
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(121) = xc00
                                      yin(121) = yc00
                                      zin(121) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   99

                                      xin(99) = xcp00
                                      yin(99) = ycp00
                                      zin(99) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  123
                                      ! i2 =  121

                                      xin(123) = xcp00*xin(121) + cp10
                                      yin(123) = ycp00*yin(121) + cp10
                                      zin(123) = zcp00*zin(121) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  121

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  145
                                      ! i3 =   97
                                      ! i4 =  121

                                      xin(145) = c10*xin(97) + xc00*xin(121)
                                      yin(145) = c10*yin(97) + yc00*yin(121)
                                      zin(145) = c10*zin(97) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  147
                                      ! i5 =  145
                                      ! i4 =  121

                                      xin(147) = xcp00*xin(145) + cp10*xin(121)
                                      yin(147) = ycp00*yin(145) + cp10*yin(121)
                                      zin(147) = zcp00*zin(145) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  145

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  169
                                      ! i3 =  121
                                      ! i4 =  145

                                      xin(169) = c10*xin(121) + xc00*xin(145)
                                      yin(169) = c10*yin(121) + yc00*yin(145)
                                      zin(169) = c10*zin(121) + zc00*zin(145)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  171
                                      ! i5 =  169
                                      ! i4 =  145

                                      xin(171) = xcp00*xin(169) + cp10*xin(145)
                                      yin(171) = ycp00*yin(169) + cp10*yin(145)
                                      zin(171) = zcp00*zin(169) + cp10*zin(145)

                                      ! ------------------

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  169

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  175
                                      ! i3 =  145
                                      ! i4 =  169

                                      xin(175) = c10*xin(145) + xc00*xin(169)
                                      yin(175) = c10*yin(145) + yc00*yin(169)
                                      zin(175) = c10*zin(145) + zc00*zin(169)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  177
                                      ! i5 =  175
                                      ! i4 =  169

                                      xin(177) = xcp00*xin(175) + cp10*xin(169)
                                      yin(177) = ycp00*yin(175) + cp10*yin(169)
                                      zin(177) = zcp00*zin(175) + cp10*zin(169)

                                      ! ------------------

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  175

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  181
                                      ! i3 =  169
                                      ! i4 =  175

                                      xin(181) = c10*xin(169) + xc00*xin(175)
                                      yin(181) = c10*yin(169) + yc00*yin(175)
                                      zin(181) = c10*zin(169) + zc00*zin(175)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  183
                                      ! i5 =  181
                                      ! i4 =  175

                                      xin(183) = xcp00*xin(181) + cp10*xin(175)
                                      yin(183) = ycp00*yin(181) + cp10*yin(175)
                                      zin(183) = zcp00*zin(181) + cp10*zin(175)

                                      ! ------------------

                                      ! i3 = i4 =  175
                                      ! i4 = i5 =  181

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  187
                                      ! i3 =  175
                                      ! i4 =  181

                                      xin(187) = c10*xin(175) + xc00*xin(181)
                                      yin(187) = c10*yin(175) + yc00*yin(181)
                                      zin(187) = c10*zin(175) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  189
                                      ! i5 =  187
                                      ! i4 =  181

                                      xin(189) = xcp00*xin(187) + cp10*xin(181)
                                      yin(189) = ycp00*yin(187) + cp10*yin(181)
                                      zin(189) = zcp00*zin(187) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   99

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  101
                                      ! i3 =   97
                                      ! i4 =   99

                                      xin(101) = cp01*xin(97) + xcp00*xin(99)
                                      yin(101) = cp01*yin(97) + ycp00*yin(99)
                                      zin(101) = cp01*zin(97) + zcp00*zin(99)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  125

                                      xin(125) = xc00*xin(101) + c01*xin(99)
                                      yin(125) = yc00*yin(101) + c01*yin(99)
                                      zin(125) = zc00*zin(101) + c01*zin(99)

                                      ! ------------------

                                      ! i3 = i4 =   99
                                      ! i4 = i5 =  101

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  102
                                      ! i3 =   99
                                      ! i4 =  101

                                      xin(102) = cp01*xin(99) + xcp00*xin(101)
                                      yin(102) = cp01*yin(99) + ycp00*yin(101)
                                      zin(102) = cp01*zin(99) + zcp00*zin(101)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  126

                                      xin(126) = xc00*xin(102) + c01*xin(101)
                                      yin(126) = yc00*yin(102) + c01*yin(101)
                                      zin(126) = zc00*zin(102) + c01*zin(101)

                                      ! ------------------

                                      ! i3 = i4 =  101
                                      ! i4 = i5 =  102

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  121

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  145

                                      xin(149) = c10*xin(101) + xc00*xin(125) + c01*xin(123)
                                      yin(149) = c10*yin(101) + yc00*yin(125) + c01*yin(123)
                                      zin(149) = c10*zin(101) + zc00*zin(125) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  145

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  169

                                      xin(173) = c10*xin(125) + xc00*xin(149) + c01*xin(147)
                                      yin(173) = c10*yin(125) + yc00*yin(149) + c01*yin(147)
                                      zin(173) = c10*zin(125) + zc00*zin(149) + c01*zin(147)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  169

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  175

                                      xin(179) = c10*xin(149) + xc00*xin(173) + c01*xin(171)
                                      yin(179) = c10*yin(149) + yc00*yin(173) + c01*yin(171)
                                      zin(179) = c10*zin(149) + zc00*zin(173) + c01*zin(171)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  175

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  181

                                      xin(185) = c10*xin(173) + xc00*xin(179) + c01*xin(177)
                                      yin(185) = c10*yin(173) + yc00*yin(179) + c01*yin(177)
                                      zin(185) = c10*zin(173) + zc00*zin(179) + c01*zin(177)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  175
                                      ! i4 = i5 =  181

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  187

                                      xin(191) = c10*xin(179) + xc00*xin(185) + c01*xin(183)
                                      yin(191) = c10*yin(179) + yc00*yin(185) + c01*yin(183)
                                      zin(191) = c10*zin(179) + zc00*zin(185) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  121

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  145

                                      xin(150) = c10*xin(102) + xc00*xin(126) + c01*xin(125)
                                      yin(150) = c10*yin(102) + yc00*yin(126) + c01*yin(125)
                                      zin(150) = c10*zin(102) + zc00*zin(126) + c01*zin(125)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  145

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  169

                                      xin(174) = c10*xin(126) + xc00*xin(150) + c01*xin(149)
                                      yin(174) = c10*yin(126) + yc00*yin(150) + c01*yin(149)
                                      zin(174) = c10*zin(126) + zc00*zin(150) + c01*zin(149)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  169

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  175

                                      xin(180) = c10*xin(150) + xc00*xin(174) + c01*xin(173)
                                      yin(180) = c10*yin(150) + yc00*yin(174) + c01*yin(173)
                                      zin(180) = c10*zin(150) + zc00*zin(174) + c01*zin(173)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  175

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  181

                                      xin(186) = c10*xin(174) + xc00*xin(180) + c01*xin(179)
                                      yin(186) = c10*yin(174) + yc00*yin(180) + c01*yin(179)
                                      zin(186) = c10*zin(174) + zc00*zin(180) + c01*zin(179)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  175
                                      ! i4 = i5 =  181

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  187

                                      xin(192) = c10*xin(180) + xc00*xin(186) + c01*xin(185)
                                      yin(192) = c10*yin(180) + yc00*yin(186) + c01*yin(185)
                                      zin(192) = c10*zin(180) + zc00*zin(186) + c01*zin(185)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  187

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  187

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  181

                                      xin(187) = xin(187) + dxij*xin(181)
                                      yin(187) = yin(187) + dyij*yin(181)
                                      zin(187) = zin(187) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  175

                                      xin(181) = xin(181) + dxij*xin(175)
                                      yin(181) = yin(181) + dyij*yin(175)
                                      zin(181) = zin(181) + dzij*zin(175)

                                      ! i3 = i4 =  175
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  169

                                      xin(175) = xin(175) + dxij*xin(169)
                                      yin(175) = yin(175) + dyij*yin(169)
                                      zin(175) = zin(175) + dzij*zin(169)

                                      ! i3 = i4 =  169
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  187

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  181

                                      xin(187) = xin(187) + dxij*xin(181)
                                      yin(187) = yin(187) + dyij*yin(181)
                                      zin(187) = zin(187) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  175

                                      xin(181) = xin(181) + dxij*xin(175)
                                      yin(181) = yin(181) + dyij*yin(175)
                                      zin(181) = zin(181) + dzij*zin(175)

                                      ! i3 = i4 =  175
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  187

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  181

                                      xin(187) = xin(187) + dxij*xin(181)
                                      yin(187) = yin(187) + dyij*yin(181)
                                      zin(187) = zin(187) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  103

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(121) + dxij*xin(97)
                                      yin(103) = yin(121) + dyij*yin(97)
                                      zin(103) = zin(121) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    2

                                      xin(127) = xin(145) + dxij*xin(121)
                                      yin(127) = yin(145) + dyij*yin(121)
                                      zin(127) = zin(145) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  151

                                      ! ni =    3

                                      xin(151) = xin(169) + dxij*xin(145)
                                      yin(151) = yin(169) + dyij*yin(145)
                                      zin(151) = zin(169) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  109

                                      ! nj =    2

                                      ! i4 = i3 =  109

                                      ! do ni = 1,    3

                                      xin(109) = xin(127) + dxij*xin(103)
                                      yin(109) = yin(127) + dyij*yin(103)
                                      zin(109) = zin(127) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    2

                                      xin(133) = xin(151) + dxij*xin(127)
                                      yin(133) = yin(151) + dyij*yin(127)
                                      zin(133) = zin(151) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! ni =    3

                                      xin(157) = xin(175) + dxij*xin(151)
                                      yin(157) = yin(175) + dyij*yin(151)
                                      zin(157) = zin(175) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  115

                                      ! nj =    3

                                      ! i4 = i3 =  115

                                      ! do ni = 1,    3

                                      xin(115) = xin(133) + dxij*xin(109)
                                      yin(115) = yin(133) + dyij*yin(109)
                                      zin(115) = zin(133) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    2

                                      xin(139) = xin(157) + dxij*xin(133)
                                      yin(139) = yin(157) + dyij*yin(133)
                                      zin(139) = zin(157) + dzij*zin(133)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! ni =    3

                                      xin(163) = xin(181) + dxij*xin(157)
                                      yin(163) = yin(181) + dyij*yin(157)
                                      zin(163) = zin(181) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  121

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  183

                                      xin(189) = xin(189) + dxij*xin(183)
                                      yin(189) = yin(189) + dyij*yin(183)
                                      zin(189) = zin(189) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  177

                                      xin(183) = xin(183) + dxij*xin(177)
                                      yin(183) = yin(183) + dyij*yin(177)
                                      zin(183) = zin(183) + dzij*zin(177)

                                      ! i3 = i4 =  177
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  171

                                      xin(177) = xin(177) + dxij*xin(171)
                                      yin(177) = yin(177) + dyij*yin(171)
                                      zin(177) = zin(177) + dzij*zin(171)

                                      ! i3 = i4 =  171
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  183

                                      xin(189) = xin(189) + dxij*xin(183)
                                      yin(189) = yin(189) + dyij*yin(183)
                                      zin(189) = zin(189) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  177

                                      xin(183) = xin(183) + dxij*xin(177)
                                      yin(183) = yin(183) + dyij*yin(177)
                                      zin(183) = zin(183) + dzij*zin(177)

                                      ! i3 = i4 =  177
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  183

                                      xin(189) = xin(189) + dxij*xin(183)
                                      yin(189) = yin(189) + dyij*yin(183)
                                      zin(189) = zin(189) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  105

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  105

                                      ! do ni = 1,    3

                                      xin(105) = xin(123) + dxij*xin(99)
                                      yin(105) = yin(123) + dyij*yin(99)
                                      zin(105) = zin(123) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    2

                                      xin(129) = xin(147) + dxij*xin(123)
                                      yin(129) = yin(147) + dyij*yin(123)
                                      zin(129) = zin(147) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  153

                                      ! ni =    3

                                      xin(153) = xin(171) + dxij*xin(147)
                                      yin(153) = yin(171) + dyij*yin(147)
                                      zin(153) = zin(171) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  111

                                      ! nj =    2

                                      ! i4 = i3 =  111

                                      ! do ni = 1,    3

                                      xin(111) = xin(129) + dxij*xin(105)
                                      yin(111) = yin(129) + dyij*yin(105)
                                      zin(111) = zin(129) + dzij*zin(105)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    2

                                      xin(135) = xin(153) + dxij*xin(129)
                                      yin(135) = yin(153) + dyij*yin(129)
                                      zin(135) = zin(153) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  159

                                      ! ni =    3

                                      xin(159) = xin(177) + dxij*xin(153)
                                      yin(159) = yin(177) + dyij*yin(153)
                                      zin(159) = zin(177) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  183

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  117

                                      ! nj =    3

                                      ! i4 = i3 =  117

                                      ! do ni = 1,    3

                                      xin(117) = xin(135) + dxij*xin(111)
                                      yin(117) = yin(135) + dyij*yin(111)
                                      zin(117) = zin(135) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    2

                                      xin(141) = xin(159) + dxij*xin(135)
                                      yin(141) = yin(159) + dyij*yin(135)
                                      zin(141) = zin(159) + dzij*zin(135)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    3

                                      xin(165) = xin(183) + dxij*xin(159)
                                      yin(165) = yin(183) + dyij*yin(159)
                                      zin(165) = zin(183) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  123

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(191) = xin(191) + dxij*xin(185)
                                      yin(191) = yin(191) + dyij*yin(185)
                                      zin(191) = zin(191) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  179

                                      xin(185) = xin(185) + dxij*xin(179)
                                      yin(185) = yin(185) + dyij*yin(179)
                                      zin(185) = zin(185) + dzij*zin(179)

                                      ! i3 = i4 =  179
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  173

                                      xin(179) = xin(179) + dxij*xin(173)
                                      yin(179) = yin(179) + dyij*yin(173)
                                      zin(179) = zin(179) + dzij*zin(173)

                                      ! i3 = i4 =  173
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(191) = xin(191) + dxij*xin(185)
                                      yin(191) = yin(191) + dyij*yin(185)
                                      zin(191) = zin(191) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  179

                                      xin(185) = xin(185) + dxij*xin(179)
                                      yin(185) = yin(185) + dyij*yin(179)
                                      zin(185) = zin(185) + dzij*zin(179)

                                      ! i3 = i4 =  179
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(191) = xin(191) + dxij*xin(185)
                                      yin(191) = yin(191) + dyij*yin(185)
                                      zin(191) = zin(191) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  107

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  107

                                      ! do ni = 1,    3

                                      xin(107) = xin(125) + dxij*xin(101)
                                      yin(107) = yin(125) + dyij*yin(101)
                                      zin(107) = zin(125) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    2

                                      xin(131) = xin(149) + dxij*xin(125)
                                      yin(131) = yin(149) + dyij*yin(125)
                                      zin(131) = zin(149) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  155

                                      ! ni =    3

                                      xin(155) = xin(173) + dxij*xin(149)
                                      yin(155) = yin(173) + dyij*yin(149)
                                      zin(155) = zin(173) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  113

                                      ! nj =    2

                                      ! i4 = i3 =  113

                                      ! do ni = 1,    3

                                      xin(113) = xin(131) + dxij*xin(107)
                                      yin(113) = yin(131) + dyij*yin(107)
                                      zin(113) = zin(131) + dzij*zin(107)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    2

                                      xin(137) = xin(155) + dxij*xin(131)
                                      yin(137) = yin(155) + dyij*yin(131)
                                      zin(137) = zin(155) + dzij*zin(131)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    3

                                      xin(161) = xin(179) + dxij*xin(155)
                                      yin(161) = yin(179) + dyij*yin(155)
                                      zin(161) = zin(179) + dzij*zin(155)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  119

                                      ! nj =    3

                                      ! i4 = i3 =  119

                                      ! do ni = 1,    3

                                      xin(119) = xin(137) + dxij*xin(113)
                                      yin(119) = yin(137) + dyij*yin(113)
                                      zin(119) = zin(137) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    2

                                      xin(143) = xin(161) + dxij*xin(137)
                                      yin(143) = yin(161) + dyij*yin(137)
                                      zin(143) = zin(161) + dzij*zin(137)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    3

                                      xin(167) = xin(185) + dxij*xin(161)
                                      yin(167) = yin(185) + dyij*yin(161)
                                      zin(167) = zin(185) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  125

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(192) = xin(192) + dxij*xin(186)
                                      yin(192) = yin(192) + dyij*yin(186)
                                      zin(192) = zin(192) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  180

                                      xin(186) = xin(186) + dxij*xin(180)
                                      yin(186) = yin(186) + dyij*yin(180)
                                      zin(186) = zin(186) + dzij*zin(180)

                                      ! i3 = i4 =  180
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  174

                                      xin(180) = xin(180) + dxij*xin(174)
                                      yin(180) = yin(180) + dyij*yin(174)
                                      zin(180) = zin(180) + dzij*zin(174)

                                      ! i3 = i4 =  174
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(192) = xin(192) + dxij*xin(186)
                                      yin(192) = yin(192) + dyij*yin(186)
                                      zin(192) = zin(192) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  180

                                      xin(186) = xin(186) + dxij*xin(180)
                                      yin(186) = yin(186) + dyij*yin(180)
                                      zin(186) = zin(186) + dzij*zin(180)

                                      ! i3 = i4 =  180
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(192) = xin(192) + dxij*xin(186)
                                      yin(192) = yin(192) + dyij*yin(186)
                                      zin(192) = zin(192) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  108

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  108

                                      ! do ni = 1,    3

                                      xin(108) = xin(126) + dxij*xin(102)
                                      yin(108) = yin(126) + dyij*yin(102)
                                      zin(108) = zin(126) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    2

                                      xin(132) = xin(150) + dxij*xin(126)
                                      yin(132) = yin(150) + dyij*yin(126)
                                      zin(132) = zin(150) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  156

                                      ! ni =    3

                                      xin(156) = xin(174) + dxij*xin(150)
                                      yin(156) = yin(174) + dyij*yin(150)
                                      zin(156) = zin(174) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  114

                                      ! nj =    2

                                      ! i4 = i3 =  114

                                      ! do ni = 1,    3

                                      xin(114) = xin(132) + dxij*xin(108)
                                      yin(114) = yin(132) + dyij*yin(108)
                                      zin(114) = zin(132) + dzij*zin(108)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    2

                                      xin(138) = xin(156) + dxij*xin(132)
                                      yin(138) = yin(156) + dyij*yin(132)
                                      zin(138) = zin(156) + dzij*zin(132)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  162

                                      ! ni =    3

                                      xin(162) = xin(180) + dxij*xin(156)
                                      yin(162) = yin(180) + dyij*yin(156)
                                      zin(162) = zin(180) + dzij*zin(156)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  186

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  120

                                      ! nj =    3

                                      ! i4 = i3 =  120

                                      ! do ni = 1,    3

                                      xin(120) = xin(138) + dxij*xin(114)
                                      yin(120) = yin(138) + dyij*yin(114)
                                      zin(120) = zin(138) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    2

                                      xin(144) = xin(162) + dxij*xin(138)
                                      yin(144) = yin(162) + dyij*yin(138)
                                      zin(144) = zin(162) + dzij*zin(138)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    3

                                      xin(168) = xin(186) + dxij*xin(162)
                                      yin(168) = yin(186) + dyij*yin(162)
                                      zin(168) = zin(186) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  126

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =   97

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  102

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  101

                                      xin(102) = xin(102) + dxkl*xin(101)
                                      yin(102) = yin(102) + dykl*yin(101)
                                      zin(102) = zin(102) + dzkl*zin(101)

                                      ! i3 = i4 =  101
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =   98

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   98

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =   99

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  103

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  104

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  104

                                      ! do nk = 1,    2

                                      xin(104) = xin(105) + dxkl*xin(103)
                                      yin(104) = yin(105) + dykl*yin(103)
                                      zin(104) = zin(105) + dzkl*zin(103)
                                      ! i4 = i4 + lang+1 =  106

                                      ! nk =    2

                                      xin(106) = xin(107) + dxkl*xin(105)
                                      yin(106) = yin(107) + dykl*yin(105)
                                      zin(106) = zin(107) + dzkl*zin(105)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  105

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  109

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  114

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  113

                                      xin(114) = xin(114) + dxkl*xin(113)
                                      yin(114) = yin(114) + dykl*yin(113)
                                      zin(114) = zin(114) + dzkl*zin(113)

                                      ! i3 = i4 =  113
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  110

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  110

                                      ! do nk = 1,    2

                                      xin(110) = xin(111) + dxkl*xin(109)
                                      yin(110) = yin(111) + dykl*yin(109)
                                      zin(110) = zin(111) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  112

                                      ! nk =    2

                                      xin(112) = xin(113) + dxkl*xin(111)
                                      yin(112) = yin(113) + dykl*yin(111)
                                      zin(112) = zin(113) + dzkl*zin(111)
                                      ! i4 = i4 + lang+1 =  114

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  111

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  115

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  120

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  119

                                      xin(120) = xin(120) + dxkl*xin(119)
                                      yin(120) = yin(120) + dykl*yin(119)
                                      zin(120) = zin(120) + dzkl*zin(119)

                                      ! i3 = i4 =  119
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  116

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  116

                                      ! do nk = 1,    2

                                      xin(116) = xin(117) + dxkl*xin(115)
                                      yin(116) = yin(117) + dykl*yin(115)
                                      zin(116) = zin(117) + dzkl*zin(115)
                                      ! i4 = i4 + lang+1 =  118

                                      ! nk =    2

                                      xin(118) = xin(119) + dxkl*xin(117)
                                      yin(118) = yin(119) + dykl*yin(117)
                                      zin(118) = zin(119) + dzkl*zin(117)
                                      ! i4 = i4 + lang+1 =  120

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  117

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  121

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  126

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  125

                                      xin(126) = xin(126) + dxkl*xin(125)
                                      yin(126) = yin(126) + dykl*yin(125)
                                      zin(126) = zin(126) + dzkl*zin(125)

                                      ! i3 = i4 =  125
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  122

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  122

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  123

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  127

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  132

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  131

                                      xin(132) = xin(132) + dxkl*xin(131)
                                      yin(132) = yin(132) + dykl*yin(131)
                                      zin(132) = zin(132) + dzkl*zin(131)

                                      ! i3 = i4 =  131
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  128

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  128

                                      ! do nk = 1,    2

                                      xin(128) = xin(129) + dxkl*xin(127)
                                      yin(128) = yin(129) + dykl*yin(127)
                                      zin(128) = zin(129) + dzkl*zin(127)
                                      ! i4 = i4 + lang+1 =  130

                                      ! nk =    2

                                      xin(130) = xin(131) + dxkl*xin(129)
                                      yin(130) = yin(131) + dykl*yin(129)
                                      zin(130) = zin(131) + dzkl*zin(129)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  129

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  133

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  138

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  137

                                      xin(138) = xin(138) + dxkl*xin(137)
                                      yin(138) = yin(138) + dykl*yin(137)
                                      zin(138) = zin(138) + dzkl*zin(137)

                                      ! i3 = i4 =  137
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  134

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  134

                                      ! do nk = 1,    2

                                      xin(134) = xin(135) + dxkl*xin(133)
                                      yin(134) = yin(135) + dykl*yin(133)
                                      zin(134) = zin(135) + dzkl*zin(133)
                                      ! i4 = i4 + lang+1 =  136

                                      ! nk =    2

                                      xin(136) = xin(137) + dxkl*xin(135)
                                      yin(136) = yin(137) + dykl*yin(135)
                                      zin(136) = zin(137) + dzkl*zin(135)
                                      ! i4 = i4 + lang+1 =  138

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  135

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  139

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  140

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  140

                                      ! do nk = 1,    2

                                      xin(140) = xin(141) + dxkl*xin(139)
                                      yin(140) = yin(141) + dykl*yin(139)
                                      zin(140) = zin(141) + dzkl*zin(139)
                                      ! i4 = i4 + lang+1 =  142

                                      ! nk =    2

                                      xin(142) = xin(143) + dxkl*xin(141)
                                      yin(142) = yin(143) + dykl*yin(141)
                                      zin(142) = zin(143) + dzkl*zin(141)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  141

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  150

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  149

                                      xin(150) = xin(150) + dxkl*xin(149)
                                      yin(150) = yin(150) + dykl*yin(149)
                                      zin(150) = zin(150) + dzkl*zin(149)

                                      ! i3 = i4 =  149
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  151

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  156

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  155

                                      xin(156) = xin(156) + dxkl*xin(155)
                                      yin(156) = yin(156) + dykl*yin(155)
                                      zin(156) = zin(156) + dzkl*zin(155)

                                      ! i3 = i4 =  155
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  152

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  152

                                      ! do nk = 1,    2

                                      xin(152) = xin(153) + dxkl*xin(151)
                                      yin(152) = yin(153) + dykl*yin(151)
                                      zin(152) = zin(153) + dzkl*zin(151)
                                      ! i4 = i4 + lang+1 =  154

                                      ! nk =    2

                                      xin(154) = xin(155) + dxkl*xin(153)
                                      yin(154) = yin(155) + dykl*yin(153)
                                      zin(154) = zin(155) + dzkl*zin(153)
                                      ! i4 = i4 + lang+1 =  156

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  153

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  157

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  162

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  161

                                      xin(162) = xin(162) + dxkl*xin(161)
                                      yin(162) = yin(162) + dykl*yin(161)
                                      zin(162) = zin(162) + dzkl*zin(161)

                                      ! i3 = i4 =  161
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  158

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  158

                                      ! do nk = 1,    2

                                      xin(158) = xin(159) + dxkl*xin(157)
                                      yin(158) = yin(159) + dykl*yin(157)
                                      zin(158) = zin(159) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  160

                                      ! nk =    2

                                      xin(160) = xin(161) + dxkl*xin(159)
                                      yin(160) = yin(161) + dykl*yin(159)
                                      zin(160) = zin(161) + dzkl*zin(159)
                                      ! i4 = i4 + lang+1 =  162

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  159

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  163

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  168

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  167

                                      xin(168) = xin(168) + dxkl*xin(167)
                                      yin(168) = yin(168) + dykl*yin(167)
                                      zin(168) = zin(168) + dzkl*zin(167)

                                      ! i3 = i4 =  167
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  164

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  164

                                      ! do nk = 1,    2

                                      xin(164) = xin(165) + dxkl*xin(163)
                                      yin(164) = yin(165) + dykl*yin(163)
                                      zin(164) = zin(165) + dzkl*zin(163)
                                      ! i4 = i4 + lang+1 =  166

                                      ! nk =    2

                                      xin(166) = xin(167) + dxkl*xin(165)
                                      yin(166) = yin(167) + dykl*yin(165)
                                      zin(166) = zin(167) + dzkl*zin(165)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  165

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  169

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  174

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  173

                                      xin(174) = xin(174) + dxkl*xin(173)
                                      yin(174) = yin(174) + dykl*yin(173)
                                      zin(174) = zin(174) + dzkl*zin(173)

                                      ! i3 = i4 =  173
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  170

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  170

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  171

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  175

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  176

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  176

                                      ! do nk = 1,    2

                                      xin(176) = xin(177) + dxkl*xin(175)
                                      yin(176) = yin(177) + dykl*yin(175)
                                      zin(176) = zin(177) + dzkl*zin(175)
                                      ! i4 = i4 + lang+1 =  178

                                      ! nk =    2

                                      xin(178) = xin(179) + dxkl*xin(177)
                                      yin(178) = yin(179) + dykl*yin(177)
                                      zin(178) = zin(179) + dzkl*zin(177)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  177

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  181

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  186

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  185

                                      xin(186) = xin(186) + dxkl*xin(185)
                                      yin(186) = yin(186) + dykl*yin(185)
                                      zin(186) = zin(186) + dzkl*zin(185)

                                      ! i3 = i4 =  185
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  182

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  182

                                      ! do nk = 1,    2

                                      xin(182) = xin(183) + dxkl*xin(181)
                                      yin(182) = yin(183) + dykl*yin(181)
                                      zin(182) = zin(183) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  184

                                      ! nk =    2

                                      xin(184) = xin(185) + dxkl*xin(183)
                                      yin(184) = yin(185) + dykl*yin(183)
                                      zin(184) = zin(185) + dzkl*zin(183)
                                      ! i4 = i4 + lang+1 =  186

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  183

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  187

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  192

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  191

                                      xin(192) = xin(192) + dxkl*xin(191)
                                      yin(192) = yin(192) + dykl*yin(191)
                                      zin(192) = zin(192) + dzkl*zin(191)

                                      ! i3 = i4 =  191
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  188

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  188

                                      ! do nk = 1,    2

                                      xin(188) = xin(189) + dxkl*xin(187)
                                      yin(188) = yin(189) + dykl*yin(187)
                                      zin(188) = zin(189) + dzkl*zin(187)
                                      ! i4 = i4 + lang+1 =  190

                                      ! nk =    2

                                      xin(190) = xin(191) + dxkl*xin(189)
                                      yin(190) = yin(191) + dykl*yin(189)
                                      zin(190) = zin(191) + dzkl*zin(189)
                                      ! i4 = i4 + lang+1 =  192

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  189

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  193

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

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

                                      ! i1 = in(1) =  193

                                      xin(193) = 1.0_dp
                                      yin(193) = 1.0_dp
                                      zin(193) = f00

                                      ! i2 = in(2) =  217
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(217) = xc00
                                      yin(217) = yc00
                                      zin(217) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  195

                                      xin(195) = xcp00
                                      yin(195) = ycp00
                                      zin(195) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  219
                                      ! i2 =  217

                                      xin(219) = xcp00*xin(217) + cp10
                                      yin(219) = ycp00*yin(217) + cp10
                                      zin(219) = zcp00*zin(217) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  217

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  241
                                      ! i3 =  193
                                      ! i4 =  217

                                      xin(241) = c10*xin(193) + xc00*xin(217)
                                      yin(241) = c10*yin(193) + yc00*yin(217)
                                      zin(241) = c10*zin(193) + zc00*zin(217)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  243
                                      ! i5 =  241
                                      ! i4 =  217

                                      xin(243) = xcp00*xin(241) + cp10*xin(217)
                                      yin(243) = ycp00*yin(241) + cp10*yin(217)
                                      zin(243) = zcp00*zin(241) + cp10*zin(217)

                                      ! ------------------

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  241

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  265
                                      ! i3 =  217
                                      ! i4 =  241

                                      xin(265) = c10*xin(217) + xc00*xin(241)
                                      yin(265) = c10*yin(217) + yc00*yin(241)
                                      zin(265) = c10*zin(217) + zc00*zin(241)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  267
                                      ! i5 =  265
                                      ! i4 =  241

                                      xin(267) = xcp00*xin(265) + cp10*xin(241)
                                      yin(267) = ycp00*yin(265) + cp10*yin(241)
                                      zin(267) = zcp00*zin(265) + cp10*zin(241)

                                      ! ------------------

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  271
                                      ! i3 =  241
                                      ! i4 =  265

                                      xin(271) = c10*xin(241) + xc00*xin(265)
                                      yin(271) = c10*yin(241) + yc00*yin(265)
                                      zin(271) = c10*zin(241) + zc00*zin(265)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  273
                                      ! i5 =  271
                                      ! i4 =  265

                                      xin(273) = xcp00*xin(271) + cp10*xin(265)
                                      yin(273) = ycp00*yin(271) + cp10*yin(265)
                                      zin(273) = zcp00*zin(271) + cp10*zin(265)

                                      ! ------------------

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  271

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  277
                                      ! i3 =  265
                                      ! i4 =  271

                                      xin(277) = c10*xin(265) + xc00*xin(271)
                                      yin(277) = c10*yin(265) + yc00*yin(271)
                                      zin(277) = c10*zin(265) + zc00*zin(271)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  279
                                      ! i5 =  277
                                      ! i4 =  271

                                      xin(279) = xcp00*xin(277) + cp10*xin(271)
                                      yin(279) = ycp00*yin(277) + cp10*yin(271)
                                      zin(279) = zcp00*zin(277) + cp10*zin(271)

                                      ! ------------------

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  277

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  283
                                      ! i3 =  271
                                      ! i4 =  277

                                      xin(283) = c10*xin(271) + xc00*xin(277)
                                      yin(283) = c10*yin(271) + yc00*yin(277)
                                      zin(283) = c10*zin(271) + zc00*zin(277)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  285
                                      ! i5 =  283
                                      ! i4 =  277

                                      xin(285) = xcp00*xin(283) + cp10*xin(277)
                                      yin(285) = ycp00*yin(283) + cp10*yin(277)
                                      zin(285) = zcp00*zin(283) + cp10*zin(277)

                                      ! ------------------

                                      ! i3 = i4 =  277
                                      ! i4 = i5 =  283

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  193
                                      ! i4 = i1+k2 =  195

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  197
                                      ! i3 =  193
                                      ! i4 =  195

                                      xin(197) = cp01*xin(193) + xcp00*xin(195)
                                      yin(197) = cp01*yin(193) + ycp00*yin(195)
                                      zin(197) = cp01*zin(193) + zcp00*zin(195)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  221

                                      xin(221) = xc00*xin(197) + c01*xin(195)
                                      yin(221) = yc00*yin(197) + c01*yin(195)
                                      zin(221) = zc00*zin(197) + c01*zin(195)

                                      ! ------------------

                                      ! i3 = i4 =  195
                                      ! i4 = i5 =  197

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  198
                                      ! i3 =  195
                                      ! i4 =  197

                                      xin(198) = cp01*xin(195) + xcp00*xin(197)
                                      yin(198) = cp01*yin(195) + ycp00*yin(197)
                                      zin(198) = cp01*zin(195) + zcp00*zin(197)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  222

                                      xin(222) = xc00*xin(198) + c01*xin(197)
                                      yin(222) = yc00*yin(198) + c01*yin(197)
                                      zin(222) = zc00*zin(198) + c01*zin(197)

                                      ! ------------------

                                      ! i3 = i4 =  197
                                      ! i4 = i5 =  198

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  217

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  241

                                      xin(245) = c10*xin(197) + xc00*xin(221) + c01*xin(219)
                                      yin(245) = c10*yin(197) + yc00*yin(221) + c01*yin(219)
                                      zin(245) = c10*zin(197) + zc00*zin(221) + c01*zin(219)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  241

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  265

                                      xin(269) = c10*xin(221) + xc00*xin(245) + c01*xin(243)
                                      yin(269) = c10*yin(221) + yc00*yin(245) + c01*yin(243)
                                      zin(269) = c10*zin(221) + zc00*zin(245) + c01*zin(243)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  271

                                      xin(275) = c10*xin(245) + xc00*xin(269) + c01*xin(267)
                                      yin(275) = c10*yin(245) + yc00*yin(269) + c01*yin(267)
                                      zin(275) = c10*zin(245) + zc00*zin(269) + c01*zin(267)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  271

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  277

                                      xin(281) = c10*xin(269) + xc00*xin(275) + c01*xin(273)
                                      yin(281) = c10*yin(269) + yc00*yin(275) + c01*yin(273)
                                      zin(281) = c10*zin(269) + zc00*zin(275) + c01*zin(273)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  277

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  283

                                      xin(287) = c10*xin(275) + xc00*xin(281) + c01*xin(279)
                                      yin(287) = c10*yin(275) + yc00*yin(281) + c01*yin(279)
                                      zin(287) = c10*zin(275) + zc00*zin(281) + c01*zin(279)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  277
                                      ! i4 = i5 =  283

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  217

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  241

                                      xin(246) = c10*xin(198) + xc00*xin(222) + c01*xin(221)
                                      yin(246) = c10*yin(198) + yc00*yin(222) + c01*yin(221)
                                      zin(246) = c10*zin(198) + zc00*zin(222) + c01*zin(221)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  241

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  265

                                      xin(270) = c10*xin(222) + xc00*xin(246) + c01*xin(245)
                                      yin(270) = c10*yin(222) + yc00*yin(246) + c01*yin(245)
                                      zin(270) = c10*zin(222) + zc00*zin(246) + c01*zin(245)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  271

                                      xin(276) = c10*xin(246) + xc00*xin(270) + c01*xin(269)
                                      yin(276) = c10*yin(246) + yc00*yin(270) + c01*yin(269)
                                      zin(276) = c10*zin(246) + zc00*zin(270) + c01*zin(269)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  271

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  277

                                      xin(282) = c10*xin(270) + xc00*xin(276) + c01*xin(275)
                                      yin(282) = c10*yin(270) + yc00*yin(276) + c01*yin(275)
                                      zin(282) = c10*zin(270) + zc00*zin(276) + c01*zin(275)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  277

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  283

                                      xin(288) = c10*xin(276) + xc00*xin(282) + c01*xin(281)
                                      yin(288) = c10*yin(276) + yc00*yin(282) + c01*yin(281)
                                      zin(288) = c10*zin(276) + zc00*zin(282) + c01*zin(281)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  277
                                      ! i4 = i5 =  283

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  283

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(283) = xin(283) + dxij*xin(277)
                                      yin(283) = yin(283) + dyij*yin(277)
                                      zin(283) = zin(283) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  271

                                      xin(277) = xin(277) + dxij*xin(271)
                                      yin(277) = yin(277) + dyij*yin(271)
                                      zin(277) = zin(277) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  265

                                      xin(271) = xin(271) + dxij*xin(265)
                                      yin(271) = yin(271) + dyij*yin(265)
                                      zin(271) = zin(271) + dzij*zin(265)

                                      ! i3 = i4 =  265
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(283) = xin(283) + dxij*xin(277)
                                      yin(283) = yin(283) + dyij*yin(277)
                                      zin(283) = zin(283) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  271

                                      xin(277) = xin(277) + dxij*xin(271)
                                      yin(277) = yin(277) + dyij*yin(271)
                                      zin(277) = zin(277) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(283) = xin(283) + dxij*xin(277)
                                      yin(283) = yin(283) + dyij*yin(277)
                                      zin(283) = zin(283) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  199

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  199

                                      ! do ni = 1,    3

                                      xin(199) = xin(217) + dxij*xin(193)
                                      yin(199) = yin(217) + dyij*yin(193)
                                      zin(199) = zin(217) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  223

                                      ! ni =    2

                                      xin(223) = xin(241) + dxij*xin(217)
                                      yin(223) = yin(241) + dyij*yin(217)
                                      zin(223) = zin(241) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                      ! ni =    3

                                      xin(247) = xin(265) + dxij*xin(241)
                                      yin(247) = yin(265) + dyij*yin(241)
                                      zin(247) = zin(265) + dzij*zin(241)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  271

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  205

                                      ! nj =    2

                                      ! i4 = i3 =  205

                                      ! do ni = 1,    3

                                      xin(205) = xin(223) + dxij*xin(199)
                                      yin(205) = yin(223) + dyij*yin(199)
                                      zin(205) = zin(223) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  229

                                      ! ni =    2

                                      xin(229) = xin(247) + dxij*xin(223)
                                      yin(229) = yin(247) + dyij*yin(223)
                                      zin(229) = zin(247) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! ni =    3

                                      xin(253) = xin(271) + dxij*xin(247)
                                      yin(253) = yin(271) + dyij*yin(247)
                                      zin(253) = zin(271) + dzij*zin(247)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  211

                                      ! nj =    3

                                      ! i4 = i3 =  211

                                      ! do ni = 1,    3

                                      xin(211) = xin(229) + dxij*xin(205)
                                      yin(211) = yin(229) + dyij*yin(205)
                                      zin(211) = zin(229) + dzij*zin(205)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    2

                                      xin(235) = xin(253) + dxij*xin(229)
                                      yin(235) = yin(253) + dyij*yin(229)
                                      zin(235) = zin(253) + dzij*zin(229)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  259

                                      ! ni =    3

                                      xin(259) = xin(277) + dxij*xin(253)
                                      yin(259) = yin(277) + dyij*yin(253)
                                      zin(259) = zin(277) + dzij*zin(253)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  217

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(285) = xin(285) + dxij*xin(279)
                                      yin(285) = yin(285) + dyij*yin(279)
                                      zin(285) = zin(285) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  273

                                      xin(279) = xin(279) + dxij*xin(273)
                                      yin(279) = yin(279) + dyij*yin(273)
                                      zin(279) = zin(279) + dzij*zin(273)

                                      ! i3 = i4 =  273
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  267

                                      xin(273) = xin(273) + dxij*xin(267)
                                      yin(273) = yin(273) + dyij*yin(267)
                                      zin(273) = zin(273) + dzij*zin(267)

                                      ! i3 = i4 =  267
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(285) = xin(285) + dxij*xin(279)
                                      yin(285) = yin(285) + dyij*yin(279)
                                      zin(285) = zin(285) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  273

                                      xin(279) = xin(279) + dxij*xin(273)
                                      yin(279) = yin(279) + dyij*yin(273)
                                      zin(279) = zin(279) + dzij*zin(273)

                                      ! i3 = i4 =  273
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(285) = xin(285) + dxij*xin(279)
                                      yin(285) = yin(285) + dyij*yin(279)
                                      zin(285) = zin(285) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  201

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  201

                                      ! do ni = 1,    3

                                      xin(201) = xin(219) + dxij*xin(195)
                                      yin(201) = yin(219) + dyij*yin(195)
                                      zin(201) = zin(219) + dzij*zin(195)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  225

                                      ! ni =    2

                                      xin(225) = xin(243) + dxij*xin(219)
                                      yin(225) = yin(243) + dyij*yin(219)
                                      zin(225) = zin(243) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  249

                                      ! ni =    3

                                      xin(249) = xin(267) + dxij*xin(243)
                                      yin(249) = yin(267) + dyij*yin(243)
                                      zin(249) = zin(267) + dzij*zin(243)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  273

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  207

                                      ! nj =    2

                                      ! i4 = i3 =  207

                                      ! do ni = 1,    3

                                      xin(207) = xin(225) + dxij*xin(201)
                                      yin(207) = yin(225) + dyij*yin(201)
                                      zin(207) = zin(225) + dzij*zin(201)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  231

                                      ! ni =    2

                                      xin(231) = xin(249) + dxij*xin(225)
                                      yin(231) = yin(249) + dyij*yin(225)
                                      zin(231) = zin(249) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  255

                                      ! ni =    3

                                      xin(255) = xin(273) + dxij*xin(249)
                                      yin(255) = yin(273) + dyij*yin(249)
                                      zin(255) = zin(273) + dzij*zin(249)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  279

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  213

                                      ! nj =    3

                                      ! i4 = i3 =  213

                                      ! do ni = 1,    3

                                      xin(213) = xin(231) + dxij*xin(207)
                                      yin(213) = yin(231) + dyij*yin(207)
                                      zin(213) = zin(231) + dzij*zin(207)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  237

                                      ! ni =    2

                                      xin(237) = xin(255) + dxij*xin(231)
                                      yin(237) = yin(255) + dyij*yin(231)
                                      zin(237) = zin(255) + dzij*zin(231)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  261

                                      ! ni =    3

                                      xin(261) = xin(279) + dxij*xin(255)
                                      yin(261) = yin(279) + dyij*yin(255)
                                      zin(261) = zin(279) + dzij*zin(255)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  285

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  219

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  281

                                      xin(287) = xin(287) + dxij*xin(281)
                                      yin(287) = yin(287) + dyij*yin(281)
                                      zin(287) = zin(287) + dzij*zin(281)

                                      ! i3 = i4 =  281
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  275

                                      xin(281) = xin(281) + dxij*xin(275)
                                      yin(281) = yin(281) + dyij*yin(275)
                                      zin(281) = zin(281) + dzij*zin(275)

                                      ! i3 = i4 =  275
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  269

                                      xin(275) = xin(275) + dxij*xin(269)
                                      yin(275) = yin(275) + dyij*yin(269)
                                      zin(275) = zin(275) + dzij*zin(269)

                                      ! i3 = i4 =  269
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  281

                                      xin(287) = xin(287) + dxij*xin(281)
                                      yin(287) = yin(287) + dyij*yin(281)
                                      zin(287) = zin(287) + dzij*zin(281)

                                      ! i3 = i4 =  281
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  275

                                      xin(281) = xin(281) + dxij*xin(275)
                                      yin(281) = yin(281) + dyij*yin(275)
                                      zin(281) = zin(281) + dzij*zin(275)

                                      ! i3 = i4 =  275
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  281

                                      xin(287) = xin(287) + dxij*xin(281)
                                      yin(287) = yin(287) + dyij*yin(281)
                                      zin(287) = zin(287) + dzij*zin(281)

                                      ! i3 = i4 =  281
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  203

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  203

                                      ! do ni = 1,    3

                                      xin(203) = xin(221) + dxij*xin(197)
                                      yin(203) = yin(221) + dyij*yin(197)
                                      zin(203) = zin(221) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  227

                                      ! ni =    2

                                      xin(227) = xin(245) + dxij*xin(221)
                                      yin(227) = yin(245) + dyij*yin(221)
                                      zin(227) = zin(245) + dzij*zin(221)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    3

                                      xin(251) = xin(269) + dxij*xin(245)
                                      yin(251) = yin(269) + dyij*yin(245)
                                      zin(251) = zin(269) + dzij*zin(245)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  275

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  209

                                      ! nj =    2

                                      ! i4 = i3 =  209

                                      ! do ni = 1,    3

                                      xin(209) = xin(227) + dxij*xin(203)
                                      yin(209) = yin(227) + dyij*yin(203)
                                      zin(209) = zin(227) + dzij*zin(203)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    2

                                      xin(233) = xin(251) + dxij*xin(227)
                                      yin(233) = yin(251) + dyij*yin(227)
                                      zin(233) = zin(251) + dzij*zin(227)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  257

                                      ! ni =    3

                                      xin(257) = xin(275) + dxij*xin(251)
                                      yin(257) = yin(275) + dyij*yin(251)
                                      zin(257) = zin(275) + dzij*zin(251)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  281

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  215

                                      ! nj =    3

                                      ! i4 = i3 =  215

                                      ! do ni = 1,    3

                                      xin(215) = xin(233) + dxij*xin(209)
                                      yin(215) = yin(233) + dyij*yin(209)
                                      zin(215) = zin(233) + dzij*zin(209)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                      ! ni =    2

                                      xin(239) = xin(257) + dxij*xin(233)
                                      yin(239) = yin(257) + dyij*yin(233)
                                      zin(239) = zin(257) + dzij*zin(233)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  263

                                      ! ni =    3

                                      xin(263) = xin(281) + dxij*xin(257)
                                      yin(263) = yin(281) + dyij*yin(257)
                                      zin(263) = zin(281) + dzij*zin(257)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  221

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  282

                                      xin(288) = xin(288) + dxij*xin(282)
                                      yin(288) = yin(288) + dyij*yin(282)
                                      zin(288) = zin(288) + dzij*zin(282)

                                      ! i3 = i4 =  282
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  276

                                      xin(282) = xin(282) + dxij*xin(276)
                                      yin(282) = yin(282) + dyij*yin(276)
                                      zin(282) = zin(282) + dzij*zin(276)

                                      ! i3 = i4 =  276
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  270

                                      xin(276) = xin(276) + dxij*xin(270)
                                      yin(276) = yin(276) + dyij*yin(270)
                                      zin(276) = zin(276) + dzij*zin(270)

                                      ! i3 = i4 =  270
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  282

                                      xin(288) = xin(288) + dxij*xin(282)
                                      yin(288) = yin(288) + dyij*yin(282)
                                      zin(288) = zin(288) + dzij*zin(282)

                                      ! i3 = i4 =  282
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  276

                                      xin(282) = xin(282) + dxij*xin(276)
                                      yin(282) = yin(282) + dyij*yin(276)
                                      zin(282) = zin(282) + dzij*zin(276)

                                      ! i3 = i4 =  276
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  282

                                      xin(288) = xin(288) + dxij*xin(282)
                                      yin(288) = yin(288) + dyij*yin(282)
                                      zin(288) = zin(288) + dzij*zin(282)

                                      ! i3 = i4 =  282
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  204

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  204

                                      ! do ni = 1,    3

                                      xin(204) = xin(222) + dxij*xin(198)
                                      yin(204) = yin(222) + dyij*yin(198)
                                      zin(204) = zin(222) + dzij*zin(198)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  228

                                      ! ni =    2

                                      xin(228) = xin(246) + dxij*xin(222)
                                      yin(228) = yin(246) + dyij*yin(222)
                                      zin(228) = zin(246) + dzij*zin(222)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  252

                                      ! ni =    3

                                      xin(252) = xin(270) + dxij*xin(246)
                                      yin(252) = yin(270) + dyij*yin(246)
                                      zin(252) = zin(270) + dzij*zin(246)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  276

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  210

                                      ! nj =    2

                                      ! i4 = i3 =  210

                                      ! do ni = 1,    3

                                      xin(210) = xin(228) + dxij*xin(204)
                                      yin(210) = yin(228) + dyij*yin(204)
                                      zin(210) = zin(228) + dzij*zin(204)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  234

                                      ! ni =    2

                                      xin(234) = xin(252) + dxij*xin(228)
                                      yin(234) = yin(252) + dyij*yin(228)
                                      zin(234) = zin(252) + dzij*zin(228)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  258

                                      ! ni =    3

                                      xin(258) = xin(276) + dxij*xin(252)
                                      yin(258) = yin(276) + dyij*yin(252)
                                      zin(258) = zin(276) + dzij*zin(252)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  282

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  216

                                      ! nj =    3

                                      ! i4 = i3 =  216

                                      ! do ni = 1,    3

                                      xin(216) = xin(234) + dxij*xin(210)
                                      yin(216) = yin(234) + dyij*yin(210)
                                      zin(216) = zin(234) + dzij*zin(210)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                      ! ni =    2

                                      xin(240) = xin(258) + dxij*xin(234)
                                      yin(240) = yin(258) + dyij*yin(234)
                                      zin(240) = zin(258) + dzij*zin(234)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  264

                                      ! ni =    3

                                      xin(264) = xin(282) + dxij*xin(258)
                                      yin(264) = yin(282) + dyij*yin(258)
                                      zin(264) = zin(282) + dzij*zin(258)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  222

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  193

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  198

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  197

                                      xin(198) = xin(198) + dxkl*xin(197)
                                      yin(198) = yin(198) + dykl*yin(197)
                                      zin(198) = zin(198) + dzkl*zin(197)

                                      ! i3 = i4 =  197
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  194

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  194

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  195

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  199

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  204

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  203

                                      xin(204) = xin(204) + dxkl*xin(203)
                                      yin(204) = yin(204) + dykl*yin(203)
                                      zin(204) = zin(204) + dzkl*zin(203)

                                      ! i3 = i4 =  203
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  200

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  200

                                      ! do nk = 1,    2

                                      xin(200) = xin(201) + dxkl*xin(199)
                                      yin(200) = yin(201) + dykl*yin(199)
                                      zin(200) = zin(201) + dzkl*zin(199)
                                      ! i4 = i4 + lang+1 =  202

                                      ! nk =    2

                                      xin(202) = xin(203) + dxkl*xin(201)
                                      yin(202) = yin(203) + dykl*yin(201)
                                      zin(202) = zin(203) + dzkl*zin(201)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  201

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  205

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  210

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  209

                                      xin(210) = xin(210) + dxkl*xin(209)
                                      yin(210) = yin(210) + dykl*yin(209)
                                      zin(210) = zin(210) + dzkl*zin(209)

                                      ! i3 = i4 =  209
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  206

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  206

                                      ! do nk = 1,    2

                                      xin(206) = xin(207) + dxkl*xin(205)
                                      yin(206) = yin(207) + dykl*yin(205)
                                      zin(206) = zin(207) + dzkl*zin(205)
                                      ! i4 = i4 + lang+1 =  208

                                      ! nk =    2

                                      xin(208) = xin(209) + dxkl*xin(207)
                                      yin(208) = yin(209) + dykl*yin(207)
                                      zin(208) = zin(209) + dzkl*zin(207)
                                      ! i4 = i4 + lang+1 =  210

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  207

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  211

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  212

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  212

                                      ! do nk = 1,    2

                                      xin(212) = xin(213) + dxkl*xin(211)
                                      yin(212) = yin(213) + dykl*yin(211)
                                      zin(212) = zin(213) + dzkl*zin(211)
                                      ! i4 = i4 + lang+1 =  214

                                      ! nk =    2

                                      xin(214) = xin(215) + dxkl*xin(213)
                                      yin(214) = yin(215) + dykl*yin(213)
                                      zin(214) = zin(215) + dzkl*zin(213)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  213

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  222

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  221

                                      xin(222) = xin(222) + dxkl*xin(221)
                                      yin(222) = yin(222) + dykl*yin(221)
                                      zin(222) = zin(222) + dzkl*zin(221)

                                      ! i3 = i4 =  221
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  223

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  228

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  227

                                      xin(228) = xin(228) + dxkl*xin(227)
                                      yin(228) = yin(228) + dykl*yin(227)
                                      zin(228) = zin(228) + dzkl*zin(227)

                                      ! i3 = i4 =  227
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  224

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  224

                                      ! do nk = 1,    2

                                      xin(224) = xin(225) + dxkl*xin(223)
                                      yin(224) = yin(225) + dykl*yin(223)
                                      zin(224) = zin(225) + dzkl*zin(223)
                                      ! i4 = i4 + lang+1 =  226

                                      ! nk =    2

                                      xin(226) = xin(227) + dxkl*xin(225)
                                      yin(226) = yin(227) + dykl*yin(225)
                                      zin(226) = zin(227) + dzkl*zin(225)
                                      ! i4 = i4 + lang+1 =  228

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  225

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  229

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  234

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  233

                                      xin(234) = xin(234) + dxkl*xin(233)
                                      yin(234) = yin(234) + dykl*yin(233)
                                      zin(234) = zin(234) + dzkl*zin(233)

                                      ! i3 = i4 =  233
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  230

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  230

                                      ! do nk = 1,    2

                                      xin(230) = xin(231) + dxkl*xin(229)
                                      yin(230) = yin(231) + dykl*yin(229)
                                      zin(230) = zin(231) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  232

                                      ! nk =    2

                                      xin(232) = xin(233) + dxkl*xin(231)
                                      yin(232) = yin(233) + dykl*yin(231)
                                      zin(232) = zin(233) + dzkl*zin(231)
                                      ! i4 = i4 + lang+1 =  234

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  231

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  235

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  240

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  239

                                      xin(240) = xin(240) + dxkl*xin(239)
                                      yin(240) = yin(240) + dykl*yin(239)
                                      zin(240) = zin(240) + dzkl*zin(239)

                                      ! i3 = i4 =  239
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  236

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  236

                                      ! do nk = 1,    2

                                      xin(236) = xin(237) + dxkl*xin(235)
                                      yin(236) = yin(237) + dykl*yin(235)
                                      zin(236) = zin(237) + dzkl*zin(235)
                                      ! i4 = i4 + lang+1 =  238

                                      ! nk =    2

                                      xin(238) = xin(239) + dxkl*xin(237)
                                      yin(238) = yin(239) + dykl*yin(237)
                                      zin(238) = zin(239) + dzkl*zin(237)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  237

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  241

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  246

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  245

                                      xin(246) = xin(246) + dxkl*xin(245)
                                      yin(246) = yin(246) + dykl*yin(245)
                                      zin(246) = zin(246) + dzkl*zin(245)

                                      ! i3 = i4 =  245
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  242

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  242

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  243

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  247

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  248

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  248

                                      ! do nk = 1,    2

                                      xin(248) = xin(249) + dxkl*xin(247)
                                      yin(248) = yin(249) + dykl*yin(247)
                                      zin(248) = zin(249) + dzkl*zin(247)
                                      ! i4 = i4 + lang+1 =  250

                                      ! nk =    2

                                      xin(250) = xin(251) + dxkl*xin(249)
                                      yin(250) = yin(251) + dykl*yin(249)
                                      zin(250) = zin(251) + dzkl*zin(249)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  249

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  253

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  258

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  257

                                      xin(258) = xin(258) + dxkl*xin(257)
                                      yin(258) = yin(258) + dykl*yin(257)
                                      zin(258) = zin(258) + dzkl*zin(257)

                                      ! i3 = i4 =  257
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  254

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  254

                                      ! do nk = 1,    2

                                      xin(254) = xin(255) + dxkl*xin(253)
                                      yin(254) = yin(255) + dykl*yin(253)
                                      zin(254) = zin(255) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  256

                                      ! nk =    2

                                      xin(256) = xin(257) + dxkl*xin(255)
                                      yin(256) = yin(257) + dykl*yin(255)
                                      zin(256) = zin(257) + dzkl*zin(255)
                                      ! i4 = i4 + lang+1 =  258

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  255

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  259

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  264

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  263

                                      xin(264) = xin(264) + dxkl*xin(263)
                                      yin(264) = yin(264) + dykl*yin(263)
                                      zin(264) = zin(264) + dzkl*zin(263)

                                      ! i3 = i4 =  263
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  260

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  260

                                      ! do nk = 1,    2

                                      xin(260) = xin(261) + dxkl*xin(259)
                                      yin(260) = yin(261) + dykl*yin(259)
                                      zin(260) = zin(261) + dzkl*zin(259)
                                      ! i4 = i4 + lang+1 =  262

                                      ! nk =    2

                                      xin(262) = xin(263) + dxkl*xin(261)
                                      yin(262) = yin(263) + dykl*yin(261)
                                      zin(262) = zin(263) + dzkl*zin(261)
                                      ! i4 = i4 + lang+1 =  264

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  261

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  265

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  265

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  270

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  269

                                      xin(270) = xin(270) + dxkl*xin(269)
                                      yin(270) = yin(270) + dykl*yin(269)
                                      zin(270) = zin(270) + dzkl*zin(269)

                                      ! i3 = i4 =  269
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  266

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  266

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  267

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  271

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  276

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  275

                                      xin(276) = xin(276) + dxkl*xin(275)
                                      yin(276) = yin(276) + dykl*yin(275)
                                      zin(276) = zin(276) + dzkl*zin(275)

                                      ! i3 = i4 =  275
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  272

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  272

                                      ! do nk = 1,    2

                                      xin(272) = xin(273) + dxkl*xin(271)
                                      yin(272) = yin(273) + dykl*yin(271)
                                      zin(272) = zin(273) + dzkl*zin(271)
                                      ! i4 = i4 + lang+1 =  274

                                      ! nk =    2

                                      xin(274) = xin(275) + dxkl*xin(273)
                                      yin(274) = yin(275) + dykl*yin(273)
                                      zin(274) = zin(275) + dzkl*zin(273)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  273

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  277

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  282

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  281

                                      xin(282) = xin(282) + dxkl*xin(281)
                                      yin(282) = yin(282) + dykl*yin(281)
                                      zin(282) = zin(282) + dzkl*zin(281)

                                      ! i3 = i4 =  281
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  278

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  278

                                      ! do nk = 1,    2

                                      xin(278) = xin(279) + dxkl*xin(277)
                                      yin(278) = yin(279) + dykl*yin(277)
                                      zin(278) = zin(279) + dzkl*zin(277)
                                      ! i4 = i4 + lang+1 =  280

                                      ! nk =    2

                                      xin(280) = xin(281) + dxkl*xin(279)
                                      yin(280) = yin(281) + dykl*yin(279)
                                      zin(280) = zin(281) + dzkl*zin(279)
                                      ! i4 = i4 + lang+1 =  282

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  279

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  283

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  284

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  284

                                      ! do nk = 1,    2

                                      xin(284) = xin(285) + dxkl*xin(283)
                                      yin(284) = yin(285) + dykl*yin(283)
                                      zin(284) = zin(285) + dzkl*zin(283)
                                      ! i4 = i4 + lang+1 =  286

                                      ! nk =    2

                                      xin(286) = xin(287) + dxkl*xin(285)
                                      yin(286) = yin(287) + dykl*yin(285)
                                      zin(286) = zin(287) + dzkl*zin(285)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  285

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  288

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

                                      ! i1 = in(1) =  289

                                      xin(289) = 1.0_dp
                                      yin(289) = 1.0_dp
                                      zin(289) = f00

                                      ! i2 = in(2) =  313
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(313) = xc00
                                      yin(313) = yc00
                                      zin(313) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  291

                                      xin(291) = xcp00
                                      yin(291) = ycp00
                                      zin(291) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  315
                                      ! i2 =  313

                                      xin(315) = xcp00*xin(313) + cp10
                                      yin(315) = ycp00*yin(313) + cp10
                                      zin(315) = zcp00*zin(313) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  337
                                      ! i3 =  289
                                      ! i4 =  313

                                      xin(337) = c10*xin(289) + xc00*xin(313)
                                      yin(337) = c10*yin(289) + yc00*yin(313)
                                      zin(337) = c10*zin(289) + zc00*zin(313)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  339
                                      ! i5 =  337
                                      ! i4 =  313

                                      xin(339) = xcp00*xin(337) + cp10*xin(313)
                                      yin(339) = ycp00*yin(337) + cp10*yin(313)
                                      zin(339) = zcp00*zin(337) + cp10*zin(313)

                                      ! ------------------

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  361
                                      ! i3 =  313
                                      ! i4 =  337

                                      xin(361) = c10*xin(313) + xc00*xin(337)
                                      yin(361) = c10*yin(313) + yc00*yin(337)
                                      zin(361) = c10*zin(313) + zc00*zin(337)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  363
                                      ! i5 =  361
                                      ! i4 =  337

                                      xin(363) = xcp00*xin(361) + cp10*xin(337)
                                      yin(363) = ycp00*yin(361) + cp10*yin(337)
                                      zin(363) = zcp00*zin(361) + cp10*zin(337)

                                      ! ------------------

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  361

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  367
                                      ! i3 =  337
                                      ! i4 =  361

                                      xin(367) = c10*xin(337) + xc00*xin(361)
                                      yin(367) = c10*yin(337) + yc00*yin(361)
                                      zin(367) = c10*zin(337) + zc00*zin(361)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  369
                                      ! i5 =  367
                                      ! i4 =  361

                                      xin(369) = xcp00*xin(367) + cp10*xin(361)
                                      yin(369) = ycp00*yin(367) + cp10*yin(361)
                                      zin(369) = zcp00*zin(367) + cp10*zin(361)

                                      ! ------------------

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  367

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  373
                                      ! i3 =  361
                                      ! i4 =  367

                                      xin(373) = c10*xin(361) + xc00*xin(367)
                                      yin(373) = c10*yin(361) + yc00*yin(367)
                                      zin(373) = c10*zin(361) + zc00*zin(367)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  375
                                      ! i5 =  373
                                      ! i4 =  367

                                      xin(375) = xcp00*xin(373) + cp10*xin(367)
                                      yin(375) = ycp00*yin(373) + cp10*yin(367)
                                      zin(375) = zcp00*zin(373) + cp10*zin(367)

                                      ! ------------------

                                      ! i3 = i4 =  367
                                      ! i4 = i5 =  373

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  379
                                      ! i3 =  367
                                      ! i4 =  373

                                      xin(379) = c10*xin(367) + xc00*xin(373)
                                      yin(379) = c10*yin(367) + yc00*yin(373)
                                      zin(379) = c10*zin(367) + zc00*zin(373)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  381
                                      ! i5 =  379
                                      ! i4 =  373

                                      xin(381) = xcp00*xin(379) + cp10*xin(373)
                                      yin(381) = ycp00*yin(379) + cp10*yin(373)
                                      zin(381) = zcp00*zin(379) + cp10*zin(373)

                                      ! ------------------

                                      ! i3 = i4 =  373
                                      ! i4 = i5 =  379

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  289
                                      ! i4 = i1+k2 =  291

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  293
                                      ! i3 =  289
                                      ! i4 =  291

                                      xin(293) = cp01*xin(289) + xcp00*xin(291)
                                      yin(293) = cp01*yin(289) + ycp00*yin(291)
                                      zin(293) = cp01*zin(289) + zcp00*zin(291)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  317

                                      xin(317) = xc00*xin(293) + c01*xin(291)
                                      yin(317) = yc00*yin(293) + c01*yin(291)
                                      zin(317) = zc00*zin(293) + c01*zin(291)

                                      ! ------------------

                                      ! i3 = i4 =  291
                                      ! i4 = i5 =  293

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  294
                                      ! i3 =  291
                                      ! i4 =  293

                                      xin(294) = cp01*xin(291) + xcp00*xin(293)
                                      yin(294) = cp01*yin(291) + ycp00*yin(293)
                                      zin(294) = cp01*zin(291) + zcp00*zin(293)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  318

                                      xin(318) = xc00*xin(294) + c01*xin(293)
                                      yin(318) = yc00*yin(294) + c01*yin(293)
                                      zin(318) = zc00*zin(294) + c01*zin(293)

                                      ! ------------------

                                      ! i3 = i4 =  293
                                      ! i4 = i5 =  294

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  337

                                      xin(341) = c10*xin(293) + xc00*xin(317) + c01*xin(315)
                                      yin(341) = c10*yin(293) + yc00*yin(317) + c01*yin(315)
                                      zin(341) = c10*zin(293) + zc00*zin(317) + c01*zin(315)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  361

                                      xin(365) = c10*xin(317) + xc00*xin(341) + c01*xin(339)
                                      yin(365) = c10*yin(317) + yc00*yin(341) + c01*yin(339)
                                      zin(365) = c10*zin(317) + zc00*zin(341) + c01*zin(339)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  361

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  367

                                      xin(371) = c10*xin(341) + xc00*xin(365) + c01*xin(363)
                                      yin(371) = c10*yin(341) + yc00*yin(365) + c01*yin(363)
                                      zin(371) = c10*zin(341) + zc00*zin(365) + c01*zin(363)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  367

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  373

                                      xin(377) = c10*xin(365) + xc00*xin(371) + c01*xin(369)
                                      yin(377) = c10*yin(365) + yc00*yin(371) + c01*yin(369)
                                      zin(377) = c10*zin(365) + zc00*zin(371) + c01*zin(369)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  367
                                      ! i4 = i5 =  373

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  379

                                      xin(383) = c10*xin(371) + xc00*xin(377) + c01*xin(375)
                                      yin(383) = c10*yin(371) + yc00*yin(377) + c01*yin(375)
                                      zin(383) = c10*zin(371) + zc00*zin(377) + c01*zin(375)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  373
                                      ! i4 = i5 =  379

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  337

                                      xin(342) = c10*xin(294) + xc00*xin(318) + c01*xin(317)
                                      yin(342) = c10*yin(294) + yc00*yin(318) + c01*yin(317)
                                      zin(342) = c10*zin(294) + zc00*zin(318) + c01*zin(317)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  361

                                      xin(366) = c10*xin(318) + xc00*xin(342) + c01*xin(341)
                                      yin(366) = c10*yin(318) + yc00*yin(342) + c01*yin(341)
                                      zin(366) = c10*zin(318) + zc00*zin(342) + c01*zin(341)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  361

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  367

                                      xin(372) = c10*xin(342) + xc00*xin(366) + c01*xin(365)
                                      yin(372) = c10*yin(342) + yc00*yin(366) + c01*yin(365)
                                      zin(372) = c10*zin(342) + zc00*zin(366) + c01*zin(365)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  361
                                      ! i4 = i5 =  367

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  373

                                      xin(378) = c10*xin(366) + xc00*xin(372) + c01*xin(371)
                                      yin(378) = c10*yin(366) + yc00*yin(372) + c01*yin(371)
                                      zin(378) = c10*zin(366) + zc00*zin(372) + c01*zin(371)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  367
                                      ! i4 = i5 =  373

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  379

                                      xin(384) = c10*xin(372) + xc00*xin(378) + c01*xin(377)
                                      yin(384) = c10*yin(372) + yc00*yin(378) + c01*yin(377)
                                      zin(384) = c10*zin(372) + zc00*zin(378) + c01*zin(377)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  373
                                      ! i4 = i5 =  379

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  379

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(379) = xin(379) + dxij*xin(373)
                                      yin(379) = yin(379) + dyij*yin(373)
                                      zin(379) = zin(379) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  367

                                      xin(373) = xin(373) + dxij*xin(367)
                                      yin(373) = yin(373) + dyij*yin(367)
                                      zin(373) = zin(373) + dzij*zin(367)

                                      ! i3 = i4 =  367
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  361

                                      xin(367) = xin(367) + dxij*xin(361)
                                      yin(367) = yin(367) + dyij*yin(361)
                                      zin(367) = zin(367) + dzij*zin(361)

                                      ! i3 = i4 =  361
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(379) = xin(379) + dxij*xin(373)
                                      yin(379) = yin(379) + dyij*yin(373)
                                      zin(379) = zin(379) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  367

                                      xin(373) = xin(373) + dxij*xin(367)
                                      yin(373) = yin(373) + dyij*yin(367)
                                      zin(373) = zin(373) + dzij*zin(367)

                                      ! i3 = i4 =  367
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  379

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  373

                                      xin(379) = xin(379) + dxij*xin(373)
                                      yin(379) = yin(379) + dyij*yin(373)
                                      zin(379) = zin(379) + dzij*zin(373)

                                      ! i3 = i4 =  373
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  295

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  295

                                      ! do ni = 1,    3

                                      xin(295) = xin(313) + dxij*xin(289)
                                      yin(295) = yin(313) + dyij*yin(289)
                                      zin(295) = zin(313) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  319

                                      ! ni =    2

                                      xin(319) = xin(337) + dxij*xin(313)
                                      yin(319) = yin(337) + dyij*yin(313)
                                      zin(319) = zin(337) + dzij*zin(313)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  343

                                      ! ni =    3

                                      xin(343) = xin(361) + dxij*xin(337)
                                      yin(343) = yin(361) + dyij*yin(337)
                                      zin(343) = zin(361) + dzij*zin(337)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  367

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  301

                                      ! nj =    2

                                      ! i4 = i3 =  301

                                      ! do ni = 1,    3

                                      xin(301) = xin(319) + dxij*xin(295)
                                      yin(301) = yin(319) + dyij*yin(295)
                                      zin(301) = zin(319) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  325

                                      ! ni =    2

                                      xin(325) = xin(343) + dxij*xin(319)
                                      yin(325) = yin(343) + dyij*yin(319)
                                      zin(325) = zin(343) + dzij*zin(319)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                      ! ni =    3

                                      xin(349) = xin(367) + dxij*xin(343)
                                      yin(349) = yin(367) + dyij*yin(343)
                                      zin(349) = zin(367) + dzij*zin(343)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  373

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  307

                                      ! nj =    3

                                      ! i4 = i3 =  307

                                      ! do ni = 1,    3

                                      xin(307) = xin(325) + dxij*xin(301)
                                      yin(307) = yin(325) + dyij*yin(301)
                                      zin(307) = zin(325) + dzij*zin(301)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  331

                                      ! ni =    2

                                      xin(331) = xin(349) + dxij*xin(325)
                                      yin(331) = yin(349) + dyij*yin(325)
                                      zin(331) = zin(349) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  355

                                      ! ni =    3

                                      xin(355) = xin(373) + dxij*xin(349)
                                      yin(355) = yin(373) + dyij*yin(349)
                                      zin(355) = zin(373) + dzij*zin(349)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  379

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  313

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(381) = xin(381) + dxij*xin(375)
                                      yin(381) = yin(381) + dyij*yin(375)
                                      zin(381) = zin(381) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  369

                                      xin(375) = xin(375) + dxij*xin(369)
                                      yin(375) = yin(375) + dyij*yin(369)
                                      zin(375) = zin(375) + dzij*zin(369)

                                      ! i3 = i4 =  369
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  363

                                      xin(369) = xin(369) + dxij*xin(363)
                                      yin(369) = yin(369) + dyij*yin(363)
                                      zin(369) = zin(369) + dzij*zin(363)

                                      ! i3 = i4 =  363
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(381) = xin(381) + dxij*xin(375)
                                      yin(381) = yin(381) + dyij*yin(375)
                                      zin(381) = zin(381) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  369

                                      xin(375) = xin(375) + dxij*xin(369)
                                      yin(375) = yin(375) + dyij*yin(369)
                                      zin(375) = zin(375) + dzij*zin(369)

                                      ! i3 = i4 =  369
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  381

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  375

                                      xin(381) = xin(381) + dxij*xin(375)
                                      yin(381) = yin(381) + dyij*yin(375)
                                      zin(381) = zin(381) + dzij*zin(375)

                                      ! i3 = i4 =  375
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  297

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  297

                                      ! do ni = 1,    3

                                      xin(297) = xin(315) + dxij*xin(291)
                                      yin(297) = yin(315) + dyij*yin(291)
                                      zin(297) = zin(315) + dzij*zin(291)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  321

                                      ! ni =    2

                                      xin(321) = xin(339) + dxij*xin(315)
                                      yin(321) = yin(339) + dyij*yin(315)
                                      zin(321) = zin(339) + dzij*zin(315)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  345

                                      ! ni =    3

                                      xin(345) = xin(363) + dxij*xin(339)
                                      yin(345) = yin(363) + dyij*yin(339)
                                      zin(345) = zin(363) + dzij*zin(339)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  369

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  303

                                      ! nj =    2

                                      ! i4 = i3 =  303

                                      ! do ni = 1,    3

                                      xin(303) = xin(321) + dxij*xin(297)
                                      yin(303) = yin(321) + dyij*yin(297)
                                      zin(303) = zin(321) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  327

                                      ! ni =    2

                                      xin(327) = xin(345) + dxij*xin(321)
                                      yin(327) = yin(345) + dyij*yin(321)
                                      zin(327) = zin(345) + dzij*zin(321)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  351

                                      ! ni =    3

                                      xin(351) = xin(369) + dxij*xin(345)
                                      yin(351) = yin(369) + dyij*yin(345)
                                      zin(351) = zin(369) + dzij*zin(345)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  375

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  309

                                      ! nj =    3

                                      ! i4 = i3 =  309

                                      ! do ni = 1,    3

                                      xin(309) = xin(327) + dxij*xin(303)
                                      yin(309) = yin(327) + dyij*yin(303)
                                      zin(309) = zin(327) + dzij*zin(303)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  333

                                      ! ni =    2

                                      xin(333) = xin(351) + dxij*xin(327)
                                      yin(333) = yin(351) + dyij*yin(327)
                                      zin(333) = zin(351) + dzij*zin(327)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  357

                                      ! ni =    3

                                      xin(357) = xin(375) + dxij*xin(351)
                                      yin(357) = yin(375) + dyij*yin(351)
                                      zin(357) = zin(375) + dzij*zin(351)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  381

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  315

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  377

                                      xin(383) = xin(383) + dxij*xin(377)
                                      yin(383) = yin(383) + dyij*yin(377)
                                      zin(383) = zin(383) + dzij*zin(377)

                                      ! i3 = i4 =  377
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  371

                                      xin(377) = xin(377) + dxij*xin(371)
                                      yin(377) = yin(377) + dyij*yin(371)
                                      zin(377) = zin(377) + dzij*zin(371)

                                      ! i3 = i4 =  371
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  365

                                      xin(371) = xin(371) + dxij*xin(365)
                                      yin(371) = yin(371) + dyij*yin(365)
                                      zin(371) = zin(371) + dzij*zin(365)

                                      ! i3 = i4 =  365
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  377

                                      xin(383) = xin(383) + dxij*xin(377)
                                      yin(383) = yin(383) + dyij*yin(377)
                                      zin(383) = zin(383) + dzij*zin(377)

                                      ! i3 = i4 =  377
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  371

                                      xin(377) = xin(377) + dxij*xin(371)
                                      yin(377) = yin(377) + dyij*yin(371)
                                      zin(377) = zin(377) + dzij*zin(371)

                                      ! i3 = i4 =  371
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  383

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  377

                                      xin(383) = xin(383) + dxij*xin(377)
                                      yin(383) = yin(383) + dyij*yin(377)
                                      zin(383) = zin(383) + dzij*zin(377)

                                      ! i3 = i4 =  377
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  299

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  299

                                      ! do ni = 1,    3

                                      xin(299) = xin(317) + dxij*xin(293)
                                      yin(299) = yin(317) + dyij*yin(293)
                                      zin(299) = zin(317) + dzij*zin(293)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  323

                                      ! ni =    2

                                      xin(323) = xin(341) + dxij*xin(317)
                                      yin(323) = yin(341) + dyij*yin(317)
                                      zin(323) = zin(341) + dzij*zin(317)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  347

                                      ! ni =    3

                                      xin(347) = xin(365) + dxij*xin(341)
                                      yin(347) = yin(365) + dyij*yin(341)
                                      zin(347) = zin(365) + dzij*zin(341)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  371

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  305

                                      ! nj =    2

                                      ! i4 = i3 =  305

                                      ! do ni = 1,    3

                                      xin(305) = xin(323) + dxij*xin(299)
                                      yin(305) = yin(323) + dyij*yin(299)
                                      zin(305) = zin(323) + dzij*zin(299)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  329

                                      ! ni =    2

                                      xin(329) = xin(347) + dxij*xin(323)
                                      yin(329) = yin(347) + dyij*yin(323)
                                      zin(329) = zin(347) + dzij*zin(323)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  353

                                      ! ni =    3

                                      xin(353) = xin(371) + dxij*xin(347)
                                      yin(353) = yin(371) + dyij*yin(347)
                                      zin(353) = zin(371) + dzij*zin(347)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  377

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  311

                                      ! nj =    3

                                      ! i4 = i3 =  311

                                      ! do ni = 1,    3

                                      xin(311) = xin(329) + dxij*xin(305)
                                      yin(311) = yin(329) + dyij*yin(305)
                                      zin(311) = zin(329) + dzij*zin(305)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  335

                                      ! ni =    2

                                      xin(335) = xin(353) + dxij*xin(329)
                                      yin(335) = yin(353) + dyij*yin(329)
                                      zin(335) = zin(353) + dzij*zin(329)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  359

                                      ! ni =    3

                                      xin(359) = xin(377) + dxij*xin(353)
                                      yin(359) = yin(377) + dyij*yin(353)
                                      zin(359) = zin(377) + dzij*zin(353)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  383

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  317

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  378

                                      xin(384) = xin(384) + dxij*xin(378)
                                      yin(384) = yin(384) + dyij*yin(378)
                                      zin(384) = zin(384) + dzij*zin(378)

                                      ! i3 = i4 =  378
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  372

                                      xin(378) = xin(378) + dxij*xin(372)
                                      yin(378) = yin(378) + dyij*yin(372)
                                      zin(378) = zin(378) + dzij*zin(372)

                                      ! i3 = i4 =  372
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  366

                                      xin(372) = xin(372) + dxij*xin(366)
                                      yin(372) = yin(372) + dyij*yin(366)
                                      zin(372) = zin(372) + dzij*zin(366)

                                      ! i3 = i4 =  366
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  378

                                      xin(384) = xin(384) + dxij*xin(378)
                                      yin(384) = yin(384) + dyij*yin(378)
                                      zin(384) = zin(384) + dzij*zin(378)

                                      ! i3 = i4 =  378
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  372

                                      xin(378) = xin(378) + dxij*xin(372)
                                      yin(378) = yin(378) + dyij*yin(372)
                                      zin(378) = zin(378) + dzij*zin(372)

                                      ! i3 = i4 =  372
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  384

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  378

                                      xin(384) = xin(384) + dxij*xin(378)
                                      yin(384) = yin(384) + dyij*yin(378)
                                      zin(384) = zin(384) + dzij*zin(378)

                                      ! i3 = i4 =  378
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  300

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  300

                                      ! do ni = 1,    3

                                      xin(300) = xin(318) + dxij*xin(294)
                                      yin(300) = yin(318) + dyij*yin(294)
                                      zin(300) = zin(318) + dzij*zin(294)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  324

                                      ! ni =    2

                                      xin(324) = xin(342) + dxij*xin(318)
                                      yin(324) = yin(342) + dyij*yin(318)
                                      zin(324) = zin(342) + dzij*zin(318)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  348

                                      ! ni =    3

                                      xin(348) = xin(366) + dxij*xin(342)
                                      yin(348) = yin(366) + dyij*yin(342)
                                      zin(348) = zin(366) + dzij*zin(342)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  372

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  306

                                      ! nj =    2

                                      ! i4 = i3 =  306

                                      ! do ni = 1,    3

                                      xin(306) = xin(324) + dxij*xin(300)
                                      yin(306) = yin(324) + dyij*yin(300)
                                      zin(306) = zin(324) + dzij*zin(300)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  330

                                      ! ni =    2

                                      xin(330) = xin(348) + dxij*xin(324)
                                      yin(330) = yin(348) + dyij*yin(324)
                                      zin(330) = zin(348) + dzij*zin(324)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  354

                                      ! ni =    3

                                      xin(354) = xin(372) + dxij*xin(348)
                                      yin(354) = yin(372) + dyij*yin(348)
                                      zin(354) = zin(372) + dzij*zin(348)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  378

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  312

                                      ! nj =    3

                                      ! i4 = i3 =  312

                                      ! do ni = 1,    3

                                      xin(312) = xin(330) + dxij*xin(306)
                                      yin(312) = yin(330) + dyij*yin(306)
                                      zin(312) = zin(330) + dzij*zin(306)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  336

                                      ! ni =    2

                                      xin(336) = xin(354) + dxij*xin(330)
                                      yin(336) = yin(354) + dyij*yin(330)
                                      zin(336) = zin(354) + dzij*zin(330)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  360

                                      ! ni =    3

                                      xin(360) = xin(378) + dxij*xin(354)
                                      yin(360) = yin(378) + dyij*yin(354)
                                      zin(360) = zin(378) + dzij*zin(354)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  384

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  318

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  289

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  294

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  293

                                      xin(294) = xin(294) + dxkl*xin(293)
                                      yin(294) = yin(294) + dykl*yin(293)
                                      zin(294) = zin(294) + dzkl*zin(293)

                                      ! i3 = i4 =  293
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  295

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  300

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  299

                                      xin(300) = xin(300) + dxkl*xin(299)
                                      yin(300) = yin(300) + dykl*yin(299)
                                      zin(300) = zin(300) + dzkl*zin(299)

                                      ! i3 = i4 =  299
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  296

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  296

                                      ! do nk = 1,    2

                                      xin(296) = xin(297) + dxkl*xin(295)
                                      yin(296) = yin(297) + dykl*yin(295)
                                      zin(296) = zin(297) + dzkl*zin(295)
                                      ! i4 = i4 + lang+1 =  298

                                      ! nk =    2

                                      xin(298) = xin(299) + dxkl*xin(297)
                                      yin(298) = yin(299) + dykl*yin(297)
                                      zin(298) = zin(299) + dzkl*zin(297)
                                      ! i4 = i4 + lang+1 =  300

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  297

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  301

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  306

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  305

                                      xin(306) = xin(306) + dxkl*xin(305)
                                      yin(306) = yin(306) + dykl*yin(305)
                                      zin(306) = zin(306) + dzkl*zin(305)

                                      ! i3 = i4 =  305
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  302

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  302

                                      ! do nk = 1,    2

                                      xin(302) = xin(303) + dxkl*xin(301)
                                      yin(302) = yin(303) + dykl*yin(301)
                                      zin(302) = zin(303) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  304

                                      ! nk =    2

                                      xin(304) = xin(305) + dxkl*xin(303)
                                      yin(304) = yin(305) + dykl*yin(303)
                                      zin(304) = zin(305) + dzkl*zin(303)
                                      ! i4 = i4 + lang+1 =  306

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  303

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  307

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  312

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  311

                                      xin(312) = xin(312) + dxkl*xin(311)
                                      yin(312) = yin(312) + dykl*yin(311)
                                      zin(312) = zin(312) + dzkl*zin(311)

                                      ! i3 = i4 =  311
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  308

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  308

                                      ! do nk = 1,    2

                                      xin(308) = xin(309) + dxkl*xin(307)
                                      yin(308) = yin(309) + dykl*yin(307)
                                      zin(308) = zin(309) + dzkl*zin(307)
                                      ! i4 = i4 + lang+1 =  310

                                      ! nk =    2

                                      xin(310) = xin(311) + dxkl*xin(309)
                                      yin(310) = yin(311) + dykl*yin(309)
                                      zin(310) = zin(311) + dzkl*zin(309)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  309

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  313

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  318

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  317

                                      xin(318) = xin(318) + dxkl*xin(317)
                                      yin(318) = yin(318) + dykl*yin(317)
                                      zin(318) = zin(318) + dzkl*zin(317)

                                      ! i3 = i4 =  317
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  314

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  314

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  315

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  319

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  324

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  323

                                      xin(324) = xin(324) + dxkl*xin(323)
                                      yin(324) = yin(324) + dykl*yin(323)
                                      zin(324) = zin(324) + dzkl*zin(323)

                                      ! i3 = i4 =  323
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  320

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  320

                                      ! do nk = 1,    2

                                      xin(320) = xin(321) + dxkl*xin(319)
                                      yin(320) = yin(321) + dykl*yin(319)
                                      zin(320) = zin(321) + dzkl*zin(319)
                                      ! i4 = i4 + lang+1 =  322

                                      ! nk =    2

                                      xin(322) = xin(323) + dxkl*xin(321)
                                      yin(322) = yin(323) + dykl*yin(321)
                                      zin(322) = zin(323) + dzkl*zin(321)
                                      ! i4 = i4 + lang+1 =  324

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  321

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  325

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  330

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  329

                                      xin(330) = xin(330) + dxkl*xin(329)
                                      yin(330) = yin(330) + dykl*yin(329)
                                      zin(330) = zin(330) + dzkl*zin(329)

                                      ! i3 = i4 =  329
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  326

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  326

                                      ! do nk = 1,    2

                                      xin(326) = xin(327) + dxkl*xin(325)
                                      yin(326) = yin(327) + dykl*yin(325)
                                      zin(326) = zin(327) + dzkl*zin(325)
                                      ! i4 = i4 + lang+1 =  328

                                      ! nk =    2

                                      xin(328) = xin(329) + dxkl*xin(327)
                                      yin(328) = yin(329) + dykl*yin(327)
                                      zin(328) = zin(329) + dzkl*zin(327)
                                      ! i4 = i4 + lang+1 =  330

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  327

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  331

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  336

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  335

                                      xin(336) = xin(336) + dxkl*xin(335)
                                      yin(336) = yin(336) + dykl*yin(335)
                                      zin(336) = zin(336) + dzkl*zin(335)

                                      ! i3 = i4 =  335
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  332

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  332

                                      ! do nk = 1,    2

                                      xin(332) = xin(333) + dxkl*xin(331)
                                      yin(332) = yin(333) + dykl*yin(331)
                                      zin(332) = zin(333) + dzkl*zin(331)
                                      ! i4 = i4 + lang+1 =  334

                                      ! nk =    2

                                      xin(334) = xin(335) + dxkl*xin(333)
                                      yin(334) = yin(335) + dykl*yin(333)
                                      zin(334) = zin(335) + dzkl*zin(333)
                                      ! i4 = i4 + lang+1 =  336

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  333

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  337

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  337

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  342

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  341

                                      xin(342) = xin(342) + dxkl*xin(341)
                                      yin(342) = yin(342) + dykl*yin(341)
                                      zin(342) = zin(342) + dzkl*zin(341)

                                      ! i3 = i4 =  341
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  338

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  338

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  339

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  343

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  348

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  347

                                      xin(348) = xin(348) + dxkl*xin(347)
                                      yin(348) = yin(348) + dykl*yin(347)
                                      zin(348) = zin(348) + dzkl*zin(347)

                                      ! i3 = i4 =  347
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  344

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  344

                                      ! do nk = 1,    2

                                      xin(344) = xin(345) + dxkl*xin(343)
                                      yin(344) = yin(345) + dykl*yin(343)
                                      zin(344) = zin(345) + dzkl*zin(343)
                                      ! i4 = i4 + lang+1 =  346

                                      ! nk =    2

                                      xin(346) = xin(347) + dxkl*xin(345)
                                      yin(346) = yin(347) + dykl*yin(345)
                                      zin(346) = zin(347) + dzkl*zin(345)
                                      ! i4 = i4 + lang+1 =  348

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  345

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  349

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  354

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  353

                                      xin(354) = xin(354) + dxkl*xin(353)
                                      yin(354) = yin(354) + dykl*yin(353)
                                      zin(354) = zin(354) + dzkl*zin(353)

                                      ! i3 = i4 =  353
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  350

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  350

                                      ! do nk = 1,    2

                                      xin(350) = xin(351) + dxkl*xin(349)
                                      yin(350) = yin(351) + dykl*yin(349)
                                      zin(350) = zin(351) + dzkl*zin(349)
                                      ! i4 = i4 + lang+1 =  352

                                      ! nk =    2

                                      xin(352) = xin(353) + dxkl*xin(351)
                                      yin(352) = yin(353) + dykl*yin(351)
                                      zin(352) = zin(353) + dzkl*zin(351)
                                      ! i4 = i4 + lang+1 =  354

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  351

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  355

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  356

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  356

                                      ! do nk = 1,    2

                                      xin(356) = xin(357) + dxkl*xin(355)
                                      yin(356) = yin(357) + dykl*yin(355)
                                      zin(356) = zin(357) + dzkl*zin(355)
                                      ! i4 = i4 + lang+1 =  358

                                      ! nk =    2

                                      xin(358) = xin(359) + dxkl*xin(357)
                                      yin(358) = yin(359) + dykl*yin(357)
                                      zin(358) = zin(359) + dzkl*zin(357)
                                      ! i4 = i4 + lang+1 =  360

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  357

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  361

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  366

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  365

                                      xin(366) = xin(366) + dxkl*xin(365)
                                      yin(366) = yin(366) + dykl*yin(365)
                                      zin(366) = zin(366) + dzkl*zin(365)

                                      ! i3 = i4 =  365
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  362

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  362

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  363

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  367

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  372

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  371

                                      xin(372) = xin(372) + dxkl*xin(371)
                                      yin(372) = yin(372) + dykl*yin(371)
                                      zin(372) = zin(372) + dzkl*zin(371)

                                      ! i3 = i4 =  371
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  368

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  368

                                      ! do nk = 1,    2

                                      xin(368) = xin(369) + dxkl*xin(367)
                                      yin(368) = yin(369) + dykl*yin(367)
                                      zin(368) = zin(369) + dzkl*zin(367)
                                      ! i4 = i4 + lang+1 =  370

                                      ! nk =    2

                                      xin(370) = xin(371) + dxkl*xin(369)
                                      yin(370) = yin(371) + dykl*yin(369)
                                      zin(370) = zin(371) + dzkl*zin(369)
                                      ! i4 = i4 + lang+1 =  372

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  369

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  373

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  378

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  377

                                      xin(378) = xin(378) + dxkl*xin(377)
                                      yin(378) = yin(378) + dykl*yin(377)
                                      zin(378) = zin(378) + dzkl*zin(377)

                                      ! i3 = i4 =  377
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  374

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  374

                                      ! do nk = 1,    2

                                      xin(374) = xin(375) + dxkl*xin(373)
                                      yin(374) = yin(375) + dykl*yin(373)
                                      zin(374) = zin(375) + dzkl*zin(373)
                                      ! i4 = i4 + lang+1 =  376

                                      ! nk =    2

                                      xin(376) = xin(377) + dxkl*xin(375)
                                      yin(376) = yin(377) + dykl*yin(375)
                                      zin(376) = zin(377) + dzkl*zin(375)
                                      ! i4 = i4 + lang+1 =  378

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  375

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  379

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  384

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  383

                                      xin(384) = xin(384) + dxkl*xin(383)
                                      yin(384) = yin(384) + dykl*yin(383)
                                      zin(384) = zin(384) + dzkl*zin(383)

                                      ! i3 = i4 =  383
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  380

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  380

                                      ! do nk = 1,    2

                                      xin(380) = xin(381) + dxkl*xin(379)
                                      yin(380) = yin(381) + dykl*yin(379)
                                      zin(380) = zin(381) + dzkl*zin(379)
                                      ! i4 = i4 + lang+1 =  382

                                      ! nk =    2

                                      xin(382) = xin(383) + dxkl*xin(381)
                                      yin(382) = yin(383) + dykl*yin(381)
                                      zin(382) = zin(383) + dzkl*zin(381)
                                      ! i4 = i4 + lang+1 =  384

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  381

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  385

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  385

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  384

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

                                      ! i1 = in(1) =  385

                                      xin(385) = 1.0_dp
                                      yin(385) = 1.0_dp
                                      zin(385) = f00

                                      ! i2 = in(2) =  409
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(409) = xc00
                                      yin(409) = yc00
                                      zin(409) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  387

                                      xin(387) = xcp00
                                      yin(387) = ycp00
                                      zin(387) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  411
                                      ! i2 =  409

                                      xin(411) = xcp00*xin(409) + cp10
                                      yin(411) = ycp00*yin(409) + cp10
                                      zin(411) = zcp00*zin(409) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  409

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  433
                                      ! i3 =  385
                                      ! i4 =  409

                                      xin(433) = c10*xin(385) + xc00*xin(409)
                                      yin(433) = c10*yin(385) + yc00*yin(409)
                                      zin(433) = c10*zin(385) + zc00*zin(409)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  435
                                      ! i5 =  433
                                      ! i4 =  409

                                      xin(435) = xcp00*xin(433) + cp10*xin(409)
                                      yin(435) = ycp00*yin(433) + cp10*yin(409)
                                      zin(435) = zcp00*zin(433) + cp10*zin(409)

                                      ! ------------------

                                      ! i3 = i4 =  409
                                      ! i4 = i5 =  433

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  457
                                      ! i3 =  409
                                      ! i4 =  433

                                      xin(457) = c10*xin(409) + xc00*xin(433)
                                      yin(457) = c10*yin(409) + yc00*yin(433)
                                      zin(457) = c10*zin(409) + zc00*zin(433)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  459
                                      ! i5 =  457
                                      ! i4 =  433

                                      xin(459) = xcp00*xin(457) + cp10*xin(433)
                                      yin(459) = ycp00*yin(457) + cp10*yin(433)
                                      zin(459) = zcp00*zin(457) + cp10*zin(433)

                                      ! ------------------

                                      ! i3 = i4 =  433
                                      ! i4 = i5 =  457

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  463
                                      ! i3 =  433
                                      ! i4 =  457

                                      xin(463) = c10*xin(433) + xc00*xin(457)
                                      yin(463) = c10*yin(433) + yc00*yin(457)
                                      zin(463) = c10*zin(433) + zc00*zin(457)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  465
                                      ! i5 =  463
                                      ! i4 =  457

                                      xin(465) = xcp00*xin(463) + cp10*xin(457)
                                      yin(465) = ycp00*yin(463) + cp10*yin(457)
                                      zin(465) = zcp00*zin(463) + cp10*zin(457)

                                      ! ------------------

                                      ! i3 = i4 =  457
                                      ! i4 = i5 =  463

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  469
                                      ! i3 =  457
                                      ! i4 =  463

                                      xin(469) = c10*xin(457) + xc00*xin(463)
                                      yin(469) = c10*yin(457) + yc00*yin(463)
                                      zin(469) = c10*zin(457) + zc00*zin(463)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  471
                                      ! i5 =  469
                                      ! i4 =  463

                                      xin(471) = xcp00*xin(469) + cp10*xin(463)
                                      yin(471) = ycp00*yin(469) + cp10*yin(463)
                                      zin(471) = zcp00*zin(469) + cp10*zin(463)

                                      ! ------------------

                                      ! i3 = i4 =  463
                                      ! i4 = i5 =  469

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  475
                                      ! i3 =  463
                                      ! i4 =  469

                                      xin(475) = c10*xin(463) + xc00*xin(469)
                                      yin(475) = c10*yin(463) + yc00*yin(469)
                                      zin(475) = c10*zin(463) + zc00*zin(469)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  477
                                      ! i5 =  475
                                      ! i4 =  469

                                      xin(477) = xcp00*xin(475) + cp10*xin(469)
                                      yin(477) = ycp00*yin(475) + cp10*yin(469)
                                      zin(477) = zcp00*zin(475) + cp10*zin(469)

                                      ! ------------------

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  475

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  385
                                      ! i4 = i1+k2 =  387

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  389
                                      ! i3 =  385
                                      ! i4 =  387

                                      xin(389) = cp01*xin(385) + xcp00*xin(387)
                                      yin(389) = cp01*yin(385) + ycp00*yin(387)
                                      zin(389) = cp01*zin(385) + zcp00*zin(387)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  413

                                      xin(413) = xc00*xin(389) + c01*xin(387)
                                      yin(413) = yc00*yin(389) + c01*yin(387)
                                      zin(413) = zc00*zin(389) + c01*zin(387)

                                      ! ------------------

                                      ! i3 = i4 =  387
                                      ! i4 = i5 =  389

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  390
                                      ! i3 =  387
                                      ! i4 =  389

                                      xin(390) = cp01*xin(387) + xcp00*xin(389)
                                      yin(390) = cp01*yin(387) + ycp00*yin(389)
                                      zin(390) = cp01*zin(387) + zcp00*zin(389)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  414

                                      xin(414) = xc00*xin(390) + c01*xin(389)
                                      yin(414) = yc00*yin(390) + c01*yin(389)
                                      zin(414) = zc00*zin(390) + c01*zin(389)

                                      ! ------------------

                                      ! i3 = i4 =  389
                                      ! i4 = i5 =  390

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  409

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  433

                                      xin(437) = c10*xin(389) + xc00*xin(413) + c01*xin(411)
                                      yin(437) = c10*yin(389) + yc00*yin(413) + c01*yin(411)
                                      zin(437) = c10*zin(389) + zc00*zin(413) + c01*zin(411)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  409
                                      ! i4 = i5 =  433

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  457

                                      xin(461) = c10*xin(413) + xc00*xin(437) + c01*xin(435)
                                      yin(461) = c10*yin(413) + yc00*yin(437) + c01*yin(435)
                                      zin(461) = c10*zin(413) + zc00*zin(437) + c01*zin(435)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  433
                                      ! i4 = i5 =  457

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  463

                                      xin(467) = c10*xin(437) + xc00*xin(461) + c01*xin(459)
                                      yin(467) = c10*yin(437) + yc00*yin(461) + c01*yin(459)
                                      zin(467) = c10*zin(437) + zc00*zin(461) + c01*zin(459)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  457
                                      ! i4 = i5 =  463

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  469

                                      xin(473) = c10*xin(461) + xc00*xin(467) + c01*xin(465)
                                      yin(473) = c10*yin(461) + yc00*yin(467) + c01*yin(465)
                                      zin(473) = c10*zin(461) + zc00*zin(467) + c01*zin(465)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  463
                                      ! i4 = i5 =  469

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  475

                                      xin(479) = c10*xin(467) + xc00*xin(473) + c01*xin(471)
                                      yin(479) = c10*yin(467) + yc00*yin(473) + c01*yin(471)
                                      zin(479) = c10*zin(467) + zc00*zin(473) + c01*zin(471)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  475

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  385
                                      ! i4 = i2 =  409

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  433

                                      xin(438) = c10*xin(390) + xc00*xin(414) + c01*xin(413)
                                      yin(438) = c10*yin(390) + yc00*yin(414) + c01*yin(413)
                                      zin(438) = c10*zin(390) + zc00*zin(414) + c01*zin(413)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  409
                                      ! i4 = i5 =  433

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  457

                                      xin(462) = c10*xin(414) + xc00*xin(438) + c01*xin(437)
                                      yin(462) = c10*yin(414) + yc00*yin(438) + c01*yin(437)
                                      zin(462) = c10*zin(414) + zc00*zin(438) + c01*zin(437)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  433
                                      ! i4 = i5 =  457

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  463

                                      xin(468) = c10*xin(438) + xc00*xin(462) + c01*xin(461)
                                      yin(468) = c10*yin(438) + yc00*yin(462) + c01*yin(461)
                                      zin(468) = c10*zin(438) + zc00*zin(462) + c01*zin(461)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  457
                                      ! i4 = i5 =  463

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  469

                                      xin(474) = c10*xin(462) + xc00*xin(468) + c01*xin(467)
                                      yin(474) = c10*yin(462) + yc00*yin(468) + c01*yin(467)
                                      zin(474) = c10*zin(462) + zc00*zin(468) + c01*zin(467)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  463
                                      ! i4 = i5 =  469

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  475

                                      xin(480) = c10*xin(468) + xc00*xin(474) + c01*xin(473)
                                      yin(480) = c10*yin(468) + yc00*yin(474) + c01*yin(473)
                                      zin(480) = c10*zin(468) + zc00*zin(474) + c01*zin(473)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  469
                                      ! i4 = i5 =  475

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  475

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  475

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  469

                                      xin(475) = xin(475) + dxij*xin(469)
                                      yin(475) = yin(475) + dyij*yin(469)
                                      zin(475) = zin(475) + dzij*zin(469)

                                      ! i3 = i4 =  469
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  463

                                      xin(469) = xin(469) + dxij*xin(463)
                                      yin(469) = yin(469) + dyij*yin(463)
                                      zin(469) = zin(469) + dzij*zin(463)

                                      ! i3 = i4 =  463
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  457

                                      xin(463) = xin(463) + dxij*xin(457)
                                      yin(463) = yin(463) + dyij*yin(457)
                                      zin(463) = zin(463) + dzij*zin(457)

                                      ! i3 = i4 =  457
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  475

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  469

                                      xin(475) = xin(475) + dxij*xin(469)
                                      yin(475) = yin(475) + dyij*yin(469)
                                      zin(475) = zin(475) + dzij*zin(469)

                                      ! i3 = i4 =  469
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  463

                                      xin(469) = xin(469) + dxij*xin(463)
                                      yin(469) = yin(469) + dyij*yin(463)
                                      zin(469) = zin(469) + dzij*zin(463)

                                      ! i3 = i4 =  463
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  475

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  469

                                      xin(475) = xin(475) + dxij*xin(469)
                                      yin(475) = yin(475) + dyij*yin(469)
                                      zin(475) = zin(475) + dzij*zin(469)

                                      ! i3 = i4 =  469
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  391

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  391

                                      ! do ni = 1,    3

                                      xin(391) = xin(409) + dxij*xin(385)
                                      yin(391) = yin(409) + dyij*yin(385)
                                      zin(391) = zin(409) + dzij*zin(385)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  415

                                      ! ni =    2

                                      xin(415) = xin(433) + dxij*xin(409)
                                      yin(415) = yin(433) + dyij*yin(409)
                                      zin(415) = zin(433) + dzij*zin(409)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  439

                                      ! ni =    3

                                      xin(439) = xin(457) + dxij*xin(433)
                                      yin(439) = yin(457) + dyij*yin(433)
                                      zin(439) = zin(457) + dzij*zin(433)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  463

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  397

                                      ! nj =    2

                                      ! i4 = i3 =  397

                                      ! do ni = 1,    3

                                      xin(397) = xin(415) + dxij*xin(391)
                                      yin(397) = yin(415) + dyij*yin(391)
                                      zin(397) = zin(415) + dzij*zin(391)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  421

                                      ! ni =    2

                                      xin(421) = xin(439) + dxij*xin(415)
                                      yin(421) = yin(439) + dyij*yin(415)
                                      zin(421) = zin(439) + dzij*zin(415)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  445

                                      ! ni =    3

                                      xin(445) = xin(463) + dxij*xin(439)
                                      yin(445) = yin(463) + dyij*yin(439)
                                      zin(445) = zin(463) + dzij*zin(439)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  469

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  403

                                      ! nj =    3

                                      ! i4 = i3 =  403

                                      ! do ni = 1,    3

                                      xin(403) = xin(421) + dxij*xin(397)
                                      yin(403) = yin(421) + dyij*yin(397)
                                      zin(403) = zin(421) + dzij*zin(397)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  427

                                      ! ni =    2

                                      xin(427) = xin(445) + dxij*xin(421)
                                      yin(427) = yin(445) + dyij*yin(421)
                                      zin(427) = zin(445) + dzij*zin(421)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  451

                                      ! ni =    3

                                      xin(451) = xin(469) + dxij*xin(445)
                                      yin(451) = yin(469) + dyij*yin(445)
                                      zin(451) = zin(469) + dzij*zin(445)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  475

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  409

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  477

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  471

                                      xin(477) = xin(477) + dxij*xin(471)
                                      yin(477) = yin(477) + dyij*yin(471)
                                      zin(477) = zin(477) + dzij*zin(471)

                                      ! i3 = i4 =  471
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  465

                                      xin(471) = xin(471) + dxij*xin(465)
                                      yin(471) = yin(471) + dyij*yin(465)
                                      zin(471) = zin(471) + dzij*zin(465)

                                      ! i3 = i4 =  465
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  459

                                      xin(465) = xin(465) + dxij*xin(459)
                                      yin(465) = yin(465) + dyij*yin(459)
                                      zin(465) = zin(465) + dzij*zin(459)

                                      ! i3 = i4 =  459
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  477

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  471

                                      xin(477) = xin(477) + dxij*xin(471)
                                      yin(477) = yin(477) + dyij*yin(471)
                                      zin(477) = zin(477) + dzij*zin(471)

                                      ! i3 = i4 =  471
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  465

                                      xin(471) = xin(471) + dxij*xin(465)
                                      yin(471) = yin(471) + dyij*yin(465)
                                      zin(471) = zin(471) + dzij*zin(465)

                                      ! i3 = i4 =  465
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  477

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  471

                                      xin(477) = xin(477) + dxij*xin(471)
                                      yin(477) = yin(477) + dyij*yin(471)
                                      zin(477) = zin(477) + dzij*zin(471)

                                      ! i3 = i4 =  471
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  393

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  393

                                      ! do ni = 1,    3

                                      xin(393) = xin(411) + dxij*xin(387)
                                      yin(393) = yin(411) + dyij*yin(387)
                                      zin(393) = zin(411) + dzij*zin(387)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  417

                                      ! ni =    2

                                      xin(417) = xin(435) + dxij*xin(411)
                                      yin(417) = yin(435) + dyij*yin(411)
                                      zin(417) = zin(435) + dzij*zin(411)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  441

                                      ! ni =    3

                                      xin(441) = xin(459) + dxij*xin(435)
                                      yin(441) = yin(459) + dyij*yin(435)
                                      zin(441) = zin(459) + dzij*zin(435)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  465

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  399

                                      ! nj =    2

                                      ! i4 = i3 =  399

                                      ! do ni = 1,    3

                                      xin(399) = xin(417) + dxij*xin(393)
                                      yin(399) = yin(417) + dyij*yin(393)
                                      zin(399) = zin(417) + dzij*zin(393)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  423

                                      ! ni =    2

                                      xin(423) = xin(441) + dxij*xin(417)
                                      yin(423) = yin(441) + dyij*yin(417)
                                      zin(423) = zin(441) + dzij*zin(417)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  447

                                      ! ni =    3

                                      xin(447) = xin(465) + dxij*xin(441)
                                      yin(447) = yin(465) + dyij*yin(441)
                                      zin(447) = zin(465) + dzij*zin(441)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  471

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  405

                                      ! nj =    3

                                      ! i4 = i3 =  405

                                      ! do ni = 1,    3

                                      xin(405) = xin(423) + dxij*xin(399)
                                      yin(405) = yin(423) + dyij*yin(399)
                                      zin(405) = zin(423) + dzij*zin(399)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  429

                                      ! ni =    2

                                      xin(429) = xin(447) + dxij*xin(423)
                                      yin(429) = yin(447) + dyij*yin(423)
                                      zin(429) = zin(447) + dzij*zin(423)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  453

                                      ! ni =    3

                                      xin(453) = xin(471) + dxij*xin(447)
                                      yin(453) = yin(471) + dyij*yin(447)
                                      zin(453) = zin(471) + dzij*zin(447)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  477

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  411

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  479

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  473

                                      xin(479) = xin(479) + dxij*xin(473)
                                      yin(479) = yin(479) + dyij*yin(473)
                                      zin(479) = zin(479) + dzij*zin(473)

                                      ! i3 = i4 =  473
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  467

                                      xin(473) = xin(473) + dxij*xin(467)
                                      yin(473) = yin(473) + dyij*yin(467)
                                      zin(473) = zin(473) + dzij*zin(467)

                                      ! i3 = i4 =  467
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  461

                                      xin(467) = xin(467) + dxij*xin(461)
                                      yin(467) = yin(467) + dyij*yin(461)
                                      zin(467) = zin(467) + dzij*zin(461)

                                      ! i3 = i4 =  461
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  479

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  473

                                      xin(479) = xin(479) + dxij*xin(473)
                                      yin(479) = yin(479) + dyij*yin(473)
                                      zin(479) = zin(479) + dzij*zin(473)

                                      ! i3 = i4 =  473
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  467

                                      xin(473) = xin(473) + dxij*xin(467)
                                      yin(473) = yin(473) + dyij*yin(467)
                                      zin(473) = zin(473) + dzij*zin(467)

                                      ! i3 = i4 =  467
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  479

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  473

                                      xin(479) = xin(479) + dxij*xin(473)
                                      yin(479) = yin(479) + dyij*yin(473)
                                      zin(479) = zin(479) + dzij*zin(473)

                                      ! i3 = i4 =  473
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  395

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  395

                                      ! do ni = 1,    3

                                      xin(395) = xin(413) + dxij*xin(389)
                                      yin(395) = yin(413) + dyij*yin(389)
                                      zin(395) = zin(413) + dzij*zin(389)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  419

                                      ! ni =    2

                                      xin(419) = xin(437) + dxij*xin(413)
                                      yin(419) = yin(437) + dyij*yin(413)
                                      zin(419) = zin(437) + dzij*zin(413)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  443

                                      ! ni =    3

                                      xin(443) = xin(461) + dxij*xin(437)
                                      yin(443) = yin(461) + dyij*yin(437)
                                      zin(443) = zin(461) + dzij*zin(437)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  467

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  401

                                      ! nj =    2

                                      ! i4 = i3 =  401

                                      ! do ni = 1,    3

                                      xin(401) = xin(419) + dxij*xin(395)
                                      yin(401) = yin(419) + dyij*yin(395)
                                      zin(401) = zin(419) + dzij*zin(395)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  425

                                      ! ni =    2

                                      xin(425) = xin(443) + dxij*xin(419)
                                      yin(425) = yin(443) + dyij*yin(419)
                                      zin(425) = zin(443) + dzij*zin(419)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  449

                                      ! ni =    3

                                      xin(449) = xin(467) + dxij*xin(443)
                                      yin(449) = yin(467) + dyij*yin(443)
                                      zin(449) = zin(467) + dzij*zin(443)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  473

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  407

                                      ! nj =    3

                                      ! i4 = i3 =  407

                                      ! do ni = 1,    3

                                      xin(407) = xin(425) + dxij*xin(401)
                                      yin(407) = yin(425) + dyij*yin(401)
                                      zin(407) = zin(425) + dzij*zin(401)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  431

                                      ! ni =    2

                                      xin(431) = xin(449) + dxij*xin(425)
                                      yin(431) = yin(449) + dyij*yin(425)
                                      zin(431) = zin(449) + dzij*zin(425)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  455

                                      ! ni =    3

                                      xin(455) = xin(473) + dxij*xin(449)
                                      yin(455) = yin(473) + dyij*yin(449)
                                      zin(455) = zin(473) + dzij*zin(449)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  479

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  413

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  480

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  474

                                      xin(480) = xin(480) + dxij*xin(474)
                                      yin(480) = yin(480) + dyij*yin(474)
                                      zin(480) = zin(480) + dzij*zin(474)

                                      ! i3 = i4 =  474
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  468

                                      xin(474) = xin(474) + dxij*xin(468)
                                      yin(474) = yin(474) + dyij*yin(468)
                                      zin(474) = zin(474) + dzij*zin(468)

                                      ! i3 = i4 =  468
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  462

                                      xin(468) = xin(468) + dxij*xin(462)
                                      yin(468) = yin(468) + dyij*yin(462)
                                      zin(468) = zin(468) + dzij*zin(462)

                                      ! i3 = i4 =  462
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  480

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  474

                                      xin(480) = xin(480) + dxij*xin(474)
                                      yin(480) = yin(480) + dyij*yin(474)
                                      zin(480) = zin(480) + dzij*zin(474)

                                      ! i3 = i4 =  474
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  468

                                      xin(474) = xin(474) + dxij*xin(468)
                                      yin(474) = yin(474) + dyij*yin(468)
                                      zin(474) = zin(474) + dzij*zin(468)

                                      ! i3 = i4 =  468
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  480

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  474

                                      xin(480) = xin(480) + dxij*xin(474)
                                      yin(480) = yin(480) + dyij*yin(474)
                                      zin(480) = zin(480) + dzij*zin(474)

                                      ! i3 = i4 =  474
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  396

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  396

                                      ! do ni = 1,    3

                                      xin(396) = xin(414) + dxij*xin(390)
                                      yin(396) = yin(414) + dyij*yin(390)
                                      zin(396) = zin(414) + dzij*zin(390)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  420

                                      ! ni =    2

                                      xin(420) = xin(438) + dxij*xin(414)
                                      yin(420) = yin(438) + dyij*yin(414)
                                      zin(420) = zin(438) + dzij*zin(414)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  444

                                      ! ni =    3

                                      xin(444) = xin(462) + dxij*xin(438)
                                      yin(444) = yin(462) + dyij*yin(438)
                                      zin(444) = zin(462) + dzij*zin(438)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  468

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  402

                                      ! nj =    2

                                      ! i4 = i3 =  402

                                      ! do ni = 1,    3

                                      xin(402) = xin(420) + dxij*xin(396)
                                      yin(402) = yin(420) + dyij*yin(396)
                                      zin(402) = zin(420) + dzij*zin(396)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  426

                                      ! ni =    2

                                      xin(426) = xin(444) + dxij*xin(420)
                                      yin(426) = yin(444) + dyij*yin(420)
                                      zin(426) = zin(444) + dzij*zin(420)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  450

                                      ! ni =    3

                                      xin(450) = xin(468) + dxij*xin(444)
                                      yin(450) = yin(468) + dyij*yin(444)
                                      zin(450) = zin(468) + dzij*zin(444)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  474

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  408

                                      ! nj =    3

                                      ! i4 = i3 =  408

                                      ! do ni = 1,    3

                                      xin(408) = xin(426) + dxij*xin(402)
                                      yin(408) = yin(426) + dyij*yin(402)
                                      zin(408) = zin(426) + dzij*zin(402)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  432

                                      ! ni =    2

                                      xin(432) = xin(450) + dxij*xin(426)
                                      yin(432) = yin(450) + dyij*yin(426)
                                      zin(432) = zin(450) + dzij*zin(426)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  456

                                      ! ni =    3

                                      xin(456) = xin(474) + dxij*xin(450)
                                      yin(456) = yin(474) + dyij*yin(450)
                                      zin(456) = zin(474) + dzij*zin(450)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  480

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  414

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  385

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  390

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  389

                                      xin(390) = xin(390) + dxkl*xin(389)
                                      yin(390) = yin(390) + dykl*yin(389)
                                      zin(390) = zin(390) + dzkl*zin(389)

                                      ! i3 = i4 =  389
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  386

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  386

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  387

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  391

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  396

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  395

                                      xin(396) = xin(396) + dxkl*xin(395)
                                      yin(396) = yin(396) + dykl*yin(395)
                                      zin(396) = zin(396) + dzkl*zin(395)

                                      ! i3 = i4 =  395
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  392

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  392

                                      ! do nk = 1,    2

                                      xin(392) = xin(393) + dxkl*xin(391)
                                      yin(392) = yin(393) + dykl*yin(391)
                                      zin(392) = zin(393) + dzkl*zin(391)
                                      ! i4 = i4 + lang+1 =  394

                                      ! nk =    2

                                      xin(394) = xin(395) + dxkl*xin(393)
                                      yin(394) = yin(395) + dykl*yin(393)
                                      zin(394) = zin(395) + dzkl*zin(393)
                                      ! i4 = i4 + lang+1 =  396

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  393

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  397

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  402

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  401

                                      xin(402) = xin(402) + dxkl*xin(401)
                                      yin(402) = yin(402) + dykl*yin(401)
                                      zin(402) = zin(402) + dzkl*zin(401)

                                      ! i3 = i4 =  401
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  398

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  398

                                      ! do nk = 1,    2

                                      xin(398) = xin(399) + dxkl*xin(397)
                                      yin(398) = yin(399) + dykl*yin(397)
                                      zin(398) = zin(399) + dzkl*zin(397)
                                      ! i4 = i4 + lang+1 =  400

                                      ! nk =    2

                                      xin(400) = xin(401) + dxkl*xin(399)
                                      yin(400) = yin(401) + dykl*yin(399)
                                      zin(400) = zin(401) + dzkl*zin(399)
                                      ! i4 = i4 + lang+1 =  402

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  399

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  403

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  408

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  407

                                      xin(408) = xin(408) + dxkl*xin(407)
                                      yin(408) = yin(408) + dykl*yin(407)
                                      zin(408) = zin(408) + dzkl*zin(407)

                                      ! i3 = i4 =  407
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  404

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  404

                                      ! do nk = 1,    2

                                      xin(404) = xin(405) + dxkl*xin(403)
                                      yin(404) = yin(405) + dykl*yin(403)
                                      zin(404) = zin(405) + dzkl*zin(403)
                                      ! i4 = i4 + lang+1 =  406

                                      ! nk =    2

                                      xin(406) = xin(407) + dxkl*xin(405)
                                      yin(406) = yin(407) + dykl*yin(405)
                                      zin(406) = zin(407) + dzkl*zin(405)
                                      ! i4 = i4 + lang+1 =  408

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  405

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  409

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  409

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  414

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  413

                                      xin(414) = xin(414) + dxkl*xin(413)
                                      yin(414) = yin(414) + dykl*yin(413)
                                      zin(414) = zin(414) + dzkl*zin(413)

                                      ! i3 = i4 =  413
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  410

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  410

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  411

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  415

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  420

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  419

                                      xin(420) = xin(420) + dxkl*xin(419)
                                      yin(420) = yin(420) + dykl*yin(419)
                                      zin(420) = zin(420) + dzkl*zin(419)

                                      ! i3 = i4 =  419
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  416

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  416

                                      ! do nk = 1,    2

                                      xin(416) = xin(417) + dxkl*xin(415)
                                      yin(416) = yin(417) + dykl*yin(415)
                                      zin(416) = zin(417) + dzkl*zin(415)
                                      ! i4 = i4 + lang+1 =  418

                                      ! nk =    2

                                      xin(418) = xin(419) + dxkl*xin(417)
                                      yin(418) = yin(419) + dykl*yin(417)
                                      zin(418) = zin(419) + dzkl*zin(417)
                                      ! i4 = i4 + lang+1 =  420

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  417

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  421

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  426

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  425

                                      xin(426) = xin(426) + dxkl*xin(425)
                                      yin(426) = yin(426) + dykl*yin(425)
                                      zin(426) = zin(426) + dzkl*zin(425)

                                      ! i3 = i4 =  425
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  422

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  422

                                      ! do nk = 1,    2

                                      xin(422) = xin(423) + dxkl*xin(421)
                                      yin(422) = yin(423) + dykl*yin(421)
                                      zin(422) = zin(423) + dzkl*zin(421)
                                      ! i4 = i4 + lang+1 =  424

                                      ! nk =    2

                                      xin(424) = xin(425) + dxkl*xin(423)
                                      yin(424) = yin(425) + dykl*yin(423)
                                      zin(424) = zin(425) + dzkl*zin(423)
                                      ! i4 = i4 + lang+1 =  426

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  423

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  427

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  432

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  431

                                      xin(432) = xin(432) + dxkl*xin(431)
                                      yin(432) = yin(432) + dykl*yin(431)
                                      zin(432) = zin(432) + dzkl*zin(431)

                                      ! i3 = i4 =  431
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  428

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  428

                                      ! do nk = 1,    2

                                      xin(428) = xin(429) + dxkl*xin(427)
                                      yin(428) = yin(429) + dykl*yin(427)
                                      zin(428) = zin(429) + dzkl*zin(427)
                                      ! i4 = i4 + lang+1 =  430

                                      ! nk =    2

                                      xin(430) = xin(431) + dxkl*xin(429)
                                      yin(430) = yin(431) + dykl*yin(429)
                                      zin(430) = zin(431) + dzkl*zin(429)
                                      ! i4 = i4 + lang+1 =  432

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  429

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  433

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  433

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  438

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  437

                                      xin(438) = xin(438) + dxkl*xin(437)
                                      yin(438) = yin(438) + dykl*yin(437)
                                      zin(438) = zin(438) + dzkl*zin(437)

                                      ! i3 = i4 =  437
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  434

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  434

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  435

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  439

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  444

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  443

                                      xin(444) = xin(444) + dxkl*xin(443)
                                      yin(444) = yin(444) + dykl*yin(443)
                                      zin(444) = zin(444) + dzkl*zin(443)

                                      ! i3 = i4 =  443
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  440

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  440

                                      ! do nk = 1,    2

                                      xin(440) = xin(441) + dxkl*xin(439)
                                      yin(440) = yin(441) + dykl*yin(439)
                                      zin(440) = zin(441) + dzkl*zin(439)
                                      ! i4 = i4 + lang+1 =  442

                                      ! nk =    2

                                      xin(442) = xin(443) + dxkl*xin(441)
                                      yin(442) = yin(443) + dykl*yin(441)
                                      zin(442) = zin(443) + dzkl*zin(441)
                                      ! i4 = i4 + lang+1 =  444

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  441

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  445

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  450

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  449

                                      xin(450) = xin(450) + dxkl*xin(449)
                                      yin(450) = yin(450) + dykl*yin(449)
                                      zin(450) = zin(450) + dzkl*zin(449)

                                      ! i3 = i4 =  449
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  446

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  446

                                      ! do nk = 1,    2

                                      xin(446) = xin(447) + dxkl*xin(445)
                                      yin(446) = yin(447) + dykl*yin(445)
                                      zin(446) = zin(447) + dzkl*zin(445)
                                      ! i4 = i4 + lang+1 =  448

                                      ! nk =    2

                                      xin(448) = xin(449) + dxkl*xin(447)
                                      yin(448) = yin(449) + dykl*yin(447)
                                      zin(448) = zin(449) + dzkl*zin(447)
                                      ! i4 = i4 + lang+1 =  450

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  447

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  451

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  456

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  455

                                      xin(456) = xin(456) + dxkl*xin(455)
                                      yin(456) = yin(456) + dykl*yin(455)
                                      zin(456) = zin(456) + dzkl*zin(455)

                                      ! i3 = i4 =  455
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  452

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  452

                                      ! do nk = 1,    2

                                      xin(452) = xin(453) + dxkl*xin(451)
                                      yin(452) = yin(453) + dykl*yin(451)
                                      zin(452) = zin(453) + dzkl*zin(451)
                                      ! i4 = i4 + lang+1 =  454

                                      ! nk =    2

                                      xin(454) = xin(455) + dxkl*xin(453)
                                      yin(454) = yin(455) + dykl*yin(453)
                                      zin(454) = zin(455) + dzkl*zin(453)
                                      ! i4 = i4 + lang+1 =  456

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  453

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  457

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  457

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  462

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  461

                                      xin(462) = xin(462) + dxkl*xin(461)
                                      yin(462) = yin(462) + dykl*yin(461)
                                      zin(462) = zin(462) + dzkl*zin(461)

                                      ! i3 = i4 =  461
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  458

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  458

                                      ! do nk = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + 1 =  459

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  463

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  468

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  467

                                      xin(468) = xin(468) + dxkl*xin(467)
                                      yin(468) = yin(468) + dykl*yin(467)
                                      zin(468) = zin(468) + dzkl*zin(467)

                                      ! i3 = i4 =  467
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  464

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  464

                                      ! do nk = 1,    2

                                      xin(464) = xin(465) + dxkl*xin(463)
                                      yin(464) = yin(465) + dykl*yin(463)
                                      zin(464) = zin(465) + dzkl*zin(463)
                                      ! i4 = i4 + lang+1 =  466

                                      ! nk =    2

                                      xin(466) = xin(467) + dxkl*xin(465)
                                      yin(466) = yin(467) + dykl*yin(465)
                                      zin(466) = zin(467) + dzkl*zin(465)
                                      ! i4 = i4 + lang+1 =  468

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  465

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  469

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  474

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  473

                                      xin(474) = xin(474) + dxkl*xin(473)
                                      yin(474) = yin(474) + dykl*yin(473)
                                      zin(474) = zin(474) + dzkl*zin(473)

                                      ! i3 = i4 =  473
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  470

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  470

                                      ! do nk = 1,    2

                                      xin(470) = xin(471) + dxkl*xin(469)
                                      yin(470) = yin(471) + dykl*yin(469)
                                      zin(470) = zin(471) + dzkl*zin(469)
                                      ! i4 = i4 + lang+1 =  472

                                      ! nk =    2

                                      xin(472) = xin(473) + dxkl*xin(471)
                                      yin(472) = yin(473) + dykl*yin(471)
                                      zin(472) = zin(473) + dzkl*zin(471)
                                      ! i4 = i4 + lang+1 =  474

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  471

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  475

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    3

                                      ! i3 = ib+i5 =  480

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  479

                                      xin(480) = xin(480) + dxkl*xin(479)
                                      yin(480) = yin(480) + dykl*yin(479)
                                      zin(480) = zin(480) + dzkl*zin(479)

                                      ! i3 = i4 =  479
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! end do

                                      ! i3 = ib + 1 =  476

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  476

                                      ! do nk = 1,    2

                                      xin(476) = xin(477) + dxkl*xin(475)
                                      yin(476) = yin(477) + dykl*yin(475)
                                      zin(476) = zin(477) + dzkl*zin(475)
                                      ! i4 = i4 + lang+1 =  478

                                      ! nk =    2

                                      xin(478) = xin(479) + dxkl*xin(477)
                                      yin(478) = yin(479) + dykl*yin(477)
                                      zin(478) = zin(479) + dzkl*zin(477)
                                      ! i4 = i4 + lang+1 =  480

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  477

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  481

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  481

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  480

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 1800! loop over all integrals

                                        l = n - 18*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d12ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 96)*yin(my + 96)*zin(mz + 96) & ! root  2
                                                        + xin(mx + 192)*yin(my + 192)*zin(mz + 192) & ! root  3
                                                        + xin(mx + 288)*yin(my + 288)*zin(mz + 288) & ! root  4
                                                        + xin(mx + 384)*yin(my + 384)*zin(mz + 384)) ! root  5

                                        j = int(n/18) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*180 ! Stride between functions in i

                                    do j = 1, maxj2

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*18 + ip ! Add stride between functions in j

                                      do k = 1, 6 ! # of cartesians in k

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
                              deallocate (n12ket)
                              deallocate (xint12ket)

                              end subroutine int3321
                              end submodule
