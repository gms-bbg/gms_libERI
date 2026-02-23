! The total angular momentum of this class is:           8
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3311_impl
contains
  module subroutine int3311(ff_pair, pp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, pp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n11ket(:)
    real(dp), allocatable :: xint33bra(:), xint11ket(:)
    integer(kind=int64) :: nffbra, nppket
    real(dp) :: scutffbra, scutppket, test
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
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(320), yin(320), zin(320)
    real(dp) :: eri_value(900)
    real(dp) :: d33bra(100), d11ket(9)
    integer(kind=int64) :: ix(10), jx(10), kx(3), lx(3)
    integer(kind=int64) :: iy(10), jy(10), ky(3), ly(3)
    integer(kind=int64) :: iz(10), jz(10), kz(3), lz(3)
    integer(kind=int64) :: in(7), in1(7), kn(3)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(9), kly(9), klz(9)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj, kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 17
    in1(3) = 33
    in1(4) = 49
    in1(5) = 53
    in1(6) = 57
    in1(7) = 61

    kn(1) = 0
    kn(2) = 2
    kn(3) = 3

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 1
    lx(2) = 0
    lx(3) = 0

    kx(1) = 2
    kx(2) = 0
    kx(3) = 0

    jx(1) = 12
    jx(2) = 0
    jx(3) = 0
    jx(4) = 8
    jx(5) = 8
    jx(6) = 4
    jx(7) = 0
    jx(8) = 4
    jx(9) = 0
    jx(10) = 4

    ix(1) = 49
    ix(2) = 1
    ix(3) = 1
    ix(4) = 33
    ix(5) = 33
    ix(6) = 17
    ix(7) = 1
    ix(8) = 17
    ix(9) = 1
    ix(10) = 17

    ! y-arrays

    ly(1) = 0
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0

    jy(1) = 0
    jy(2) = 12
    jy(3) = 0
    jy(4) = 4
    jy(5) = 0
    jy(6) = 8
    jy(7) = 8
    jy(8) = 0
    jy(9) = 4
    jy(10) = 4

    iy(1) = 1
    iy(2) = 49
    iy(3) = 1
    iy(4) = 17
    iy(5) = 1
    iy(6) = 33
    iy(7) = 33
    iy(8) = 1
    iy(9) = 17
    iy(10) = 17

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 12
    jz(4) = 0
    jz(5) = 4
    jz(6) = 0
    jz(7) = 4
    jz(8) = 8
    jz(9) = 8
    jz(10) = 4

    iz(1) = 1
    iz(2) = 1
    iz(3) = 49
    iz(4) = 1
    iz(5) = 17
    iz(6) = 1
    iz(7) = 17
    iz(8) = 33
    iz(9) = 33
    iz(10) = 17

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 61
    ijx(2) = 49
    ijx(3) = 49
    ijx(4) = 57
    ijx(5) = 57
    ijx(6) = 53
    ijx(7) = 49
    ijx(8) = 53
    ijx(9) = 49
    ijx(10) = 53
    ijx(11) = 13
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 9
    ijx(15) = 9
    ijx(16) = 5
    ijx(17) = 1
    ijx(18) = 5
    ijx(19) = 1
    ijx(20) = 5
    ijx(21) = 13
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 9
    ijx(25) = 9
    ijx(26) = 5
    ijx(27) = 1
    ijx(28) = 5
    ijx(29) = 1
    ijx(30) = 5
    ijx(31) = 45
    ijx(32) = 33
    ijx(33) = 33
    ijx(34) = 41
    ijx(35) = 41
    ijx(36) = 37
    ijx(37) = 33
    ijx(38) = 37
    ijx(39) = 33
    ijx(40) = 37
    ijx(41) = 45
    ijx(42) = 33
    ijx(43) = 33
    ijx(44) = 41
    ijx(45) = 41
    ijx(46) = 37
    ijx(47) = 33
    ijx(48) = 37
    ijx(49) = 33
    ijx(50) = 37
    ijx(51) = 29
    ijx(52) = 17
    ijx(53) = 17
    ijx(54) = 25
    ijx(55) = 25
    ijx(56) = 21
    ijx(57) = 17
    ijx(58) = 21
    ijx(59) = 17
    ijx(60) = 21
    ijx(61) = 13
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 9
    ijx(65) = 9
    ijx(66) = 5
    ijx(67) = 1
    ijx(68) = 5
    ijx(69) = 1
    ijx(70) = 5
    ijx(71) = 29
    ijx(72) = 17
    ijx(73) = 17
    ijx(74) = 25
    ijx(75) = 25
    ijx(76) = 21
    ijx(77) = 17
    ijx(78) = 21
    ijx(79) = 17
    ijx(80) = 21
    ijx(81) = 13
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 9
    ijx(85) = 9
    ijx(86) = 5
    ijx(87) = 1
    ijx(88) = 5
    ijx(89) = 1
    ijx(90) = 5
    ijx(91) = 29
    ijx(92) = 17
    ijx(93) = 17
    ijx(94) = 25
    ijx(95) = 25
    ijx(96) = 21
    ijx(97) = 17
    ijx(98) = 21
    ijx(99) = 17
    ijx(100) = 21

    ijy(1) = 1
    ijy(2) = 13
    ijy(3) = 1
    ijy(4) = 5
    ijy(5) = 1
    ijy(6) = 9
    ijy(7) = 9
    ijy(8) = 1
    ijy(9) = 5
    ijy(10) = 5
    ijy(11) = 49
    ijy(12) = 61
    ijy(13) = 49
    ijy(14) = 53
    ijy(15) = 49
    ijy(16) = 57
    ijy(17) = 57
    ijy(18) = 49
    ijy(19) = 53
    ijy(20) = 53
    ijy(21) = 1
    ijy(22) = 13
    ijy(23) = 1
    ijy(24) = 5
    ijy(25) = 1
    ijy(26) = 9
    ijy(27) = 9
    ijy(28) = 1
    ijy(29) = 5
    ijy(30) = 5
    ijy(31) = 17
    ijy(32) = 29
    ijy(33) = 17
    ijy(34) = 21
    ijy(35) = 17
    ijy(36) = 25
    ijy(37) = 25
    ijy(38) = 17
    ijy(39) = 21
    ijy(40) = 21
    ijy(41) = 1
    ijy(42) = 13
    ijy(43) = 1
    ijy(44) = 5
    ijy(45) = 1
    ijy(46) = 9
    ijy(47) = 9
    ijy(48) = 1
    ijy(49) = 5
    ijy(50) = 5
    ijy(51) = 33
    ijy(52) = 45
    ijy(53) = 33
    ijy(54) = 37
    ijy(55) = 33
    ijy(56) = 41
    ijy(57) = 41
    ijy(58) = 33
    ijy(59) = 37
    ijy(60) = 37
    ijy(61) = 33
    ijy(62) = 45
    ijy(63) = 33
    ijy(64) = 37
    ijy(65) = 33
    ijy(66) = 41
    ijy(67) = 41
    ijy(68) = 33
    ijy(69) = 37
    ijy(70) = 37
    ijy(71) = 1
    ijy(72) = 13
    ijy(73) = 1
    ijy(74) = 5
    ijy(75) = 1
    ijy(76) = 9
    ijy(77) = 9
    ijy(78) = 1
    ijy(79) = 5
    ijy(80) = 5
    ijy(81) = 17
    ijy(82) = 29
    ijy(83) = 17
    ijy(84) = 21
    ijy(85) = 17
    ijy(86) = 25
    ijy(87) = 25
    ijy(88) = 17
    ijy(89) = 21
    ijy(90) = 21
    ijy(91) = 17
    ijy(92) = 29
    ijy(93) = 17
    ijy(94) = 21
    ijy(95) = 17
    ijy(96) = 25
    ijy(97) = 25
    ijy(98) = 17
    ijy(99) = 21
    ijy(100) = 21

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 13
    ijz(4) = 1
    ijz(5) = 5
    ijz(6) = 1
    ijz(7) = 5
    ijz(8) = 9
    ijz(9) = 9
    ijz(10) = 5
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 13
    ijz(14) = 1
    ijz(15) = 5
    ijz(16) = 1
    ijz(17) = 5
    ijz(18) = 9
    ijz(19) = 9
    ijz(20) = 5
    ijz(21) = 49
    ijz(22) = 49
    ijz(23) = 61
    ijz(24) = 49
    ijz(25) = 53
    ijz(26) = 49
    ijz(27) = 53
    ijz(28) = 57
    ijz(29) = 57
    ijz(30) = 53
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 13
    ijz(34) = 1
    ijz(35) = 5
    ijz(36) = 1
    ijz(37) = 5
    ijz(38) = 9
    ijz(39) = 9
    ijz(40) = 5
    ijz(41) = 17
    ijz(42) = 17
    ijz(43) = 29
    ijz(44) = 17
    ijz(45) = 21
    ijz(46) = 17
    ijz(47) = 21
    ijz(48) = 25
    ijz(49) = 25
    ijz(50) = 21
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 13
    ijz(54) = 1
    ijz(55) = 5
    ijz(56) = 1
    ijz(57) = 5
    ijz(58) = 9
    ijz(59) = 9
    ijz(60) = 5
    ijz(61) = 17
    ijz(62) = 17
    ijz(63) = 29
    ijz(64) = 17
    ijz(65) = 21
    ijz(66) = 17
    ijz(67) = 21
    ijz(68) = 25
    ijz(69) = 25
    ijz(70) = 21
    ijz(71) = 33
    ijz(72) = 33
    ijz(73) = 45
    ijz(74) = 33
    ijz(75) = 37
    ijz(76) = 33
    ijz(77) = 37
    ijz(78) = 41
    ijz(79) = 41
    ijz(80) = 37
    ijz(81) = 33
    ijz(82) = 33
    ijz(83) = 45
    ijz(84) = 33
    ijz(85) = 37
    ijz(86) = 33
    ijz(87) = 37
    ijz(88) = 41
    ijz(89) = 41
    ijz(90) = 37
    ijz(91) = 17
    ijz(92) = 17
    ijz(93) = 29
    ijz(94) = 17
    ijz(95) = 21
    ijz(96) = 17
    ijz(97) = 21
    ijz(98) = 25
    ijz(99) = 25
    ijz(100) = 21

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 3
    klx(2) = 2
    klx(3) = 2
    klx(4) = 1
    klx(5) = 0
    klx(6) = 0
    klx(7) = 1
    klx(8) = 0
    klx(9) = 0

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0
    kly(4) = 2
    kly(5) = 3
    kly(6) = 2
    kly(7) = 0
    kly(8) = 1
    kly(9) = 0

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1
    klz(4) = 0
    klz(5) = 0
    klz(6) = 1
    klz(7) = 2
    klz(8) = 2
    klz(9) = 3

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11ket(res%n_p_shl*(res%n_p_shl + 1)/2))

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

    scutppket = cutoff_schwarz/maxval(pp_pair%xints)
    nppket = 0
    do ij = 1, res%n_p_shl*(res%n_p_shl + 1)/2
      if (pp_pair%xints(ij) .ge. scutppket) then
        nppket = nppket + 1
        xint11ket(nppket) = pp_pair%xints(ij)
        n11ket(nppket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*nppket) .le. nchunksize_int64) nchunksize_int64 = nffbra*nppket
    ntile = int(nffbra*nppket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*nppket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint11ket, n11ket, ff_pair, pp_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d11ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
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

              test = xint33bra(ij_tmp)*xint11ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n11ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_p_shl(ksh_tmp)
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

                  t_expon_cd = pp_pair%t_expon_ab(pp_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = pp_pair%expon_a(pp_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = pp_pair%expon_b(pp_pair%pair_loc(kl) + ket_loop)
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

                  d11ket(1) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(2) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(3) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(4) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(5) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(6) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(7) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(8) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(9) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2

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

                                      ! i2 = in(2) =   17
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(17) = xc00
                                      yin(17) = yc00
                                      zin(17) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   19
                                      ! i2 =   17

                                      xin(19) = xcp00*xin(17) + cp10
                                      yin(19) = ycp00*yin(17) + cp10
                                      zin(19) = zcp00*zin(17) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   17

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   33
                                      ! i3 =    1
                                      ! i4 =   17

                                      xin(33) = c10*xin(1) + xc00*xin(17)
                                      yin(33) = c10*yin(1) + yc00*yin(17)
                                      zin(33) = c10*zin(1) + zc00*zin(17)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   35
                                      ! i5 =   33
                                      ! i4 =   17

                                      xin(35) = xcp00*xin(33) + cp10*xin(17)
                                      yin(35) = ycp00*yin(33) + cp10*yin(17)
                                      zin(35) = zcp00*zin(33) + cp10*zin(17)

                                      ! ------------------

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   33

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =   17
                                      ! i4 =   33

                                      xin(49) = c10*xin(17) + xc00*xin(33)
                                      yin(49) = c10*yin(17) + yc00*yin(33)
                                      zin(49) = c10*zin(17) + zc00*zin(33)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   51
                                      ! i5 =   49
                                      ! i4 =   33

                                      xin(51) = xcp00*xin(49) + cp10*xin(33)
                                      yin(51) = ycp00*yin(49) + cp10*yin(33)
                                      zin(51) = zcp00*zin(49) + cp10*zin(33)

                                      ! ------------------

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   49

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   53
                                      ! i3 =   33
                                      ! i4 =   49

                                      xin(53) = c10*xin(33) + xc00*xin(49)
                                      yin(53) = c10*yin(33) + yc00*yin(49)
                                      zin(53) = c10*zin(33) + zc00*zin(49)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   55
                                      ! i5 =   53
                                      ! i4 =   49

                                      xin(55) = xcp00*xin(53) + cp10*xin(49)
                                      yin(55) = ycp00*yin(53) + cp10*yin(49)
                                      zin(55) = zcp00*zin(53) + cp10*zin(49)

                                      ! ------------------

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   53

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   57
                                      ! i3 =   49
                                      ! i4 =   53

                                      xin(57) = c10*xin(49) + xc00*xin(53)
                                      yin(57) = c10*yin(49) + yc00*yin(53)
                                      zin(57) = c10*zin(49) + zc00*zin(53)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   59
                                      ! i5 =   57
                                      ! i4 =   53

                                      xin(59) = xcp00*xin(57) + cp10*xin(53)
                                      yin(59) = ycp00*yin(57) + cp10*yin(53)
                                      zin(59) = zcp00*zin(57) + cp10*zin(53)

                                      ! ------------------

                                      ! i3 = i4 =   53
                                      ! i4 = i5 =   57

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   53
                                      ! i4 =   57

                                      xin(61) = c10*xin(53) + xc00*xin(57)
                                      yin(61) = c10*yin(53) + yc00*yin(57)
                                      zin(61) = c10*zin(53) + zc00*zin(57)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   63
                                      ! i5 =   61
                                      ! i4 =   57

                                      xin(63) = xcp00*xin(61) + cp10*xin(57)
                                      yin(63) = ycp00*yin(61) + cp10*yin(57)
                                      zin(63) = zcp00*zin(61) + cp10*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    3

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    4
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(4) = cp01*xin(1) + xcp00*xin(3)
                                      yin(4) = cp01*yin(1) + ycp00*yin(3)
                                      zin(4) = cp01*zin(1) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   20

                                      xin(20) = xc00*xin(4) + c01*xin(3)
                                      yin(20) = yc00*yin(4) + c01*yin(3)
                                      zin(20) = zc00*zin(4) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    4

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   17

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   33

                                      xin(36) = c10*xin(4) + xc00*xin(20) + c01*xin(19)
                                      yin(36) = c10*yin(4) + yc00*yin(20) + c01*yin(19)
                                      zin(36) = c10*zin(4) + zc00*zin(20) + c01*zin(19)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   33

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   49

                                      xin(52) = c10*xin(20) + xc00*xin(36) + c01*xin(35)
                                      yin(52) = c10*yin(20) + yc00*yin(36) + c01*yin(35)
                                      zin(52) = c10*zin(20) + zc00*zin(36) + c01*zin(35)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   49

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   53

                                      xin(56) = c10*xin(36) + xc00*xin(52) + c01*xin(51)
                                      yin(56) = c10*yin(36) + yc00*yin(52) + c01*yin(51)
                                      zin(56) = c10*zin(36) + zc00*zin(52) + c01*zin(51)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   53

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   57

                                      xin(60) = c10*xin(52) + xc00*xin(56) + c01*xin(55)
                                      yin(60) = c10*yin(52) + yc00*yin(56) + c01*yin(55)
                                      zin(60) = c10*zin(52) + zc00*zin(56) + c01*zin(55)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   53
                                      ! i4 = i5 =   57

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   61

                                      xin(64) = c10*xin(56) + xc00*xin(60) + c01*xin(59)
                                      yin(64) = c10*yin(56) + yc00*yin(60) + c01*yin(59)
                                      zin(64) = c10*zin(56) + zc00*zin(60) + c01*zin(59)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   61

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   61

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(61) = xin(61) + dxij*xin(57)
                                      yin(61) = yin(61) + dyij*yin(57)
                                      zin(61) = zin(61) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   53

                                      xin(57) = xin(57) + dxij*xin(53)
                                      yin(57) = yin(57) + dyij*yin(53)
                                      zin(57) = zin(57) + dzij*zin(53)

                                      ! i3 = i4 =   53
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   49

                                      xin(53) = xin(53) + dxij*xin(49)
                                      yin(53) = yin(53) + dyij*yin(49)
                                      zin(53) = zin(53) + dzij*zin(49)

                                      ! i3 = i4 =   49
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   61

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(61) = xin(61) + dxij*xin(57)
                                      yin(61) = yin(61) + dyij*yin(57)
                                      zin(61) = zin(61) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   53

                                      xin(57) = xin(57) + dxij*xin(53)
                                      yin(57) = yin(57) + dyij*yin(53)
                                      zin(57) = zin(57) + dzij*zin(53)

                                      ! i3 = i4 =   53
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   61

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(61) = xin(61) + dxij*xin(57)
                                      yin(61) = yin(61) + dyij*yin(57)
                                      zin(61) = zin(61) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(17) + dxij*xin(1)
                                      yin(5) = yin(17) + dyij*yin(1)
                                      zin(5) = zin(17) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    2

                                      xin(21) = xin(33) + dxij*xin(17)
                                      yin(21) = yin(33) + dyij*yin(17)
                                      zin(21) = zin(33) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! ni =    3

                                      xin(37) = xin(49) + dxij*xin(33)
                                      yin(37) = yin(49) + dyij*yin(33)
                                      zin(37) = zin(49) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(21) + dxij*xin(5)
                                      yin(9) = yin(21) + dyij*yin(5)
                                      zin(9) = zin(21) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! ni =    2

                                      xin(25) = xin(37) + dxij*xin(21)
                                      yin(25) = yin(37) + dyij*yin(21)
                                      zin(25) = zin(37) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    3

                                      xin(41) = xin(53) + dxij*xin(37)
                                      yin(41) = yin(53) + dyij*yin(37)
                                      zin(41) = zin(53) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    3

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(25) + dxij*xin(9)
                                      yin(13) = yin(25) + dyij*yin(9)
                                      zin(13) = zin(25) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    2

                                      xin(29) = xin(41) + dxij*xin(25)
                                      yin(29) = yin(41) + dyij*yin(25)
                                      zin(29) = zin(41) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    3

                                      xin(45) = xin(57) + dxij*xin(41)
                                      yin(45) = yin(57) + dyij*yin(41)
                                      zin(45) = zin(57) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   63

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(63) = xin(63) + dxij*xin(59)
                                      yin(63) = yin(63) + dyij*yin(59)
                                      zin(63) = zin(63) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   55

                                      xin(59) = xin(59) + dxij*xin(55)
                                      yin(59) = yin(59) + dyij*yin(55)
                                      zin(59) = zin(59) + dzij*zin(55)

                                      ! i3 = i4 =   55
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   51

                                      xin(55) = xin(55) + dxij*xin(51)
                                      yin(55) = yin(55) + dyij*yin(51)
                                      zin(55) = zin(55) + dzij*zin(51)

                                      ! i3 = i4 =   51
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   63

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(63) = xin(63) + dxij*xin(59)
                                      yin(63) = yin(63) + dyij*yin(59)
                                      zin(63) = zin(63) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   55

                                      xin(59) = xin(59) + dxij*xin(55)
                                      yin(59) = yin(59) + dyij*yin(55)
                                      zin(59) = zin(59) + dzij*zin(55)

                                      ! i3 = i4 =   55
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   63

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(63) = xin(63) + dxij*xin(59)
                                      yin(63) = yin(63) + dyij*yin(59)
                                      zin(63) = zin(63) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(19) + dxij*xin(3)
                                      yin(7) = yin(19) + dyij*yin(3)
                                      zin(7) = zin(19) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    2

                                      xin(23) = xin(35) + dxij*xin(19)
                                      yin(23) = yin(35) + dyij*yin(19)
                                      zin(23) = zin(35) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    3

                                      xin(39) = xin(51) + dxij*xin(35)
                                      yin(39) = yin(51) + dyij*yin(35)
                                      zin(39) = zin(51) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(23) + dxij*xin(7)
                                      yin(11) = yin(23) + dyij*yin(7)
                                      zin(11) = zin(23) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   27

                                      ! ni =    2

                                      xin(27) = xin(39) + dxij*xin(23)
                                      yin(27) = yin(39) + dyij*yin(23)
                                      zin(27) = zin(39) + dzij*zin(23)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    3

                                      xin(43) = xin(55) + dxij*xin(39)
                                      yin(43) = yin(55) + dyij*yin(39)
                                      zin(43) = zin(55) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    3

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    3

                                      xin(15) = xin(27) + dxij*xin(11)
                                      yin(15) = yin(27) + dyij*yin(11)
                                      zin(15) = zin(27) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    2

                                      xin(31) = xin(43) + dxij*xin(27)
                                      yin(31) = yin(43) + dyij*yin(27)
                                      zin(31) = zin(43) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    3

                                      xin(47) = xin(59) + dxij*xin(43)
                                      yin(47) = yin(59) + dyij*yin(43)
                                      zin(47) = zin(59) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   60

                                      xin(64) = xin(64) + dxij*xin(60)
                                      yin(64) = yin(64) + dyij*yin(60)
                                      zin(64) = zin(64) + dzij*zin(60)

                                      ! i3 = i4 =   60
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   56

                                      xin(60) = xin(60) + dxij*xin(56)
                                      yin(60) = yin(60) + dyij*yin(56)
                                      zin(60) = zin(60) + dzij*zin(56)

                                      ! i3 = i4 =   56
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   52

                                      xin(56) = xin(56) + dxij*xin(52)
                                      yin(56) = yin(56) + dyij*yin(52)
                                      zin(56) = zin(56) + dzij*zin(52)

                                      ! i3 = i4 =   52
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   60

                                      xin(64) = xin(64) + dxij*xin(60)
                                      yin(64) = yin(64) + dyij*yin(60)
                                      zin(64) = zin(64) + dzij*zin(60)

                                      ! i3 = i4 =   60
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   56

                                      xin(60) = xin(60) + dxij*xin(56)
                                      yin(60) = yin(60) + dyij*yin(56)
                                      zin(60) = zin(60) + dzij*zin(56)

                                      ! i3 = i4 =   56
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   60

                                      xin(64) = xin(64) + dxij*xin(60)
                                      yin(64) = yin(64) + dyij*yin(60)
                                      zin(64) = zin(64) + dzij*zin(60)

                                      ! i3 = i4 =   60
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

                                      xin(8) = xin(20) + dxij*xin(4)
                                      yin(8) = yin(20) + dyij*yin(4)
                                      zin(8) = zin(20) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    2

                                      xin(24) = xin(36) + dxij*xin(20)
                                      yin(24) = yin(36) + dyij*yin(20)
                                      zin(24) = zin(36) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    3

                                      xin(40) = xin(52) + dxij*xin(36)
                                      yin(40) = yin(52) + dyij*yin(36)
                                      zin(40) = zin(52) + dzij*zin(36)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   56

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(24) + dxij*xin(8)
                                      yin(12) = yin(24) + dyij*yin(8)
                                      zin(12) = zin(24) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   28

                                      ! ni =    2

                                      xin(28) = xin(40) + dxij*xin(24)
                                      yin(28) = yin(40) + dyij*yin(24)
                                      zin(28) = zin(40) + dzij*zin(24)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    3

                                      xin(44) = xin(56) + dxij*xin(40)
                                      yin(44) = yin(56) + dyij*yin(40)
                                      zin(44) = zin(56) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   16

                                      ! nj =    3

                                      ! i4 = i3 =   16

                                      ! do ni = 1,    3

                                      xin(16) = xin(28) + dxij*xin(12)
                                      yin(16) = yin(28) + dyij*yin(12)
                                      zin(16) = zin(28) + dzij*zin(12)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    2

                                      xin(32) = xin(44) + dxij*xin(28)
                                      yin(32) = yin(44) + dyij*yin(28)
                                      zin(32) = zin(44) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    3

                                      xin(48) = xin(60) + dxij*xin(44)
                                      yin(48) = yin(60) + dyij*yin(44)
                                      zin(48) = zin(60) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   20

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =    4

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    3

                                      xin(4) = xin(4) + dxkl*xin(3)
                                      yin(4) = yin(4) + dykl*yin(3)
                                      zin(4) = zin(4) + dzkl*zin(3)

                                      ! i3 = i4 =    3
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    1

                                      xin(2) = xin(3) + dxkl*xin(1)
                                      yin(2) = yin(3) + dykl*yin(1)
                                      zin(2) = zin(3) + dzkl*zin(1)
                                      ! i4 = i4 + lang+1 =    4

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    5

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =    8

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =    6

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    6

                                      ! do nk = 1,    1

                                      xin(6) = xin(7) + dxkl*xin(5)
                                      yin(6) = yin(7) + dykl*yin(5)
                                      zin(6) = zin(7) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =    7

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    9

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   12

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   11

                                      xin(12) = xin(12) + dxkl*xin(11)
                                      yin(12) = yin(12) + dykl*yin(11)
                                      zin(12) = zin(12) + dzkl*zin(11)

                                      ! i3 = i4 =   11
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   10

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   10

                                      ! do nk = 1,    1

                                      xin(10) = xin(11) + dxkl*xin(9)
                                      yin(10) = yin(11) + dykl*yin(9)
                                      zin(10) = zin(11) + dzkl*zin(9)
                                      ! i4 = i4 + lang+1 =   12

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   11

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   13

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   16

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   15

                                      xin(16) = xin(16) + dxkl*xin(15)
                                      yin(16) = yin(16) + dykl*yin(15)
                                      zin(16) = zin(16) + dzkl*zin(15)

                                      ! i3 = i4 =   15
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   14

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   14

                                      ! do nk = 1,    1

                                      xin(14) = xin(15) + dxkl*xin(13)
                                      yin(14) = yin(15) + dykl*yin(13)
                                      zin(14) = zin(15) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   16

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   15

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   17

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   20

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   19

                                      xin(20) = xin(20) + dxkl*xin(19)
                                      yin(20) = yin(20) + dykl*yin(19)
                                      zin(20) = zin(20) + dzkl*zin(19)

                                      ! i3 = i4 =   19
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   18

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   18

                                      ! do nk = 1,    1

                                      xin(18) = xin(19) + dxkl*xin(17)
                                      yin(18) = yin(19) + dykl*yin(17)
                                      zin(18) = zin(19) + dzkl*zin(17)
                                      ! i4 = i4 + lang+1 =   20

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   19

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   21

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   24

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   23

                                      xin(24) = xin(24) + dxkl*xin(23)
                                      yin(24) = yin(24) + dykl*yin(23)
                                      zin(24) = zin(24) + dzkl*zin(23)

                                      ! i3 = i4 =   23
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   22

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   22

                                      ! do nk = 1,    1

                                      xin(22) = xin(23) + dxkl*xin(21)
                                      yin(22) = yin(23) + dykl*yin(21)
                                      zin(22) = zin(23) + dzkl*zin(21)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   23

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   25

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   28

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   27

                                      xin(28) = xin(28) + dxkl*xin(27)
                                      yin(28) = yin(28) + dykl*yin(27)
                                      zin(28) = zin(28) + dzkl*zin(27)

                                      ! i3 = i4 =   27
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   26

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   26

                                      ! do nk = 1,    1

                                      xin(26) = xin(27) + dxkl*xin(25)
                                      yin(26) = yin(27) + dykl*yin(25)
                                      zin(26) = zin(27) + dzkl*zin(25)
                                      ! i4 = i4 + lang+1 =   28

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   27

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   29

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   32

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   31

                                      xin(32) = xin(32) + dxkl*xin(31)
                                      yin(32) = yin(32) + dykl*yin(31)
                                      zin(32) = zin(32) + dzkl*zin(31)

                                      ! i3 = i4 =   31
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   30

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   30

                                      ! do nk = 1,    1

                                      xin(30) = xin(31) + dxkl*xin(29)
                                      yin(30) = yin(31) + dykl*yin(29)
                                      zin(30) = zin(31) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   31

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   33

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   34

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   34

                                      ! do nk = 1,    1

                                      xin(34) = xin(35) + dxkl*xin(33)
                                      yin(34) = yin(35) + dykl*yin(33)
                                      zin(34) = zin(35) + dzkl*zin(33)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   35

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   37

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   40

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   39

                                      xin(40) = xin(40) + dxkl*xin(39)
                                      yin(40) = yin(40) + dykl*yin(39)
                                      zin(40) = zin(40) + dzkl*zin(39)

                                      ! i3 = i4 =   39
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   38

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   38

                                      ! do nk = 1,    1

                                      xin(38) = xin(39) + dxkl*xin(37)
                                      yin(38) = yin(39) + dykl*yin(37)
                                      zin(38) = zin(39) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   40

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   39

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   41

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   44

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   43

                                      xin(44) = xin(44) + dxkl*xin(43)
                                      yin(44) = yin(44) + dykl*yin(43)
                                      zin(44) = zin(44) + dzkl*zin(43)

                                      ! i3 = i4 =   43
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   42

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   42

                                      ! do nk = 1,    1

                                      xin(42) = xin(43) + dxkl*xin(41)
                                      yin(42) = yin(43) + dykl*yin(41)
                                      zin(42) = zin(43) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   43

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   45

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   48

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   47

                                      xin(48) = xin(48) + dxkl*xin(47)
                                      yin(48) = yin(48) + dykl*yin(47)
                                      zin(48) = zin(48) + dzkl*zin(47)

                                      ! i3 = i4 =   47
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   46

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   46

                                      ! do nk = 1,    1

                                      xin(46) = xin(47) + dxkl*xin(45)
                                      yin(46) = yin(47) + dykl*yin(45)
                                      zin(46) = zin(47) + dzkl*zin(45)
                                      ! i4 = i4 + lang+1 =   48

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   47

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   49

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   52

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   51

                                      xin(52) = xin(52) + dxkl*xin(51)
                                      yin(52) = yin(52) + dykl*yin(51)
                                      zin(52) = zin(52) + dzkl*zin(51)

                                      ! i3 = i4 =   51
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   50

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   50

                                      ! do nk = 1,    1

                                      xin(50) = xin(51) + dxkl*xin(49)
                                      yin(50) = yin(51) + dykl*yin(49)
                                      zin(50) = zin(51) + dzkl*zin(49)
                                      ! i4 = i4 + lang+1 =   52

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   51

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   53

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   56

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   55

                                      xin(56) = xin(56) + dxkl*xin(55)
                                      yin(56) = yin(56) + dykl*yin(55)
                                      zin(56) = zin(56) + dzkl*zin(55)

                                      ! i3 = i4 =   55
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   54

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   54

                                      ! do nk = 1,    1

                                      xin(54) = xin(55) + dxkl*xin(53)
                                      yin(54) = yin(55) + dykl*yin(53)
                                      zin(54) = zin(55) + dzkl*zin(53)
                                      ! i4 = i4 + lang+1 =   56

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   55

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   57

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   60

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   59

                                      xin(60) = xin(60) + dxkl*xin(59)
                                      yin(60) = yin(60) + dykl*yin(59)
                                      zin(60) = zin(60) + dzkl*zin(59)

                                      ! i3 = i4 =   59
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   58

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   58

                                      ! do nk = 1,    1

                                      xin(58) = xin(59) + dxkl*xin(57)
                                      yin(58) = yin(59) + dykl*yin(57)
                                      zin(58) = zin(59) + dzkl*zin(57)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   59

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   61

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   64

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   63

                                      xin(64) = xin(64) + dxkl*xin(63)
                                      yin(64) = yin(64) + dykl*yin(63)
                                      zin(64) = zin(64) + dzkl*zin(63)

                                      ! i3 = i4 =   63
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   62

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   62

                                      ! do nk = 1,    1

                                      xin(62) = xin(63) + dxkl*xin(61)
                                      yin(62) = yin(63) + dykl*yin(61)
                                      zin(62) = zin(63) + dzkl*zin(61)
                                      ! i4 = i4 + lang+1 =   64

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   63

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   65

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   64

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

                                      ! i1 = in(1) =   65

                                      xin(65) = 1.0_dp
                                      yin(65) = 1.0_dp
                                      zin(65) = f00

                                      ! i2 = in(2) =   81
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(81) = xc00
                                      yin(81) = yc00
                                      zin(81) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   67

                                      xin(67) = xcp00
                                      yin(67) = ycp00
                                      zin(67) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   83
                                      ! i2 =   81

                                      xin(83) = xcp00*xin(81) + cp10
                                      yin(83) = ycp00*yin(81) + cp10
                                      zin(83) = zcp00*zin(81) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   65
                                      ! i4 = i2 =   81

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   97
                                      ! i3 =   65
                                      ! i4 =   81

                                      xin(97) = c10*xin(65) + xc00*xin(81)
                                      yin(97) = c10*yin(65) + yc00*yin(81)
                                      zin(97) = c10*zin(65) + zc00*zin(81)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   99
                                      ! i5 =   97
                                      ! i4 =   81

                                      xin(99) = xcp00*xin(97) + cp10*xin(81)
                                      yin(99) = ycp00*yin(97) + cp10*yin(81)
                                      zin(99) = zcp00*zin(97) + cp10*zin(81)

                                      ! ------------------

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   97

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  113
                                      ! i3 =   81
                                      ! i4 =   97

                                      xin(113) = c10*xin(81) + xc00*xin(97)
                                      yin(113) = c10*yin(81) + yc00*yin(97)
                                      zin(113) = c10*zin(81) + zc00*zin(97)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  115
                                      ! i5 =  113
                                      ! i4 =   97

                                      xin(115) = xcp00*xin(113) + cp10*xin(97)
                                      yin(115) = ycp00*yin(113) + cp10*yin(97)
                                      zin(115) = zcp00*zin(113) + cp10*zin(97)

                                      ! ------------------

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  113

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  117
                                      ! i3 =   97
                                      ! i4 =  113

                                      xin(117) = c10*xin(97) + xc00*xin(113)
                                      yin(117) = c10*yin(97) + yc00*yin(113)
                                      zin(117) = c10*zin(97) + zc00*zin(113)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  119
                                      ! i5 =  117
                                      ! i4 =  113

                                      xin(119) = xcp00*xin(117) + cp10*xin(113)
                                      yin(119) = ycp00*yin(117) + cp10*yin(113)
                                      zin(119) = zcp00*zin(117) + cp10*zin(113)

                                      ! ------------------

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  117

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =  113
                                      ! i4 =  117

                                      xin(121) = c10*xin(113) + xc00*xin(117)
                                      yin(121) = c10*yin(113) + yc00*yin(117)
                                      zin(121) = c10*zin(113) + zc00*zin(117)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  123
                                      ! i5 =  121
                                      ! i4 =  117

                                      xin(123) = xcp00*xin(121) + cp10*xin(117)
                                      yin(123) = ycp00*yin(121) + cp10*yin(117)
                                      zin(123) = zcp00*zin(121) + cp10*zin(117)

                                      ! ------------------

                                      ! i3 = i4 =  117
                                      ! i4 = i5 =  121

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  125
                                      ! i3 =  117
                                      ! i4 =  121

                                      xin(125) = c10*xin(117) + xc00*xin(121)
                                      yin(125) = c10*yin(117) + yc00*yin(121)
                                      zin(125) = c10*zin(117) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  127
                                      ! i5 =  125
                                      ! i4 =  121

                                      xin(127) = xcp00*xin(125) + cp10*xin(121)
                                      yin(127) = ycp00*yin(125) + cp10*yin(121)
                                      zin(127) = zcp00*zin(125) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   65
                                      ! i4 = i1+k2 =   67

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   68
                                      ! i3 =   65
                                      ! i4 =   67

                                      xin(68) = cp01*xin(65) + xcp00*xin(67)
                                      yin(68) = cp01*yin(65) + ycp00*yin(67)
                                      zin(68) = cp01*zin(65) + zcp00*zin(67)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   84

                                      xin(84) = xc00*xin(68) + c01*xin(67)
                                      yin(84) = yc00*yin(68) + c01*yin(67)
                                      zin(84) = zc00*zin(68) + c01*zin(67)

                                      ! ------------------

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   68

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   65
                                      ! i4 = i2 =   81

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   97

                                      xin(100) = c10*xin(68) + xc00*xin(84) + c01*xin(83)
                                      yin(100) = c10*yin(68) + yc00*yin(84) + c01*yin(83)
                                      zin(100) = c10*zin(68) + zc00*zin(84) + c01*zin(83)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   97

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  113

                                      xin(116) = c10*xin(84) + xc00*xin(100) + c01*xin(99)
                                      yin(116) = c10*yin(84) + yc00*yin(100) + c01*yin(99)
                                      zin(116) = c10*zin(84) + zc00*zin(100) + c01*zin(99)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  113

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  117

                                      xin(120) = c10*xin(100) + xc00*xin(116) + c01*xin(115)
                                      yin(120) = c10*yin(100) + yc00*yin(116) + c01*yin(115)
                                      zin(120) = c10*zin(100) + zc00*zin(116) + c01*zin(115)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  117

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  121

                                      xin(124) = c10*xin(116) + xc00*xin(120) + c01*xin(119)
                                      yin(124) = c10*yin(116) + yc00*yin(120) + c01*yin(119)
                                      zin(124) = c10*zin(116) + zc00*zin(120) + c01*zin(119)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  117
                                      ! i4 = i5 =  121

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  125

                                      xin(128) = c10*xin(120) + xc00*xin(124) + c01*xin(123)
                                      yin(128) = c10*yin(120) + yc00*yin(124) + c01*yin(123)
                                      zin(128) = c10*zin(120) + zc00*zin(124) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  125

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  121

                                      xin(125) = xin(125) + dxij*xin(121)
                                      yin(125) = yin(125) + dyij*yin(121)
                                      zin(125) = zin(125) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  117

                                      xin(121) = xin(121) + dxij*xin(117)
                                      yin(121) = yin(121) + dyij*yin(117)
                                      zin(121) = zin(121) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  113

                                      xin(117) = xin(117) + dxij*xin(113)
                                      yin(117) = yin(117) + dyij*yin(113)
                                      zin(117) = zin(117) + dzij*zin(113)

                                      ! i3 = i4 =  113
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  121

                                      xin(125) = xin(125) + dxij*xin(121)
                                      yin(125) = yin(125) + dyij*yin(121)
                                      zin(125) = zin(125) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  117

                                      xin(121) = xin(121) + dxij*xin(117)
                                      yin(121) = yin(121) + dyij*yin(117)
                                      zin(121) = zin(121) + dzij*zin(117)

                                      ! i3 = i4 =  117
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  121

                                      xin(125) = xin(125) + dxij*xin(121)
                                      yin(125) = yin(125) + dyij*yin(121)
                                      zin(125) = zin(125) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   69

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   69

                                      ! do ni = 1,    3

                                      xin(69) = xin(81) + dxij*xin(65)
                                      yin(69) = yin(81) + dyij*yin(65)
                                      zin(69) = zin(81) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    2

                                      xin(85) = xin(97) + dxij*xin(81)
                                      yin(85) = yin(97) + dyij*yin(81)
                                      zin(85) = zin(97) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  101

                                      ! ni =    3

                                      xin(101) = xin(113) + dxij*xin(97)
                                      yin(101) = yin(113) + dyij*yin(97)
                                      zin(101) = zin(113) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   73

                                      ! nj =    2

                                      ! i4 = i3 =   73

                                      ! do ni = 1,    3

                                      xin(73) = xin(85) + dxij*xin(69)
                                      yin(73) = yin(85) + dyij*yin(69)
                                      zin(73) = zin(85) + dzij*zin(69)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    2

                                      xin(89) = xin(101) + dxij*xin(85)
                                      yin(89) = yin(101) + dyij*yin(85)
                                      zin(89) = zin(101) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    3

                                      xin(105) = xin(117) + dxij*xin(101)
                                      yin(105) = yin(117) + dyij*yin(101)
                                      zin(105) = zin(117) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   77

                                      ! nj =    3

                                      ! i4 = i3 =   77

                                      ! do ni = 1,    3

                                      xin(77) = xin(89) + dxij*xin(73)
                                      yin(77) = yin(89) + dyij*yin(73)
                                      zin(77) = zin(89) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    2

                                      xin(93) = xin(105) + dxij*xin(89)
                                      yin(93) = yin(105) + dyij*yin(89)
                                      zin(93) = zin(105) + dzij*zin(89)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! ni =    3

                                      xin(109) = xin(121) + dxij*xin(105)
                                      yin(109) = yin(121) + dyij*yin(105)
                                      zin(109) = zin(121) + dzij*zin(105)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   81

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  123

                                      xin(127) = xin(127) + dxij*xin(123)
                                      yin(127) = yin(127) + dyij*yin(123)
                                      zin(127) = zin(127) + dzij*zin(123)

                                      ! i3 = i4 =  123
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  119

                                      xin(123) = xin(123) + dxij*xin(119)
                                      yin(123) = yin(123) + dyij*yin(119)
                                      zin(123) = zin(123) + dzij*zin(119)

                                      ! i3 = i4 =  119
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  115

                                      xin(119) = xin(119) + dxij*xin(115)
                                      yin(119) = yin(119) + dyij*yin(115)
                                      zin(119) = zin(119) + dzij*zin(115)

                                      ! i3 = i4 =  115
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  123

                                      xin(127) = xin(127) + dxij*xin(123)
                                      yin(127) = yin(127) + dyij*yin(123)
                                      zin(127) = zin(127) + dzij*zin(123)

                                      ! i3 = i4 =  123
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  119

                                      xin(123) = xin(123) + dxij*xin(119)
                                      yin(123) = yin(123) + dyij*yin(119)
                                      zin(123) = zin(123) + dzij*zin(119)

                                      ! i3 = i4 =  119
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  123

                                      xin(127) = xin(127) + dxij*xin(123)
                                      yin(127) = yin(127) + dyij*yin(123)
                                      zin(127) = zin(127) + dzij*zin(123)

                                      ! i3 = i4 =  123
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   71

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   71

                                      ! do ni = 1,    3

                                      xin(71) = xin(83) + dxij*xin(67)
                                      yin(71) = yin(83) + dyij*yin(67)
                                      zin(71) = zin(83) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    2

                                      xin(87) = xin(99) + dxij*xin(83)
                                      yin(87) = yin(99) + dyij*yin(83)
                                      zin(87) = zin(99) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    3

                                      xin(103) = xin(115) + dxij*xin(99)
                                      yin(103) = yin(115) + dyij*yin(99)
                                      zin(103) = zin(115) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   75

                                      ! nj =    2

                                      ! i4 = i3 =   75

                                      ! do ni = 1,    3

                                      xin(75) = xin(87) + dxij*xin(71)
                                      yin(75) = yin(87) + dyij*yin(71)
                                      zin(75) = zin(87) + dzij*zin(71)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    2

                                      xin(91) = xin(103) + dxij*xin(87)
                                      yin(91) = yin(103) + dyij*yin(87)
                                      zin(91) = zin(103) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    3

                                      xin(107) = xin(119) + dxij*xin(103)
                                      yin(107) = yin(119) + dyij*yin(103)
                                      zin(107) = zin(119) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  123

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   79

                                      ! nj =    3

                                      ! i4 = i3 =   79

                                      ! do ni = 1,    3

                                      xin(79) = xin(91) + dxij*xin(75)
                                      yin(79) = yin(91) + dyij*yin(75)
                                      zin(79) = zin(91) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    2

                                      xin(95) = xin(107) + dxij*xin(91)
                                      yin(95) = yin(107) + dyij*yin(91)
                                      zin(95) = zin(107) + dzij*zin(91)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  111

                                      ! ni =    3

                                      xin(111) = xin(123) + dxij*xin(107)
                                      yin(111) = yin(123) + dyij*yin(107)
                                      zin(111) = zin(123) + dzij*zin(107)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   83

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  124

                                      xin(128) = xin(128) + dxij*xin(124)
                                      yin(128) = yin(128) + dyij*yin(124)
                                      zin(128) = zin(128) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  120

                                      xin(124) = xin(124) + dxij*xin(120)
                                      yin(124) = yin(124) + dyij*yin(120)
                                      zin(124) = zin(124) + dzij*zin(120)

                                      ! i3 = i4 =  120
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  116

                                      xin(120) = xin(120) + dxij*xin(116)
                                      yin(120) = yin(120) + dyij*yin(116)
                                      zin(120) = zin(120) + dzij*zin(116)

                                      ! i3 = i4 =  116
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  124

                                      xin(128) = xin(128) + dxij*xin(124)
                                      yin(128) = yin(128) + dyij*yin(124)
                                      zin(128) = zin(128) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  120

                                      xin(124) = xin(124) + dxij*xin(120)
                                      yin(124) = yin(124) + dyij*yin(120)
                                      zin(124) = zin(124) + dzij*zin(120)

                                      ! i3 = i4 =  120
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  124

                                      xin(128) = xin(128) + dxij*xin(124)
                                      yin(128) = yin(128) + dyij*yin(124)
                                      zin(128) = zin(128) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   72

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   72

                                      ! do ni = 1,    3

                                      xin(72) = xin(84) + dxij*xin(68)
                                      yin(72) = yin(84) + dyij*yin(68)
                                      zin(72) = zin(84) + dzij*zin(68)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    2

                                      xin(88) = xin(100) + dxij*xin(84)
                                      yin(88) = yin(100) + dyij*yin(84)
                                      zin(88) = zin(100) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  104

                                      ! ni =    3

                                      xin(104) = xin(116) + dxij*xin(100)
                                      yin(104) = yin(116) + dyij*yin(100)
                                      zin(104) = zin(116) + dzij*zin(100)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   76

                                      ! nj =    2

                                      ! i4 = i3 =   76

                                      ! do ni = 1,    3

                                      xin(76) = xin(88) + dxij*xin(72)
                                      yin(76) = yin(88) + dyij*yin(72)
                                      zin(76) = zin(88) + dzij*zin(72)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   92

                                      ! ni =    2

                                      xin(92) = xin(104) + dxij*xin(88)
                                      yin(92) = yin(104) + dyij*yin(88)
                                      zin(92) = zin(104) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    3

                                      xin(108) = xin(120) + dxij*xin(104)
                                      yin(108) = yin(120) + dyij*yin(104)
                                      zin(108) = zin(120) + dzij*zin(104)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  124

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   80

                                      ! nj =    3

                                      ! i4 = i3 =   80

                                      ! do ni = 1,    3

                                      xin(80) = xin(92) + dxij*xin(76)
                                      yin(80) = yin(92) + dyij*yin(76)
                                      zin(80) = zin(92) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    2

                                      xin(96) = xin(108) + dxij*xin(92)
                                      yin(96) = yin(108) + dyij*yin(92)
                                      zin(96) = zin(108) + dzij*zin(92)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  112

                                      ! ni =    3

                                      xin(112) = xin(124) + dxij*xin(108)
                                      yin(112) = yin(124) + dyij*yin(108)
                                      zin(112) = zin(124) + dzij*zin(108)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   84

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   65

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   68

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   67

                                      xin(68) = xin(68) + dxkl*xin(67)
                                      yin(68) = yin(68) + dykl*yin(67)
                                      zin(68) = zin(68) + dzkl*zin(67)

                                      ! i3 = i4 =   67
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   66

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    1

                                      xin(66) = xin(67) + dxkl*xin(65)
                                      yin(66) = yin(67) + dykl*yin(65)
                                      zin(66) = zin(67) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   69

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   70

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   70

                                      ! do nk = 1,    1

                                      xin(70) = xin(71) + dxkl*xin(69)
                                      yin(70) = yin(71) + dykl*yin(69)
                                      zin(70) = zin(71) + dzkl*zin(69)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   71

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   76

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   75

                                      xin(76) = xin(76) + dxkl*xin(75)
                                      yin(76) = yin(76) + dykl*yin(75)
                                      zin(76) = zin(76) + dzkl*zin(75)

                                      ! i3 = i4 =   75
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    1

                                      xin(74) = xin(75) + dxkl*xin(73)
                                      yin(74) = yin(75) + dykl*yin(73)
                                      zin(74) = zin(75) + dzkl*zin(73)
                                      ! i4 = i4 + lang+1 =   76

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   77

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   80

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   78

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   78

                                      ! do nk = 1,    1

                                      xin(78) = xin(79) + dxkl*xin(77)
                                      yin(78) = yin(79) + dykl*yin(77)
                                      zin(78) = zin(79) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   79

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   81

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   84

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   83

                                      xin(84) = xin(84) + dxkl*xin(83)
                                      yin(84) = yin(84) + dykl*yin(83)
                                      zin(84) = zin(84) + dzkl*zin(83)

                                      ! i3 = i4 =   83
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   82

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   82

                                      ! do nk = 1,    1

                                      xin(82) = xin(83) + dxkl*xin(81)
                                      yin(82) = yin(83) + dykl*yin(81)
                                      zin(82) = zin(83) + dzkl*zin(81)
                                      ! i4 = i4 + lang+1 =   84

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   83

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   85

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   88

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   87

                                      xin(88) = xin(88) + dxkl*xin(87)
                                      yin(88) = yin(88) + dykl*yin(87)
                                      zin(88) = zin(88) + dzkl*zin(87)

                                      ! i3 = i4 =   87
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   86

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   86

                                      ! do nk = 1,    1

                                      xin(86) = xin(87) + dxkl*xin(85)
                                      yin(86) = yin(87) + dykl*yin(85)
                                      zin(86) = zin(87) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   88

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   87

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   89

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   92

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   91

                                      xin(92) = xin(92) + dxkl*xin(91)
                                      yin(92) = yin(92) + dykl*yin(91)
                                      zin(92) = zin(92) + dzkl*zin(91)

                                      ! i3 = i4 =   91
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   90

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   90

                                      ! do nk = 1,    1

                                      xin(90) = xin(91) + dxkl*xin(89)
                                      yin(90) = yin(91) + dykl*yin(89)
                                      zin(90) = zin(91) + dzkl*zin(89)
                                      ! i4 = i4 + lang+1 =   92

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   91

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   93

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   96

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   95

                                      xin(96) = xin(96) + dxkl*xin(95)
                                      yin(96) = yin(96) + dykl*yin(95)
                                      zin(96) = zin(96) + dzkl*zin(95)

                                      ! i3 = i4 =   95
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   94

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   94

                                      ! do nk = 1,    1

                                      xin(94) = xin(95) + dxkl*xin(93)
                                      yin(94) = yin(95) + dykl*yin(93)
                                      zin(94) = zin(95) + dzkl*zin(93)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   95

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   97

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  100

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   99

                                      xin(100) = xin(100) + dxkl*xin(99)
                                      yin(100) = yin(100) + dykl*yin(99)
                                      zin(100) = zin(100) + dzkl*zin(99)

                                      ! i3 = i4 =   99
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   98

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   98

                                      ! do nk = 1,    1

                                      xin(98) = xin(99) + dxkl*xin(97)
                                      yin(98) = yin(99) + dykl*yin(97)
                                      zin(98) = zin(99) + dzkl*zin(97)
                                      ! i4 = i4 + lang+1 =  100

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   99

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  101

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  104

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  103

                                      xin(104) = xin(104) + dxkl*xin(103)
                                      yin(104) = yin(104) + dykl*yin(103)
                                      zin(104) = zin(104) + dzkl*zin(103)

                                      ! i3 = i4 =  103
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  102

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  102

                                      ! do nk = 1,    1

                                      xin(102) = xin(103) + dxkl*xin(101)
                                      yin(102) = yin(103) + dykl*yin(101)
                                      zin(102) = zin(103) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  103

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  105

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  106

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  106

                                      ! do nk = 1,    1

                                      xin(106) = xin(107) + dxkl*xin(105)
                                      yin(106) = yin(107) + dykl*yin(105)
                                      zin(106) = zin(107) + dzkl*zin(105)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  107

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  109

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  112

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  111

                                      xin(112) = xin(112) + dxkl*xin(111)
                                      yin(112) = yin(112) + dykl*yin(111)
                                      zin(112) = zin(112) + dzkl*zin(111)

                                      ! i3 = i4 =  111
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  110

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  110

                                      ! do nk = 1,    1

                                      xin(110) = xin(111) + dxkl*xin(109)
                                      yin(110) = yin(111) + dykl*yin(109)
                                      zin(110) = zin(111) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  112

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  111

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  113

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  116

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  115

                                      xin(116) = xin(116) + dxkl*xin(115)
                                      yin(116) = yin(116) + dykl*yin(115)
                                      zin(116) = zin(116) + dzkl*zin(115)

                                      ! i3 = i4 =  115
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  114

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  114

                                      ! do nk = 1,    1

                                      xin(114) = xin(115) + dxkl*xin(113)
                                      yin(114) = yin(115) + dykl*yin(113)
                                      zin(114) = zin(115) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  115

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  117

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  120

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  119

                                      xin(120) = xin(120) + dxkl*xin(119)
                                      yin(120) = yin(120) + dykl*yin(119)
                                      zin(120) = zin(120) + dzkl*zin(119)

                                      ! i3 = i4 =  119
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  118

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  118

                                      ! do nk = 1,    1

                                      xin(118) = xin(119) + dxkl*xin(117)
                                      yin(118) = yin(119) + dykl*yin(117)
                                      zin(118) = zin(119) + dzkl*zin(117)
                                      ! i4 = i4 + lang+1 =  120

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  119

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  121

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  124

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  123

                                      xin(124) = xin(124) + dxkl*xin(123)
                                      yin(124) = yin(124) + dykl*yin(123)
                                      zin(124) = zin(124) + dzkl*zin(123)

                                      ! i3 = i4 =  123
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  122

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  122

                                      ! do nk = 1,    1

                                      xin(122) = xin(123) + dxkl*xin(121)
                                      yin(122) = yin(123) + dykl*yin(121)
                                      zin(122) = zin(123) + dzkl*zin(121)
                                      ! i4 = i4 + lang+1 =  124

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  123

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  125

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  128

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  127

                                      xin(128) = xin(128) + dxkl*xin(127)
                                      yin(128) = yin(128) + dykl*yin(127)
                                      zin(128) = zin(128) + dzkl*zin(127)

                                      ! i3 = i4 =  127
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  126

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  126

                                      ! do nk = 1,    1

                                      xin(126) = xin(127) + dxkl*xin(125)
                                      yin(126) = yin(127) + dykl*yin(125)
                                      zin(126) = zin(127) + dzkl*zin(125)
                                      ! i4 = i4 + lang+1 =  128

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  127

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  129

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  128

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

                                      ! i1 = in(1) =  129

                                      xin(129) = 1.0_dp
                                      yin(129) = 1.0_dp
                                      zin(129) = f00

                                      ! i2 = in(2) =  145
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(145) = xc00
                                      yin(145) = yc00
                                      zin(145) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  131

                                      xin(131) = xcp00
                                      yin(131) = ycp00
                                      zin(131) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  147
                                      ! i2 =  145

                                      xin(147) = xcp00*xin(145) + cp10
                                      yin(147) = ycp00*yin(145) + cp10
                                      zin(147) = zcp00*zin(145) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  145

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  161
                                      ! i3 =  129
                                      ! i4 =  145

                                      xin(161) = c10*xin(129) + xc00*xin(145)
                                      yin(161) = c10*yin(129) + yc00*yin(145)
                                      zin(161) = c10*zin(129) + zc00*zin(145)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  163
                                      ! i5 =  161
                                      ! i4 =  145

                                      xin(163) = xcp00*xin(161) + cp10*xin(145)
                                      yin(163) = ycp00*yin(161) + cp10*yin(145)
                                      zin(163) = zcp00*zin(161) + cp10*zin(145)

                                      ! ------------------

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  161

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  177
                                      ! i3 =  145
                                      ! i4 =  161

                                      xin(177) = c10*xin(145) + xc00*xin(161)
                                      yin(177) = c10*yin(145) + yc00*yin(161)
                                      zin(177) = c10*zin(145) + zc00*zin(161)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  179
                                      ! i5 =  177
                                      ! i4 =  161

                                      xin(179) = xcp00*xin(177) + cp10*xin(161)
                                      yin(179) = ycp00*yin(177) + cp10*yin(161)
                                      zin(179) = zcp00*zin(177) + cp10*zin(161)

                                      ! ------------------

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  177

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  181
                                      ! i3 =  161
                                      ! i4 =  177

                                      xin(181) = c10*xin(161) + xc00*xin(177)
                                      yin(181) = c10*yin(161) + yc00*yin(177)
                                      zin(181) = c10*zin(161) + zc00*zin(177)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  183
                                      ! i5 =  181
                                      ! i4 =  177

                                      xin(183) = xcp00*xin(181) + cp10*xin(177)
                                      yin(183) = ycp00*yin(181) + cp10*yin(177)
                                      zin(183) = zcp00*zin(181) + cp10*zin(177)

                                      ! ------------------

                                      ! i3 = i4 =  177
                                      ! i4 = i5 =  181

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  185
                                      ! i3 =  177
                                      ! i4 =  181

                                      xin(185) = c10*xin(177) + xc00*xin(181)
                                      yin(185) = c10*yin(177) + yc00*yin(181)
                                      zin(185) = c10*zin(177) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  187
                                      ! i5 =  185
                                      ! i4 =  181

                                      xin(187) = xcp00*xin(185) + cp10*xin(181)
                                      yin(187) = ycp00*yin(185) + cp10*yin(181)
                                      zin(187) = zcp00*zin(185) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  189
                                      ! i3 =  181
                                      ! i4 =  185

                                      xin(189) = c10*xin(181) + xc00*xin(185)
                                      yin(189) = c10*yin(181) + yc00*yin(185)
                                      zin(189) = c10*zin(181) + zc00*zin(185)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  191
                                      ! i5 =  189
                                      ! i4 =  185

                                      xin(191) = xcp00*xin(189) + cp10*xin(185)
                                      yin(191) = ycp00*yin(189) + cp10*yin(185)
                                      zin(191) = zcp00*zin(189) + cp10*zin(185)

                                      ! ------------------

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  129
                                      ! i4 = i1+k2 =  131

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  132
                                      ! i3 =  129
                                      ! i4 =  131

                                      xin(132) = cp01*xin(129) + xcp00*xin(131)
                                      yin(132) = cp01*yin(129) + ycp00*yin(131)
                                      zin(132) = cp01*zin(129) + zcp00*zin(131)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  148

                                      xin(148) = xc00*xin(132) + c01*xin(131)
                                      yin(148) = yc00*yin(132) + c01*yin(131)
                                      zin(148) = zc00*zin(132) + c01*zin(131)

                                      ! ------------------

                                      ! i3 = i4 =  131
                                      ! i4 = i5 =  132

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  145

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  161

                                      xin(164) = c10*xin(132) + xc00*xin(148) + c01*xin(147)
                                      yin(164) = c10*yin(132) + yc00*yin(148) + c01*yin(147)
                                      zin(164) = c10*zin(132) + zc00*zin(148) + c01*zin(147)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  161

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  177

                                      xin(180) = c10*xin(148) + xc00*xin(164) + c01*xin(163)
                                      yin(180) = c10*yin(148) + yc00*yin(164) + c01*yin(163)
                                      zin(180) = c10*zin(148) + zc00*zin(164) + c01*zin(163)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  177

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  181

                                      xin(184) = c10*xin(164) + xc00*xin(180) + c01*xin(179)
                                      yin(184) = c10*yin(164) + yc00*yin(180) + c01*yin(179)
                                      zin(184) = c10*zin(164) + zc00*zin(180) + c01*zin(179)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  177
                                      ! i4 = i5 =  181

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  185

                                      xin(188) = c10*xin(180) + xc00*xin(184) + c01*xin(183)
                                      yin(188) = c10*yin(180) + yc00*yin(184) + c01*yin(183)
                                      zin(188) = c10*zin(180) + zc00*zin(184) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  189

                                      xin(192) = c10*xin(184) + xc00*xin(188) + c01*xin(187)
                                      yin(192) = c10*yin(184) + yc00*yin(188) + c01*yin(187)
                                      zin(192) = c10*zin(184) + zc00*zin(188) + c01*zin(187)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  189

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  181

                                      xin(185) = xin(185) + dxij*xin(181)
                                      yin(185) = yin(185) + dyij*yin(181)
                                      zin(185) = zin(185) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  177

                                      xin(181) = xin(181) + dxij*xin(177)
                                      yin(181) = yin(181) + dyij*yin(177)
                                      zin(181) = zin(181) + dzij*zin(177)

                                      ! i3 = i4 =  177
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  181

                                      xin(185) = xin(185) + dxij*xin(181)
                                      yin(185) = yin(185) + dyij*yin(181)
                                      zin(185) = zin(185) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  133

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  133

                                      ! do ni = 1,    3

                                      xin(133) = xin(145) + dxij*xin(129)
                                      yin(133) = yin(145) + dyij*yin(129)
                                      zin(133) = zin(145) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  149

                                      ! ni =    2

                                      xin(149) = xin(161) + dxij*xin(145)
                                      yin(149) = yin(161) + dyij*yin(145)
                                      zin(149) = zin(161) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    3

                                      xin(165) = xin(177) + dxij*xin(161)
                                      yin(165) = yin(177) + dyij*yin(161)
                                      zin(165) = zin(177) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  137

                                      ! nj =    2

                                      ! i4 = i3 =  137

                                      ! do ni = 1,    3

                                      xin(137) = xin(149) + dxij*xin(133)
                                      yin(137) = yin(149) + dyij*yin(133)
                                      zin(137) = zin(149) + dzij*zin(133)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  153

                                      ! ni =    2

                                      xin(153) = xin(165) + dxij*xin(149)
                                      yin(153) = yin(165) + dyij*yin(149)
                                      zin(153) = zin(165) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! ni =    3

                                      xin(169) = xin(181) + dxij*xin(165)
                                      yin(169) = yin(181) + dyij*yin(165)
                                      zin(169) = zin(181) + dzij*zin(165)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  141

                                      ! nj =    3

                                      ! i4 = i3 =  141

                                      ! do ni = 1,    3

                                      xin(141) = xin(153) + dxij*xin(137)
                                      yin(141) = yin(153) + dyij*yin(137)
                                      zin(141) = zin(153) + dzij*zin(137)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! ni =    2

                                      xin(157) = xin(169) + dxij*xin(153)
                                      yin(157) = yin(169) + dyij*yin(153)
                                      zin(157) = zin(169) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    3

                                      xin(173) = xin(185) + dxij*xin(169)
                                      yin(173) = yin(185) + dyij*yin(169)
                                      zin(173) = zin(185) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  145

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  183

                                      xin(187) = xin(187) + dxij*xin(183)
                                      yin(187) = yin(187) + dyij*yin(183)
                                      zin(187) = zin(187) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  179

                                      xin(183) = xin(183) + dxij*xin(179)
                                      yin(183) = yin(183) + dyij*yin(179)
                                      zin(183) = zin(183) + dzij*zin(179)

                                      ! i3 = i4 =  179
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  183

                                      xin(187) = xin(187) + dxij*xin(183)
                                      yin(187) = yin(187) + dyij*yin(183)
                                      zin(187) = zin(187) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  135

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  135

                                      ! do ni = 1,    3

                                      xin(135) = xin(147) + dxij*xin(131)
                                      yin(135) = yin(147) + dyij*yin(131)
                                      zin(135) = zin(147) + dzij*zin(131)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  151

                                      ! ni =    2

                                      xin(151) = xin(163) + dxij*xin(147)
                                      yin(151) = yin(163) + dyij*yin(147)
                                      zin(151) = zin(163) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    3

                                      xin(167) = xin(179) + dxij*xin(163)
                                      yin(167) = yin(179) + dyij*yin(163)
                                      zin(167) = zin(179) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  183

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  139

                                      ! nj =    2

                                      ! i4 = i3 =  139

                                      ! do ni = 1,    3

                                      xin(139) = xin(151) + dxij*xin(135)
                                      yin(139) = yin(151) + dyij*yin(135)
                                      zin(139) = zin(151) + dzij*zin(135)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  155

                                      ! ni =    2

                                      xin(155) = xin(167) + dxij*xin(151)
                                      yin(155) = yin(167) + dyij*yin(151)
                                      zin(155) = zin(167) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  171

                                      ! ni =    3

                                      xin(171) = xin(183) + dxij*xin(167)
                                      yin(171) = yin(183) + dyij*yin(167)
                                      zin(171) = zin(183) + dzij*zin(167)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  143

                                      ! nj =    3

                                      ! i4 = i3 =  143

                                      ! do ni = 1,    3

                                      xin(143) = xin(155) + dxij*xin(139)
                                      yin(143) = yin(155) + dyij*yin(139)
                                      zin(143) = zin(155) + dzij*zin(139)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  159

                                      ! ni =    2

                                      xin(159) = xin(171) + dxij*xin(155)
                                      yin(159) = yin(171) + dyij*yin(155)
                                      zin(159) = zin(171) + dzij*zin(155)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(187) + dxij*xin(171)
                                      yin(175) = yin(187) + dyij*yin(171)
                                      zin(175) = zin(187) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  147

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  184

                                      xin(188) = xin(188) + dxij*xin(184)
                                      yin(188) = yin(188) + dyij*yin(184)
                                      zin(188) = zin(188) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  180

                                      xin(184) = xin(184) + dxij*xin(180)
                                      yin(184) = yin(184) + dyij*yin(180)
                                      zin(184) = zin(184) + dzij*zin(180)

                                      ! i3 = i4 =  180
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  184

                                      xin(188) = xin(188) + dxij*xin(184)
                                      yin(188) = yin(188) + dyij*yin(184)
                                      zin(188) = zin(188) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  136

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  136

                                      ! do ni = 1,    3

                                      xin(136) = xin(148) + dxij*xin(132)
                                      yin(136) = yin(148) + dyij*yin(132)
                                      zin(136) = zin(148) + dzij*zin(132)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  152

                                      ! ni =    2

                                      xin(152) = xin(164) + dxij*xin(148)
                                      yin(152) = yin(164) + dyij*yin(148)
                                      zin(152) = zin(164) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    3

                                      xin(168) = xin(180) + dxij*xin(164)
                                      yin(168) = yin(180) + dyij*yin(164)
                                      zin(168) = zin(180) + dzij*zin(164)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  184

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  140

                                      ! nj =    2

                                      ! i4 = i3 =  140

                                      ! do ni = 1,    3

                                      xin(140) = xin(152) + dxij*xin(136)
                                      yin(140) = yin(152) + dyij*yin(136)
                                      zin(140) = zin(152) + dzij*zin(136)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  156

                                      ! ni =    2

                                      xin(156) = xin(168) + dxij*xin(152)
                                      yin(156) = yin(168) + dyij*yin(152)
                                      zin(156) = zin(168) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  172

                                      ! ni =    3

                                      xin(172) = xin(184) + dxij*xin(168)
                                      yin(172) = yin(184) + dyij*yin(168)
                                      zin(172) = zin(184) + dzij*zin(168)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  188

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  144

                                      ! nj =    3

                                      ! i4 = i3 =  144

                                      ! do ni = 1,    3

                                      xin(144) = xin(156) + dxij*xin(140)
                                      yin(144) = yin(156) + dyij*yin(140)
                                      zin(144) = zin(156) + dzij*zin(140)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  160

                                      ! ni =    2

                                      xin(160) = xin(172) + dxij*xin(156)
                                      yin(160) = yin(172) + dyij*yin(156)
                                      zin(160) = zin(172) + dzij*zin(156)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  176

                                      ! ni =    3

                                      xin(176) = xin(188) + dxij*xin(172)
                                      yin(176) = yin(188) + dyij*yin(172)
                                      zin(176) = zin(188) + dzij*zin(172)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  148

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =  129

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  132

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  131

                                      xin(132) = xin(132) + dxkl*xin(131)
                                      yin(132) = yin(132) + dykl*yin(131)
                                      zin(132) = zin(132) + dzkl*zin(131)

                                      ! i3 = i4 =  131
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  130

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  130

                                      ! do nk = 1,    1

                                      xin(130) = xin(131) + dxkl*xin(129)
                                      yin(130) = yin(131) + dykl*yin(129)
                                      zin(130) = zin(131) + dzkl*zin(129)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  131

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  133

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  136

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  135

                                      xin(136) = xin(136) + dxkl*xin(135)
                                      yin(136) = yin(136) + dykl*yin(135)
                                      zin(136) = zin(136) + dzkl*zin(135)

                                      ! i3 = i4 =  135
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  134

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  134

                                      ! do nk = 1,    1

                                      xin(134) = xin(135) + dxkl*xin(133)
                                      yin(134) = yin(135) + dykl*yin(133)
                                      zin(134) = zin(135) + dzkl*zin(133)
                                      ! i4 = i4 + lang+1 =  136

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  135

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  137

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  140

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  139

                                      xin(140) = xin(140) + dxkl*xin(139)
                                      yin(140) = yin(140) + dykl*yin(139)
                                      zin(140) = zin(140) + dzkl*zin(139)

                                      ! i3 = i4 =  139
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  138

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  138

                                      ! do nk = 1,    1

                                      xin(138) = xin(139) + dxkl*xin(137)
                                      yin(138) = yin(139) + dykl*yin(137)
                                      zin(138) = zin(139) + dzkl*zin(137)
                                      ! i4 = i4 + lang+1 =  140

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  139

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  141

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  142

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  142

                                      ! do nk = 1,    1

                                      xin(142) = xin(143) + dxkl*xin(141)
                                      yin(142) = yin(143) + dykl*yin(141)
                                      zin(142) = zin(143) + dzkl*zin(141)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  143

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  148

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  147

                                      xin(148) = xin(148) + dxkl*xin(147)
                                      yin(148) = yin(148) + dykl*yin(147)
                                      zin(148) = zin(148) + dzkl*zin(147)

                                      ! i3 = i4 =  147
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    1

                                      xin(146) = xin(147) + dxkl*xin(145)
                                      yin(146) = yin(147) + dykl*yin(145)
                                      zin(146) = zin(147) + dzkl*zin(145)
                                      ! i4 = i4 + lang+1 =  148

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  149

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  152

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  151

                                      xin(152) = xin(152) + dxkl*xin(151)
                                      yin(152) = yin(152) + dykl*yin(151)
                                      zin(152) = zin(152) + dzkl*zin(151)

                                      ! i3 = i4 =  151
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  150

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  150

                                      ! do nk = 1,    1

                                      xin(150) = xin(151) + dxkl*xin(149)
                                      yin(150) = yin(151) + dykl*yin(149)
                                      zin(150) = zin(151) + dzkl*zin(149)
                                      ! i4 = i4 + lang+1 =  152

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  151

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  153

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  156

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  155

                                      xin(156) = xin(156) + dxkl*xin(155)
                                      yin(156) = yin(156) + dykl*yin(155)
                                      zin(156) = zin(156) + dzkl*zin(155)

                                      ! i3 = i4 =  155
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  154

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  154

                                      ! do nk = 1,    1

                                      xin(154) = xin(155) + dxkl*xin(153)
                                      yin(154) = yin(155) + dykl*yin(153)
                                      zin(154) = zin(155) + dzkl*zin(153)
                                      ! i4 = i4 + lang+1 =  156

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  155

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  157

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  160

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  159

                                      xin(160) = xin(160) + dxkl*xin(159)
                                      yin(160) = yin(160) + dykl*yin(159)
                                      zin(160) = zin(160) + dzkl*zin(159)

                                      ! i3 = i4 =  159
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  158

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  158

                                      ! do nk = 1,    1

                                      xin(158) = xin(159) + dxkl*xin(157)
                                      yin(158) = yin(159) + dykl*yin(157)
                                      zin(158) = zin(159) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  160

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  159

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  161

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  164

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  163

                                      xin(164) = xin(164) + dxkl*xin(163)
                                      yin(164) = yin(164) + dykl*yin(163)
                                      zin(164) = zin(164) + dzkl*zin(163)

                                      ! i3 = i4 =  163
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  162

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  162

                                      ! do nk = 1,    1

                                      xin(162) = xin(163) + dxkl*xin(161)
                                      yin(162) = yin(163) + dykl*yin(161)
                                      zin(162) = zin(163) + dzkl*zin(161)
                                      ! i4 = i4 + lang+1 =  164

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  163

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  165

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  168

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  167

                                      xin(168) = xin(168) + dxkl*xin(167)
                                      yin(168) = yin(168) + dykl*yin(167)
                                      zin(168) = zin(168) + dzkl*zin(167)

                                      ! i3 = i4 =  167
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  166

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  166

                                      ! do nk = 1,    1

                                      xin(166) = xin(167) + dxkl*xin(165)
                                      yin(166) = yin(167) + dykl*yin(165)
                                      zin(166) = zin(167) + dzkl*zin(165)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  167

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  169

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  172

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  171

                                      xin(172) = xin(172) + dxkl*xin(171)
                                      yin(172) = yin(172) + dykl*yin(171)
                                      zin(172) = zin(172) + dzkl*zin(171)

                                      ! i3 = i4 =  171
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  170

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  170

                                      ! do nk = 1,    1

                                      xin(170) = xin(171) + dxkl*xin(169)
                                      yin(170) = yin(171) + dykl*yin(169)
                                      zin(170) = zin(171) + dzkl*zin(169)
                                      ! i4 = i4 + lang+1 =  172

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  171

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  173

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  176

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  175

                                      xin(176) = xin(176) + dxkl*xin(175)
                                      yin(176) = yin(176) + dykl*yin(175)
                                      zin(176) = zin(176) + dzkl*zin(175)

                                      ! i3 = i4 =  175
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  174

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  174

                                      ! do nk = 1,    1

                                      xin(174) = xin(175) + dxkl*xin(173)
                                      yin(174) = yin(175) + dykl*yin(173)
                                      zin(174) = zin(175) + dzkl*zin(173)
                                      ! i4 = i4 + lang+1 =  176

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  175

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  177

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  178

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  178

                                      ! do nk = 1,    1

                                      xin(178) = xin(179) + dxkl*xin(177)
                                      yin(178) = yin(179) + dykl*yin(177)
                                      zin(178) = zin(179) + dzkl*zin(177)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  179

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  181

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  184

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  183

                                      xin(184) = xin(184) + dxkl*xin(183)
                                      yin(184) = yin(184) + dykl*yin(183)
                                      zin(184) = zin(184) + dzkl*zin(183)

                                      ! i3 = i4 =  183
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  182

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  182

                                      ! do nk = 1,    1

                                      xin(182) = xin(183) + dxkl*xin(181)
                                      yin(182) = yin(183) + dykl*yin(181)
                                      zin(182) = zin(183) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  184

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  183

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  185

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  188

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  187

                                      xin(188) = xin(188) + dxkl*xin(187)
                                      yin(188) = yin(188) + dykl*yin(187)
                                      zin(188) = zin(188) + dzkl*zin(187)

                                      ! i3 = i4 =  187
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  186

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  186

                                      ! do nk = 1,    1

                                      xin(186) = xin(187) + dxkl*xin(185)
                                      yin(186) = yin(187) + dykl*yin(185)
                                      zin(186) = zin(187) + dzkl*zin(185)
                                      ! i4 = i4 + lang+1 =  188

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  187

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  189

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  192

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  191

                                      xin(192) = xin(192) + dxkl*xin(191)
                                      yin(192) = yin(192) + dykl*yin(191)
                                      zin(192) = zin(192) + dzkl*zin(191)

                                      ! i3 = i4 =  191
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  190

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  190

                                      ! do nk = 1,    1

                                      xin(190) = xin(191) + dxkl*xin(189)
                                      yin(190) = yin(191) + dykl*yin(189)
                                      zin(190) = zin(191) + dzkl*zin(189)
                                      ! i4 = i4 + lang+1 =  192

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  191

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  193

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

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

                                      ! i1 = in(1) =  193

                                      xin(193) = 1.0_dp
                                      yin(193) = 1.0_dp
                                      zin(193) = f00

                                      ! i2 = in(2) =  209
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(209) = xc00
                                      yin(209) = yc00
                                      zin(209) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  195

                                      xin(195) = xcp00
                                      yin(195) = ycp00
                                      zin(195) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  211
                                      ! i2 =  209

                                      xin(211) = xcp00*xin(209) + cp10
                                      yin(211) = ycp00*yin(209) + cp10
                                      zin(211) = zcp00*zin(209) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  209

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  225
                                      ! i3 =  193
                                      ! i4 =  209

                                      xin(225) = c10*xin(193) + xc00*xin(209)
                                      yin(225) = c10*yin(193) + yc00*yin(209)
                                      zin(225) = c10*zin(193) + zc00*zin(209)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  227
                                      ! i5 =  225
                                      ! i4 =  209

                                      xin(227) = xcp00*xin(225) + cp10*xin(209)
                                      yin(227) = ycp00*yin(225) + cp10*yin(209)
                                      zin(227) = zcp00*zin(225) + cp10*zin(209)

                                      ! ------------------

                                      ! i3 = i4 =  209
                                      ! i4 = i5 =  225

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  241
                                      ! i3 =  209
                                      ! i4 =  225

                                      xin(241) = c10*xin(209) + xc00*xin(225)
                                      yin(241) = c10*yin(209) + yc00*yin(225)
                                      zin(241) = c10*zin(209) + zc00*zin(225)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  243
                                      ! i5 =  241
                                      ! i4 =  225

                                      xin(243) = xcp00*xin(241) + cp10*xin(225)
                                      yin(243) = ycp00*yin(241) + cp10*yin(225)
                                      zin(243) = zcp00*zin(241) + cp10*zin(225)

                                      ! ------------------

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  241

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  245
                                      ! i3 =  225
                                      ! i4 =  241

                                      xin(245) = c10*xin(225) + xc00*xin(241)
                                      yin(245) = c10*yin(225) + yc00*yin(241)
                                      zin(245) = c10*zin(225) + zc00*zin(241)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  247
                                      ! i5 =  245
                                      ! i4 =  241

                                      xin(247) = xcp00*xin(245) + cp10*xin(241)
                                      yin(247) = ycp00*yin(245) + cp10*yin(241)
                                      zin(247) = zcp00*zin(245) + cp10*zin(241)

                                      ! ------------------

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  245

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  249
                                      ! i3 =  241
                                      ! i4 =  245

                                      xin(249) = c10*xin(241) + xc00*xin(245)
                                      yin(249) = c10*yin(241) + yc00*yin(245)
                                      zin(249) = c10*zin(241) + zc00*zin(245)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  251
                                      ! i5 =  249
                                      ! i4 =  245

                                      xin(251) = xcp00*xin(249) + cp10*xin(245)
                                      yin(251) = ycp00*yin(249) + cp10*yin(245)
                                      zin(251) = zcp00*zin(249) + cp10*zin(245)

                                      ! ------------------

                                      ! i3 = i4 =  245
                                      ! i4 = i5 =  249

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  253
                                      ! i3 =  245
                                      ! i4 =  249

                                      xin(253) = c10*xin(245) + xc00*xin(249)
                                      yin(253) = c10*yin(245) + yc00*yin(249)
                                      zin(253) = c10*zin(245) + zc00*zin(249)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  255
                                      ! i5 =  253
                                      ! i4 =  249

                                      xin(255) = xcp00*xin(253) + cp10*xin(249)
                                      yin(255) = ycp00*yin(253) + cp10*yin(249)
                                      zin(255) = zcp00*zin(253) + cp10*zin(249)

                                      ! ------------------

                                      ! i3 = i4 =  249
                                      ! i4 = i5 =  253

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  193
                                      ! i4 = i1+k2 =  195

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  196
                                      ! i3 =  193
                                      ! i4 =  195

                                      xin(196) = cp01*xin(193) + xcp00*xin(195)
                                      yin(196) = cp01*yin(193) + ycp00*yin(195)
                                      zin(196) = cp01*zin(193) + zcp00*zin(195)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  212

                                      xin(212) = xc00*xin(196) + c01*xin(195)
                                      yin(212) = yc00*yin(196) + c01*yin(195)
                                      zin(212) = zc00*zin(196) + c01*zin(195)

                                      ! ------------------

                                      ! i3 = i4 =  195
                                      ! i4 = i5 =  196

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  209

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  225

                                      xin(228) = c10*xin(196) + xc00*xin(212) + c01*xin(211)
                                      yin(228) = c10*yin(196) + yc00*yin(212) + c01*yin(211)
                                      zin(228) = c10*zin(196) + zc00*zin(212) + c01*zin(211)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  209
                                      ! i4 = i5 =  225

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  241

                                      xin(244) = c10*xin(212) + xc00*xin(228) + c01*xin(227)
                                      yin(244) = c10*yin(212) + yc00*yin(228) + c01*yin(227)
                                      zin(244) = c10*zin(212) + zc00*zin(228) + c01*zin(227)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  241

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  245

                                      xin(248) = c10*xin(228) + xc00*xin(244) + c01*xin(243)
                                      yin(248) = c10*yin(228) + yc00*yin(244) + c01*yin(243)
                                      zin(248) = c10*zin(228) + zc00*zin(244) + c01*zin(243)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  245

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  249

                                      xin(252) = c10*xin(244) + xc00*xin(248) + c01*xin(247)
                                      yin(252) = c10*yin(244) + yc00*yin(248) + c01*yin(247)
                                      zin(252) = c10*zin(244) + zc00*zin(248) + c01*zin(247)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  245
                                      ! i4 = i5 =  249

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  253

                                      xin(256) = c10*xin(248) + xc00*xin(252) + c01*xin(251)
                                      yin(256) = c10*yin(248) + yc00*yin(252) + c01*yin(251)
                                      zin(256) = c10*zin(248) + zc00*zin(252) + c01*zin(251)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  249
                                      ! i4 = i5 =  253

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  253

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  249

                                      xin(253) = xin(253) + dxij*xin(249)
                                      yin(253) = yin(253) + dyij*yin(249)
                                      zin(253) = zin(253) + dzij*zin(249)

                                      ! i3 = i4 =  249
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  245

                                      xin(249) = xin(249) + dxij*xin(245)
                                      yin(249) = yin(249) + dyij*yin(245)
                                      zin(249) = zin(249) + dzij*zin(245)

                                      ! i3 = i4 =  245
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  241

                                      xin(245) = xin(245) + dxij*xin(241)
                                      yin(245) = yin(245) + dyij*yin(241)
                                      zin(245) = zin(245) + dzij*zin(241)

                                      ! i3 = i4 =  241
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  249

                                      xin(253) = xin(253) + dxij*xin(249)
                                      yin(253) = yin(253) + dyij*yin(249)
                                      zin(253) = zin(253) + dzij*zin(249)

                                      ! i3 = i4 =  249
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  245

                                      xin(249) = xin(249) + dxij*xin(245)
                                      yin(249) = yin(249) + dyij*yin(245)
                                      zin(249) = zin(249) + dzij*zin(245)

                                      ! i3 = i4 =  245
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  253

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  249

                                      xin(253) = xin(253) + dxij*xin(249)
                                      yin(253) = yin(253) + dyij*yin(249)
                                      zin(253) = zin(253) + dzij*zin(249)

                                      ! i3 = i4 =  249
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  197

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  197

                                      ! do ni = 1,    3

                                      xin(197) = xin(209) + dxij*xin(193)
                                      yin(197) = yin(209) + dyij*yin(193)
                                      zin(197) = zin(209) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    2

                                      xin(213) = xin(225) + dxij*xin(209)
                                      yin(213) = yin(225) + dyij*yin(209)
                                      zin(213) = zin(225) + dzij*zin(209)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  229

                                      ! ni =    3

                                      xin(229) = xin(241) + dxij*xin(225)
                                      yin(229) = yin(241) + dyij*yin(225)
                                      zin(229) = zin(241) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  245

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  201

                                      ! nj =    2

                                      ! i4 = i3 =  201

                                      ! do ni = 1,    3

                                      xin(201) = xin(213) + dxij*xin(197)
                                      yin(201) = yin(213) + dyij*yin(197)
                                      zin(201) = zin(213) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! ni =    2

                                      xin(217) = xin(229) + dxij*xin(213)
                                      yin(217) = yin(229) + dyij*yin(213)
                                      zin(217) = zin(229) + dzij*zin(213)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    3

                                      xin(233) = xin(245) + dxij*xin(229)
                                      yin(233) = yin(245) + dyij*yin(229)
                                      zin(233) = zin(245) + dzij*zin(229)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  249

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  205

                                      ! nj =    3

                                      ! i4 = i3 =  205

                                      ! do ni = 1,    3

                                      xin(205) = xin(217) + dxij*xin(201)
                                      yin(205) = yin(217) + dyij*yin(201)
                                      zin(205) = zin(217) + dzij*zin(201)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  221

                                      ! ni =    2

                                      xin(221) = xin(233) + dxij*xin(217)
                                      yin(221) = yin(233) + dyij*yin(217)
                                      zin(221) = zin(233) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  237

                                      ! ni =    3

                                      xin(237) = xin(249) + dxij*xin(233)
                                      yin(237) = yin(249) + dyij*yin(233)
                                      zin(237) = zin(249) + dzij*zin(233)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  209

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  251

                                      xin(255) = xin(255) + dxij*xin(251)
                                      yin(255) = yin(255) + dyij*yin(251)
                                      zin(255) = zin(255) + dzij*zin(251)

                                      ! i3 = i4 =  251
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  247

                                      xin(251) = xin(251) + dxij*xin(247)
                                      yin(251) = yin(251) + dyij*yin(247)
                                      zin(251) = zin(251) + dzij*zin(247)

                                      ! i3 = i4 =  247
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  243

                                      xin(247) = xin(247) + dxij*xin(243)
                                      yin(247) = yin(247) + dyij*yin(243)
                                      zin(247) = zin(247) + dzij*zin(243)

                                      ! i3 = i4 =  243
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  251

                                      xin(255) = xin(255) + dxij*xin(251)
                                      yin(255) = yin(255) + dyij*yin(251)
                                      zin(255) = zin(255) + dzij*zin(251)

                                      ! i3 = i4 =  251
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  247

                                      xin(251) = xin(251) + dxij*xin(247)
                                      yin(251) = yin(251) + dyij*yin(247)
                                      zin(251) = zin(251) + dzij*zin(247)

                                      ! i3 = i4 =  247
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  255

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  251

                                      xin(255) = xin(255) + dxij*xin(251)
                                      yin(255) = yin(255) + dyij*yin(251)
                                      zin(255) = zin(255) + dzij*zin(251)

                                      ! i3 = i4 =  251
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  199

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  199

                                      ! do ni = 1,    3

                                      xin(199) = xin(211) + dxij*xin(195)
                                      yin(199) = yin(211) + dyij*yin(195)
                                      zin(199) = zin(211) + dzij*zin(195)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    2

                                      xin(215) = xin(227) + dxij*xin(211)
                                      yin(215) = yin(227) + dyij*yin(211)
                                      zin(215) = zin(227) + dzij*zin(211)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  231

                                      ! ni =    3

                                      xin(231) = xin(243) + dxij*xin(227)
                                      yin(231) = yin(243) + dyij*yin(227)
                                      zin(231) = zin(243) + dzij*zin(227)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  203

                                      ! nj =    2

                                      ! i4 = i3 =  203

                                      ! do ni = 1,    3

                                      xin(203) = xin(215) + dxij*xin(199)
                                      yin(203) = yin(215) + dyij*yin(199)
                                      zin(203) = zin(215) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  219

                                      ! ni =    2

                                      xin(219) = xin(231) + dxij*xin(215)
                                      yin(219) = yin(231) + dyij*yin(215)
                                      zin(219) = zin(231) + dzij*zin(215)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    3

                                      xin(235) = xin(247) + dxij*xin(231)
                                      yin(235) = yin(247) + dyij*yin(231)
                                      zin(235) = zin(247) + dzij*zin(231)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  207

                                      ! nj =    3

                                      ! i4 = i3 =  207

                                      ! do ni = 1,    3

                                      xin(207) = xin(219) + dxij*xin(203)
                                      yin(207) = yin(219) + dyij*yin(203)
                                      zin(207) = zin(219) + dzij*zin(203)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  223

                                      ! ni =    2

                                      xin(223) = xin(235) + dxij*xin(219)
                                      yin(223) = yin(235) + dyij*yin(219)
                                      zin(223) = zin(235) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                      ! ni =    3

                                      xin(239) = xin(251) + dxij*xin(235)
                                      yin(239) = yin(251) + dyij*yin(235)
                                      zin(239) = zin(251) + dzij*zin(235)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  255

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  211

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  252

                                      xin(256) = xin(256) + dxij*xin(252)
                                      yin(256) = yin(256) + dyij*yin(252)
                                      zin(256) = zin(256) + dzij*zin(252)

                                      ! i3 = i4 =  252
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  248

                                      xin(252) = xin(252) + dxij*xin(248)
                                      yin(252) = yin(252) + dyij*yin(248)
                                      zin(252) = zin(252) + dzij*zin(248)

                                      ! i3 = i4 =  248
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  244

                                      xin(248) = xin(248) + dxij*xin(244)
                                      yin(248) = yin(248) + dyij*yin(244)
                                      zin(248) = zin(248) + dzij*zin(244)

                                      ! i3 = i4 =  244
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  252

                                      xin(256) = xin(256) + dxij*xin(252)
                                      yin(256) = yin(256) + dyij*yin(252)
                                      zin(256) = zin(256) + dzij*zin(252)

                                      ! i3 = i4 =  252
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  248

                                      xin(252) = xin(252) + dxij*xin(248)
                                      yin(252) = yin(252) + dyij*yin(248)
                                      zin(252) = zin(252) + dzij*zin(248)

                                      ! i3 = i4 =  248
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  256

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  252

                                      xin(256) = xin(256) + dxij*xin(252)
                                      yin(256) = yin(256) + dyij*yin(252)
                                      zin(256) = zin(256) + dzij*zin(252)

                                      ! i3 = i4 =  252
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  200

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  200

                                      ! do ni = 1,    3

                                      xin(200) = xin(212) + dxij*xin(196)
                                      yin(200) = yin(212) + dyij*yin(196)
                                      zin(200) = zin(212) + dzij*zin(196)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    2

                                      xin(216) = xin(228) + dxij*xin(212)
                                      yin(216) = yin(228) + dyij*yin(212)
                                      zin(216) = zin(228) + dzij*zin(212)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  232

                                      ! ni =    3

                                      xin(232) = xin(244) + dxij*xin(228)
                                      yin(232) = yin(244) + dyij*yin(228)
                                      zin(232) = zin(244) + dzij*zin(228)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  248

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  204

                                      ! nj =    2

                                      ! i4 = i3 =  204

                                      ! do ni = 1,    3

                                      xin(204) = xin(216) + dxij*xin(200)
                                      yin(204) = yin(216) + dyij*yin(200)
                                      zin(204) = zin(216) + dzij*zin(200)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  220

                                      ! ni =    2

                                      xin(220) = xin(232) + dxij*xin(216)
                                      yin(220) = yin(232) + dyij*yin(216)
                                      zin(220) = zin(232) + dzij*zin(216)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  236

                                      ! ni =    3

                                      xin(236) = xin(248) + dxij*xin(232)
                                      yin(236) = yin(248) + dyij*yin(232)
                                      zin(236) = zin(248) + dzij*zin(232)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  252

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  208

                                      ! nj =    3

                                      ! i4 = i3 =  208

                                      ! do ni = 1,    3

                                      xin(208) = xin(220) + dxij*xin(204)
                                      yin(208) = yin(220) + dyij*yin(204)
                                      zin(208) = zin(220) + dzij*zin(204)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  224

                                      ! ni =    2

                                      xin(224) = xin(236) + dxij*xin(220)
                                      yin(224) = yin(236) + dyij*yin(220)
                                      zin(224) = zin(236) + dzij*zin(220)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                      ! ni =    3

                                      xin(240) = xin(252) + dxij*xin(236)
                                      yin(240) = yin(252) + dyij*yin(236)
                                      zin(240) = zin(252) + dzij*zin(236)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  256

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  212

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =  193

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  196

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  195

                                      xin(196) = xin(196) + dxkl*xin(195)
                                      yin(196) = yin(196) + dykl*yin(195)
                                      zin(196) = zin(196) + dzkl*zin(195)

                                      ! i3 = i4 =  195
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  194

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  194

                                      ! do nk = 1,    1

                                      xin(194) = xin(195) + dxkl*xin(193)
                                      yin(194) = yin(195) + dykl*yin(193)
                                      zin(194) = zin(195) + dzkl*zin(193)
                                      ! i4 = i4 + lang+1 =  196

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  195

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  197

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  200

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  199

                                      xin(200) = xin(200) + dxkl*xin(199)
                                      yin(200) = yin(200) + dykl*yin(199)
                                      zin(200) = zin(200) + dzkl*zin(199)

                                      ! i3 = i4 =  199
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  198

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  198

                                      ! do nk = 1,    1

                                      xin(198) = xin(199) + dxkl*xin(197)
                                      yin(198) = yin(199) + dykl*yin(197)
                                      zin(198) = zin(199) + dzkl*zin(197)
                                      ! i4 = i4 + lang+1 =  200

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  199

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  201

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  204

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  203

                                      xin(204) = xin(204) + dxkl*xin(203)
                                      yin(204) = yin(204) + dykl*yin(203)
                                      zin(204) = zin(204) + dzkl*zin(203)

                                      ! i3 = i4 =  203
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  202

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  202

                                      ! do nk = 1,    1

                                      xin(202) = xin(203) + dxkl*xin(201)
                                      yin(202) = yin(203) + dykl*yin(201)
                                      zin(202) = zin(203) + dzkl*zin(201)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  203

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  205

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  208

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  207

                                      xin(208) = xin(208) + dxkl*xin(207)
                                      yin(208) = yin(208) + dykl*yin(207)
                                      zin(208) = zin(208) + dzkl*zin(207)

                                      ! i3 = i4 =  207
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  206

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  206

                                      ! do nk = 1,    1

                                      xin(206) = xin(207) + dxkl*xin(205)
                                      yin(206) = yin(207) + dykl*yin(205)
                                      zin(206) = zin(207) + dzkl*zin(205)
                                      ! i4 = i4 + lang+1 =  208

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  207

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  209

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  212

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  211

                                      xin(212) = xin(212) + dxkl*xin(211)
                                      yin(212) = yin(212) + dykl*yin(211)
                                      zin(212) = zin(212) + dzkl*zin(211)

                                      ! i3 = i4 =  211
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  210

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  210

                                      ! do nk = 1,    1

                                      xin(210) = xin(211) + dxkl*xin(209)
                                      yin(210) = yin(211) + dykl*yin(209)
                                      zin(210) = zin(211) + dzkl*zin(209)
                                      ! i4 = i4 + lang+1 =  212

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  211

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  213

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  214

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  214

                                      ! do nk = 1,    1

                                      xin(214) = xin(215) + dxkl*xin(213)
                                      yin(214) = yin(215) + dykl*yin(213)
                                      zin(214) = zin(215) + dzkl*zin(213)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  215

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  220

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  219

                                      xin(220) = xin(220) + dxkl*xin(219)
                                      yin(220) = yin(220) + dykl*yin(219)
                                      zin(220) = zin(220) + dzkl*zin(219)

                                      ! i3 = i4 =  219
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    1

                                      xin(218) = xin(219) + dxkl*xin(217)
                                      yin(218) = yin(219) + dykl*yin(217)
                                      zin(218) = zin(219) + dzkl*zin(217)
                                      ! i4 = i4 + lang+1 =  220

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  221

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  224

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  223

                                      xin(224) = xin(224) + dxkl*xin(223)
                                      yin(224) = yin(224) + dykl*yin(223)
                                      zin(224) = zin(224) + dzkl*zin(223)

                                      ! i3 = i4 =  223
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  222

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  222

                                      ! do nk = 1,    1

                                      xin(222) = xin(223) + dxkl*xin(221)
                                      yin(222) = yin(223) + dykl*yin(221)
                                      zin(222) = zin(223) + dzkl*zin(221)
                                      ! i4 = i4 + lang+1 =  224

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  223

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  225

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  225

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  228

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  227

                                      xin(228) = xin(228) + dxkl*xin(227)
                                      yin(228) = yin(228) + dykl*yin(227)
                                      zin(228) = zin(228) + dzkl*zin(227)

                                      ! i3 = i4 =  227
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  226

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  226

                                      ! do nk = 1,    1

                                      xin(226) = xin(227) + dxkl*xin(225)
                                      yin(226) = yin(227) + dykl*yin(225)
                                      zin(226) = zin(227) + dzkl*zin(225)
                                      ! i4 = i4 + lang+1 =  228

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  227

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  229

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  232

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  231

                                      xin(232) = xin(232) + dxkl*xin(231)
                                      yin(232) = yin(232) + dykl*yin(231)
                                      zin(232) = zin(232) + dzkl*zin(231)

                                      ! i3 = i4 =  231
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  230

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  230

                                      ! do nk = 1,    1

                                      xin(230) = xin(231) + dxkl*xin(229)
                                      yin(230) = yin(231) + dykl*yin(229)
                                      zin(230) = zin(231) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  232

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  231

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  233

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  236

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  235

                                      xin(236) = xin(236) + dxkl*xin(235)
                                      yin(236) = yin(236) + dykl*yin(235)
                                      zin(236) = zin(236) + dzkl*zin(235)

                                      ! i3 = i4 =  235
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  234

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  234

                                      ! do nk = 1,    1

                                      xin(234) = xin(235) + dxkl*xin(233)
                                      yin(234) = yin(235) + dykl*yin(233)
                                      zin(234) = zin(235) + dzkl*zin(233)
                                      ! i4 = i4 + lang+1 =  236

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  235

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  237

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  240

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  239

                                      xin(240) = xin(240) + dxkl*xin(239)
                                      yin(240) = yin(240) + dykl*yin(239)
                                      zin(240) = zin(240) + dzkl*zin(239)

                                      ! i3 = i4 =  239
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  238

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  238

                                      ! do nk = 1,    1

                                      xin(238) = xin(239) + dxkl*xin(237)
                                      yin(238) = yin(239) + dykl*yin(237)
                                      zin(238) = zin(239) + dzkl*zin(237)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  239

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  241

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  244

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  243

                                      xin(244) = xin(244) + dxkl*xin(243)
                                      yin(244) = yin(244) + dykl*yin(243)
                                      zin(244) = zin(244) + dzkl*zin(243)

                                      ! i3 = i4 =  243
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  242

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  242

                                      ! do nk = 1,    1

                                      xin(242) = xin(243) + dxkl*xin(241)
                                      yin(242) = yin(243) + dykl*yin(241)
                                      zin(242) = zin(243) + dzkl*zin(241)
                                      ! i4 = i4 + lang+1 =  244

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  243

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  245

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  248

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  247

                                      xin(248) = xin(248) + dxkl*xin(247)
                                      yin(248) = yin(248) + dykl*yin(247)
                                      zin(248) = zin(248) + dzkl*zin(247)

                                      ! i3 = i4 =  247
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  246

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  246

                                      ! do nk = 1,    1

                                      xin(246) = xin(247) + dxkl*xin(245)
                                      yin(246) = yin(247) + dykl*yin(245)
                                      zin(246) = zin(247) + dzkl*zin(245)
                                      ! i4 = i4 + lang+1 =  248

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  247

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  249

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  250

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  250

                                      ! do nk = 1,    1

                                      xin(250) = xin(251) + dxkl*xin(249)
                                      yin(250) = yin(251) + dykl*yin(249)
                                      zin(250) = zin(251) + dzkl*zin(249)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  251

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  253

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  256

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  255

                                      xin(256) = xin(256) + dxkl*xin(255)
                                      yin(256) = yin(256) + dykl*yin(255)
                                      zin(256) = zin(256) + dzkl*zin(255)

                                      ! i3 = i4 =  255
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  254

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  254

                                      ! do nk = 1,    1

                                      xin(254) = xin(255) + dxkl*xin(253)
                                      yin(254) = yin(255) + dykl*yin(253)
                                      zin(254) = zin(255) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  256

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  255

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  257

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  257

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  256

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

                                      ! i1 = in(1) =  257

                                      xin(257) = 1.0_dp
                                      yin(257) = 1.0_dp
                                      zin(257) = f00

                                      ! i2 = in(2) =  273
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(273) = xc00
                                      yin(273) = yc00
                                      zin(273) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  259

                                      xin(259) = xcp00
                                      yin(259) = ycp00
                                      zin(259) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  275
                                      ! i2 =  273

                                      xin(275) = xcp00*xin(273) + cp10
                                      yin(275) = ycp00*yin(273) + cp10
                                      zin(275) = zcp00*zin(273) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  273

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  289
                                      ! i3 =  257
                                      ! i4 =  273

                                      xin(289) = c10*xin(257) + xc00*xin(273)
                                      yin(289) = c10*yin(257) + yc00*yin(273)
                                      zin(289) = c10*zin(257) + zc00*zin(273)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  291
                                      ! i5 =  289
                                      ! i4 =  273

                                      xin(291) = xcp00*xin(289) + cp10*xin(273)
                                      yin(291) = ycp00*yin(289) + cp10*yin(273)
                                      zin(291) = zcp00*zin(289) + cp10*zin(273)

                                      ! ------------------

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  289

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  305
                                      ! i3 =  273
                                      ! i4 =  289

                                      xin(305) = c10*xin(273) + xc00*xin(289)
                                      yin(305) = c10*yin(273) + yc00*yin(289)
                                      zin(305) = c10*zin(273) + zc00*zin(289)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  307
                                      ! i5 =  305
                                      ! i4 =  289

                                      xin(307) = xcp00*xin(305) + cp10*xin(289)
                                      yin(307) = ycp00*yin(305) + cp10*yin(289)
                                      zin(307) = zcp00*zin(305) + cp10*zin(289)

                                      ! ------------------

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  305

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  309
                                      ! i3 =  289
                                      ! i4 =  305

                                      xin(309) = c10*xin(289) + xc00*xin(305)
                                      yin(309) = c10*yin(289) + yc00*yin(305)
                                      zin(309) = c10*zin(289) + zc00*zin(305)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  311
                                      ! i5 =  309
                                      ! i4 =  305

                                      xin(311) = xcp00*xin(309) + cp10*xin(305)
                                      yin(311) = ycp00*yin(309) + cp10*yin(305)
                                      zin(311) = zcp00*zin(309) + cp10*zin(305)

                                      ! ------------------

                                      ! i3 = i4 =  305
                                      ! i4 = i5 =  309

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  313
                                      ! i3 =  305
                                      ! i4 =  309

                                      xin(313) = c10*xin(305) + xc00*xin(309)
                                      yin(313) = c10*yin(305) + yc00*yin(309)
                                      zin(313) = c10*zin(305) + zc00*zin(309)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  315
                                      ! i5 =  313
                                      ! i4 =  309

                                      xin(315) = xcp00*xin(313) + cp10*xin(309)
                                      yin(315) = ycp00*yin(313) + cp10*yin(309)
                                      zin(315) = zcp00*zin(313) + cp10*zin(309)

                                      ! ------------------

                                      ! i3 = i4 =  309
                                      ! i4 = i5 =  313

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  317
                                      ! i3 =  309
                                      ! i4 =  313

                                      xin(317) = c10*xin(309) + xc00*xin(313)
                                      yin(317) = c10*yin(309) + yc00*yin(313)
                                      zin(317) = c10*zin(309) + zc00*zin(313)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  319
                                      ! i5 =  317
                                      ! i4 =  313

                                      xin(319) = xcp00*xin(317) + cp10*xin(313)
                                      yin(319) = ycp00*yin(317) + cp10*yin(313)
                                      zin(319) = zcp00*zin(317) + cp10*zin(313)

                                      ! ------------------

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  317

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  257
                                      ! i4 = i1+k2 =  259

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  260
                                      ! i3 =  257
                                      ! i4 =  259

                                      xin(260) = cp01*xin(257) + xcp00*xin(259)
                                      yin(260) = cp01*yin(257) + ycp00*yin(259)
                                      zin(260) = cp01*zin(257) + zcp00*zin(259)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  276

                                      xin(276) = xc00*xin(260) + c01*xin(259)
                                      yin(276) = yc00*yin(260) + c01*yin(259)
                                      zin(276) = zc00*zin(260) + c01*zin(259)

                                      ! ------------------

                                      ! i3 = i4 =  259
                                      ! i4 = i5 =  260

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  273

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  289

                                      xin(292) = c10*xin(260) + xc00*xin(276) + c01*xin(275)
                                      yin(292) = c10*yin(260) + yc00*yin(276) + c01*yin(275)
                                      zin(292) = c10*zin(260) + zc00*zin(276) + c01*zin(275)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  289

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  305

                                      xin(308) = c10*xin(276) + xc00*xin(292) + c01*xin(291)
                                      yin(308) = c10*yin(276) + yc00*yin(292) + c01*yin(291)
                                      zin(308) = c10*zin(276) + zc00*zin(292) + c01*zin(291)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  305

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  309

                                      xin(312) = c10*xin(292) + xc00*xin(308) + c01*xin(307)
                                      yin(312) = c10*yin(292) + yc00*yin(308) + c01*yin(307)
                                      zin(312) = c10*zin(292) + zc00*zin(308) + c01*zin(307)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  305
                                      ! i4 = i5 =  309

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  313

                                      xin(316) = c10*xin(308) + xc00*xin(312) + c01*xin(311)
                                      yin(316) = c10*yin(308) + yc00*yin(312) + c01*yin(311)
                                      zin(316) = c10*zin(308) + zc00*zin(312) + c01*zin(311)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  309
                                      ! i4 = i5 =  313

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  317

                                      xin(320) = c10*xin(312) + xc00*xin(316) + c01*xin(315)
                                      yin(320) = c10*yin(312) + yc00*yin(316) + c01*yin(315)
                                      zin(320) = c10*zin(312) + zc00*zin(316) + c01*zin(315)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  317

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  317

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  317

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  313

                                      xin(317) = xin(317) + dxij*xin(313)
                                      yin(317) = yin(317) + dyij*yin(313)
                                      zin(317) = zin(317) + dzij*zin(313)

                                      ! i3 = i4 =  313
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  309

                                      xin(313) = xin(313) + dxij*xin(309)
                                      yin(313) = yin(313) + dyij*yin(309)
                                      zin(313) = zin(313) + dzij*zin(309)

                                      ! i3 = i4 =  309
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  305

                                      xin(309) = xin(309) + dxij*xin(305)
                                      yin(309) = yin(309) + dyij*yin(305)
                                      zin(309) = zin(309) + dzij*zin(305)

                                      ! i3 = i4 =  305
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  317

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  313

                                      xin(317) = xin(317) + dxij*xin(313)
                                      yin(317) = yin(317) + dyij*yin(313)
                                      zin(317) = zin(317) + dzij*zin(313)

                                      ! i3 = i4 =  313
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  309

                                      xin(313) = xin(313) + dxij*xin(309)
                                      yin(313) = yin(313) + dyij*yin(309)
                                      zin(313) = zin(313) + dzij*zin(309)

                                      ! i3 = i4 =  309
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  317

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  313

                                      xin(317) = xin(317) + dxij*xin(313)
                                      yin(317) = yin(317) + dyij*yin(313)
                                      zin(317) = zin(317) + dzij*zin(313)

                                      ! i3 = i4 =  313
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  261

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  261

                                      ! do ni = 1,    3

                                      xin(261) = xin(273) + dxij*xin(257)
                                      yin(261) = yin(273) + dyij*yin(257)
                                      zin(261) = zin(273) + dzij*zin(257)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                      ! ni =    2

                                      xin(277) = xin(289) + dxij*xin(273)
                                      yin(277) = yin(289) + dyij*yin(273)
                                      zin(277) = zin(289) + dzij*zin(273)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  293

                                      ! ni =    3

                                      xin(293) = xin(305) + dxij*xin(289)
                                      yin(293) = yin(305) + dyij*yin(289)
                                      zin(293) = zin(305) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  309

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  265

                                      ! nj =    2

                                      ! i4 = i3 =  265

                                      ! do ni = 1,    3

                                      xin(265) = xin(277) + dxij*xin(261)
                                      yin(265) = yin(277) + dyij*yin(261)
                                      zin(265) = zin(277) + dzij*zin(261)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  281

                                      ! ni =    2

                                      xin(281) = xin(293) + dxij*xin(277)
                                      yin(281) = yin(293) + dyij*yin(277)
                                      zin(281) = zin(293) + dzij*zin(277)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  297

                                      ! ni =    3

                                      xin(297) = xin(309) + dxij*xin(293)
                                      yin(297) = yin(309) + dyij*yin(293)
                                      zin(297) = zin(309) + dzij*zin(293)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  269

                                      ! nj =    3

                                      ! i4 = i3 =  269

                                      ! do ni = 1,    3

                                      xin(269) = xin(281) + dxij*xin(265)
                                      yin(269) = yin(281) + dyij*yin(265)
                                      zin(269) = zin(281) + dzij*zin(265)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  285

                                      ! ni =    2

                                      xin(285) = xin(297) + dxij*xin(281)
                                      yin(285) = yin(297) + dyij*yin(281)
                                      zin(285) = zin(297) + dzij*zin(281)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  301

                                      ! ni =    3

                                      xin(301) = xin(313) + dxij*xin(297)
                                      yin(301) = yin(313) + dyij*yin(297)
                                      zin(301) = zin(313) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  317

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  273

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  319

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  315

                                      xin(319) = xin(319) + dxij*xin(315)
                                      yin(319) = yin(319) + dyij*yin(315)
                                      zin(319) = zin(319) + dzij*zin(315)

                                      ! i3 = i4 =  315
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  311

                                      xin(315) = xin(315) + dxij*xin(311)
                                      yin(315) = yin(315) + dyij*yin(311)
                                      zin(315) = zin(315) + dzij*zin(311)

                                      ! i3 = i4 =  311
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  307

                                      xin(311) = xin(311) + dxij*xin(307)
                                      yin(311) = yin(311) + dyij*yin(307)
                                      zin(311) = zin(311) + dzij*zin(307)

                                      ! i3 = i4 =  307
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  319

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  315

                                      xin(319) = xin(319) + dxij*xin(315)
                                      yin(319) = yin(319) + dyij*yin(315)
                                      zin(319) = zin(319) + dzij*zin(315)

                                      ! i3 = i4 =  315
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  311

                                      xin(315) = xin(315) + dxij*xin(311)
                                      yin(315) = yin(315) + dyij*yin(311)
                                      zin(315) = zin(315) + dzij*zin(311)

                                      ! i3 = i4 =  311
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  319

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  315

                                      xin(319) = xin(319) + dxij*xin(315)
                                      yin(319) = yin(319) + dyij*yin(315)
                                      zin(319) = zin(319) + dzij*zin(315)

                                      ! i3 = i4 =  315
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  263

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  263

                                      ! do ni = 1,    3

                                      xin(263) = xin(275) + dxij*xin(259)
                                      yin(263) = yin(275) + dyij*yin(259)
                                      zin(263) = zin(275) + dzij*zin(259)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  279

                                      ! ni =    2

                                      xin(279) = xin(291) + dxij*xin(275)
                                      yin(279) = yin(291) + dyij*yin(275)
                                      zin(279) = zin(291) + dzij*zin(275)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  295

                                      ! ni =    3

                                      xin(295) = xin(307) + dxij*xin(291)
                                      yin(295) = yin(307) + dyij*yin(291)
                                      zin(295) = zin(307) + dzij*zin(291)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  311

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  267

                                      ! nj =    2

                                      ! i4 = i3 =  267

                                      ! do ni = 1,    3

                                      xin(267) = xin(279) + dxij*xin(263)
                                      yin(267) = yin(279) + dyij*yin(263)
                                      zin(267) = zin(279) + dzij*zin(263)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    2

                                      xin(283) = xin(295) + dxij*xin(279)
                                      yin(283) = yin(295) + dyij*yin(279)
                                      zin(283) = zin(295) + dzij*zin(279)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  299

                                      ! ni =    3

                                      xin(299) = xin(311) + dxij*xin(295)
                                      yin(299) = yin(311) + dyij*yin(295)
                                      zin(299) = zin(311) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  315

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  271

                                      ! nj =    3

                                      ! i4 = i3 =  271

                                      ! do ni = 1,    3

                                      xin(271) = xin(283) + dxij*xin(267)
                                      yin(271) = yin(283) + dyij*yin(267)
                                      zin(271) = zin(283) + dzij*zin(267)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    2

                                      xin(287) = xin(299) + dxij*xin(283)
                                      yin(287) = yin(299) + dyij*yin(283)
                                      zin(287) = zin(299) + dzij*zin(283)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  303

                                      ! ni =    3

                                      xin(303) = xin(315) + dxij*xin(299)
                                      yin(303) = yin(315) + dyij*yin(299)
                                      zin(303) = zin(315) + dzij*zin(299)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  319

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  275

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  320

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  316

                                      xin(320) = xin(320) + dxij*xin(316)
                                      yin(320) = yin(320) + dyij*yin(316)
                                      zin(320) = zin(320) + dzij*zin(316)

                                      ! i3 = i4 =  316
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  312

                                      xin(316) = xin(316) + dxij*xin(312)
                                      yin(316) = yin(316) + dyij*yin(312)
                                      zin(316) = zin(316) + dzij*zin(312)

                                      ! i3 = i4 =  312
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  308

                                      xin(312) = xin(312) + dxij*xin(308)
                                      yin(312) = yin(312) + dyij*yin(308)
                                      zin(312) = zin(312) + dzij*zin(308)

                                      ! i3 = i4 =  308
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  320

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  316

                                      xin(320) = xin(320) + dxij*xin(316)
                                      yin(320) = yin(320) + dyij*yin(316)
                                      zin(320) = zin(320) + dzij*zin(316)

                                      ! i3 = i4 =  316
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  312

                                      xin(316) = xin(316) + dxij*xin(312)
                                      yin(316) = yin(316) + dyij*yin(312)
                                      zin(316) = zin(316) + dzij*zin(312)

                                      ! i3 = i4 =  312
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  320

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  316

                                      xin(320) = xin(320) + dxij*xin(316)
                                      yin(320) = yin(320) + dyij*yin(316)
                                      zin(320) = zin(320) + dzij*zin(316)

                                      ! i3 = i4 =  316
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  264

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  264

                                      ! do ni = 1,    3

                                      xin(264) = xin(276) + dxij*xin(260)
                                      yin(264) = yin(276) + dyij*yin(260)
                                      zin(264) = zin(276) + dzij*zin(260)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  280

                                      ! ni =    2

                                      xin(280) = xin(292) + dxij*xin(276)
                                      yin(280) = yin(292) + dyij*yin(276)
                                      zin(280) = zin(292) + dzij*zin(276)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  296

                                      ! ni =    3

                                      xin(296) = xin(308) + dxij*xin(292)
                                      yin(296) = yin(308) + dyij*yin(292)
                                      zin(296) = zin(308) + dzij*zin(292)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  312

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  268

                                      ! nj =    2

                                      ! i4 = i3 =  268

                                      ! do ni = 1,    3

                                      xin(268) = xin(280) + dxij*xin(264)
                                      yin(268) = yin(280) + dyij*yin(264)
                                      zin(268) = zin(280) + dzij*zin(264)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  284

                                      ! ni =    2

                                      xin(284) = xin(296) + dxij*xin(280)
                                      yin(284) = yin(296) + dyij*yin(280)
                                      zin(284) = zin(296) + dzij*zin(280)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  300

                                      ! ni =    3

                                      xin(300) = xin(312) + dxij*xin(296)
                                      yin(300) = yin(312) + dyij*yin(296)
                                      zin(300) = zin(312) + dzij*zin(296)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  316

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  272

                                      ! nj =    3

                                      ! i4 = i3 =  272

                                      ! do ni = 1,    3

                                      xin(272) = xin(284) + dxij*xin(268)
                                      yin(272) = yin(284) + dyij*yin(268)
                                      zin(272) = zin(284) + dzij*zin(268)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    2

                                      xin(288) = xin(300) + dxij*xin(284)
                                      yin(288) = yin(300) + dyij*yin(284)
                                      zin(288) = zin(300) + dzij*zin(284)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  304

                                      ! ni =    3

                                      xin(304) = xin(316) + dxij*xin(300)
                                      yin(304) = yin(316) + dyij*yin(300)
                                      zin(304) = zin(316) + dzij*zin(300)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  320

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  276

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =  257

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  260

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  259

                                      xin(260) = xin(260) + dxkl*xin(259)
                                      yin(260) = yin(260) + dykl*yin(259)
                                      zin(260) = zin(260) + dzkl*zin(259)

                                      ! i3 = i4 =  259
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  258

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  258

                                      ! do nk = 1,    1

                                      xin(258) = xin(259) + dxkl*xin(257)
                                      yin(258) = yin(259) + dykl*yin(257)
                                      zin(258) = zin(259) + dzkl*zin(257)
                                      ! i4 = i4 + lang+1 =  260

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  259

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  261

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  264

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  263

                                      xin(264) = xin(264) + dxkl*xin(263)
                                      yin(264) = yin(264) + dykl*yin(263)
                                      zin(264) = zin(264) + dzkl*zin(263)

                                      ! i3 = i4 =  263
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  262

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  262

                                      ! do nk = 1,    1

                                      xin(262) = xin(263) + dxkl*xin(261)
                                      yin(262) = yin(263) + dykl*yin(261)
                                      zin(262) = zin(263) + dzkl*zin(261)
                                      ! i4 = i4 + lang+1 =  264

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  263

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  265

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  268

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  267

                                      xin(268) = xin(268) + dxkl*xin(267)
                                      yin(268) = yin(268) + dykl*yin(267)
                                      zin(268) = zin(268) + dzkl*zin(267)

                                      ! i3 = i4 =  267
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  266

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  266

                                      ! do nk = 1,    1

                                      xin(266) = xin(267) + dxkl*xin(265)
                                      yin(266) = yin(267) + dykl*yin(265)
                                      zin(266) = zin(267) + dzkl*zin(265)
                                      ! i4 = i4 + lang+1 =  268

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  267

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  269

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  272

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  271

                                      xin(272) = xin(272) + dxkl*xin(271)
                                      yin(272) = yin(272) + dykl*yin(271)
                                      zin(272) = zin(272) + dzkl*zin(271)

                                      ! i3 = i4 =  271
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  270

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  270

                                      ! do nk = 1,    1

                                      xin(270) = xin(271) + dxkl*xin(269)
                                      yin(270) = yin(271) + dykl*yin(269)
                                      zin(270) = zin(271) + dzkl*zin(269)
                                      ! i4 = i4 + lang+1 =  272

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  271

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  273

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  273

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  276

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  275

                                      xin(276) = xin(276) + dxkl*xin(275)
                                      yin(276) = yin(276) + dykl*yin(275)
                                      zin(276) = zin(276) + dzkl*zin(275)

                                      ! i3 = i4 =  275
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  274

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  274

                                      ! do nk = 1,    1

                                      xin(274) = xin(275) + dxkl*xin(273)
                                      yin(274) = yin(275) + dykl*yin(273)
                                      zin(274) = zin(275) + dzkl*zin(273)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  275

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  277

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  280

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  279

                                      xin(280) = xin(280) + dxkl*xin(279)
                                      yin(280) = yin(280) + dykl*yin(279)
                                      zin(280) = zin(280) + dzkl*zin(279)

                                      ! i3 = i4 =  279
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  278

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  278

                                      ! do nk = 1,    1

                                      xin(278) = xin(279) + dxkl*xin(277)
                                      yin(278) = yin(279) + dykl*yin(277)
                                      zin(278) = zin(279) + dzkl*zin(277)
                                      ! i4 = i4 + lang+1 =  280

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  279

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  281

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  284

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  283

                                      xin(284) = xin(284) + dxkl*xin(283)
                                      yin(284) = yin(284) + dykl*yin(283)
                                      zin(284) = zin(284) + dzkl*zin(283)

                                      ! i3 = i4 =  283
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  282

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  282

                                      ! do nk = 1,    1

                                      xin(282) = xin(283) + dxkl*xin(281)
                                      yin(282) = yin(283) + dykl*yin(281)
                                      zin(282) = zin(283) + dzkl*zin(281)
                                      ! i4 = i4 + lang+1 =  284

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  283

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  285

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  286

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  286

                                      ! do nk = 1,    1

                                      xin(286) = xin(287) + dxkl*xin(285)
                                      yin(286) = yin(287) + dykl*yin(285)
                                      zin(286) = zin(287) + dzkl*zin(285)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  287

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  292

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  291

                                      xin(292) = xin(292) + dxkl*xin(291)
                                      yin(292) = yin(292) + dykl*yin(291)
                                      zin(292) = zin(292) + dzkl*zin(291)

                                      ! i3 = i4 =  291
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    1

                                      xin(290) = xin(291) + dxkl*xin(289)
                                      yin(290) = yin(291) + dykl*yin(289)
                                      zin(290) = zin(291) + dzkl*zin(289)
                                      ! i4 = i4 + lang+1 =  292

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  293

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  296

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  295

                                      xin(296) = xin(296) + dxkl*xin(295)
                                      yin(296) = yin(296) + dykl*yin(295)
                                      zin(296) = zin(296) + dzkl*zin(295)

                                      ! i3 = i4 =  295
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  294

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  294

                                      ! do nk = 1,    1

                                      xin(294) = xin(295) + dxkl*xin(293)
                                      yin(294) = yin(295) + dykl*yin(293)
                                      zin(294) = zin(295) + dzkl*zin(293)
                                      ! i4 = i4 + lang+1 =  296

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  295

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  297

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  300

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  299

                                      xin(300) = xin(300) + dxkl*xin(299)
                                      yin(300) = yin(300) + dykl*yin(299)
                                      zin(300) = zin(300) + dzkl*zin(299)

                                      ! i3 = i4 =  299
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  298

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  298

                                      ! do nk = 1,    1

                                      xin(298) = xin(299) + dxkl*xin(297)
                                      yin(298) = yin(299) + dykl*yin(297)
                                      zin(298) = zin(299) + dzkl*zin(297)
                                      ! i4 = i4 + lang+1 =  300

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  299

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  301

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  304

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  303

                                      xin(304) = xin(304) + dxkl*xin(303)
                                      yin(304) = yin(304) + dykl*yin(303)
                                      zin(304) = zin(304) + dzkl*zin(303)

                                      ! i3 = i4 =  303
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  302

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  302

                                      ! do nk = 1,    1

                                      xin(302) = xin(303) + dxkl*xin(301)
                                      yin(302) = yin(303) + dykl*yin(301)
                                      zin(302) = zin(303) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  304

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  303

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  305

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  305

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  308

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  307

                                      xin(308) = xin(308) + dxkl*xin(307)
                                      yin(308) = yin(308) + dykl*yin(307)
                                      zin(308) = zin(308) + dzkl*zin(307)

                                      ! i3 = i4 =  307
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  306

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  306

                                      ! do nk = 1,    1

                                      xin(306) = xin(307) + dxkl*xin(305)
                                      yin(306) = yin(307) + dykl*yin(305)
                                      zin(306) = zin(307) + dzkl*zin(305)
                                      ! i4 = i4 + lang+1 =  308

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  307

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  309

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  312

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  311

                                      xin(312) = xin(312) + dxkl*xin(311)
                                      yin(312) = yin(312) + dykl*yin(311)
                                      zin(312) = zin(312) + dzkl*zin(311)

                                      ! i3 = i4 =  311
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  310

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  310

                                      ! do nk = 1,    1

                                      xin(310) = xin(311) + dxkl*xin(309)
                                      yin(310) = yin(311) + dykl*yin(309)
                                      zin(310) = zin(311) + dzkl*zin(309)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  311

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  313

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  316

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  315

                                      xin(316) = xin(316) + dxkl*xin(315)
                                      yin(316) = yin(316) + dykl*yin(315)
                                      zin(316) = zin(316) + dzkl*zin(315)

                                      ! i3 = i4 =  315
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  314

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  314

                                      ! do nk = 1,    1

                                      xin(314) = xin(315) + dxkl*xin(313)
                                      yin(314) = yin(315) + dykl*yin(313)
                                      zin(314) = zin(315) + dzkl*zin(313)
                                      ! i4 = i4 + lang+1 =  316

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  315

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  317

                                      ! nj = nj + 1 =    3

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  320

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  319

                                      xin(320) = xin(320) + dxkl*xin(319)
                                      yin(320) = yin(320) + dykl*yin(319)
                                      zin(320) = zin(320) + dzkl*zin(319)

                                      ! i3 = i4 =  319
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  318

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  318

                                      ! do nk = 1,    1

                                      xin(318) = xin(319) + dxkl*xin(317)
                                      yin(318) = yin(319) + dykl*yin(317)
                                      zin(318) = zin(319) + dzkl*zin(317)
                                      ! i4 = i4 + lang+1 =  320

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  319

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  321

                                      ! nj = nj + 1 =    4

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  321

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  320

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 900! loop over all integrals

                                        l = n - 9*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d11ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 64)*yin(my + 64)*zin(mz + 64) & ! root  2
                                                        + xin(mx + 128)*yin(my + 128)*zin(mz + 128) & ! root  3
                                                        + xin(mx + 192)*yin(my + 192)*zin(mz + 192) & ! root  4
                                                        + xin(mx + 256)*yin(my + 256)*zin(mz + 256)) ! root  5

                                        j = int(n/9) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxj2 = 10
                                  iandj = ish .eq. jsh
                                  maxl = 3
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    if (iandj) maxj2 = i

                                    ii1 = i + loci
                                    ip = (i - 1)*90 ! Stride between functions in i

                                    do j = 1, maxj2

                                      maxl2 = maxl

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*9 + ip ! Add stride between functions in j

                                      do k = 1, 3 ! # of cartesians in k

                                        if (kandl) maxl2 = k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*3 + ijp ! Add stride between functions in k

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
                              deallocate (n11ket)
                              deallocate (xint11ket)

                              end subroutine int3311
                              end submodule
