! The total angular momentum of this class is:           9
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3330_impl
contains
  module subroutine int3330(ff_pair, sf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n03ket(:)
    real(dp), allocatable :: xint33bra(:), xint03ket(:)
    integer(kind=int64) :: nffbra, nsfket
    real(dp) :: scutffbra, scutsfket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxj2
    integer(kind=int64) :: n, i1, i3, i4, i5, k3, k4, nn, nm, km, nj, ni, nl, nk
    real(dp) :: cp10, c10, cp01, c01
    integer(kind=int64) :: nx, ny, nz, mx, my, mz
    integer(kind=int64) :: bra_loop, ket_loop, ijtop, kltop
    real(dp) :: t_expon_ab, t_expon_cd, t_inverse_expon_ab, t_inverse_expon_cd
    real(dp) :: t_expon_abcd_inverse, rho, expe, dum, rab, rcd
    real(dp) :: brrk, akxk, akyk, akzk, t_expon_c, t_expon_d, t_expon_a, t_expon_b
    real(dp) :: xa, ya, za, axak, ayak, azak, axai, ayai, azai, bbrrk, xb, yb, zb, bxbk
    real(dp) :: bybk, bzbk, bxbi, bybi, bzbi, xx, c1x, c2x, c3x, c4x, c1y, c2y, c3y, c4y
    real(dp) :: c1z, c2z, c3z, c4z, f00, u2, duminv, dm2inv, bp01, b00, b10, xcp00, xc00
    real(dp) :: ycp00, zcp00, zc00, yc00, dij, dxij, dyij, dzij
    real(dp) :: buff(9)
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(320), yin(320), zin(320)
    real(dp) :: eri_value(1000)
    real(dp) :: d33bra(100), d03ket(10)
    integer(kind=int64) :: ix(10), jx(10), kx(10), lx(1)
    integer(kind=int64) :: iy(10), jy(10), ky(10), ly(1)
    integer(kind=int64) :: iz(10), jz(10), kz(10), lz(1)
    integer(kind=int64) :: in(7), in1(7), kn(4)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(10), kly(10), klz(10)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

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
    kn(2) = 1
    kn(3) = 2
    kn(4) = 3

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 3
    kx(2) = 0
    kx(3) = 0
    kx(4) = 2
    kx(5) = 2
    kx(6) = 1
    kx(7) = 0
    kx(8) = 1
    kx(9) = 0
    kx(10) = 1

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

    ky(1) = 0
    ky(2) = 3
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 2
    ky(7) = 2
    ky(8) = 0
    ky(9) = 1
    ky(10) = 1

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

    kz(1) = 0
    kz(2) = 0
    kz(3) = 3
    kz(4) = 0
    kz(5) = 1
    kz(6) = 0
    kz(7) = 1
    kz(8) = 2
    kz(9) = 2
    kz(10) = 1

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
    klx(2) = 0
    klx(3) = 0
    klx(4) = 2
    klx(5) = 2
    klx(6) = 1
    klx(7) = 0
    klx(8) = 1
    klx(9) = 0
    klx(10) = 1

    kly(1) = 0
    kly(2) = 3
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 2
    kly(7) = 2
    kly(8) = 0
    kly(9) = 1
    kly(10) = 1

    klz(1) = 0
    klz(2) = 0
    klz(3) = 3
    klz(4) = 0
    klz(5) = 1
    klz(6) = 0
    klz(7) = 1
    klz(8) = 2
    klz(9) = 2
    klz(10) = 1

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

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

    scutsfket = cutoff_schwarz/maxval(sf_pair%xints)
    nsfket = 0
    do ij = 1, res%n_s_shl*res%n_f_shl
      if (sf_pair%xints(ij) .ge. scutsfket) then
        nsfket = nsfket + 1
        xint03ket(nsfket) = sf_pair%xints(ij)
        n03ket(nsfket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*nsfket) .le. nchunksize_int64) nchunksize_int64 = nffbra*nsfket
    ntile = int(nffbra*nsfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*nsfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint03ket, n03ket, ff_pair, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d03ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d33bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,iandj)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nffbra) + 1
              kl_tmp = (iquart - 1)/nffbra + 1

              test = xint33bra(ij_tmp)*xint03ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n03ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_f_shl(ksh_tmp)
                lsh = res%i_s_shl(lsh_tmp)

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

                ! Prepare to begin looping

                eri_value = 0.0_dp

                ket_loop = 0

                ! --- Start looping over primitives ---

                do k = 1, kltop

                  ket_loop = ket_loop + 1

                  t_expon_cd = sf_pair%t_expon_ab(sf_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = sf_pair%expon_b(sf_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = sf_pair%expon_a(sf_pair%pair_loc(kl) + ket_loop)
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

                  d03ket(1) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(2) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(3) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(4) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(5) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(6) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(7) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(8) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(9) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(10) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(17) = xc00
                                      yin(17) = yc00
                                      zin(17) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   18
                                      ! i2 =   17

                                      xin(18) = xcp00*xin(17) + cp10
                                      yin(18) = ycp00*yin(17) + cp10
                                      zin(18) = zcp00*zin(17) + cp10*f00

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

                                      ! i3 = i5 + k2 =   34
                                      ! i5 =   33
                                      ! i4 =   17

                                      xin(34) = xcp00*xin(33) + cp10*xin(17)
                                      yin(34) = ycp00*yin(33) + cp10*yin(17)
                                      zin(34) = zcp00*zin(33) + cp10*zin(17)

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

                                      ! i3 = i5 + k2 =   50
                                      ! i5 =   49
                                      ! i4 =   33

                                      xin(50) = xcp00*xin(49) + cp10*xin(33)
                                      yin(50) = ycp00*yin(49) + cp10*yin(33)
                                      zin(50) = zcp00*zin(49) + cp10*zin(33)

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

                                      ! i3 = i5 + k2 =   54
                                      ! i5 =   53
                                      ! i4 =   49

                                      xin(54) = xcp00*xin(53) + cp10*xin(49)
                                      yin(54) = ycp00*yin(53) + cp10*yin(49)
                                      zin(54) = zcp00*zin(53) + cp10*zin(49)

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

                                      ! i3 = i5 + k2 =   58
                                      ! i5 =   57
                                      ! i4 =   53

                                      xin(58) = xcp00*xin(57) + cp10*xin(53)
                                      yin(58) = ycp00*yin(57) + cp10*yin(53)
                                      zin(58) = zcp00*zin(57) + cp10*zin(53)

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

                                      ! i3 = i5 + k2 =   62
                                      ! i5 =   61
                                      ! i4 =   57

                                      xin(62) = xcp00*xin(61) + cp10*xin(57)
                                      yin(62) = ycp00*yin(61) + cp10*yin(57)
                                      zin(62) = zcp00*zin(61) + cp10*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    2

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    3
                                      ! i3 =    1
                                      ! i4 =    2

                                      xin(3) = cp01*xin(1) + xcp00*xin(2)
                                      yin(3) = cp01*yin(1) + ycp00*yin(2)
                                      zin(3) = cp01*zin(1) + zcp00*zin(2)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   19

                                      xin(19) = xc00*xin(3) + c01*xin(2)
                                      yin(19) = yc00*yin(3) + c01*yin(2)
                                      zin(19) = zc00*zin(3) + c01*zin(2)

                                      ! ------------------

                                      ! i3 = i4 =    2
                                      ! i4 = i5 =    3

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    4
                                      ! i3 =    2
                                      ! i4 =    3

                                      xin(4) = cp01*xin(2) + xcp00*xin(3)
                                      yin(4) = cp01*yin(2) + ycp00*yin(3)
                                      zin(4) = cp01*zin(2) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   20

                                      xin(20) = xc00*xin(4) + c01*xin(3)
                                      yin(20) = yc00*yin(4) + c01*yin(3)
                                      zin(20) = zc00*zin(4) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    4

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   17

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   33

                                      xin(35) = c10*xin(3) + xc00*xin(19) + c01*xin(18)
                                      yin(35) = c10*yin(3) + yc00*yin(19) + c01*yin(18)
                                      zin(35) = c10*zin(3) + zc00*zin(19) + c01*zin(18)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   33

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   49

                                      xin(51) = c10*xin(19) + xc00*xin(35) + c01*xin(34)
                                      yin(51) = c10*yin(19) + yc00*yin(35) + c01*yin(34)
                                      zin(51) = c10*zin(19) + zc00*zin(35) + c01*zin(34)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   49

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   53

                                      xin(55) = c10*xin(35) + xc00*xin(51) + c01*xin(50)
                                      yin(55) = c10*yin(35) + yc00*yin(51) + c01*yin(50)
                                      zin(55) = c10*zin(35) + zc00*zin(51) + c01*zin(50)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   53

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   57

                                      xin(59) = c10*xin(51) + xc00*xin(55) + c01*xin(54)
                                      yin(59) = c10*yin(51) + yc00*yin(55) + c01*yin(54)
                                      zin(59) = c10*zin(51) + zc00*zin(55) + c01*zin(54)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   53
                                      ! i4 = i5 =   57

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   61

                                      xin(63) = c10*xin(55) + xc00*xin(59) + c01*xin(58)
                                      yin(63) = c10*yin(55) + yc00*yin(59) + c01*yin(58)
                                      zin(63) = c10*zin(55) + zc00*zin(59) + c01*zin(58)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   62

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   58

                                      xin(62) = xin(62) + dxij*xin(58)
                                      yin(62) = yin(62) + dyij*yin(58)
                                      zin(62) = zin(62) + dzij*zin(58)

                                      ! i3 = i4 =   58
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   54

                                      xin(58) = xin(58) + dxij*xin(54)
                                      yin(58) = yin(58) + dyij*yin(54)
                                      zin(58) = zin(58) + dzij*zin(54)

                                      ! i3 = i4 =   54
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   50

                                      xin(54) = xin(54) + dxij*xin(50)
                                      yin(54) = yin(54) + dyij*yin(50)
                                      zin(54) = zin(54) + dzij*zin(50)

                                      ! i3 = i4 =   50
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   62

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   58

                                      xin(62) = xin(62) + dxij*xin(58)
                                      yin(62) = yin(62) + dyij*yin(58)
                                      zin(62) = zin(62) + dzij*zin(58)

                                      ! i3 = i4 =   58
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   54

                                      xin(58) = xin(58) + dxij*xin(54)
                                      yin(58) = yin(58) + dyij*yin(54)
                                      zin(58) = zin(58) + dzij*zin(54)

                                      ! i3 = i4 =   54
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   62

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   58

                                      xin(62) = xin(62) + dxij*xin(58)
                                      yin(62) = yin(62) + dyij*yin(58)
                                      zin(62) = zin(62) + dzij*zin(58)

                                      ! i3 = i4 =   58
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    3

                                      xin(6) = xin(18) + dxij*xin(2)
                                      yin(6) = yin(18) + dyij*yin(2)
                                      zin(6) = zin(18) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    2

                                      xin(22) = xin(34) + dxij*xin(18)
                                      yin(22) = yin(34) + dyij*yin(18)
                                      zin(22) = zin(34) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   38

                                      ! ni =    3

                                      xin(38) = xin(50) + dxij*xin(34)
                                      yin(38) = yin(50) + dyij*yin(34)
                                      zin(38) = zin(50) + dzij*zin(34)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    2

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(22) + dxij*xin(6)
                                      yin(10) = yin(22) + dyij*yin(6)
                                      zin(10) = zin(22) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   26

                                      ! ni =    2

                                      xin(26) = xin(38) + dxij*xin(22)
                                      yin(26) = yin(38) + dyij*yin(22)
                                      zin(26) = zin(38) + dzij*zin(22)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    3

                                      xin(42) = xin(54) + dxij*xin(38)
                                      yin(42) = yin(54) + dyij*yin(38)
                                      zin(42) = zin(54) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   14

                                      ! nj =    3

                                      ! i4 = i3 =   14

                                      ! do ni = 1,    3

                                      xin(14) = xin(26) + dxij*xin(10)
                                      yin(14) = yin(26) + dyij*yin(10)
                                      zin(14) = zin(26) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    2

                                      xin(30) = xin(42) + dxij*xin(26)
                                      yin(30) = yin(42) + dyij*yin(26)
                                      zin(30) = zin(42) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    3

                                      xin(46) = xin(58) + dxij*xin(42)
                                      yin(46) = yin(58) + dyij*yin(42)
                                      zin(46) = zin(58) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   62

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   18

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

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

                                      ! nm = nm + 1 =    3

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

                                      ! nm = nm + 1 =    4

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(81) = xc00
                                      yin(81) = yc00
                                      zin(81) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   66

                                      xin(66) = xcp00
                                      yin(66) = ycp00
                                      zin(66) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   82
                                      ! i2 =   81

                                      xin(82) = xcp00*xin(81) + cp10
                                      yin(82) = ycp00*yin(81) + cp10
                                      zin(82) = zcp00*zin(81) + cp10*f00

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

                                      ! i3 = i5 + k2 =   98
                                      ! i5 =   97
                                      ! i4 =   81

                                      xin(98) = xcp00*xin(97) + cp10*xin(81)
                                      yin(98) = ycp00*yin(97) + cp10*yin(81)
                                      zin(98) = zcp00*zin(97) + cp10*zin(81)

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

                                      ! i3 = i5 + k2 =  114
                                      ! i5 =  113
                                      ! i4 =   97

                                      xin(114) = xcp00*xin(113) + cp10*xin(97)
                                      yin(114) = ycp00*yin(113) + cp10*yin(97)
                                      zin(114) = zcp00*zin(113) + cp10*zin(97)

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

                                      ! i3 = i5 + k2 =  118
                                      ! i5 =  117
                                      ! i4 =  113

                                      xin(118) = xcp00*xin(117) + cp10*xin(113)
                                      yin(118) = ycp00*yin(117) + cp10*yin(113)
                                      zin(118) = zcp00*zin(117) + cp10*zin(113)

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

                                      ! i3 = i5 + k2 =  122
                                      ! i5 =  121
                                      ! i4 =  117

                                      xin(122) = xcp00*xin(121) + cp10*xin(117)
                                      yin(122) = ycp00*yin(121) + cp10*yin(117)
                                      zin(122) = zcp00*zin(121) + cp10*zin(117)

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

                                      ! i3 = i5 + k2 =  126
                                      ! i5 =  125
                                      ! i4 =  121

                                      xin(126) = xcp00*xin(125) + cp10*xin(121)
                                      yin(126) = ycp00*yin(125) + cp10*yin(121)
                                      zin(126) = zcp00*zin(125) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   65
                                      ! i4 = i1+k2 =   66

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   67
                                      ! i3 =   65
                                      ! i4 =   66

                                      xin(67) = cp01*xin(65) + xcp00*xin(66)
                                      yin(67) = cp01*yin(65) + ycp00*yin(66)
                                      zin(67) = cp01*zin(65) + zcp00*zin(66)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   83

                                      xin(83) = xc00*xin(67) + c01*xin(66)
                                      yin(83) = yc00*yin(67) + c01*yin(66)
                                      zin(83) = zc00*zin(67) + c01*zin(66)

                                      ! ------------------

                                      ! i3 = i4 =   66
                                      ! i4 = i5 =   67

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   68
                                      ! i3 =   66
                                      ! i4 =   67

                                      xin(68) = cp01*xin(66) + xcp00*xin(67)
                                      yin(68) = cp01*yin(66) + ycp00*yin(67)
                                      zin(68) = cp01*zin(66) + zcp00*zin(67)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   84

                                      xin(84) = xc00*xin(68) + c01*xin(67)
                                      yin(84) = yc00*yin(68) + c01*yin(67)
                                      zin(84) = zc00*zin(68) + c01*zin(67)

                                      ! ------------------

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   68

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   65
                                      ! i4 = i2 =   81

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   97

                                      xin(99) = c10*xin(67) + xc00*xin(83) + c01*xin(82)
                                      yin(99) = c10*yin(67) + yc00*yin(83) + c01*yin(82)
                                      zin(99) = c10*zin(67) + zc00*zin(83) + c01*zin(82)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   97

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  113

                                      xin(115) = c10*xin(83) + xc00*xin(99) + c01*xin(98)
                                      yin(115) = c10*yin(83) + yc00*yin(99) + c01*yin(98)
                                      zin(115) = c10*zin(83) + zc00*zin(99) + c01*zin(98)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  113

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  117

                                      xin(119) = c10*xin(99) + xc00*xin(115) + c01*xin(114)
                                      yin(119) = c10*yin(99) + yc00*yin(115) + c01*yin(114)
                                      zin(119) = c10*zin(99) + zc00*zin(115) + c01*zin(114)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  117

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  121

                                      xin(123) = c10*xin(115) + xc00*xin(119) + c01*xin(118)
                                      yin(123) = c10*yin(115) + yc00*yin(119) + c01*yin(118)
                                      zin(123) = c10*zin(115) + zc00*zin(119) + c01*zin(118)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  117
                                      ! i4 = i5 =  121

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  125

                                      xin(127) = c10*xin(119) + xc00*xin(123) + c01*xin(122)
                                      yin(127) = c10*yin(119) + yc00*yin(123) + c01*yin(122)
                                      zin(127) = c10*zin(119) + zc00*zin(123) + c01*zin(122)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  126

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  122

                                      xin(126) = xin(126) + dxij*xin(122)
                                      yin(126) = yin(126) + dyij*yin(122)
                                      zin(126) = zin(126) + dzij*zin(122)

                                      ! i3 = i4 =  122
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  118

                                      xin(122) = xin(122) + dxij*xin(118)
                                      yin(122) = yin(122) + dyij*yin(118)
                                      zin(122) = zin(122) + dzij*zin(118)

                                      ! i3 = i4 =  118
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  114

                                      xin(118) = xin(118) + dxij*xin(114)
                                      yin(118) = yin(118) + dyij*yin(114)
                                      zin(118) = zin(118) + dzij*zin(114)

                                      ! i3 = i4 =  114
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  126

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  122

                                      xin(126) = xin(126) + dxij*xin(122)
                                      yin(126) = yin(126) + dyij*yin(122)
                                      zin(126) = zin(126) + dzij*zin(122)

                                      ! i3 = i4 =  122
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  118

                                      xin(122) = xin(122) + dxij*xin(118)
                                      yin(122) = yin(122) + dyij*yin(118)
                                      zin(122) = zin(122) + dzij*zin(118)

                                      ! i3 = i4 =  118
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  126

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  122

                                      xin(126) = xin(126) + dxij*xin(122)
                                      yin(126) = yin(126) + dyij*yin(122)
                                      zin(126) = zin(126) + dzij*zin(122)

                                      ! i3 = i4 =  122
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   70

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   70

                                      ! do ni = 1,    3

                                      xin(70) = xin(82) + dxij*xin(66)
                                      yin(70) = yin(82) + dyij*yin(66)
                                      zin(70) = zin(82) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   86

                                      ! ni =    2

                                      xin(86) = xin(98) + dxij*xin(82)
                                      yin(86) = yin(98) + dyij*yin(82)
                                      zin(86) = zin(98) + dzij*zin(82)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  102

                                      ! ni =    3

                                      xin(102) = xin(114) + dxij*xin(98)
                                      yin(102) = yin(114) + dyij*yin(98)
                                      zin(102) = zin(114) + dzij*zin(98)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  118

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   74

                                      ! nj =    2

                                      ! i4 = i3 =   74

                                      ! do ni = 1,    3

                                      xin(74) = xin(86) + dxij*xin(70)
                                      yin(74) = yin(86) + dyij*yin(70)
                                      zin(74) = zin(86) + dzij*zin(70)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    2

                                      xin(90) = xin(102) + dxij*xin(86)
                                      yin(90) = yin(102) + dyij*yin(86)
                                      zin(90) = zin(102) + dzij*zin(86)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                      ! ni =    3

                                      xin(106) = xin(118) + dxij*xin(102)
                                      yin(106) = yin(118) + dyij*yin(102)
                                      zin(106) = zin(118) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  122

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   78

                                      ! nj =    3

                                      ! i4 = i3 =   78

                                      ! do ni = 1,    3

                                      xin(78) = xin(90) + dxij*xin(74)
                                      yin(78) = yin(90) + dyij*yin(74)
                                      zin(78) = zin(90) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    2

                                      xin(94) = xin(106) + dxij*xin(90)
                                      yin(94) = yin(106) + dyij*yin(90)
                                      zin(94) = zin(106) + dzij*zin(90)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  110

                                      ! ni =    3

                                      xin(110) = xin(122) + dxij*xin(106)
                                      yin(110) = yin(122) + dyij*yin(106)
                                      zin(110) = zin(122) + dzij*zin(106)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   82

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

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

                                      ! nm = nm + 1 =    3

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

                                      ! nm = nm + 1 =    4

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(145) = xc00
                                      yin(145) = yc00
                                      zin(145) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  130

                                      xin(130) = xcp00
                                      yin(130) = ycp00
                                      zin(130) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  146
                                      ! i2 =  145

                                      xin(146) = xcp00*xin(145) + cp10
                                      yin(146) = ycp00*yin(145) + cp10
                                      zin(146) = zcp00*zin(145) + cp10*f00

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

                                      ! i3 = i5 + k2 =  162
                                      ! i5 =  161
                                      ! i4 =  145

                                      xin(162) = xcp00*xin(161) + cp10*xin(145)
                                      yin(162) = ycp00*yin(161) + cp10*yin(145)
                                      zin(162) = zcp00*zin(161) + cp10*zin(145)

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

                                      ! i3 = i5 + k2 =  178
                                      ! i5 =  177
                                      ! i4 =  161

                                      xin(178) = xcp00*xin(177) + cp10*xin(161)
                                      yin(178) = ycp00*yin(177) + cp10*yin(161)
                                      zin(178) = zcp00*zin(177) + cp10*zin(161)

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

                                      ! i3 = i5 + k2 =  182
                                      ! i5 =  181
                                      ! i4 =  177

                                      xin(182) = xcp00*xin(181) + cp10*xin(177)
                                      yin(182) = ycp00*yin(181) + cp10*yin(177)
                                      zin(182) = zcp00*zin(181) + cp10*zin(177)

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

                                      ! i3 = i5 + k2 =  186
                                      ! i5 =  185
                                      ! i4 =  181

                                      xin(186) = xcp00*xin(185) + cp10*xin(181)
                                      yin(186) = ycp00*yin(185) + cp10*yin(181)
                                      zin(186) = zcp00*zin(185) + cp10*zin(181)

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

                                      ! i3 = i5 + k2 =  190
                                      ! i5 =  189
                                      ! i4 =  185

                                      xin(190) = xcp00*xin(189) + cp10*xin(185)
                                      yin(190) = ycp00*yin(189) + cp10*yin(185)
                                      zin(190) = zcp00*zin(189) + cp10*zin(185)

                                      ! ------------------

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  129
                                      ! i4 = i1+k2 =  130

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  131
                                      ! i3 =  129
                                      ! i4 =  130

                                      xin(131) = cp01*xin(129) + xcp00*xin(130)
                                      yin(131) = cp01*yin(129) + ycp00*yin(130)
                                      zin(131) = cp01*zin(129) + zcp00*zin(130)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  147

                                      xin(147) = xc00*xin(131) + c01*xin(130)
                                      yin(147) = yc00*yin(131) + c01*yin(130)
                                      zin(147) = zc00*zin(131) + c01*zin(130)

                                      ! ------------------

                                      ! i3 = i4 =  130
                                      ! i4 = i5 =  131

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  132
                                      ! i3 =  130
                                      ! i4 =  131

                                      xin(132) = cp01*xin(130) + xcp00*xin(131)
                                      yin(132) = cp01*yin(130) + ycp00*yin(131)
                                      zin(132) = cp01*zin(130) + zcp00*zin(131)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  148

                                      xin(148) = xc00*xin(132) + c01*xin(131)
                                      yin(148) = yc00*yin(132) + c01*yin(131)
                                      zin(148) = zc00*zin(132) + c01*zin(131)

                                      ! ------------------

                                      ! i3 = i4 =  131
                                      ! i4 = i5 =  132

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  129
                                      ! i4 = i2 =  145

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  161

                                      xin(163) = c10*xin(131) + xc00*xin(147) + c01*xin(146)
                                      yin(163) = c10*yin(131) + yc00*yin(147) + c01*yin(146)
                                      zin(163) = c10*zin(131) + zc00*zin(147) + c01*zin(146)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  161

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  177

                                      xin(179) = c10*xin(147) + xc00*xin(163) + c01*xin(162)
                                      yin(179) = c10*yin(147) + yc00*yin(163) + c01*yin(162)
                                      zin(179) = c10*zin(147) + zc00*zin(163) + c01*zin(162)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  161
                                      ! i4 = i5 =  177

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  181

                                      xin(183) = c10*xin(163) + xc00*xin(179) + c01*xin(178)
                                      yin(183) = c10*yin(163) + yc00*yin(179) + c01*yin(178)
                                      zin(183) = c10*zin(163) + zc00*zin(179) + c01*zin(178)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  177
                                      ! i4 = i5 =  181

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  185

                                      xin(187) = c10*xin(179) + xc00*xin(183) + c01*xin(182)
                                      yin(187) = c10*yin(179) + yc00*yin(183) + c01*yin(182)
                                      zin(187) = c10*zin(179) + zc00*zin(183) + c01*zin(182)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  189

                                      xin(191) = c10*xin(183) + xc00*xin(187) + c01*xin(186)
                                      yin(191) = c10*yin(183) + yc00*yin(187) + c01*yin(186)
                                      zin(191) = c10*zin(183) + zc00*zin(187) + c01*zin(186)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(190) = xin(190) + dxij*xin(186)
                                      yin(190) = yin(190) + dyij*yin(186)
                                      zin(190) = zin(190) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  182

                                      xin(186) = xin(186) + dxij*xin(182)
                                      yin(186) = yin(186) + dyij*yin(182)
                                      zin(186) = zin(186) + dzij*zin(182)

                                      ! i3 = i4 =  182
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  178

                                      xin(182) = xin(182) + dxij*xin(178)
                                      yin(182) = yin(182) + dyij*yin(178)
                                      zin(182) = zin(182) + dzij*zin(178)

                                      ! i3 = i4 =  178
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(190) = xin(190) + dxij*xin(186)
                                      yin(190) = yin(190) + dyij*yin(186)
                                      zin(190) = zin(190) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  182

                                      xin(186) = xin(186) + dxij*xin(182)
                                      yin(186) = yin(186) + dyij*yin(182)
                                      zin(186) = zin(186) + dzij*zin(182)

                                      ! i3 = i4 =  182
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(190) = xin(190) + dxij*xin(186)
                                      yin(190) = yin(190) + dyij*yin(186)
                                      zin(190) = zin(190) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  134

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  134

                                      ! do ni = 1,    3

                                      xin(134) = xin(146) + dxij*xin(130)
                                      yin(134) = yin(146) + dyij*yin(130)
                                      zin(134) = zin(146) + dzij*zin(130)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  150

                                      ! ni =    2

                                      xin(150) = xin(162) + dxij*xin(146)
                                      yin(150) = yin(162) + dyij*yin(146)
                                      zin(150) = zin(162) + dzij*zin(146)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  166

                                      ! ni =    3

                                      xin(166) = xin(178) + dxij*xin(162)
                                      yin(166) = yin(178) + dyij*yin(162)
                                      zin(166) = zin(178) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  182

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  138

                                      ! nj =    2

                                      ! i4 = i3 =  138

                                      ! do ni = 1,    3

                                      xin(138) = xin(150) + dxij*xin(134)
                                      yin(138) = yin(150) + dyij*yin(134)
                                      zin(138) = zin(150) + dzij*zin(134)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  154

                                      ! ni =    2

                                      xin(154) = xin(166) + dxij*xin(150)
                                      yin(154) = yin(166) + dyij*yin(150)
                                      zin(154) = zin(166) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  170

                                      ! ni =    3

                                      xin(170) = xin(182) + dxij*xin(166)
                                      yin(170) = yin(182) + dyij*yin(166)
                                      zin(170) = zin(182) + dzij*zin(166)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  186

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  142

                                      ! nj =    3

                                      ! i4 = i3 =  142

                                      ! do ni = 1,    3

                                      xin(142) = xin(154) + dxij*xin(138)
                                      yin(142) = yin(154) + dyij*yin(138)
                                      zin(142) = zin(154) + dzij*zin(138)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  158

                                      ! ni =    2

                                      xin(158) = xin(170) + dxij*xin(154)
                                      yin(158) = yin(170) + dyij*yin(154)
                                      zin(158) = zin(170) + dzij*zin(154)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  174

                                      ! ni =    3

                                      xin(174) = xin(186) + dxij*xin(170)
                                      yin(174) = yin(186) + dyij*yin(170)
                                      zin(174) = zin(186) + dzij*zin(170)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  190

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  146

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

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

                                      ! nm = nm + 1 =    3

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

                                      ! nm = nm + 1 =    4

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(209) = xc00
                                      yin(209) = yc00
                                      zin(209) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  194

                                      xin(194) = xcp00
                                      yin(194) = ycp00
                                      zin(194) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  210
                                      ! i2 =  209

                                      xin(210) = xcp00*xin(209) + cp10
                                      yin(210) = ycp00*yin(209) + cp10
                                      zin(210) = zcp00*zin(209) + cp10*f00

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

                                      ! i3 = i5 + k2 =  226
                                      ! i5 =  225
                                      ! i4 =  209

                                      xin(226) = xcp00*xin(225) + cp10*xin(209)
                                      yin(226) = ycp00*yin(225) + cp10*yin(209)
                                      zin(226) = zcp00*zin(225) + cp10*zin(209)

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

                                      ! i3 = i5 + k2 =  242
                                      ! i5 =  241
                                      ! i4 =  225

                                      xin(242) = xcp00*xin(241) + cp10*xin(225)
                                      yin(242) = ycp00*yin(241) + cp10*yin(225)
                                      zin(242) = zcp00*zin(241) + cp10*zin(225)

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

                                      ! i3 = i5 + k2 =  246
                                      ! i5 =  245
                                      ! i4 =  241

                                      xin(246) = xcp00*xin(245) + cp10*xin(241)
                                      yin(246) = ycp00*yin(245) + cp10*yin(241)
                                      zin(246) = zcp00*zin(245) + cp10*zin(241)

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

                                      ! i3 = i5 + k2 =  250
                                      ! i5 =  249
                                      ! i4 =  245

                                      xin(250) = xcp00*xin(249) + cp10*xin(245)
                                      yin(250) = ycp00*yin(249) + cp10*yin(245)
                                      zin(250) = zcp00*zin(249) + cp10*zin(245)

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

                                      ! i3 = i5 + k2 =  254
                                      ! i5 =  253
                                      ! i4 =  249

                                      xin(254) = xcp00*xin(253) + cp10*xin(249)
                                      yin(254) = ycp00*yin(253) + cp10*yin(249)
                                      zin(254) = zcp00*zin(253) + cp10*zin(249)

                                      ! ------------------

                                      ! i3 = i4 =  249
                                      ! i4 = i5 =  253

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  193
                                      ! i4 = i1+k2 =  194

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  195
                                      ! i3 =  193
                                      ! i4 =  194

                                      xin(195) = cp01*xin(193) + xcp00*xin(194)
                                      yin(195) = cp01*yin(193) + ycp00*yin(194)
                                      zin(195) = cp01*zin(193) + zcp00*zin(194)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  211

                                      xin(211) = xc00*xin(195) + c01*xin(194)
                                      yin(211) = yc00*yin(195) + c01*yin(194)
                                      zin(211) = zc00*zin(195) + c01*zin(194)

                                      ! ------------------

                                      ! i3 = i4 =  194
                                      ! i4 = i5 =  195

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  196
                                      ! i3 =  194
                                      ! i4 =  195

                                      xin(196) = cp01*xin(194) + xcp00*xin(195)
                                      yin(196) = cp01*yin(194) + ycp00*yin(195)
                                      zin(196) = cp01*zin(194) + zcp00*zin(195)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  212

                                      xin(212) = xc00*xin(196) + c01*xin(195)
                                      yin(212) = yc00*yin(196) + c01*yin(195)
                                      zin(212) = zc00*zin(196) + c01*zin(195)

                                      ! ------------------

                                      ! i3 = i4 =  195
                                      ! i4 = i5 =  196

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  209

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  225

                                      xin(227) = c10*xin(195) + xc00*xin(211) + c01*xin(210)
                                      yin(227) = c10*yin(195) + yc00*yin(211) + c01*yin(210)
                                      zin(227) = c10*zin(195) + zc00*zin(211) + c01*zin(210)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  209
                                      ! i4 = i5 =  225

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  241

                                      xin(243) = c10*xin(211) + xc00*xin(227) + c01*xin(226)
                                      yin(243) = c10*yin(211) + yc00*yin(227) + c01*yin(226)
                                      zin(243) = c10*zin(211) + zc00*zin(227) + c01*zin(226)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  225
                                      ! i4 = i5 =  241

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  245

                                      xin(247) = c10*xin(227) + xc00*xin(243) + c01*xin(242)
                                      yin(247) = c10*yin(227) + yc00*yin(243) + c01*yin(242)
                                      zin(247) = c10*zin(227) + zc00*zin(243) + c01*zin(242)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  245

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  249

                                      xin(251) = c10*xin(243) + xc00*xin(247) + c01*xin(246)
                                      yin(251) = c10*yin(243) + yc00*yin(247) + c01*yin(246)
                                      zin(251) = c10*zin(243) + zc00*zin(247) + c01*zin(246)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  245
                                      ! i4 = i5 =  249

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  253

                                      xin(255) = c10*xin(247) + xc00*xin(251) + c01*xin(250)
                                      yin(255) = c10*yin(247) + yc00*yin(251) + c01*yin(250)
                                      zin(255) = c10*zin(247) + zc00*zin(251) + c01*zin(250)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  249
                                      ! i4 = i5 =  253

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  254

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  250

                                      xin(254) = xin(254) + dxij*xin(250)
                                      yin(254) = yin(254) + dyij*yin(250)
                                      zin(254) = zin(254) + dzij*zin(250)

                                      ! i3 = i4 =  250
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  246

                                      xin(250) = xin(250) + dxij*xin(246)
                                      yin(250) = yin(250) + dyij*yin(246)
                                      zin(250) = zin(250) + dzij*zin(246)

                                      ! i3 = i4 =  246
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  242

                                      xin(246) = xin(246) + dxij*xin(242)
                                      yin(246) = yin(246) + dyij*yin(242)
                                      zin(246) = zin(246) + dzij*zin(242)

                                      ! i3 = i4 =  242
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  254

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  250

                                      xin(254) = xin(254) + dxij*xin(250)
                                      yin(254) = yin(254) + dyij*yin(250)
                                      zin(254) = zin(254) + dzij*zin(250)

                                      ! i3 = i4 =  250
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  246

                                      xin(250) = xin(250) + dxij*xin(246)
                                      yin(250) = yin(250) + dyij*yin(246)
                                      zin(250) = zin(250) + dzij*zin(246)

                                      ! i3 = i4 =  246
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  254

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  250

                                      xin(254) = xin(254) + dxij*xin(250)
                                      yin(254) = yin(254) + dyij*yin(250)
                                      zin(254) = zin(254) + dzij*zin(250)

                                      ! i3 = i4 =  250
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  198

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  198

                                      ! do ni = 1,    3

                                      xin(198) = xin(210) + dxij*xin(194)
                                      yin(198) = yin(210) + dyij*yin(194)
                                      zin(198) = zin(210) + dzij*zin(194)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                      ! ni =    2

                                      xin(214) = xin(226) + dxij*xin(210)
                                      yin(214) = yin(226) + dyij*yin(210)
                                      zin(214) = zin(226) + dzij*zin(210)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  230

                                      ! ni =    3

                                      xin(230) = xin(242) + dxij*xin(226)
                                      yin(230) = yin(242) + dyij*yin(226)
                                      zin(230) = zin(242) + dzij*zin(226)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  246

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  202

                                      ! nj =    2

                                      ! i4 = i3 =  202

                                      ! do ni = 1,    3

                                      xin(202) = xin(214) + dxij*xin(198)
                                      yin(202) = yin(214) + dyij*yin(198)
                                      zin(202) = zin(214) + dzij*zin(198)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  218

                                      ! ni =    2

                                      xin(218) = xin(230) + dxij*xin(214)
                                      yin(218) = yin(230) + dyij*yin(214)
                                      zin(218) = zin(230) + dzij*zin(214)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  234

                                      ! ni =    3

                                      xin(234) = xin(246) + dxij*xin(230)
                                      yin(234) = yin(246) + dyij*yin(230)
                                      zin(234) = zin(246) + dzij*zin(230)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  250

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  206

                                      ! nj =    3

                                      ! i4 = i3 =  206

                                      ! do ni = 1,    3

                                      xin(206) = xin(218) + dxij*xin(202)
                                      yin(206) = yin(218) + dyij*yin(202)
                                      zin(206) = zin(218) + dzij*zin(202)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  222

                                      ! ni =    2

                                      xin(222) = xin(234) + dxij*xin(218)
                                      yin(222) = yin(234) + dyij*yin(218)
                                      zin(222) = zin(234) + dzij*zin(218)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  238

                                      ! ni =    3

                                      xin(238) = xin(250) + dxij*xin(234)
                                      yin(238) = yin(250) + dyij*yin(234)
                                      zin(238) = zin(250) + dzij*zin(234)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  254

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  210

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

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

                                      ! nm = nm + 1 =    3

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

                                      ! nm = nm + 1 =    4

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(273) = xc00
                                      yin(273) = yc00
                                      zin(273) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  258

                                      xin(258) = xcp00
                                      yin(258) = ycp00
                                      zin(258) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  274
                                      ! i2 =  273

                                      xin(274) = xcp00*xin(273) + cp10
                                      yin(274) = ycp00*yin(273) + cp10
                                      zin(274) = zcp00*zin(273) + cp10*f00

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

                                      ! i3 = i5 + k2 =  290
                                      ! i5 =  289
                                      ! i4 =  273

                                      xin(290) = xcp00*xin(289) + cp10*xin(273)
                                      yin(290) = ycp00*yin(289) + cp10*yin(273)
                                      zin(290) = zcp00*zin(289) + cp10*zin(273)

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

                                      ! i3 = i5 + k2 =  306
                                      ! i5 =  305
                                      ! i4 =  289

                                      xin(306) = xcp00*xin(305) + cp10*xin(289)
                                      yin(306) = ycp00*yin(305) + cp10*yin(289)
                                      zin(306) = zcp00*zin(305) + cp10*zin(289)

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

                                      ! i3 = i5 + k2 =  310
                                      ! i5 =  309
                                      ! i4 =  305

                                      xin(310) = xcp00*xin(309) + cp10*xin(305)
                                      yin(310) = ycp00*yin(309) + cp10*yin(305)
                                      zin(310) = zcp00*zin(309) + cp10*zin(305)

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

                                      ! i3 = i5 + k2 =  314
                                      ! i5 =  313
                                      ! i4 =  309

                                      xin(314) = xcp00*xin(313) + cp10*xin(309)
                                      yin(314) = ycp00*yin(313) + cp10*yin(309)
                                      zin(314) = zcp00*zin(313) + cp10*zin(309)

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

                                      ! i3 = i5 + k2 =  318
                                      ! i5 =  317
                                      ! i4 =  313

                                      xin(318) = xcp00*xin(317) + cp10*xin(313)
                                      yin(318) = ycp00*yin(317) + cp10*yin(313)
                                      zin(318) = zcp00*zin(317) + cp10*zin(313)

                                      ! ------------------

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  317

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  257
                                      ! i4 = i1+k2 =  258

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  259
                                      ! i3 =  257
                                      ! i4 =  258

                                      xin(259) = cp01*xin(257) + xcp00*xin(258)
                                      yin(259) = cp01*yin(257) + ycp00*yin(258)
                                      zin(259) = cp01*zin(257) + zcp00*zin(258)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  275

                                      xin(275) = xc00*xin(259) + c01*xin(258)
                                      yin(275) = yc00*yin(259) + c01*yin(258)
                                      zin(275) = zc00*zin(259) + c01*zin(258)

                                      ! ------------------

                                      ! i3 = i4 =  258
                                      ! i4 = i5 =  259

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  260
                                      ! i3 =  258
                                      ! i4 =  259

                                      xin(260) = cp01*xin(258) + xcp00*xin(259)
                                      yin(260) = cp01*yin(258) + ycp00*yin(259)
                                      zin(260) = cp01*zin(258) + zcp00*zin(259)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  276

                                      xin(276) = xc00*xin(260) + c01*xin(259)
                                      yin(276) = yc00*yin(260) + c01*yin(259)
                                      zin(276) = zc00*zin(260) + c01*zin(259)

                                      ! ------------------

                                      ! i3 = i4 =  259
                                      ! i4 = i5 =  260

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  257
                                      ! i4 = i2 =  273

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  289

                                      xin(291) = c10*xin(259) + xc00*xin(275) + c01*xin(274)
                                      yin(291) = c10*yin(259) + yc00*yin(275) + c01*yin(274)
                                      zin(291) = c10*zin(259) + zc00*zin(275) + c01*zin(274)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  289

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  305

                                      xin(307) = c10*xin(275) + xc00*xin(291) + c01*xin(290)
                                      yin(307) = c10*yin(275) + yc00*yin(291) + c01*yin(290)
                                      zin(307) = c10*zin(275) + zc00*zin(291) + c01*zin(290)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  289
                                      ! i4 = i5 =  305

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  309

                                      xin(311) = c10*xin(291) + xc00*xin(307) + c01*xin(306)
                                      yin(311) = c10*yin(291) + yc00*yin(307) + c01*yin(306)
                                      zin(311) = c10*zin(291) + zc00*zin(307) + c01*zin(306)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  305
                                      ! i4 = i5 =  309

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  313

                                      xin(315) = c10*xin(307) + xc00*xin(311) + c01*xin(310)
                                      yin(315) = c10*yin(307) + yc00*yin(311) + c01*yin(310)
                                      zin(315) = c10*zin(307) + zc00*zin(311) + c01*zin(310)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  309
                                      ! i4 = i5 =  313

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  317

                                      xin(319) = c10*xin(311) + xc00*xin(315) + c01*xin(314)
                                      yin(319) = c10*yin(311) + yc00*yin(315) + c01*yin(314)
                                      zin(319) = c10*zin(311) + zc00*zin(315) + c01*zin(314)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  317

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  318

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  314

                                      xin(318) = xin(318) + dxij*xin(314)
                                      yin(318) = yin(318) + dyij*yin(314)
                                      zin(318) = zin(318) + dzij*zin(314)

                                      ! i3 = i4 =  314
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  310

                                      xin(314) = xin(314) + dxij*xin(310)
                                      yin(314) = yin(314) + dyij*yin(310)
                                      zin(314) = zin(314) + dzij*zin(310)

                                      ! i3 = i4 =  310
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  306

                                      xin(310) = xin(310) + dxij*xin(306)
                                      yin(310) = yin(310) + dyij*yin(306)
                                      zin(310) = zin(310) + dzij*zin(306)

                                      ! i3 = i4 =  306
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  318

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  314

                                      xin(318) = xin(318) + dxij*xin(314)
                                      yin(318) = yin(318) + dyij*yin(314)
                                      zin(318) = zin(318) + dzij*zin(314)

                                      ! i3 = i4 =  314
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  310

                                      xin(314) = xin(314) + dxij*xin(310)
                                      yin(314) = yin(314) + dyij*yin(310)
                                      zin(314) = zin(314) + dzij*zin(310)

                                      ! i3 = i4 =  310
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  318

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  314

                                      xin(318) = xin(318) + dxij*xin(314)
                                      yin(318) = yin(318) + dyij*yin(314)
                                      zin(318) = zin(318) + dzij*zin(314)

                                      ! i3 = i4 =  314
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  262

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  262

                                      ! do ni = 1,    3

                                      xin(262) = xin(274) + dxij*xin(258)
                                      yin(262) = yin(274) + dyij*yin(258)
                                      zin(262) = zin(274) + dzij*zin(258)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  278

                                      ! ni =    2

                                      xin(278) = xin(290) + dxij*xin(274)
                                      yin(278) = yin(290) + dyij*yin(274)
                                      zin(278) = zin(290) + dzij*zin(274)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  294

                                      ! ni =    3

                                      xin(294) = xin(306) + dxij*xin(290)
                                      yin(294) = yin(306) + dyij*yin(290)
                                      zin(294) = zin(306) + dzij*zin(290)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  310

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  266

                                      ! nj =    2

                                      ! i4 = i3 =  266

                                      ! do ni = 1,    3

                                      xin(266) = xin(278) + dxij*xin(262)
                                      yin(266) = yin(278) + dyij*yin(262)
                                      zin(266) = zin(278) + dzij*zin(262)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  282

                                      ! ni =    2

                                      xin(282) = xin(294) + dxij*xin(278)
                                      yin(282) = yin(294) + dyij*yin(278)
                                      zin(282) = zin(294) + dzij*zin(278)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  298

                                      ! ni =    3

                                      xin(298) = xin(310) + dxij*xin(294)
                                      yin(298) = yin(310) + dyij*yin(294)
                                      zin(298) = zin(310) + dzij*zin(294)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  314

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  270

                                      ! nj =    3

                                      ! i4 = i3 =  270

                                      ! do ni = 1,    3

                                      xin(270) = xin(282) + dxij*xin(266)
                                      yin(270) = yin(282) + dyij*yin(266)
                                      zin(270) = zin(282) + dzij*zin(266)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  286

                                      ! ni =    2

                                      xin(286) = xin(298) + dxij*xin(282)
                                      yin(286) = yin(298) + dyij*yin(282)
                                      zin(286) = zin(298) + dzij*zin(282)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  302

                                      ! ni =    3

                                      xin(302) = xin(314) + dxij*xin(298)
                                      yin(302) = yin(314) + dyij*yin(298)
                                      zin(302) = zin(314) + dzij*zin(298)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  318

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  274

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

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

                                      ! nm = nm + 1 =    3

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

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  320

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 1000! loop over all integrals

                                        l = n - 10*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d03ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 64)*yin(my + 64)*zin(mz + 64) & ! root  2
                                                        + xin(mx + 128)*yin(my + 128)*zin(mz + 128) & ! root  3
                                                        + xin(mx + 192)*yin(my + 192)*zin(mz + 192) & ! root  4
                                                        + xin(mx + 256)*yin(my + 256)*zin(mz + 256)) ! root  5

                                        j = int(n/10) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*100 ! Stride between functions in i

                                    do j = 1, maxj2

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*10 + ip ! Add stride between functions in j

                                      do k = 1, 10 ! # of cartesians in k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*1 + ijp ! Add stride between functions in k

                                        do l = 1, 1! # of cartesians in l

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
                              deallocate (n03ket)
                              deallocate (xint03ket)

                              end subroutine int3330
                              end submodule
