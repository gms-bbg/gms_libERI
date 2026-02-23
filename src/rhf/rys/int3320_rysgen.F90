! The total angular momentum of this class is:           8
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3320_impl
contains
  module subroutine int3320(ff_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n02ket(:)
    real(dp), allocatable :: xint33bra(:), xint02ket(:)
    integer(kind=int64) :: nffbra, nsdket
    real(dp) :: scutffbra, scutsdket, test
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
    real(dp) :: xin(240), yin(240), zin(240)
    real(dp) :: eri_value(600)
    real(dp) :: d33bra(100), d02ket(6)
    integer(kind=int64) :: ix(10), jx(10), kx(6), lx(1)
    integer(kind=int64) :: iy(10), jy(10), ky(6), ly(1)
    integer(kind=int64) :: iz(10), jz(10), kz(6), lz(1)
    integer(kind=int64) :: in(7), in1(7), kn(3)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(6), kly(6), klz(6)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 37
    in1(5) = 40
    in1(6) = 43
    in1(7) = 46

    kn(1) = 0
    kn(2) = 1
    kn(3) = 2

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 2
    kx(2) = 0
    kx(3) = 0
    kx(4) = 1
    kx(5) = 1
    kx(6) = 0

    jx(1) = 9
    jx(2) = 0
    jx(3) = 0
    jx(4) = 6
    jx(5) = 6
    jx(6) = 3
    jx(7) = 0
    jx(8) = 3
    jx(9) = 0
    jx(10) = 3

    ix(1) = 37
    ix(2) = 1
    ix(3) = 1
    ix(4) = 25
    ix(5) = 25
    ix(6) = 13
    ix(7) = 1
    ix(8) = 13
    ix(9) = 1
    ix(10) = 13

    ! y-arrays

    ly(1) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 1

    jy(1) = 0
    jy(2) = 9
    jy(3) = 0
    jy(4) = 3
    jy(5) = 0
    jy(6) = 6
    jy(7) = 6
    jy(8) = 0
    jy(9) = 3
    jy(10) = 3

    iy(1) = 1
    iy(2) = 37
    iy(3) = 1
    iy(4) = 13
    iy(5) = 1
    iy(6) = 25
    iy(7) = 25
    iy(8) = 1
    iy(9) = 13
    iy(10) = 13

    ! z-arrays

    lz(1) = 0

    kz(1) = 0
    kz(2) = 0
    kz(3) = 2
    kz(4) = 0
    kz(5) = 1
    kz(6) = 1

    jz(1) = 0
    jz(2) = 0
    jz(3) = 9
    jz(4) = 0
    jz(5) = 3
    jz(6) = 0
    jz(7) = 3
    jz(8) = 6
    jz(9) = 6
    jz(10) = 3

    iz(1) = 1
    iz(2) = 1
    iz(3) = 37
    iz(4) = 1
    iz(5) = 13
    iz(6) = 1
    iz(7) = 13
    iz(8) = 25
    iz(9) = 25
    iz(10) = 13

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 46
    ijx(2) = 37
    ijx(3) = 37
    ijx(4) = 43
    ijx(5) = 43
    ijx(6) = 40
    ijx(7) = 37
    ijx(8) = 40
    ijx(9) = 37
    ijx(10) = 40
    ijx(11) = 10
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 7
    ijx(15) = 7
    ijx(16) = 4
    ijx(17) = 1
    ijx(18) = 4
    ijx(19) = 1
    ijx(20) = 4
    ijx(21) = 10
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 7
    ijx(25) = 7
    ijx(26) = 4
    ijx(27) = 1
    ijx(28) = 4
    ijx(29) = 1
    ijx(30) = 4
    ijx(31) = 34
    ijx(32) = 25
    ijx(33) = 25
    ijx(34) = 31
    ijx(35) = 31
    ijx(36) = 28
    ijx(37) = 25
    ijx(38) = 28
    ijx(39) = 25
    ijx(40) = 28
    ijx(41) = 34
    ijx(42) = 25
    ijx(43) = 25
    ijx(44) = 31
    ijx(45) = 31
    ijx(46) = 28
    ijx(47) = 25
    ijx(48) = 28
    ijx(49) = 25
    ijx(50) = 28
    ijx(51) = 22
    ijx(52) = 13
    ijx(53) = 13
    ijx(54) = 19
    ijx(55) = 19
    ijx(56) = 16
    ijx(57) = 13
    ijx(58) = 16
    ijx(59) = 13
    ijx(60) = 16
    ijx(61) = 10
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 7
    ijx(65) = 7
    ijx(66) = 4
    ijx(67) = 1
    ijx(68) = 4
    ijx(69) = 1
    ijx(70) = 4
    ijx(71) = 22
    ijx(72) = 13
    ijx(73) = 13
    ijx(74) = 19
    ijx(75) = 19
    ijx(76) = 16
    ijx(77) = 13
    ijx(78) = 16
    ijx(79) = 13
    ijx(80) = 16
    ijx(81) = 10
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 7
    ijx(85) = 7
    ijx(86) = 4
    ijx(87) = 1
    ijx(88) = 4
    ijx(89) = 1
    ijx(90) = 4
    ijx(91) = 22
    ijx(92) = 13
    ijx(93) = 13
    ijx(94) = 19
    ijx(95) = 19
    ijx(96) = 16
    ijx(97) = 13
    ijx(98) = 16
    ijx(99) = 13
    ijx(100) = 16

    ijy(1) = 1
    ijy(2) = 10
    ijy(3) = 1
    ijy(4) = 4
    ijy(5) = 1
    ijy(6) = 7
    ijy(7) = 7
    ijy(8) = 1
    ijy(9) = 4
    ijy(10) = 4
    ijy(11) = 37
    ijy(12) = 46
    ijy(13) = 37
    ijy(14) = 40
    ijy(15) = 37
    ijy(16) = 43
    ijy(17) = 43
    ijy(18) = 37
    ijy(19) = 40
    ijy(20) = 40
    ijy(21) = 1
    ijy(22) = 10
    ijy(23) = 1
    ijy(24) = 4
    ijy(25) = 1
    ijy(26) = 7
    ijy(27) = 7
    ijy(28) = 1
    ijy(29) = 4
    ijy(30) = 4
    ijy(31) = 13
    ijy(32) = 22
    ijy(33) = 13
    ijy(34) = 16
    ijy(35) = 13
    ijy(36) = 19
    ijy(37) = 19
    ijy(38) = 13
    ijy(39) = 16
    ijy(40) = 16
    ijy(41) = 1
    ijy(42) = 10
    ijy(43) = 1
    ijy(44) = 4
    ijy(45) = 1
    ijy(46) = 7
    ijy(47) = 7
    ijy(48) = 1
    ijy(49) = 4
    ijy(50) = 4
    ijy(51) = 25
    ijy(52) = 34
    ijy(53) = 25
    ijy(54) = 28
    ijy(55) = 25
    ijy(56) = 31
    ijy(57) = 31
    ijy(58) = 25
    ijy(59) = 28
    ijy(60) = 28
    ijy(61) = 25
    ijy(62) = 34
    ijy(63) = 25
    ijy(64) = 28
    ijy(65) = 25
    ijy(66) = 31
    ijy(67) = 31
    ijy(68) = 25
    ijy(69) = 28
    ijy(70) = 28
    ijy(71) = 1
    ijy(72) = 10
    ijy(73) = 1
    ijy(74) = 4
    ijy(75) = 1
    ijy(76) = 7
    ijy(77) = 7
    ijy(78) = 1
    ijy(79) = 4
    ijy(80) = 4
    ijy(81) = 13
    ijy(82) = 22
    ijy(83) = 13
    ijy(84) = 16
    ijy(85) = 13
    ijy(86) = 19
    ijy(87) = 19
    ijy(88) = 13
    ijy(89) = 16
    ijy(90) = 16
    ijy(91) = 13
    ijy(92) = 22
    ijy(93) = 13
    ijy(94) = 16
    ijy(95) = 13
    ijy(96) = 19
    ijy(97) = 19
    ijy(98) = 13
    ijy(99) = 16
    ijy(100) = 16

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 10
    ijz(4) = 1
    ijz(5) = 4
    ijz(6) = 1
    ijz(7) = 4
    ijz(8) = 7
    ijz(9) = 7
    ijz(10) = 4
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 10
    ijz(14) = 1
    ijz(15) = 4
    ijz(16) = 1
    ijz(17) = 4
    ijz(18) = 7
    ijz(19) = 7
    ijz(20) = 4
    ijz(21) = 37
    ijz(22) = 37
    ijz(23) = 46
    ijz(24) = 37
    ijz(25) = 40
    ijz(26) = 37
    ijz(27) = 40
    ijz(28) = 43
    ijz(29) = 43
    ijz(30) = 40
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 10
    ijz(34) = 1
    ijz(35) = 4
    ijz(36) = 1
    ijz(37) = 4
    ijz(38) = 7
    ijz(39) = 7
    ijz(40) = 4
    ijz(41) = 13
    ijz(42) = 13
    ijz(43) = 22
    ijz(44) = 13
    ijz(45) = 16
    ijz(46) = 13
    ijz(47) = 16
    ijz(48) = 19
    ijz(49) = 19
    ijz(50) = 16
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 10
    ijz(54) = 1
    ijz(55) = 4
    ijz(56) = 1
    ijz(57) = 4
    ijz(58) = 7
    ijz(59) = 7
    ijz(60) = 4
    ijz(61) = 13
    ijz(62) = 13
    ijz(63) = 22
    ijz(64) = 13
    ijz(65) = 16
    ijz(66) = 13
    ijz(67) = 16
    ijz(68) = 19
    ijz(69) = 19
    ijz(70) = 16
    ijz(71) = 25
    ijz(72) = 25
    ijz(73) = 34
    ijz(74) = 25
    ijz(75) = 28
    ijz(76) = 25
    ijz(77) = 28
    ijz(78) = 31
    ijz(79) = 31
    ijz(80) = 28
    ijz(81) = 25
    ijz(82) = 25
    ijz(83) = 34
    ijz(84) = 25
    ijz(85) = 28
    ijz(86) = 25
    ijz(87) = 28
    ijz(88) = 31
    ijz(89) = 31
    ijz(90) = 28
    ijz(91) = 13
    ijz(92) = 13
    ijz(93) = 22
    ijz(94) = 13
    ijz(95) = 16
    ijz(96) = 13
    ijz(97) = 16
    ijz(98) = 19
    ijz(99) = 19
    ijz(100) = 16

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 2
    klx(2) = 0
    klx(3) = 0
    klx(4) = 1
    klx(5) = 1
    klx(6) = 0

    kly(1) = 0
    kly(2) = 2
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 1

    klz(1) = 0
    klz(2) = 0
    klz(3) = 2
    klz(4) = 0
    klz(5) = 1
    klz(6) = 1

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))

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

    scutsdket = cutoff_schwarz/maxval(sd_pair%xints)
    nsdket = 0
    do ij = 1, res%n_s_shl*res%n_d_shl
      if (sd_pair%xints(ij) .ge. scutsdket) then
        nsdket = nsdket + 1
        xint02ket(nsdket) = sd_pair%xints(ij)
        n02ket(nsdket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = nffbra*nsdket
    ntile = int(nffbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*nsdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint02ket, n02ket, ff_pair, sd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d02ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
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

              test = xint33bra(ij_tmp)*xint02ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n02ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_d_shl(ksh_tmp)
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

                  t_expon_cd = sd_pair%t_expon_ab(sd_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = sd_pair%expon_b(sd_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = sd_pair%expon_a(sd_pair%pair_loc(kl) + ket_loop)
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

                  d02ket(1) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d02ket(2) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d02ket(3) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d02ket(4) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d02ket(5) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d02ket(6) = sd_pair%d_coeff_alt(sd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3

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

                                      ! i2 = in(2) =   13
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(13) = xc00
                                      yin(13) = yc00
                                      zin(13) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   14
                                      ! i2 =   13

                                      xin(14) = xcp00*xin(13) + cp10
                                      yin(14) = ycp00*yin(13) + cp10
                                      zin(14) = zcp00*zin(13) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   25
                                      ! i3 =    1
                                      ! i4 =   13

                                      xin(25) = c10*xin(1) + xc00*xin(13)
                                      yin(25) = c10*yin(1) + yc00*yin(13)
                                      zin(25) = c10*zin(1) + zc00*zin(13)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   26
                                      ! i5 =   25
                                      ! i4 =   13

                                      xin(26) = xcp00*xin(25) + cp10*xin(13)
                                      yin(26) = ycp00*yin(25) + cp10*yin(13)
                                      zin(26) = zcp00*zin(25) + cp10*zin(13)

                                      ! ------------------

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   37
                                      ! i3 =   13
                                      ! i4 =   25

                                      xin(37) = c10*xin(13) + xc00*xin(25)
                                      yin(37) = c10*yin(13) + yc00*yin(25)
                                      zin(37) = c10*zin(13) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   38
                                      ! i5 =   37
                                      ! i4 =   25

                                      xin(38) = xcp00*xin(37) + cp10*xin(25)
                                      yin(38) = ycp00*yin(37) + cp10*yin(25)
                                      zin(38) = zcp00*zin(37) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   40
                                      ! i3 =   25
                                      ! i4 =   37

                                      xin(40) = c10*xin(25) + xc00*xin(37)
                                      yin(40) = c10*yin(25) + yc00*yin(37)
                                      zin(40) = c10*zin(25) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   41
                                      ! i5 =   40
                                      ! i4 =   37

                                      xin(41) = xcp00*xin(40) + cp10*xin(37)
                                      yin(41) = ycp00*yin(40) + cp10*yin(37)
                                      zin(41) = zcp00*zin(40) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   40

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   43
                                      ! i3 =   37
                                      ! i4 =   40

                                      xin(43) = c10*xin(37) + xc00*xin(40)
                                      yin(43) = c10*yin(37) + yc00*yin(40)
                                      zin(43) = c10*zin(37) + zc00*zin(40)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   44
                                      ! i5 =   43
                                      ! i4 =   40

                                      xin(44) = xcp00*xin(43) + cp10*xin(40)
                                      yin(44) = ycp00*yin(43) + cp10*yin(40)
                                      zin(44) = zcp00*zin(43) + cp10*zin(40)

                                      ! ------------------

                                      ! i3 = i4 =   40
                                      ! i4 = i5 =   43

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   46
                                      ! i3 =   40
                                      ! i4 =   43

                                      xin(46) = c10*xin(40) + xc00*xin(43)
                                      yin(46) = c10*yin(40) + yc00*yin(43)
                                      zin(46) = c10*zin(40) + zc00*zin(43)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   47
                                      ! i5 =   46
                                      ! i4 =   43

                                      xin(47) = xcp00*xin(46) + cp10*xin(43)
                                      yin(47) = ycp00*yin(46) + cp10*yin(43)
                                      zin(47) = zcp00*zin(46) + cp10*zin(43)

                                      ! ------------------

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   46

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    2

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    3
                                      ! i3 =    1
                                      ! i4 =    2

                                      xin(3) = cp01*xin(1) + xcp00*xin(2)
                                      yin(3) = cp01*yin(1) + ycp00*yin(2)
                                      zin(3) = cp01*zin(1) + zcp00*zin(2)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   15

                                      xin(15) = xc00*xin(3) + c01*xin(2)
                                      yin(15) = yc00*yin(3) + c01*yin(2)
                                      zin(15) = zc00*zin(3) + c01*zin(2)

                                      ! ------------------

                                      ! i3 = i4 =    2
                                      ! i4 = i5 =    3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   25

                                      xin(27) = c10*xin(3) + xc00*xin(15) + c01*xin(14)
                                      yin(27) = c10*yin(3) + yc00*yin(15) + c01*yin(14)
                                      zin(27) = c10*zin(3) + zc00*zin(15) + c01*zin(14)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(39) = c10*xin(15) + xc00*xin(27) + c01*xin(26)
                                      yin(39) = c10*yin(15) + yc00*yin(27) + c01*yin(26)
                                      zin(39) = c10*zin(15) + zc00*zin(27) + c01*zin(26)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   40

                                      xin(42) = c10*xin(27) + xc00*xin(39) + c01*xin(38)
                                      yin(42) = c10*yin(27) + yc00*yin(39) + c01*yin(38)
                                      zin(42) = c10*zin(27) + zc00*zin(39) + c01*zin(38)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   40

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   43

                                      xin(45) = c10*xin(39) + xc00*xin(42) + c01*xin(41)
                                      yin(45) = c10*yin(39) + yc00*yin(42) + c01*yin(41)
                                      zin(45) = c10*zin(39) + zc00*zin(42) + c01*zin(41)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   40
                                      ! i4 = i5 =   43

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   46

                                      xin(48) = c10*xin(42) + xc00*xin(45) + c01*xin(44)
                                      yin(48) = c10*yin(42) + yc00*yin(45) + c01*yin(44)
                                      zin(48) = c10*zin(42) + zc00*zin(45) + c01*zin(44)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   46

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   46

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(46) = xin(46) + dxij*xin(43)
                                      yin(46) = yin(46) + dyij*yin(43)
                                      zin(46) = zin(46) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   40

                                      xin(43) = xin(43) + dxij*xin(40)
                                      yin(43) = yin(43) + dyij*yin(40)
                                      zin(43) = zin(43) + dzij*zin(40)

                                      ! i3 = i4 =   40
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   37

                                      xin(40) = xin(40) + dxij*xin(37)
                                      yin(40) = yin(40) + dyij*yin(37)
                                      zin(40) = zin(40) + dzij*zin(37)

                                      ! i3 = i4 =   37
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(46) = xin(46) + dxij*xin(43)
                                      yin(46) = yin(46) + dyij*yin(43)
                                      zin(46) = zin(46) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   40

                                      xin(43) = xin(43) + dxij*xin(40)
                                      yin(43) = yin(43) + dyij*yin(40)
                                      zin(43) = zin(43) + dzij*zin(40)

                                      ! i3 = i4 =   40
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(46) = xin(46) + dxij*xin(43)
                                      yin(46) = yin(46) + dyij*yin(43)
                                      zin(46) = zin(46) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    4

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    4

                                      ! do ni = 1,    3

                                      xin(4) = xin(13) + dxij*xin(1)
                                      yin(4) = yin(13) + dyij*yin(1)
                                      zin(4) = zin(13) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    2

                                      xin(16) = xin(25) + dxij*xin(13)
                                      yin(16) = yin(25) + dyij*yin(13)
                                      zin(16) = zin(25) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   28

                                      ! ni =    3

                                      xin(28) = xin(37) + dxij*xin(25)
                                      yin(28) = yin(37) + dyij*yin(25)
                                      zin(28) = zin(37) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    7

                                      ! nj =    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(16) + dxij*xin(4)
                                      yin(7) = yin(16) + dyij*yin(4)
                                      zin(7) = zin(16) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   19

                                      ! ni =    2

                                      xin(19) = xin(28) + dxij*xin(16)
                                      yin(19) = yin(28) + dyij*yin(16)
                                      zin(19) = zin(28) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    3

                                      xin(31) = xin(40) + dxij*xin(28)
                                      yin(31) = yin(40) + dyij*yin(28)
                                      zin(31) = zin(40) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    3

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(19) + dxij*xin(7)
                                      yin(10) = yin(19) + dyij*yin(7)
                                      zin(10) = zin(19) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    2

                                      xin(22) = xin(31) + dxij*xin(19)
                                      yin(22) = yin(31) + dyij*yin(19)
                                      zin(22) = zin(31) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    3

                                      xin(34) = xin(43) + dxij*xin(31)
                                      yin(34) = yin(43) + dyij*yin(31)
                                      zin(34) = zin(43) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(47) = xin(47) + dxij*xin(44)
                                      yin(47) = yin(47) + dyij*yin(44)
                                      zin(47) = zin(47) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   41

                                      xin(44) = xin(44) + dxij*xin(41)
                                      yin(44) = yin(44) + dyij*yin(41)
                                      zin(44) = zin(44) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   38

                                      xin(41) = xin(41) + dxij*xin(38)
                                      yin(41) = yin(41) + dyij*yin(38)
                                      zin(41) = zin(41) + dzij*zin(38)

                                      ! i3 = i4 =   38
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(47) = xin(47) + dxij*xin(44)
                                      yin(47) = yin(47) + dyij*yin(44)
                                      zin(47) = zin(47) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   41

                                      xin(44) = xin(44) + dxij*xin(41)
                                      yin(44) = yin(44) + dyij*yin(41)
                                      zin(44) = zin(44) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(47) = xin(47) + dxij*xin(44)
                                      yin(47) = yin(47) + dyij*yin(44)
                                      zin(47) = zin(47) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(14) + dxij*xin(2)
                                      yin(5) = yin(14) + dyij*yin(2)
                                      zin(5) = zin(14) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    2

                                      xin(17) = xin(26) + dxij*xin(14)
                                      yin(17) = yin(26) + dyij*yin(14)
                                      zin(17) = zin(26) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    3

                                      xin(29) = xin(38) + dxij*xin(26)
                                      yin(29) = yin(38) + dyij*yin(26)
                                      zin(29) = zin(38) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    8

                                      ! nj =    2

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

                                      xin(8) = xin(17) + dxij*xin(5)
                                      yin(8) = yin(17) + dyij*yin(5)
                                      zin(8) = zin(17) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   20

                                      ! ni =    2

                                      xin(20) = xin(29) + dxij*xin(17)
                                      yin(20) = yin(29) + dyij*yin(17)
                                      zin(20) = zin(29) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    3

                                      xin(32) = xin(41) + dxij*xin(29)
                                      yin(32) = yin(41) + dyij*yin(29)
                                      zin(32) = zin(41) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    3

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(20) + dxij*xin(8)
                                      yin(11) = yin(20) + dyij*yin(8)
                                      zin(11) = zin(20) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    2

                                      xin(23) = xin(32) + dxij*xin(20)
                                      yin(23) = yin(32) + dyij*yin(20)
                                      zin(23) = zin(32) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    3

                                      xin(35) = xin(44) + dxij*xin(32)
                                      yin(35) = yin(44) + dyij*yin(32)
                                      zin(35) = zin(44) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   14

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(48) = xin(48) + dxij*xin(45)
                                      yin(48) = yin(48) + dyij*yin(45)
                                      zin(48) = zin(48) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   42

                                      xin(45) = xin(45) + dxij*xin(42)
                                      yin(45) = yin(45) + dyij*yin(42)
                                      zin(45) = zin(45) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   39

                                      xin(42) = xin(42) + dxij*xin(39)
                                      yin(42) = yin(42) + dyij*yin(39)
                                      zin(42) = zin(42) + dzij*zin(39)

                                      ! i3 = i4 =   39
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(48) = xin(48) + dxij*xin(45)
                                      yin(48) = yin(48) + dyij*yin(45)
                                      zin(48) = zin(48) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   42

                                      xin(45) = xin(45) + dxij*xin(42)
                                      yin(45) = yin(45) + dyij*yin(42)
                                      zin(45) = zin(45) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(48) = xin(48) + dxij*xin(45)
                                      yin(48) = yin(48) + dyij*yin(45)
                                      zin(48) = zin(48) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    3

                                      xin(6) = xin(15) + dxij*xin(3)
                                      yin(6) = yin(15) + dyij*yin(3)
                                      zin(6) = zin(15) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    2

                                      xin(18) = xin(27) + dxij*xin(15)
                                      yin(18) = yin(27) + dyij*yin(15)
                                      zin(18) = zin(27) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    3

                                      xin(30) = xin(39) + dxij*xin(27)
                                      yin(30) = yin(39) + dyij*yin(27)
                                      zin(30) = zin(39) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(18) + dxij*xin(6)
                                      yin(9) = yin(18) + dyij*yin(6)
                                      zin(9) = zin(18) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    2

                                      xin(21) = xin(30) + dxij*xin(18)
                                      yin(21) = yin(30) + dyij*yin(18)
                                      zin(21) = zin(30) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    3

                                      xin(33) = xin(42) + dxij*xin(30)
                                      yin(33) = yin(42) + dyij*yin(30)
                                      zin(33) = zin(42) + dzij*zin(30)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    3

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(21) + dxij*xin(9)
                                      yin(12) = yin(21) + dyij*yin(9)
                                      zin(12) = zin(21) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    2

                                      xin(24) = xin(33) + dxij*xin(21)
                                      yin(24) = yin(33) + dyij*yin(21)
                                      zin(24) = zin(33) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    3

                                      xin(36) = xin(45) + dxij*xin(33)
                                      yin(36) = yin(45) + dyij*yin(33)
                                      zin(36) = zin(45) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   48

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

                                      ! i1 = in(1) =   49

                                      xin(49) = 1.0_dp
                                      yin(49) = 1.0_dp
                                      zin(49) = f00

                                      ! i2 = in(2) =   61
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(61) = xc00
                                      yin(61) = yc00
                                      zin(61) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   50

                                      xin(50) = xcp00
                                      yin(50) = ycp00
                                      zin(50) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   62
                                      ! i2 =   61

                                      xin(62) = xcp00*xin(61) + cp10
                                      yin(62) = ycp00*yin(61) + cp10
                                      zin(62) = zcp00*zin(61) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   73
                                      ! i3 =   49
                                      ! i4 =   61

                                      xin(73) = c10*xin(49) + xc00*xin(61)
                                      yin(73) = c10*yin(49) + yc00*yin(61)
                                      zin(73) = c10*zin(49) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   74
                                      ! i5 =   73
                                      ! i4 =   61

                                      xin(74) = xcp00*xin(73) + cp10*xin(61)
                                      yin(74) = ycp00*yin(73) + cp10*yin(61)
                                      zin(74) = zcp00*zin(73) + cp10*zin(61)

                                      ! ------------------

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   85
                                      ! i3 =   61
                                      ! i4 =   73

                                      xin(85) = c10*xin(61) + xc00*xin(73)
                                      yin(85) = c10*yin(61) + yc00*yin(73)
                                      zin(85) = c10*zin(61) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   86
                                      ! i5 =   85
                                      ! i4 =   73

                                      xin(86) = xcp00*xin(85) + cp10*xin(73)
                                      yin(86) = ycp00*yin(85) + cp10*yin(73)
                                      zin(86) = zcp00*zin(85) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   88
                                      ! i3 =   73
                                      ! i4 =   85

                                      xin(88) = c10*xin(73) + xc00*xin(85)
                                      yin(88) = c10*yin(73) + yc00*yin(85)
                                      zin(88) = c10*zin(73) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   89
                                      ! i5 =   88
                                      ! i4 =   85

                                      xin(89) = xcp00*xin(88) + cp10*xin(85)
                                      yin(89) = ycp00*yin(88) + cp10*yin(85)
                                      zin(89) = zcp00*zin(88) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   88

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   85
                                      ! i4 =   88

                                      xin(91) = c10*xin(85) + xc00*xin(88)
                                      yin(91) = c10*yin(85) + yc00*yin(88)
                                      zin(91) = c10*zin(85) + zc00*zin(88)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   92
                                      ! i5 =   91
                                      ! i4 =   88

                                      xin(92) = xcp00*xin(91) + cp10*xin(88)
                                      yin(92) = ycp00*yin(91) + cp10*yin(88)
                                      zin(92) = zcp00*zin(91) + cp10*zin(88)

                                      ! ------------------

                                      ! i3 = i4 =   88
                                      ! i4 = i5 =   91

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   94
                                      ! i3 =   88
                                      ! i4 =   91

                                      xin(94) = c10*xin(88) + xc00*xin(91)
                                      yin(94) = c10*yin(88) + yc00*yin(91)
                                      zin(94) = c10*zin(88) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   95
                                      ! i5 =   94
                                      ! i4 =   91

                                      xin(95) = xcp00*xin(94) + cp10*xin(91)
                                      yin(95) = ycp00*yin(94) + cp10*yin(91)
                                      zin(95) = zcp00*zin(94) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   94

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   49
                                      ! i4 = i1+k2 =   50

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   51
                                      ! i3 =   49
                                      ! i4 =   50

                                      xin(51) = cp01*xin(49) + xcp00*xin(50)
                                      yin(51) = cp01*yin(49) + ycp00*yin(50)
                                      zin(51) = cp01*zin(49) + zcp00*zin(50)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   63

                                      xin(63) = xc00*xin(51) + c01*xin(50)
                                      yin(63) = yc00*yin(51) + c01*yin(50)
                                      zin(63) = zc00*zin(51) + c01*zin(50)

                                      ! ------------------

                                      ! i3 = i4 =   50
                                      ! i4 = i5 =   51

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =   73

                                      xin(75) = c10*xin(51) + xc00*xin(63) + c01*xin(62)
                                      yin(75) = c10*yin(51) + yc00*yin(63) + c01*yin(62)
                                      zin(75) = c10*zin(51) + zc00*zin(63) + c01*zin(62)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(87) = c10*xin(63) + xc00*xin(75) + c01*xin(74)
                                      yin(87) = c10*yin(63) + yc00*yin(75) + c01*yin(74)
                                      zin(87) = c10*zin(63) + zc00*zin(75) + c01*zin(74)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   88

                                      xin(90) = c10*xin(75) + xc00*xin(87) + c01*xin(86)
                                      yin(90) = c10*yin(75) + yc00*yin(87) + c01*yin(86)
                                      zin(90) = c10*zin(75) + zc00*zin(87) + c01*zin(86)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   88

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   91

                                      xin(93) = c10*xin(87) + xc00*xin(90) + c01*xin(89)
                                      yin(93) = c10*yin(87) + yc00*yin(90) + c01*yin(89)
                                      zin(93) = c10*zin(87) + zc00*zin(90) + c01*zin(89)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   88
                                      ! i4 = i5 =   91

                                      ! nn =    6

                                      ! i5 = in(nn+1) =   94

                                      xin(96) = c10*xin(90) + xc00*xin(93) + c01*xin(92)
                                      yin(96) = c10*yin(90) + yc00*yin(93) + c01*yin(92)
                                      zin(96) = c10*zin(90) + zc00*zin(93) + c01*zin(92)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   94

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   94

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(94) = xin(94) + dxij*xin(91)
                                      yin(94) = yin(94) + dyij*yin(91)
                                      zin(94) = zin(94) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   88

                                      xin(91) = xin(91) + dxij*xin(88)
                                      yin(91) = yin(91) + dyij*yin(88)
                                      zin(91) = zin(91) + dzij*zin(88)

                                      ! i3 = i4 =   88
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   85

                                      xin(88) = xin(88) + dxij*xin(85)
                                      yin(88) = yin(88) + dyij*yin(85)
                                      zin(88) = zin(88) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(94) = xin(94) + dxij*xin(91)
                                      yin(94) = yin(94) + dyij*yin(91)
                                      zin(94) = zin(94) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   88

                                      xin(91) = xin(91) + dxij*xin(88)
                                      yin(91) = yin(91) + dyij*yin(88)
                                      zin(91) = zin(91) + dzij*zin(88)

                                      ! i3 = i4 =   88
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(94) = xin(94) + dxij*xin(91)
                                      yin(94) = yin(94) + dyij*yin(91)
                                      zin(94) = zin(94) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   52

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   52

                                      ! do ni = 1,    3

                                      xin(52) = xin(61) + dxij*xin(49)
                                      yin(52) = yin(61) + dyij*yin(49)
                                      zin(52) = zin(61) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    2

                                      xin(64) = xin(73) + dxij*xin(61)
                                      yin(64) = yin(73) + dyij*yin(61)
                                      zin(64) = zin(73) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   76

                                      ! ni =    3

                                      xin(76) = xin(85) + dxij*xin(73)
                                      yin(76) = yin(85) + dyij*yin(73)
                                      zin(76) = zin(85) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   55

                                      ! nj =    2

                                      ! i4 = i3 =   55

                                      ! do ni = 1,    3

                                      xin(55) = xin(64) + dxij*xin(52)
                                      yin(55) = yin(64) + dyij*yin(52)
                                      zin(55) = zin(64) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    2

                                      xin(67) = xin(76) + dxij*xin(64)
                                      yin(67) = yin(76) + dyij*yin(64)
                                      zin(67) = zin(76) + dzij*zin(64)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(88) + dxij*xin(76)
                                      yin(79) = yin(88) + dyij*yin(76)
                                      zin(79) = zin(88) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   58

                                      ! nj =    3

                                      ! i4 = i3 =   58

                                      ! do ni = 1,    3

                                      xin(58) = xin(67) + dxij*xin(55)
                                      yin(58) = yin(67) + dyij*yin(55)
                                      zin(58) = zin(67) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    2

                                      xin(70) = xin(79) + dxij*xin(67)
                                      yin(70) = yin(79) + dyij*yin(67)
                                      zin(70) = zin(79) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   82

                                      ! ni =    3

                                      xin(82) = xin(91) + dxij*xin(79)
                                      yin(82) = yin(91) + dyij*yin(79)
                                      zin(82) = zin(91) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   61

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(95) = xin(95) + dxij*xin(92)
                                      yin(95) = yin(95) + dyij*yin(92)
                                      zin(95) = zin(95) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   89

                                      xin(92) = xin(92) + dxij*xin(89)
                                      yin(92) = yin(92) + dyij*yin(89)
                                      zin(92) = zin(92) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   86

                                      xin(89) = xin(89) + dxij*xin(86)
                                      yin(89) = yin(89) + dyij*yin(86)
                                      zin(89) = zin(89) + dzij*zin(86)

                                      ! i3 = i4 =   86
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(95) = xin(95) + dxij*xin(92)
                                      yin(95) = yin(95) + dyij*yin(92)
                                      zin(95) = zin(95) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   89

                                      xin(92) = xin(92) + dxij*xin(89)
                                      yin(92) = yin(92) + dyij*yin(89)
                                      zin(92) = zin(92) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(95) = xin(95) + dxij*xin(92)
                                      yin(95) = yin(95) + dyij*yin(92)
                                      zin(95) = zin(95) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   53

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   53

                                      ! do ni = 1,    3

                                      xin(53) = xin(62) + dxij*xin(50)
                                      yin(53) = yin(62) + dyij*yin(50)
                                      zin(53) = zin(62) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    2

                                      xin(65) = xin(74) + dxij*xin(62)
                                      yin(65) = yin(74) + dyij*yin(62)
                                      zin(65) = zin(74) + dzij*zin(62)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   77

                                      ! ni =    3

                                      xin(77) = xin(86) + dxij*xin(74)
                                      yin(77) = yin(86) + dyij*yin(74)
                                      zin(77) = zin(86) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   56

                                      ! nj =    2

                                      ! i4 = i3 =   56

                                      ! do ni = 1,    3

                                      xin(56) = xin(65) + dxij*xin(53)
                                      yin(56) = yin(65) + dyij*yin(53)
                                      zin(56) = zin(65) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   68

                                      ! ni =    2

                                      xin(68) = xin(77) + dxij*xin(65)
                                      yin(68) = yin(77) + dyij*yin(65)
                                      zin(68) = zin(77) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    3

                                      xin(80) = xin(89) + dxij*xin(77)
                                      yin(80) = yin(89) + dyij*yin(77)
                                      zin(80) = zin(89) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   92

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   59

                                      ! nj =    3

                                      ! i4 = i3 =   59

                                      ! do ni = 1,    3

                                      xin(59) = xin(68) + dxij*xin(56)
                                      yin(59) = yin(68) + dyij*yin(56)
                                      zin(59) = zin(68) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    2

                                      xin(71) = xin(80) + dxij*xin(68)
                                      yin(71) = yin(80) + dyij*yin(68)
                                      zin(71) = zin(80) + dzij*zin(68)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    3

                                      xin(83) = xin(92) + dxij*xin(80)
                                      yin(83) = yin(92) + dyij*yin(80)
                                      zin(83) = zin(92) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   62

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   93

                                      xin(96) = xin(96) + dxij*xin(93)
                                      yin(96) = yin(96) + dyij*yin(93)
                                      zin(96) = zin(96) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   90

                                      xin(93) = xin(93) + dxij*xin(90)
                                      yin(93) = yin(93) + dyij*yin(90)
                                      zin(93) = zin(93) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   87

                                      xin(90) = xin(90) + dxij*xin(87)
                                      yin(90) = yin(90) + dyij*yin(87)
                                      zin(90) = zin(90) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   93

                                      xin(96) = xin(96) + dxij*xin(93)
                                      yin(96) = yin(96) + dyij*yin(93)
                                      zin(96) = zin(96) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   90

                                      xin(93) = xin(93) + dxij*xin(90)
                                      yin(93) = yin(93) + dyij*yin(90)
                                      zin(93) = zin(93) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   93

                                      xin(96) = xin(96) + dxij*xin(93)
                                      yin(96) = yin(96) + dyij*yin(93)
                                      zin(96) = zin(96) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   54

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   54

                                      ! do ni = 1,    3

                                      xin(54) = xin(63) + dxij*xin(51)
                                      yin(54) = yin(63) + dyij*yin(51)
                                      zin(54) = zin(63) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    2

                                      xin(66) = xin(75) + dxij*xin(63)
                                      yin(66) = yin(75) + dyij*yin(63)
                                      zin(66) = zin(75) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   78

                                      ! ni =    3

                                      xin(78) = xin(87) + dxij*xin(75)
                                      yin(78) = yin(87) + dyij*yin(75)
                                      zin(78) = zin(87) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   57

                                      ! nj =    2

                                      ! i4 = i3 =   57

                                      ! do ni = 1,    3

                                      xin(57) = xin(66) + dxij*xin(54)
                                      yin(57) = yin(66) + dyij*yin(54)
                                      zin(57) = zin(66) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    2

                                      xin(69) = xin(78) + dxij*xin(66)
                                      yin(69) = yin(78) + dyij*yin(66)
                                      zin(69) = zin(78) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(90) + dxij*xin(78)
                                      yin(81) = yin(90) + dyij*yin(78)
                                      zin(81) = zin(90) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   60

                                      ! nj =    3

                                      ! i4 = i3 =   60

                                      ! do ni = 1,    3

                                      xin(60) = xin(69) + dxij*xin(57)
                                      yin(60) = yin(69) + dyij*yin(57)
                                      zin(60) = zin(69) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    2

                                      xin(72) = xin(81) + dxij*xin(69)
                                      yin(72) = yin(81) + dyij*yin(69)
                                      zin(72) = zin(81) + dzij*zin(69)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    3

                                      xin(84) = xin(93) + dxij*xin(81)
                                      yin(84) = yin(93) + dyij*yin(81)
                                      zin(84) = zin(93) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   63

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

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

                                      ! i1 = in(1) =   97

                                      xin(97) = 1.0_dp
                                      yin(97) = 1.0_dp
                                      zin(97) = f00

                                      ! i2 = in(2) =  109
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(109) = xc00
                                      yin(109) = yc00
                                      zin(109) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   98

                                      xin(98) = xcp00
                                      yin(98) = ycp00
                                      zin(98) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  110
                                      ! i2 =  109

                                      xin(110) = xcp00*xin(109) + cp10
                                      yin(110) = ycp00*yin(109) + cp10
                                      zin(110) = zcp00*zin(109) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =   97
                                      ! i4 =  109

                                      xin(121) = c10*xin(97) + xc00*xin(109)
                                      yin(121) = c10*yin(97) + yc00*yin(109)
                                      zin(121) = c10*zin(97) + zc00*zin(109)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  122
                                      ! i5 =  121
                                      ! i4 =  109

                                      xin(122) = xcp00*xin(121) + cp10*xin(109)
                                      yin(122) = ycp00*yin(121) + cp10*yin(109)
                                      zin(122) = zcp00*zin(121) + cp10*zin(109)

                                      ! ------------------

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  133
                                      ! i3 =  109
                                      ! i4 =  121

                                      xin(133) = c10*xin(109) + xc00*xin(121)
                                      yin(133) = c10*yin(109) + yc00*yin(121)
                                      zin(133) = c10*zin(109) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  134
                                      ! i5 =  133
                                      ! i4 =  121

                                      xin(134) = xcp00*xin(133) + cp10*xin(121)
                                      yin(134) = ycp00*yin(133) + cp10*yin(121)
                                      zin(134) = zcp00*zin(133) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  136
                                      ! i3 =  121
                                      ! i4 =  133

                                      xin(136) = c10*xin(121) + xc00*xin(133)
                                      yin(136) = c10*yin(121) + yc00*yin(133)
                                      zin(136) = c10*zin(121) + zc00*zin(133)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  137
                                      ! i5 =  136
                                      ! i4 =  133

                                      xin(137) = xcp00*xin(136) + cp10*xin(133)
                                      yin(137) = ycp00*yin(136) + cp10*yin(133)
                                      zin(137) = zcp00*zin(136) + cp10*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  136

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  139
                                      ! i3 =  133
                                      ! i4 =  136

                                      xin(139) = c10*xin(133) + xc00*xin(136)
                                      yin(139) = c10*yin(133) + yc00*yin(136)
                                      zin(139) = c10*zin(133) + zc00*zin(136)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  140
                                      ! i5 =  139
                                      ! i4 =  136

                                      xin(140) = xcp00*xin(139) + cp10*xin(136)
                                      yin(140) = ycp00*yin(139) + cp10*yin(136)
                                      zin(140) = zcp00*zin(139) + cp10*zin(136)

                                      ! ------------------

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  139

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  142
                                      ! i3 =  136
                                      ! i4 =  139

                                      xin(142) = c10*xin(136) + xc00*xin(139)
                                      yin(142) = c10*yin(136) + yc00*yin(139)
                                      zin(142) = c10*zin(136) + zc00*zin(139)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  143
                                      ! i5 =  142
                                      ! i4 =  139

                                      xin(143) = xcp00*xin(142) + cp10*xin(139)
                                      yin(143) = ycp00*yin(142) + cp10*yin(139)
                                      zin(143) = zcp00*zin(142) + cp10*zin(139)

                                      ! ------------------

                                      ! i3 = i4 =  139
                                      ! i4 = i5 =  142

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   98

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   99
                                      ! i3 =   97
                                      ! i4 =   98

                                      xin(99) = cp01*xin(97) + xcp00*xin(98)
                                      yin(99) = cp01*yin(97) + ycp00*yin(98)
                                      zin(99) = cp01*zin(97) + zcp00*zin(98)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  111

                                      xin(111) = xc00*xin(99) + c01*xin(98)
                                      yin(111) = yc00*yin(99) + c01*yin(98)
                                      zin(111) = zc00*zin(99) + c01*zin(98)

                                      ! ------------------

                                      ! i3 = i4 =   98
                                      ! i4 = i5 =   99

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  121

                                      xin(123) = c10*xin(99) + xc00*xin(111) + c01*xin(110)
                                      yin(123) = c10*yin(99) + yc00*yin(111) + c01*yin(110)
                                      zin(123) = c10*zin(99) + zc00*zin(111) + c01*zin(110)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(135) = c10*xin(111) + xc00*xin(123) + c01*xin(122)
                                      yin(135) = c10*yin(111) + yc00*yin(123) + c01*yin(122)
                                      zin(135) = c10*zin(111) + zc00*zin(123) + c01*zin(122)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  136

                                      xin(138) = c10*xin(123) + xc00*xin(135) + c01*xin(134)
                                      yin(138) = c10*yin(123) + yc00*yin(135) + c01*yin(134)
                                      zin(138) = c10*zin(123) + zc00*zin(135) + c01*zin(134)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  136

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  139

                                      xin(141) = c10*xin(135) + xc00*xin(138) + c01*xin(137)
                                      yin(141) = c10*yin(135) + yc00*yin(138) + c01*yin(137)
                                      zin(141) = c10*zin(135) + zc00*zin(138) + c01*zin(137)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  139

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  142

                                      xin(144) = c10*xin(138) + xc00*xin(141) + c01*xin(140)
                                      yin(144) = c10*yin(138) + yc00*yin(141) + c01*yin(140)
                                      zin(144) = c10*zin(138) + zc00*zin(141) + c01*zin(140)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  139
                                      ! i4 = i5 =  142

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  142

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(142) = xin(142) + dxij*xin(139)
                                      yin(142) = yin(142) + dyij*yin(139)
                                      zin(142) = zin(142) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  136

                                      xin(139) = xin(139) + dxij*xin(136)
                                      yin(139) = yin(139) + dyij*yin(136)
                                      zin(139) = zin(139) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  133

                                      xin(136) = xin(136) + dxij*xin(133)
                                      yin(136) = yin(136) + dyij*yin(133)
                                      zin(136) = zin(136) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(142) = xin(142) + dxij*xin(139)
                                      yin(142) = yin(142) + dyij*yin(139)
                                      zin(142) = zin(142) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  136

                                      xin(139) = xin(139) + dxij*xin(136)
                                      yin(139) = yin(139) + dyij*yin(136)
                                      zin(139) = zin(139) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(142) = xin(142) + dxij*xin(139)
                                      yin(142) = yin(142) + dyij*yin(139)
                                      zin(142) = zin(142) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  100

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  100

                                      ! do ni = 1,    3

                                      xin(100) = xin(109) + dxij*xin(97)
                                      yin(100) = yin(109) + dyij*yin(97)
                                      zin(100) = zin(109) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  112

                                      ! ni =    2

                                      xin(112) = xin(121) + dxij*xin(109)
                                      yin(112) = yin(121) + dyij*yin(109)
                                      zin(112) = zin(121) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  124

                                      ! ni =    3

                                      xin(124) = xin(133) + dxij*xin(121)
                                      yin(124) = yin(133) + dyij*yin(121)
                                      zin(124) = zin(133) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  136

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  103

                                      ! nj =    2

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(112) + dxij*xin(100)
                                      yin(103) = yin(112) + dyij*yin(100)
                                      zin(103) = zin(112) + dzij*zin(100)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    2

                                      xin(115) = xin(124) + dxij*xin(112)
                                      yin(115) = yin(124) + dyij*yin(112)
                                      zin(115) = zin(124) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    3

                                      xin(127) = xin(136) + dxij*xin(124)
                                      yin(127) = yin(136) + dyij*yin(124)
                                      zin(127) = zin(136) + dzij*zin(124)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  106

                                      ! nj =    3

                                      ! i4 = i3 =  106

                                      ! do ni = 1,    3

                                      xin(106) = xin(115) + dxij*xin(103)
                                      yin(106) = yin(115) + dyij*yin(103)
                                      zin(106) = zin(115) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  118

                                      ! ni =    2

                                      xin(118) = xin(127) + dxij*xin(115)
                                      yin(118) = yin(127) + dyij*yin(115)
                                      zin(118) = zin(127) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                      ! ni =    3

                                      xin(130) = xin(139) + dxij*xin(127)
                                      yin(130) = yin(139) + dyij*yin(127)
                                      zin(130) = zin(139) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  109

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(143) = xin(143) + dxij*xin(140)
                                      yin(143) = yin(143) + dyij*yin(140)
                                      zin(143) = zin(143) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  137

                                      xin(140) = xin(140) + dxij*xin(137)
                                      yin(140) = yin(140) + dyij*yin(137)
                                      zin(140) = zin(140) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  134

                                      xin(137) = xin(137) + dxij*xin(134)
                                      yin(137) = yin(137) + dyij*yin(134)
                                      zin(137) = zin(137) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(143) = xin(143) + dxij*xin(140)
                                      yin(143) = yin(143) + dyij*yin(140)
                                      zin(143) = zin(143) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  137

                                      xin(140) = xin(140) + dxij*xin(137)
                                      yin(140) = yin(140) + dyij*yin(137)
                                      zin(140) = zin(140) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(143) = xin(143) + dxij*xin(140)
                                      yin(143) = yin(143) + dyij*yin(140)
                                      zin(143) = zin(143) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  101

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  101

                                      ! do ni = 1,    3

                                      xin(101) = xin(110) + dxij*xin(98)
                                      yin(101) = yin(110) + dyij*yin(98)
                                      zin(101) = zin(110) + dzij*zin(98)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! ni =    2

                                      xin(113) = xin(122) + dxij*xin(110)
                                      yin(113) = yin(122) + dyij*yin(110)
                                      zin(113) = zin(122) + dzij*zin(110)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    3

                                      xin(125) = xin(134) + dxij*xin(122)
                                      yin(125) = yin(134) + dyij*yin(122)
                                      zin(125) = zin(134) + dzij*zin(122)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  104

                                      ! nj =    2

                                      ! i4 = i3 =  104

                                      ! do ni = 1,    3

                                      xin(104) = xin(113) + dxij*xin(101)
                                      yin(104) = yin(113) + dyij*yin(101)
                                      zin(104) = zin(113) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  116

                                      ! ni =    2

                                      xin(116) = xin(125) + dxij*xin(113)
                                      yin(116) = yin(125) + dyij*yin(113)
                                      zin(116) = zin(125) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    3

                                      xin(128) = xin(137) + dxij*xin(125)
                                      yin(128) = yin(137) + dyij*yin(125)
                                      zin(128) = zin(137) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  140

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  107

                                      ! nj =    3

                                      ! i4 = i3 =  107

                                      ! do ni = 1,    3

                                      xin(107) = xin(116) + dxij*xin(104)
                                      yin(107) = yin(116) + dyij*yin(104)
                                      zin(107) = zin(116) + dzij*zin(104)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    2

                                      xin(119) = xin(128) + dxij*xin(116)
                                      yin(119) = yin(128) + dyij*yin(116)
                                      zin(119) = zin(128) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      xin(131) = xin(140) + dxij*xin(128)
                                      yin(131) = yin(140) + dyij*yin(128)
                                      zin(131) = zin(140) + dzij*zin(128)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  110

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  141

                                      xin(144) = xin(144) + dxij*xin(141)
                                      yin(144) = yin(144) + dyij*yin(141)
                                      zin(144) = zin(144) + dzij*zin(141)

                                      ! i3 = i4 =  141
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  138

                                      xin(141) = xin(141) + dxij*xin(138)
                                      yin(141) = yin(141) + dyij*yin(138)
                                      zin(141) = zin(141) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  135

                                      xin(138) = xin(138) + dxij*xin(135)
                                      yin(138) = yin(138) + dyij*yin(135)
                                      zin(138) = zin(138) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  141

                                      xin(144) = xin(144) + dxij*xin(141)
                                      yin(144) = yin(144) + dyij*yin(141)
                                      zin(144) = zin(144) + dzij*zin(141)

                                      ! i3 = i4 =  141
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  138

                                      xin(141) = xin(141) + dxij*xin(138)
                                      yin(141) = yin(141) + dyij*yin(138)
                                      zin(141) = zin(141) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  141

                                      xin(144) = xin(144) + dxij*xin(141)
                                      yin(144) = yin(144) + dyij*yin(141)
                                      zin(144) = zin(144) + dzij*zin(141)

                                      ! i3 = i4 =  141
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  102

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  102

                                      ! do ni = 1,    3

                                      xin(102) = xin(111) + dxij*xin(99)
                                      yin(102) = yin(111) + dyij*yin(99)
                                      zin(102) = zin(111) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  114

                                      ! ni =    2

                                      xin(114) = xin(123) + dxij*xin(111)
                                      yin(114) = yin(123) + dyij*yin(111)
                                      zin(114) = zin(123) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    3

                                      xin(126) = xin(135) + dxij*xin(123)
                                      yin(126) = yin(135) + dyij*yin(123)
                                      zin(126) = zin(135) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  105

                                      ! nj =    2

                                      ! i4 = i3 =  105

                                      ! do ni = 1,    3

                                      xin(105) = xin(114) + dxij*xin(102)
                                      yin(105) = yin(114) + dyij*yin(102)
                                      zin(105) = zin(114) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    2

                                      xin(117) = xin(126) + dxij*xin(114)
                                      yin(117) = yin(126) + dyij*yin(114)
                                      zin(117) = zin(126) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    3

                                      xin(129) = xin(138) + dxij*xin(126)
                                      yin(129) = yin(138) + dyij*yin(126)
                                      zin(129) = zin(138) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  108

                                      ! nj =    3

                                      ! i4 = i3 =  108

                                      ! do ni = 1,    3

                                      xin(108) = xin(117) + dxij*xin(105)
                                      yin(108) = yin(117) + dyij*yin(105)
                                      zin(108) = zin(117) + dzij*zin(105)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    2

                                      xin(120) = xin(129) + dxij*xin(117)
                                      yin(120) = yin(129) + dyij*yin(117)
                                      zin(120) = zin(129) + dzij*zin(117)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    3

                                      xin(132) = xin(141) + dxij*xin(129)
                                      yin(132) = yin(141) + dyij*yin(129)
                                      zin(132) = zin(141) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  111

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

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

                                      ! i1 = in(1) =  145

                                      xin(145) = 1.0_dp
                                      yin(145) = 1.0_dp
                                      zin(145) = f00

                                      ! i2 = in(2) =  157
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(157) = xc00
                                      yin(157) = yc00
                                      zin(157) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  146

                                      xin(146) = xcp00
                                      yin(146) = ycp00
                                      zin(146) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  158
                                      ! i2 =  157

                                      xin(158) = xcp00*xin(157) + cp10
                                      yin(158) = ycp00*yin(157) + cp10
                                      zin(158) = zcp00*zin(157) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  169
                                      ! i3 =  145
                                      ! i4 =  157

                                      xin(169) = c10*xin(145) + xc00*xin(157)
                                      yin(169) = c10*yin(145) + yc00*yin(157)
                                      zin(169) = c10*zin(145) + zc00*zin(157)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  170
                                      ! i5 =  169
                                      ! i4 =  157

                                      xin(170) = xcp00*xin(169) + cp10*xin(157)
                                      yin(170) = ycp00*yin(169) + cp10*yin(157)
                                      zin(170) = zcp00*zin(169) + cp10*zin(157)

                                      ! ------------------

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  181
                                      ! i3 =  157
                                      ! i4 =  169

                                      xin(181) = c10*xin(157) + xc00*xin(169)
                                      yin(181) = c10*yin(157) + yc00*yin(169)
                                      zin(181) = c10*zin(157) + zc00*zin(169)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  182
                                      ! i5 =  181
                                      ! i4 =  169

                                      xin(182) = xcp00*xin(181) + cp10*xin(169)
                                      yin(182) = ycp00*yin(181) + cp10*yin(169)
                                      zin(182) = zcp00*zin(181) + cp10*zin(169)

                                      ! ------------------

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  184
                                      ! i3 =  169
                                      ! i4 =  181

                                      xin(184) = c10*xin(169) + xc00*xin(181)
                                      yin(184) = c10*yin(169) + yc00*yin(181)
                                      zin(184) = c10*zin(169) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  185
                                      ! i5 =  184
                                      ! i4 =  181

                                      xin(185) = xcp00*xin(184) + cp10*xin(181)
                                      yin(185) = ycp00*yin(184) + cp10*yin(181)
                                      zin(185) = zcp00*zin(184) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  184

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  187
                                      ! i3 =  181
                                      ! i4 =  184

                                      xin(187) = c10*xin(181) + xc00*xin(184)
                                      yin(187) = c10*yin(181) + yc00*yin(184)
                                      zin(187) = c10*zin(181) + zc00*zin(184)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  188
                                      ! i5 =  187
                                      ! i4 =  184

                                      xin(188) = xcp00*xin(187) + cp10*xin(184)
                                      yin(188) = ycp00*yin(187) + cp10*yin(184)
                                      zin(188) = zcp00*zin(187) + cp10*zin(184)

                                      ! ------------------

                                      ! i3 = i4 =  184
                                      ! i4 = i5 =  187

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  190
                                      ! i3 =  184
                                      ! i4 =  187

                                      xin(190) = c10*xin(184) + xc00*xin(187)
                                      yin(190) = c10*yin(184) + yc00*yin(187)
                                      zin(190) = c10*zin(184) + zc00*zin(187)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  191
                                      ! i5 =  190
                                      ! i4 =  187

                                      xin(191) = xcp00*xin(190) + cp10*xin(187)
                                      yin(191) = ycp00*yin(190) + cp10*yin(187)
                                      zin(191) = zcp00*zin(190) + cp10*zin(187)

                                      ! ------------------

                                      ! i3 = i4 =  187
                                      ! i4 = i5 =  190

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  146

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  147
                                      ! i3 =  145
                                      ! i4 =  146

                                      xin(147) = cp01*xin(145) + xcp00*xin(146)
                                      yin(147) = cp01*yin(145) + ycp00*yin(146)
                                      zin(147) = cp01*zin(145) + zcp00*zin(146)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  159

                                      xin(159) = xc00*xin(147) + c01*xin(146)
                                      yin(159) = yc00*yin(147) + c01*yin(146)
                                      zin(159) = zc00*zin(147) + c01*zin(146)

                                      ! ------------------

                                      ! i3 = i4 =  146
                                      ! i4 = i5 =  147

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  169

                                      xin(171) = c10*xin(147) + xc00*xin(159) + c01*xin(158)
                                      yin(171) = c10*yin(147) + yc00*yin(159) + c01*yin(158)
                                      zin(171) = c10*zin(147) + zc00*zin(159) + c01*zin(158)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(183) = c10*xin(159) + xc00*xin(171) + c01*xin(170)
                                      yin(183) = c10*yin(159) + yc00*yin(171) + c01*yin(170)
                                      zin(183) = c10*zin(159) + zc00*zin(171) + c01*zin(170)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  184

                                      xin(186) = c10*xin(171) + xc00*xin(183) + c01*xin(182)
                                      yin(186) = c10*yin(171) + yc00*yin(183) + c01*yin(182)
                                      zin(186) = c10*zin(171) + zc00*zin(183) + c01*zin(182)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  184

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  187

                                      xin(189) = c10*xin(183) + xc00*xin(186) + c01*xin(185)
                                      yin(189) = c10*yin(183) + yc00*yin(186) + c01*yin(185)
                                      zin(189) = c10*zin(183) + zc00*zin(186) + c01*zin(185)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  184
                                      ! i4 = i5 =  187

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  190

                                      xin(192) = c10*xin(186) + xc00*xin(189) + c01*xin(188)
                                      yin(192) = c10*yin(186) + yc00*yin(189) + c01*yin(188)
                                      zin(192) = c10*zin(186) + zc00*zin(189) + c01*zin(188)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  187
                                      ! i4 = i5 =  190

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  190

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(190) = xin(190) + dxij*xin(187)
                                      yin(190) = yin(190) + dyij*yin(187)
                                      zin(190) = zin(190) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  184

                                      xin(187) = xin(187) + dxij*xin(184)
                                      yin(187) = yin(187) + dyij*yin(184)
                                      zin(187) = zin(187) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  181

                                      xin(184) = xin(184) + dxij*xin(181)
                                      yin(184) = yin(184) + dyij*yin(181)
                                      zin(184) = zin(184) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(190) = xin(190) + dxij*xin(187)
                                      yin(190) = yin(190) + dyij*yin(187)
                                      zin(190) = zin(190) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  184

                                      xin(187) = xin(187) + dxij*xin(184)
                                      yin(187) = yin(187) + dyij*yin(184)
                                      zin(187) = zin(187) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(190) = xin(190) + dxij*xin(187)
                                      yin(190) = yin(190) + dyij*yin(187)
                                      zin(190) = zin(190) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  148

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  148

                                      ! do ni = 1,    3

                                      xin(148) = xin(157) + dxij*xin(145)
                                      yin(148) = yin(157) + dyij*yin(145)
                                      zin(148) = zin(157) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  160

                                      ! ni =    2

                                      xin(160) = xin(169) + dxij*xin(157)
                                      yin(160) = yin(169) + dyij*yin(157)
                                      zin(160) = zin(169) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  172

                                      ! ni =    3

                                      xin(172) = xin(181) + dxij*xin(169)
                                      yin(172) = yin(181) + dyij*yin(169)
                                      zin(172) = zin(181) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  184

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  151

                                      ! nj =    2

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(160) + dxij*xin(148)
                                      yin(151) = yin(160) + dyij*yin(148)
                                      zin(151) = zin(160) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! ni =    2

                                      xin(163) = xin(172) + dxij*xin(160)
                                      yin(163) = yin(172) + dyij*yin(160)
                                      zin(163) = zin(172) + dzij*zin(160)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(184) + dxij*xin(172)
                                      yin(175) = yin(184) + dyij*yin(172)
                                      zin(175) = zin(184) + dzij*zin(172)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  154

                                      ! nj =    3

                                      ! i4 = i3 =  154

                                      ! do ni = 1,    3

                                      xin(154) = xin(163) + dxij*xin(151)
                                      yin(154) = yin(163) + dyij*yin(151)
                                      zin(154) = zin(163) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  166

                                      ! ni =    2

                                      xin(166) = xin(175) + dxij*xin(163)
                                      yin(166) = yin(175) + dyij*yin(163)
                                      zin(166) = zin(175) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  178

                                      ! ni =    3

                                      xin(178) = xin(187) + dxij*xin(175)
                                      yin(178) = yin(187) + dyij*yin(175)
                                      zin(178) = zin(187) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  190

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(191) = xin(191) + dxij*xin(188)
                                      yin(191) = yin(191) + dyij*yin(188)
                                      zin(191) = zin(191) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  185

                                      xin(188) = xin(188) + dxij*xin(185)
                                      yin(188) = yin(188) + dyij*yin(185)
                                      zin(188) = zin(188) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  182

                                      xin(185) = xin(185) + dxij*xin(182)
                                      yin(185) = yin(185) + dyij*yin(182)
                                      zin(185) = zin(185) + dzij*zin(182)

                                      ! i3 = i4 =  182
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(191) = xin(191) + dxij*xin(188)
                                      yin(191) = yin(191) + dyij*yin(188)
                                      zin(191) = zin(191) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  185

                                      xin(188) = xin(188) + dxij*xin(185)
                                      yin(188) = yin(188) + dyij*yin(185)
                                      zin(188) = zin(188) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(191) = xin(191) + dxij*xin(188)
                                      yin(191) = yin(191) + dyij*yin(188)
                                      zin(191) = zin(191) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  149

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  149

                                      ! do ni = 1,    3

                                      xin(149) = xin(158) + dxij*xin(146)
                                      yin(149) = yin(158) + dyij*yin(146)
                                      zin(149) = zin(158) + dzij*zin(146)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    2

                                      xin(161) = xin(170) + dxij*xin(158)
                                      yin(161) = yin(170) + dyij*yin(158)
                                      zin(161) = zin(170) + dzij*zin(158)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    3

                                      xin(173) = xin(182) + dxij*xin(170)
                                      yin(173) = yin(182) + dyij*yin(170)
                                      zin(173) = zin(182) + dzij*zin(170)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  152

                                      ! nj =    2

                                      ! i4 = i3 =  152

                                      ! do ni = 1,    3

                                      xin(152) = xin(161) + dxij*xin(149)
                                      yin(152) = yin(161) + dyij*yin(149)
                                      zin(152) = zin(161) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  164

                                      ! ni =    2

                                      xin(164) = xin(173) + dxij*xin(161)
                                      yin(164) = yin(173) + dyij*yin(161)
                                      zin(164) = zin(173) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  176

                                      ! ni =    3

                                      xin(176) = xin(185) + dxij*xin(173)
                                      yin(176) = yin(185) + dyij*yin(173)
                                      zin(176) = zin(185) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  188

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  155

                                      ! nj =    3

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(164) + dxij*xin(152)
                                      yin(155) = yin(164) + dyij*yin(152)
                                      zin(155) = zin(164) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    2

                                      xin(167) = xin(176) + dxij*xin(164)
                                      yin(167) = yin(176) + dyij*yin(164)
                                      zin(167) = zin(176) + dzij*zin(164)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    3

                                      xin(179) = xin(188) + dxij*xin(176)
                                      yin(179) = yin(188) + dyij*yin(176)
                                      zin(179) = zin(188) + dzij*zin(176)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  158

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  189

                                      xin(192) = xin(192) + dxij*xin(189)
                                      yin(192) = yin(192) + dyij*yin(189)
                                      zin(192) = zin(192) + dzij*zin(189)

                                      ! i3 = i4 =  189
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  186

                                      xin(189) = xin(189) + dxij*xin(186)
                                      yin(189) = yin(189) + dyij*yin(186)
                                      zin(189) = zin(189) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  183

                                      xin(186) = xin(186) + dxij*xin(183)
                                      yin(186) = yin(186) + dyij*yin(183)
                                      zin(186) = zin(186) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  189

                                      xin(192) = xin(192) + dxij*xin(189)
                                      yin(192) = yin(192) + dyij*yin(189)
                                      zin(192) = zin(192) + dzij*zin(189)

                                      ! i3 = i4 =  189
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  186

                                      xin(189) = xin(189) + dxij*xin(186)
                                      yin(189) = yin(189) + dyij*yin(186)
                                      zin(189) = zin(189) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  189

                                      xin(192) = xin(192) + dxij*xin(189)
                                      yin(192) = yin(192) + dyij*yin(189)
                                      zin(192) = zin(192) + dzij*zin(189)

                                      ! i3 = i4 =  189
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  150

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  150

                                      ! do ni = 1,    3

                                      xin(150) = xin(159) + dxij*xin(147)
                                      yin(150) = yin(159) + dyij*yin(147)
                                      zin(150) = zin(159) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  162

                                      ! ni =    2

                                      xin(162) = xin(171) + dxij*xin(159)
                                      yin(162) = yin(171) + dyij*yin(159)
                                      zin(162) = zin(171) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  174

                                      ! ni =    3

                                      xin(174) = xin(183) + dxij*xin(171)
                                      yin(174) = yin(183) + dyij*yin(171)
                                      zin(174) = zin(183) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  186

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  153

                                      ! nj =    2

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(162) + dxij*xin(150)
                                      yin(153) = yin(162) + dyij*yin(150)
                                      zin(153) = zin(162) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    2

                                      xin(165) = xin(174) + dxij*xin(162)
                                      yin(165) = yin(174) + dyij*yin(162)
                                      zin(165) = zin(174) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    3

                                      xin(177) = xin(186) + dxij*xin(174)
                                      yin(177) = yin(186) + dyij*yin(174)
                                      zin(177) = zin(186) + dzij*zin(174)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  156

                                      ! nj =    3

                                      ! i4 = i3 =  156

                                      ! do ni = 1,    3

                                      xin(156) = xin(165) + dxij*xin(153)
                                      yin(156) = yin(165) + dyij*yin(153)
                                      zin(156) = zin(165) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    2

                                      xin(168) = xin(177) + dxij*xin(165)
                                      yin(168) = yin(177) + dyij*yin(165)
                                      zin(168) = zin(177) + dzij*zin(165)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    3

                                      xin(180) = xin(189) + dxij*xin(177)
                                      yin(180) = yin(189) + dyij*yin(177)
                                      zin(180) = zin(189) + dzij*zin(177)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

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

                                      ! i1 = in(1) =  193

                                      xin(193) = 1.0_dp
                                      yin(193) = 1.0_dp
                                      zin(193) = f00

                                      ! i2 = in(2) =  205
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(205) = xc00
                                      yin(205) = yc00
                                      zin(205) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  194

                                      xin(194) = xcp00
                                      yin(194) = ycp00
                                      zin(194) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  206
                                      ! i2 =  205

                                      xin(206) = xcp00*xin(205) + cp10
                                      yin(206) = ycp00*yin(205) + cp10
                                      zin(206) = zcp00*zin(205) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  205

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  217
                                      ! i3 =  193
                                      ! i4 =  205

                                      xin(217) = c10*xin(193) + xc00*xin(205)
                                      yin(217) = c10*yin(193) + yc00*yin(205)
                                      zin(217) = c10*zin(193) + zc00*zin(205)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  218
                                      ! i5 =  217
                                      ! i4 =  205

                                      xin(218) = xcp00*xin(217) + cp10*xin(205)
                                      yin(218) = ycp00*yin(217) + cp10*yin(205)
                                      zin(218) = zcp00*zin(217) + cp10*zin(205)

                                      ! ------------------

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  217

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  229
                                      ! i3 =  205
                                      ! i4 =  217

                                      xin(229) = c10*xin(205) + xc00*xin(217)
                                      yin(229) = c10*yin(205) + yc00*yin(217)
                                      zin(229) = c10*zin(205) + zc00*zin(217)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  230
                                      ! i5 =  229
                                      ! i4 =  217

                                      xin(230) = xcp00*xin(229) + cp10*xin(217)
                                      yin(230) = ycp00*yin(229) + cp10*yin(217)
                                      zin(230) = zcp00*zin(229) + cp10*zin(217)

                                      ! ------------------

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  229

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  232
                                      ! i3 =  217
                                      ! i4 =  229

                                      xin(232) = c10*xin(217) + xc00*xin(229)
                                      yin(232) = c10*yin(217) + yc00*yin(229)
                                      zin(232) = c10*zin(217) + zc00*zin(229)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  233
                                      ! i5 =  232
                                      ! i4 =  229

                                      xin(233) = xcp00*xin(232) + cp10*xin(229)
                                      yin(233) = ycp00*yin(232) + cp10*yin(229)
                                      zin(233) = zcp00*zin(232) + cp10*zin(229)

                                      ! ------------------

                                      ! i3 = i4 =  229
                                      ! i4 = i5 =  232

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  235
                                      ! i3 =  229
                                      ! i4 =  232

                                      xin(235) = c10*xin(229) + xc00*xin(232)
                                      yin(235) = c10*yin(229) + yc00*yin(232)
                                      zin(235) = c10*zin(229) + zc00*zin(232)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  236
                                      ! i5 =  235
                                      ! i4 =  232

                                      xin(236) = xcp00*xin(235) + cp10*xin(232)
                                      yin(236) = ycp00*yin(235) + cp10*yin(232)
                                      zin(236) = zcp00*zin(235) + cp10*zin(232)

                                      ! ------------------

                                      ! i3 = i4 =  232
                                      ! i4 = i5 =  235

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  238
                                      ! i3 =  232
                                      ! i4 =  235

                                      xin(238) = c10*xin(232) + xc00*xin(235)
                                      yin(238) = c10*yin(232) + yc00*yin(235)
                                      zin(238) = c10*zin(232) + zc00*zin(235)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  239
                                      ! i5 =  238
                                      ! i4 =  235

                                      xin(239) = xcp00*xin(238) + cp10*xin(235)
                                      yin(239) = ycp00*yin(238) + cp10*yin(235)
                                      zin(239) = zcp00*zin(238) + cp10*zin(235)

                                      ! ------------------

                                      ! i3 = i4 =  235
                                      ! i4 = i5 =  238

                                      ! n =    7

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  193
                                      ! i4 = i1+k2 =  194

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  195
                                      ! i3 =  193
                                      ! i4 =  194

                                      xin(195) = cp01*xin(193) + xcp00*xin(194)
                                      yin(195) = cp01*yin(193) + ycp00*yin(194)
                                      zin(195) = cp01*zin(193) + zcp00*zin(194)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  207

                                      xin(207) = xc00*xin(195) + c01*xin(194)
                                      yin(207) = yc00*yin(195) + c01*yin(194)
                                      zin(207) = zc00*zin(195) + c01*zin(194)

                                      ! ------------------

                                      ! i3 = i4 =  194
                                      ! i4 = i5 =  195

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  205

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    6

                                      ! i5 = in(nn+1) =  217

                                      xin(219) = c10*xin(195) + xc00*xin(207) + c01*xin(206)
                                      yin(219) = c10*yin(195) + yc00*yin(207) + c01*yin(206)
                                      zin(219) = c10*zin(195) + zc00*zin(207) + c01*zin(206)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  229

                                      xin(231) = c10*xin(207) + xc00*xin(219) + c01*xin(218)
                                      yin(231) = c10*yin(207) + yc00*yin(219) + c01*yin(218)
                                      zin(231) = c10*zin(207) + zc00*zin(219) + c01*zin(218)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  229

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  232

                                      xin(234) = c10*xin(219) + xc00*xin(231) + c01*xin(230)
                                      yin(234) = c10*yin(219) + yc00*yin(231) + c01*yin(230)
                                      zin(234) = c10*zin(219) + zc00*zin(231) + c01*zin(230)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  229
                                      ! i4 = i5 =  232

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  235

                                      xin(237) = c10*xin(231) + xc00*xin(234) + c01*xin(233)
                                      yin(237) = c10*yin(231) + yc00*yin(234) + c01*yin(233)
                                      zin(237) = c10*zin(231) + zc00*zin(234) + c01*zin(233)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  232
                                      ! i4 = i5 =  235

                                      ! nn =    6

                                      ! i5 = in(nn+1) =  238

                                      xin(240) = c10*xin(234) + xc00*xin(237) + c01*xin(236)
                                      yin(240) = c10*yin(234) + yc00*yin(237) + c01*yin(236)
                                      zin(240) = c10*zin(234) + zc00*zin(237) + c01*zin(236)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  235
                                      ! i4 = i5 =  238

                                      ! nn =    7

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  238

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  238

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  235

                                      xin(238) = xin(238) + dxij*xin(235)
                                      yin(238) = yin(238) + dyij*yin(235)
                                      zin(238) = zin(238) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  232

                                      xin(235) = xin(235) + dxij*xin(232)
                                      yin(235) = yin(235) + dyij*yin(232)
                                      zin(235) = zin(235) + dzij*zin(232)

                                      ! i3 = i4 =  232
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  229

                                      xin(232) = xin(232) + dxij*xin(229)
                                      yin(232) = yin(232) + dyij*yin(229)
                                      zin(232) = zin(232) + dzij*zin(229)

                                      ! i3 = i4 =  229
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  238

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  235

                                      xin(238) = xin(238) + dxij*xin(235)
                                      yin(238) = yin(238) + dyij*yin(235)
                                      zin(238) = zin(238) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  232

                                      xin(235) = xin(235) + dxij*xin(232)
                                      yin(235) = yin(235) + dyij*yin(232)
                                      zin(235) = zin(235) + dzij*zin(232)

                                      ! i3 = i4 =  232
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  238

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  235

                                      xin(238) = xin(238) + dxij*xin(235)
                                      yin(238) = yin(238) + dyij*yin(235)
                                      zin(238) = zin(238) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  196

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  196

                                      ! do ni = 1,    3

                                      xin(196) = xin(205) + dxij*xin(193)
                                      yin(196) = yin(205) + dyij*yin(193)
                                      zin(196) = zin(205) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                      ! ni =    2

                                      xin(208) = xin(217) + dxij*xin(205)
                                      yin(208) = yin(217) + dyij*yin(205)
                                      zin(208) = zin(217) + dzij*zin(205)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  220

                                      ! ni =    3

                                      xin(220) = xin(229) + dxij*xin(217)
                                      yin(220) = yin(229) + dyij*yin(217)
                                      zin(220) = zin(229) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  232

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  199

                                      ! nj =    2

                                      ! i4 = i3 =  199

                                      ! do ni = 1,    3

                                      xin(199) = xin(208) + dxij*xin(196)
                                      yin(199) = yin(208) + dyij*yin(196)
                                      zin(199) = zin(208) + dzij*zin(196)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    2

                                      xin(211) = xin(220) + dxij*xin(208)
                                      yin(211) = yin(220) + dyij*yin(208)
                                      zin(211) = zin(220) + dzij*zin(208)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  223

                                      ! ni =    3

                                      xin(223) = xin(232) + dxij*xin(220)
                                      yin(223) = yin(232) + dyij*yin(220)
                                      zin(223) = zin(232) + dzij*zin(220)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  202

                                      ! nj =    3

                                      ! i4 = i3 =  202

                                      ! do ni = 1,    3

                                      xin(202) = xin(211) + dxij*xin(199)
                                      yin(202) = yin(211) + dyij*yin(199)
                                      zin(202) = zin(211) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                      ! ni =    2

                                      xin(214) = xin(223) + dxij*xin(211)
                                      yin(214) = yin(223) + dyij*yin(211)
                                      zin(214) = zin(223) + dzij*zin(211)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  226

                                      ! ni =    3

                                      xin(226) = xin(235) + dxij*xin(223)
                                      yin(226) = yin(235) + dyij*yin(223)
                                      zin(226) = zin(235) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  238

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  205

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  239

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  236

                                      xin(239) = xin(239) + dxij*xin(236)
                                      yin(239) = yin(239) + dyij*yin(236)
                                      zin(239) = zin(239) + dzij*zin(236)

                                      ! i3 = i4 =  236
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  233

                                      xin(236) = xin(236) + dxij*xin(233)
                                      yin(236) = yin(236) + dyij*yin(233)
                                      zin(236) = zin(236) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  230

                                      xin(233) = xin(233) + dxij*xin(230)
                                      yin(233) = yin(233) + dyij*yin(230)
                                      zin(233) = zin(233) + dzij*zin(230)

                                      ! i3 = i4 =  230
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  239

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  236

                                      xin(239) = xin(239) + dxij*xin(236)
                                      yin(239) = yin(239) + dyij*yin(236)
                                      zin(239) = zin(239) + dzij*zin(236)

                                      ! i3 = i4 =  236
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  233

                                      xin(236) = xin(236) + dxij*xin(233)
                                      yin(236) = yin(236) + dyij*yin(233)
                                      zin(236) = zin(236) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  239

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  236

                                      xin(239) = xin(239) + dxij*xin(236)
                                      yin(239) = yin(239) + dyij*yin(236)
                                      zin(239) = zin(239) + dzij*zin(236)

                                      ! i3 = i4 =  236
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  197

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  197

                                      ! do ni = 1,    3

                                      xin(197) = xin(206) + dxij*xin(194)
                                      yin(197) = yin(206) + dyij*yin(194)
                                      zin(197) = zin(206) + dzij*zin(194)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    2

                                      xin(209) = xin(218) + dxij*xin(206)
                                      yin(209) = yin(218) + dyij*yin(206)
                                      zin(209) = zin(218) + dzij*zin(206)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  221

                                      ! ni =    3

                                      xin(221) = xin(230) + dxij*xin(218)
                                      yin(221) = yin(230) + dyij*yin(218)
                                      zin(221) = zin(230) + dzij*zin(218)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  200

                                      ! nj =    2

                                      ! i4 = i3 =  200

                                      ! do ni = 1,    3

                                      xin(200) = xin(209) + dxij*xin(197)
                                      yin(200) = yin(209) + dyij*yin(197)
                                      zin(200) = zin(209) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  212

                                      ! ni =    2

                                      xin(212) = xin(221) + dxij*xin(209)
                                      yin(212) = yin(221) + dyij*yin(209)
                                      zin(212) = zin(221) + dzij*zin(209)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  224

                                      ! ni =    3

                                      xin(224) = xin(233) + dxij*xin(221)
                                      yin(224) = yin(233) + dyij*yin(221)
                                      zin(224) = zin(233) + dzij*zin(221)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  236

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  203

                                      ! nj =    3

                                      ! i4 = i3 =  203

                                      ! do ni = 1,    3

                                      xin(203) = xin(212) + dxij*xin(200)
                                      yin(203) = yin(212) + dyij*yin(200)
                                      zin(203) = zin(212) + dzij*zin(200)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    2

                                      xin(215) = xin(224) + dxij*xin(212)
                                      yin(215) = yin(224) + dyij*yin(212)
                                      zin(215) = zin(224) + dzij*zin(212)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  227

                                      ! ni =    3

                                      xin(227) = xin(236) + dxij*xin(224)
                                      yin(227) = yin(236) + dyij*yin(224)
                                      zin(227) = zin(236) + dzij*zin(224)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  206

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  240

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  237

                                      xin(240) = xin(240) + dxij*xin(237)
                                      yin(240) = yin(240) + dyij*yin(237)
                                      zin(240) = zin(240) + dzij*zin(237)

                                      ! i3 = i4 =  237
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  234

                                      xin(237) = xin(237) + dxij*xin(234)
                                      yin(237) = yin(237) + dyij*yin(234)
                                      zin(237) = zin(237) + dzij*zin(234)

                                      ! i3 = i4 =  234
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  231

                                      xin(234) = xin(234) + dxij*xin(231)
                                      yin(234) = yin(234) + dyij*yin(231)
                                      zin(234) = zin(234) + dzij*zin(231)

                                      ! i3 = i4 =  231
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  240

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  237

                                      xin(240) = xin(240) + dxij*xin(237)
                                      yin(240) = yin(240) + dyij*yin(237)
                                      zin(240) = zin(240) + dzij*zin(237)

                                      ! i3 = i4 =  237
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =  234

                                      xin(237) = xin(237) + dxij*xin(234)
                                      yin(237) = yin(237) + dyij*yin(234)
                                      zin(237) = zin(237) + dzij*zin(234)

                                      ! i3 = i4 =  234
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =  240

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  237

                                      xin(240) = xin(240) + dxij*xin(237)
                                      yin(240) = yin(240) + dyij*yin(237)
                                      zin(240) = zin(240) + dzij*zin(237)

                                      ! i3 = i4 =  237
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  198

                                      ! do nj = 1,    3

                                      ! i4 = i3 =  198

                                      ! do ni = 1,    3

                                      xin(198) = xin(207) + dxij*xin(195)
                                      yin(198) = yin(207) + dyij*yin(195)
                                      zin(198) = zin(207) + dzij*zin(195)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  210

                                      ! ni =    2

                                      xin(210) = xin(219) + dxij*xin(207)
                                      yin(210) = yin(219) + dyij*yin(207)
                                      zin(210) = zin(219) + dzij*zin(207)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  222

                                      ! ni =    3

                                      xin(222) = xin(231) + dxij*xin(219)
                                      yin(222) = yin(231) + dyij*yin(219)
                                      zin(222) = zin(231) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  234

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  201

                                      ! nj =    2

                                      ! i4 = i3 =  201

                                      ! do ni = 1,    3

                                      xin(201) = xin(210) + dxij*xin(198)
                                      yin(201) = yin(210) + dyij*yin(198)
                                      zin(201) = zin(210) + dzij*zin(198)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    2

                                      xin(213) = xin(222) + dxij*xin(210)
                                      yin(213) = yin(222) + dyij*yin(210)
                                      zin(213) = zin(222) + dzij*zin(210)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  225

                                      ! ni =    3

                                      xin(225) = xin(234) + dxij*xin(222)
                                      yin(225) = yin(234) + dyij*yin(222)
                                      zin(225) = zin(234) + dzij*zin(222)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  237

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  204

                                      ! nj =    3

                                      ! i4 = i3 =  204

                                      ! do ni = 1,    3

                                      xin(204) = xin(213) + dxij*xin(201)
                                      yin(204) = yin(213) + dyij*yin(201)
                                      zin(204) = zin(213) + dzij*zin(201)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    2

                                      xin(216) = xin(225) + dxij*xin(213)
                                      yin(216) = yin(225) + dyij*yin(213)
                                      zin(216) = zin(225) + dzij*zin(213)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  228

                                      ! ni =    3

                                      xin(228) = xin(237) + dxij*xin(225)
                                      yin(228) = yin(237) + dyij*yin(225)
                                      zin(228) = zin(237) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  207

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  240

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 600! loop over all integrals

                                        l = n - 6*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d33bra(j)*d02ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 48)*yin(my + 48)*zin(mz + 48) & ! root  2
                                                        + xin(mx + 96)*yin(my + 96)*zin(mz + 96) & ! root  3
                                                        + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  4
                                                        + xin(mx + 192)*yin(my + 192)*zin(mz + 192)) ! root  5

                                        j = int(n/6) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*60 ! Stride between functions in i

                                    do j = 1, maxj2

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*6 + ip ! Add stride between functions in j

                                      do k = 1, 6 ! # of cartesians in k

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
                              deallocate (n02ket)
                              deallocate (xint02ket)

                              end subroutine int3320
                              end submodule
