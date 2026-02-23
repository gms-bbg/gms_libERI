! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3300_impl
contains
  module subroutine int3300(ff_pair, ss_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ff_pair, ss_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n33bra(:), n00ket(:)
    real(dp), allocatable :: xint33bra(:), xint00ket(:)
    integer(kind=int64) :: nffbra, nssket
    real(dp) :: scutffbra, scutssket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxj2, maxl, maxl2
    integer(kind=int64) :: n, i1, i3, i4, i5, nm, nn, km, nj, ni
    real(dp) :: cp10, c10
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(64), yin(64), zin(64)
    real(dp) :: eri_value(100)
    real(dp) :: d33bra(100), d00ket(1)
    integer(kind=int64) :: ix(10), jx(10), kx(1), lx(1)
    integer(kind=int64) :: iy(10), jy(10), ky(1), ly(1)
    integer(kind=int64) :: iz(10), jz(10), kz(1), lz(1)
    integer(kind=int64) :: in(7), in1(7), kn(1)
    integer(kind=int64) :: ijx(100), ijy(100), ijz(100)
    integer(kind=int64) :: klx(1), kly(1), klz(1)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj, kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 5
    in1(3) = 9
    in1(4) = 13
    in1(5) = 14
    in1(6) = 15
    in1(7) = 16

    kn(1) = 0

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 0

    jx(1) = 3
    jx(2) = 0
    jx(3) = 0
    jx(4) = 2
    jx(5) = 2
    jx(6) = 1
    jx(7) = 0
    jx(8) = 1
    jx(9) = 0
    jx(10) = 1

    ix(1) = 13
    ix(2) = 1
    ix(3) = 1
    ix(4) = 9
    ix(5) = 9
    ix(6) = 5
    ix(7) = 1
    ix(8) = 5
    ix(9) = 1
    ix(10) = 5

    ! y-arrays

    ly(1) = 0

    ky(1) = 0

    jy(1) = 0
    jy(2) = 3
    jy(3) = 0
    jy(4) = 1
    jy(5) = 0
    jy(6) = 2
    jy(7) = 2
    jy(8) = 0
    jy(9) = 1
    jy(10) = 1

    iy(1) = 1
    iy(2) = 13
    iy(3) = 1
    iy(4) = 5
    iy(5) = 1
    iy(6) = 9
    iy(7) = 9
    iy(8) = 1
    iy(9) = 5
    iy(10) = 5

    ! z-arrays

    lz(1) = 0

    kz(1) = 0

    jz(1) = 0
    jz(2) = 0
    jz(3) = 3
    jz(4) = 0
    jz(5) = 1
    jz(6) = 0
    jz(7) = 1
    jz(8) = 2
    jz(9) = 2
    jz(10) = 1

    iz(1) = 1
    iz(2) = 1
    iz(3) = 13
    iz(4) = 1
    iz(5) = 5
    iz(6) = 1
    iz(7) = 5
    iz(8) = 9
    iz(9) = 9
    iz(10) = 5

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 16
    ijx(2) = 13
    ijx(3) = 13
    ijx(4) = 15
    ijx(5) = 15
    ijx(6) = 14
    ijx(7) = 13
    ijx(8) = 14
    ijx(9) = 13
    ijx(10) = 14
    ijx(11) = 4
    ijx(12) = 1
    ijx(13) = 1
    ijx(14) = 3
    ijx(15) = 3
    ijx(16) = 2
    ijx(17) = 1
    ijx(18) = 2
    ijx(19) = 1
    ijx(20) = 2
    ijx(21) = 4
    ijx(22) = 1
    ijx(23) = 1
    ijx(24) = 3
    ijx(25) = 3
    ijx(26) = 2
    ijx(27) = 1
    ijx(28) = 2
    ijx(29) = 1
    ijx(30) = 2
    ijx(31) = 12
    ijx(32) = 9
    ijx(33) = 9
    ijx(34) = 11
    ijx(35) = 11
    ijx(36) = 10
    ijx(37) = 9
    ijx(38) = 10
    ijx(39) = 9
    ijx(40) = 10
    ijx(41) = 12
    ijx(42) = 9
    ijx(43) = 9
    ijx(44) = 11
    ijx(45) = 11
    ijx(46) = 10
    ijx(47) = 9
    ijx(48) = 10
    ijx(49) = 9
    ijx(50) = 10
    ijx(51) = 8
    ijx(52) = 5
    ijx(53) = 5
    ijx(54) = 7
    ijx(55) = 7
    ijx(56) = 6
    ijx(57) = 5
    ijx(58) = 6
    ijx(59) = 5
    ijx(60) = 6
    ijx(61) = 4
    ijx(62) = 1
    ijx(63) = 1
    ijx(64) = 3
    ijx(65) = 3
    ijx(66) = 2
    ijx(67) = 1
    ijx(68) = 2
    ijx(69) = 1
    ijx(70) = 2
    ijx(71) = 8
    ijx(72) = 5
    ijx(73) = 5
    ijx(74) = 7
    ijx(75) = 7
    ijx(76) = 6
    ijx(77) = 5
    ijx(78) = 6
    ijx(79) = 5
    ijx(80) = 6
    ijx(81) = 4
    ijx(82) = 1
    ijx(83) = 1
    ijx(84) = 3
    ijx(85) = 3
    ijx(86) = 2
    ijx(87) = 1
    ijx(88) = 2
    ijx(89) = 1
    ijx(90) = 2
    ijx(91) = 8
    ijx(92) = 5
    ijx(93) = 5
    ijx(94) = 7
    ijx(95) = 7
    ijx(96) = 6
    ijx(97) = 5
    ijx(98) = 6
    ijx(99) = 5
    ijx(100) = 6

    ijy(1) = 1
    ijy(2) = 4
    ijy(3) = 1
    ijy(4) = 2
    ijy(5) = 1
    ijy(6) = 3
    ijy(7) = 3
    ijy(8) = 1
    ijy(9) = 2
    ijy(10) = 2
    ijy(11) = 13
    ijy(12) = 16
    ijy(13) = 13
    ijy(14) = 14
    ijy(15) = 13
    ijy(16) = 15
    ijy(17) = 15
    ijy(18) = 13
    ijy(19) = 14
    ijy(20) = 14
    ijy(21) = 1
    ijy(22) = 4
    ijy(23) = 1
    ijy(24) = 2
    ijy(25) = 1
    ijy(26) = 3
    ijy(27) = 3
    ijy(28) = 1
    ijy(29) = 2
    ijy(30) = 2
    ijy(31) = 5
    ijy(32) = 8
    ijy(33) = 5
    ijy(34) = 6
    ijy(35) = 5
    ijy(36) = 7
    ijy(37) = 7
    ijy(38) = 5
    ijy(39) = 6
    ijy(40) = 6
    ijy(41) = 1
    ijy(42) = 4
    ijy(43) = 1
    ijy(44) = 2
    ijy(45) = 1
    ijy(46) = 3
    ijy(47) = 3
    ijy(48) = 1
    ijy(49) = 2
    ijy(50) = 2
    ijy(51) = 9
    ijy(52) = 12
    ijy(53) = 9
    ijy(54) = 10
    ijy(55) = 9
    ijy(56) = 11
    ijy(57) = 11
    ijy(58) = 9
    ijy(59) = 10
    ijy(60) = 10
    ijy(61) = 9
    ijy(62) = 12
    ijy(63) = 9
    ijy(64) = 10
    ijy(65) = 9
    ijy(66) = 11
    ijy(67) = 11
    ijy(68) = 9
    ijy(69) = 10
    ijy(70) = 10
    ijy(71) = 1
    ijy(72) = 4
    ijy(73) = 1
    ijy(74) = 2
    ijy(75) = 1
    ijy(76) = 3
    ijy(77) = 3
    ijy(78) = 1
    ijy(79) = 2
    ijy(80) = 2
    ijy(81) = 5
    ijy(82) = 8
    ijy(83) = 5
    ijy(84) = 6
    ijy(85) = 5
    ijy(86) = 7
    ijy(87) = 7
    ijy(88) = 5
    ijy(89) = 6
    ijy(90) = 6
    ijy(91) = 5
    ijy(92) = 8
    ijy(93) = 5
    ijy(94) = 6
    ijy(95) = 5
    ijy(96) = 7
    ijy(97) = 7
    ijy(98) = 5
    ijy(99) = 6
    ijy(100) = 6

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 4
    ijz(4) = 1
    ijz(5) = 2
    ijz(6) = 1
    ijz(7) = 2
    ijz(8) = 3
    ijz(9) = 3
    ijz(10) = 2
    ijz(11) = 1
    ijz(12) = 1
    ijz(13) = 4
    ijz(14) = 1
    ijz(15) = 2
    ijz(16) = 1
    ijz(17) = 2
    ijz(18) = 3
    ijz(19) = 3
    ijz(20) = 2
    ijz(21) = 13
    ijz(22) = 13
    ijz(23) = 16
    ijz(24) = 13
    ijz(25) = 14
    ijz(26) = 13
    ijz(27) = 14
    ijz(28) = 15
    ijz(29) = 15
    ijz(30) = 14
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 4
    ijz(34) = 1
    ijz(35) = 2
    ijz(36) = 1
    ijz(37) = 2
    ijz(38) = 3
    ijz(39) = 3
    ijz(40) = 2
    ijz(41) = 5
    ijz(42) = 5
    ijz(43) = 8
    ijz(44) = 5
    ijz(45) = 6
    ijz(46) = 5
    ijz(47) = 6
    ijz(48) = 7
    ijz(49) = 7
    ijz(50) = 6
    ijz(51) = 1
    ijz(52) = 1
    ijz(53) = 4
    ijz(54) = 1
    ijz(55) = 2
    ijz(56) = 1
    ijz(57) = 2
    ijz(58) = 3
    ijz(59) = 3
    ijz(60) = 2
    ijz(61) = 5
    ijz(62) = 5
    ijz(63) = 8
    ijz(64) = 5
    ijz(65) = 6
    ijz(66) = 5
    ijz(67) = 6
    ijz(68) = 7
    ijz(69) = 7
    ijz(70) = 6
    ijz(71) = 9
    ijz(72) = 9
    ijz(73) = 12
    ijz(74) = 9
    ijz(75) = 10
    ijz(76) = 9
    ijz(77) = 10
    ijz(78) = 11
    ijz(79) = 11
    ijz(80) = 10
    ijz(81) = 9
    ijz(82) = 9
    ijz(83) = 12
    ijz(84) = 9
    ijz(85) = 10
    ijz(86) = 9
    ijz(87) = 10
    ijz(88) = 11
    ijz(89) = 11
    ijz(90) = 10
    ijz(91) = 5
    ijz(92) = 5
    ijz(93) = 8
    ijz(94) = 5
    ijz(95) = 6
    ijz(96) = 5
    ijz(97) = 6
    ijz(98) = 7
    ijz(99) = 7
    ijz(100) = 6

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 0

    kly(1) = 0

    klz(1) = 0

    allocate (n33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (xint33bra(res%n_f_shl*(res%n_f_shl + 1)/2))
    allocate (n00ket(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00ket(res%n_s_shl*(res%n_s_shl + 1)/2))

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

    scutssket = cutoff_schwarz/maxval(ss_pair%xints)
    nssket = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
      if (ss_pair%xints(ij) .ge. scutssket) then
        nssket = nssket + 1
        xint00ket(nssket) = ss_pair%xints(ij)
        n00ket(nssket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nffbra*nssket) .le. nchunksize_int64) nchunksize_int64 = nffbra*nssket
    ntile = int(nffbra*nssket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nffbra*nssket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nffbra, xint33bra, n33bra, xint00ket, n00ket, ff_pair, ss_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d00ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d33bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,nm,nn,km,nj,ni) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,maxl,maxl2,iandj,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nffbra) + 1
              kl_tmp = (iquart - 1)/nffbra + 1

              test = xint33bra(ij_tmp)*xint00ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n33bra(ij_tmp)
                kl = n00ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_f_shl(jsh_tmp)
                ksh = res%i_s_shl(ksh_tmp)
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

                  t_expon_cd = ss_pair%t_expon_ab(ss_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = ss_pair%expon_a(ss_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = ss_pair%expon_b(ss_pair%pair_loc(kl) + ket_loop)
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

                  d00ket(1) = ss_pair%d_coeff_alt(ss_pair%pair_loc(kl) + ket_loop)*twopi_5_2

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

                    if (xx .ge. 49.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.1453035215033170D+00
                      rts(2) = factr*0.1339097288126361D+01
                      rts(3) = factr*0.3926963501358284D+01
                      rts(4) = factr*0.8588635689012030D+01

                      wts(1) = factw*0.6611470125582407D+00
                      wts(2) = factw*0.2078023258148924D+00
                      wts(3) = factw*0.1707798300741349D-01
                      wts(4) = factw*0.1996040722113677D-03

                    else ! "regular" evaluation

                      rgrid(1) = 0.2412610298614493D-05
                      rgrid(2) = 0.6668254930138793D-04
                      rgrid(3) = 0.3995628201530626D-03
                      rgrid(4) = 0.1361608249860346D-02
                      rgrid(5) = 0.3448006938362446D-02
                      rgrid(6) = 0.7261957338041792D-02
                      rgrid(7) = 0.1348183025995716D-01
                      rgrid(8) = 0.2282358087416106D-01
                      rgrid(9) = 0.3600009444917834D-01
                      rgrid(10) = 0.5367929502127710D-01
                      rgrid(11) = 0.7644291300781368D-01
                      rgrid(12) = 0.1047477930875140D+00
                      rgrid(13) = 0.1388915279581130D+00
                      rgrid(14) = 0.1789840307741864D+00
                      rgrid(15) = 0.2249264163663510D+00
                      rgrid(16) = 0.2763982589216690D+00
                      rgrid(17) = 0.3328539443827702D+00
                      rgrid(18) = 0.3935284541260027D+00
                      rgrid(19) = 0.4574525186183916D+00
                      rgrid(20) = 0.5234766825459025D+00
                      rgrid(21) = 0.5903034431632973D+00
                      rgrid(22) = 0.6565262774384205D+00
                      rgrid(23) = 0.7206740756674773D+00
                      rgrid(24) = 0.7812592623647830D+00
                      rgrid(25) = 0.8368277197208105D+00
                      rgrid(26) = 0.8860085427304147D+00
                      rgrid(27) = 0.9275616556791347D+00
                      rgrid(28) = 0.9604214277884603D+00
                      rgrid(29) = 0.9837348058290485D+00
                      rgrid(30) = 0.9968958966849477D+00

                      wgrid(1) = 0.3984096248083391D-02*exp(-xx*0.2412610298614493D-05)
                      wgrid(2) = 0.9233234155545691D-02*exp(-xx*0.6668254930138793D-04)
                      wgrid(3) = 0.1439235394166157D-01*exp(-xx*0.3995628201530626D-03)
                      wgrid(4) = 0.1939959628481344D-01*exp(-xx*0.1361608249860346D-02)
                      wgrid(5) = 0.2420133641529743D-01*exp(-xx*0.3448006938362446D-02)
                      wgrid(6) = 0.2874657810880918D-01*exp(-xx*0.7261957338041792D-02)
                      wgrid(7) = 0.3298711494109012D-01*exp(-xx*0.1348183025995716D-01)
                      wgrid(8) = 0.3687798736885270D-01*exp(-xx*0.2282358087416106D-01)
                      wgrid(9) = 0.4037794761470943D-01*exp(-xx*0.3600009444917834D-01)
                      wgrid(10) = 0.4344989360054208D-01*exp(-xx*0.5367929502127710D-01)
                      wgrid(11) = 0.4606126111889302D-01*exp(-xx*0.7644291300781368D-01)
                      wgrid(12) = 0.4818436858732206D-01*exp(-xx*0.1047477930875140D+00)
                      wgrid(13) = 0.4979671029339763D-01*exp(-xx*0.1388915279581130D+00)
                      wgrid(14) = 0.5088119487420216D-01*exp(-xx*0.1789840307741864D+00)
                      wgrid(15) = 0.5142632644677943D-01*exp(-xx*0.2249264163663510D+00)
                      wgrid(16) = 0.5142632644677995D-01*exp(-xx*0.2763982589216690D+00)
                      wgrid(17) = 0.5088119487420249D-01*exp(-xx*0.3328539443827702D+00)
                      wgrid(18) = 0.4979671029339827D-01*exp(-xx*0.3935284541260027D+00)
                      wgrid(19) = 0.4818436858732192D-01*exp(-xx*0.4574525186183916D+00)
                      wgrid(20) = 0.4606126111889324D-01*exp(-xx*0.5234766825459025D+00)
                      wgrid(21) = 0.4344989360054120D-01*exp(-xx*0.5903034431632973D+00)
                      wgrid(22) = 0.4037794761471003D-01*exp(-xx*0.6565262774384205D+00)
                      wgrid(23) = 0.3687798736885173D-01*exp(-xx*0.7206740756674773D+00)
                      wgrid(24) = 0.3298711494109043D-01*exp(-xx*0.7812592623647830D+00)
                      wgrid(25) = 0.2874657810880946D-01*exp(-xx*0.8368277197208105D+00)
                      wgrid(26) = 0.2420133641529735D-01*exp(-xx*0.8860085427304147D+00)
                      wgrid(27) = 0.1939959628481361D-01*exp(-xx*0.9275616556791347D+00)
                      wgrid(28) = 0.1439235394166174D-01*exp(-xx*0.9604214277884603D+00)
                      wgrid(29) = 0.9233234155545632D-02*exp(-xx*0.9837348058290485D+00)
                      wgrid(30) = 0.3984096248083247D-02*exp(-xx*0.9968958966849477D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 30
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 30
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 3

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 30

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
                        wrk(4) = 0.0D+00
                        do 100 kk = 2, 4

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 4

                            jj = 0

105                         do 110 m = l, 4
                              if (m .eq. 4) go to 120
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

                                do 300 ii = 2, 4

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 4
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

                                    do 310 kk = 1, 4
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 4
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

                                      ! i2 = in(2) =    5
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(5) = xc00
                                      yin(5) = yc00
                                      zin(5) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    5

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    9
                                      ! i3 =    1
                                      ! i4 =    5

                                      xin(9) = c10*xin(1) + xc00*xin(5)
                                      yin(9) = c10*yin(1) + yc00*yin(5)
                                      zin(9) = c10*zin(1) + zc00*zin(5)

                                      ! i3 = i4 =    5
                                      ! i4 = i5 =    9

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   13
                                      ! i3 =    5
                                      ! i4 =    9

                                      xin(13) = c10*xin(5) + xc00*xin(9)
                                      yin(13) = c10*yin(5) + yc00*yin(9)
                                      zin(13) = c10*zin(5) + zc00*zin(9)

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   13

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   14
                                      ! i3 =    9
                                      ! i4 =   13

                                      xin(14) = c10*xin(9) + xc00*xin(13)
                                      yin(14) = c10*yin(9) + yc00*yin(13)
                                      zin(14) = c10*zin(9) + zc00*zin(13)

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   14

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   15
                                      ! i3 =   13
                                      ! i4 =   14

                                      xin(15) = c10*xin(13) + xc00*xin(14)
                                      yin(15) = c10*yin(13) + yc00*yin(14)
                                      zin(15) = c10*zin(13) + zc00*zin(14)

                                      ! i3 = i4 =   14
                                      ! i4 = i5 =   15

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   16
                                      ! i3 =   14
                                      ! i4 =   15

                                      xin(16) = c10*xin(14) + xc00*xin(15)
                                      yin(16) = c10*yin(14) + yc00*yin(15)
                                      zin(16) = c10*zin(14) + zc00*zin(15)

                                      ! i3 = i4 =   15
                                      ! i4 = i5 =   16

                                      ! n =    7

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   16

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   16

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   15

                                      xin(16) = xin(16) + dxij*xin(15)
                                      yin(16) = yin(16) + dyij*yin(15)
                                      zin(16) = zin(16) + dzij*zin(15)

                                      ! i3 = i4 =   15
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   14

                                      xin(15) = xin(15) + dxij*xin(14)
                                      yin(15) = yin(15) + dyij*yin(14)
                                      zin(15) = zin(15) + dzij*zin(14)

                                      ! i3 = i4 =   14
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   13

                                      xin(14) = xin(14) + dxij*xin(13)
                                      yin(14) = yin(14) + dyij*yin(13)
                                      zin(14) = zin(14) + dzij*zin(13)

                                      ! i3 = i4 =   13
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   16

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   15

                                      xin(16) = xin(16) + dxij*xin(15)
                                      yin(16) = yin(16) + dyij*yin(15)
                                      zin(16) = zin(16) + dzij*zin(15)

                                      ! i3 = i4 =   15
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   14

                                      xin(15) = xin(15) + dxij*xin(14)
                                      yin(15) = yin(15) + dyij*yin(14)
                                      zin(15) = zin(15) + dzij*zin(14)

                                      ! i3 = i4 =   14
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   16

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   15

                                      xin(16) = xin(16) + dxij*xin(15)
                                      yin(16) = yin(16) + dyij*yin(15)
                                      zin(16) = zin(16) + dzij*zin(15)

                                      ! i3 = i4 =   15
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    2

                                      ! do nj = 1,    3

                                      ! i4 = i3 =    2

                                      ! do ni = 1,    3

                                      xin(2) = xin(5) + dxij*xin(1)
                                      yin(2) = yin(5) + dyij*yin(1)
                                      zin(2) = zin(5) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    6

                                      ! ni =    2

                                      xin(6) = xin(9) + dxij*xin(5)
                                      yin(6) = yin(9) + dyij*yin(5)
                                      zin(6) = zin(9) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   10

                                      ! ni =    3

                                      xin(10) = xin(13) + dxij*xin(9)
                                      yin(10) = yin(13) + dyij*yin(9)
                                      zin(10) = zin(13) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   14

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    3

                                      ! nj =    2

                                      ! i4 = i3 =    3

                                      ! do ni = 1,    3

                                      xin(3) = xin(6) + dxij*xin(2)
                                      yin(3) = yin(6) + dyij*yin(2)
                                      zin(3) = zin(6) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    7

                                      ! ni =    2

                                      xin(7) = xin(10) + dxij*xin(6)
                                      yin(7) = yin(10) + dyij*yin(6)
                                      zin(7) = zin(10) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   11

                                      ! ni =    3

                                      xin(11) = xin(14) + dxij*xin(10)
                                      yin(11) = yin(14) + dyij*yin(10)
                                      zin(11) = zin(14) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   15

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    4

                                      ! nj =    3

                                      ! i4 = i3 =    4

                                      ! do ni = 1,    3

                                      xin(4) = xin(7) + dxij*xin(3)
                                      yin(4) = yin(7) + dyij*yin(3)
                                      zin(4) = zin(7) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    8

                                      ! ni =    2

                                      xin(8) = xin(11) + dxij*xin(7)
                                      yin(8) = yin(11) + dyij*yin(7)
                                      zin(8) = zin(11) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   12

                                      ! ni =    3

                                      xin(12) = xin(15) + dxij*xin(11)
                                      yin(12) = yin(15) + dyij*yin(11)
                                      zin(12) = zin(15) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    5

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   16

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

                                      ! i1 = in(1) =   17

                                      xin(17) = 1.0_dp
                                      yin(17) = 1.0_dp
                                      zin(17) = f00

                                      ! i2 = in(2) =   21
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(21) = xc00
                                      yin(21) = yc00
                                      zin(21) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   17
                                      ! i4 = i2 =   21

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   25
                                      ! i3 =   17
                                      ! i4 =   21

                                      xin(25) = c10*xin(17) + xc00*xin(21)
                                      yin(25) = c10*yin(17) + yc00*yin(21)
                                      zin(25) = c10*zin(17) + zc00*zin(21)

                                      ! i3 = i4 =   21
                                      ! i4 = i5 =   25

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   29
                                      ! i3 =   21
                                      ! i4 =   25

                                      xin(29) = c10*xin(21) + xc00*xin(25)
                                      yin(29) = c10*yin(21) + yc00*yin(25)
                                      zin(29) = c10*zin(21) + zc00*zin(25)

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   30
                                      ! i3 =   25
                                      ! i4 =   29

                                      xin(30) = c10*xin(25) + xc00*xin(29)
                                      yin(30) = c10*yin(25) + yc00*yin(29)
                                      zin(30) = c10*zin(25) + zc00*zin(29)

                                      ! i3 = i4 =   29
                                      ! i4 = i5 =   30

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   31
                                      ! i3 =   29
                                      ! i4 =   30

                                      xin(31) = c10*xin(29) + xc00*xin(30)
                                      yin(31) = c10*yin(29) + yc00*yin(30)
                                      zin(31) = c10*zin(29) + zc00*zin(30)

                                      ! i3 = i4 =   30
                                      ! i4 = i5 =   31

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   32
                                      ! i3 =   30
                                      ! i4 =   31

                                      xin(32) = c10*xin(30) + xc00*xin(31)
                                      yin(32) = c10*yin(30) + yc00*yin(31)
                                      zin(32) = c10*zin(30) + zc00*zin(31)

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   32

                                      ! n =    7

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   32

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   32

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(32) = xin(32) + dxij*xin(31)
                                      yin(32) = yin(32) + dyij*yin(31)
                                      zin(32) = zin(32) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   30

                                      xin(31) = xin(31) + dxij*xin(30)
                                      yin(31) = yin(31) + dyij*yin(30)
                                      zin(31) = zin(31) + dzij*zin(30)

                                      ! i3 = i4 =   30
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   29

                                      xin(30) = xin(30) + dxij*xin(29)
                                      yin(30) = yin(30) + dyij*yin(29)
                                      zin(30) = zin(30) + dzij*zin(29)

                                      ! i3 = i4 =   29
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   32

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(32) = xin(32) + dxij*xin(31)
                                      yin(32) = yin(32) + dyij*yin(31)
                                      zin(32) = zin(32) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   30

                                      xin(31) = xin(31) + dxij*xin(30)
                                      yin(31) = yin(31) + dyij*yin(30)
                                      zin(31) = zin(31) + dzij*zin(30)

                                      ! i3 = i4 =   30
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   32

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(32) = xin(32) + dxij*xin(31)
                                      yin(32) = yin(32) + dyij*yin(31)
                                      zin(32) = zin(32) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   18

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(21) + dxij*xin(17)
                                      yin(18) = yin(21) + dyij*yin(17)
                                      zin(18) = zin(21) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    2

                                      xin(22) = xin(25) + dxij*xin(21)
                                      yin(22) = yin(25) + dyij*yin(21)
                                      zin(22) = zin(25) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   26

                                      ! ni =    3

                                      xin(26) = xin(29) + dxij*xin(25)
                                      yin(26) = yin(29) + dyij*yin(25)
                                      zin(26) = zin(29) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    3

                                      xin(19) = xin(22) + dxij*xin(18)
                                      yin(19) = yin(22) + dyij*yin(18)
                                      zin(19) = zin(22) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    2

                                      xin(23) = xin(26) + dxij*xin(22)
                                      yin(23) = yin(26) + dyij*yin(22)
                                      zin(23) = zin(26) + dzij*zin(22)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   27

                                      ! ni =    3

                                      xin(27) = xin(30) + dxij*xin(26)
                                      yin(27) = yin(30) + dyij*yin(26)
                                      zin(27) = zin(30) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   20

                                      ! nj =    3

                                      ! i4 = i3 =   20

                                      ! do ni = 1,    3

                                      xin(20) = xin(23) + dxij*xin(19)
                                      yin(20) = yin(23) + dyij*yin(19)
                                      zin(20) = zin(23) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    2

                                      xin(24) = xin(27) + dxij*xin(23)
                                      yin(24) = yin(27) + dyij*yin(23)
                                      zin(24) = zin(27) + dzij*zin(23)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   28

                                      ! ni =    3

                                      xin(28) = xin(31) + dxij*xin(27)
                                      yin(28) = yin(31) + dyij*yin(27)
                                      zin(28) = zin(31) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   32

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

                                      ! i1 = in(1) =   33

                                      xin(33) = 1.0_dp
                                      yin(33) = 1.0_dp
                                      zin(33) = f00

                                      ! i2 = in(2) =   37
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(37) = xc00
                                      yin(37) = yc00
                                      zin(37) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   33
                                      ! i4 = i2 =   37

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   41
                                      ! i3 =   33
                                      ! i4 =   37

                                      xin(41) = c10*xin(33) + xc00*xin(37)
                                      yin(41) = c10*yin(33) + yc00*yin(37)
                                      zin(41) = c10*zin(33) + zc00*zin(37)

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   45
                                      ! i3 =   37
                                      ! i4 =   41

                                      xin(45) = c10*xin(37) + xc00*xin(41)
                                      yin(45) = c10*yin(37) + yc00*yin(41)
                                      zin(45) = c10*zin(37) + zc00*zin(41)

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   46
                                      ! i3 =   41
                                      ! i4 =   45

                                      xin(46) = c10*xin(41) + xc00*xin(45)
                                      yin(46) = c10*yin(41) + yc00*yin(45)
                                      zin(46) = c10*zin(41) + zc00*zin(45)

                                      ! i3 = i4 =   45
                                      ! i4 = i5 =   46

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   47
                                      ! i3 =   45
                                      ! i4 =   46

                                      xin(47) = c10*xin(45) + xc00*xin(46)
                                      yin(47) = c10*yin(45) + yc00*yin(46)
                                      zin(47) = c10*zin(45) + zc00*zin(46)

                                      ! i3 = i4 =   46
                                      ! i4 = i5 =   47

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   48
                                      ! i3 =   46
                                      ! i4 =   47

                                      xin(48) = c10*xin(46) + xc00*xin(47)
                                      yin(48) = c10*yin(46) + yc00*yin(47)
                                      zin(48) = c10*zin(46) + zc00*zin(47)

                                      ! i3 = i4 =   47
                                      ! i4 = i5 =   48

                                      ! n =    7

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   48

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   47

                                      xin(48) = xin(48) + dxij*xin(47)
                                      yin(48) = yin(48) + dyij*yin(47)
                                      zin(48) = zin(48) + dzij*zin(47)

                                      ! i3 = i4 =   47
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   46

                                      xin(47) = xin(47) + dxij*xin(46)
                                      yin(47) = yin(47) + dyij*yin(46)
                                      zin(47) = zin(47) + dzij*zin(46)

                                      ! i3 = i4 =   46
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   45

                                      xin(46) = xin(46) + dxij*xin(45)
                                      yin(46) = yin(46) + dyij*yin(45)
                                      zin(46) = zin(46) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   47

                                      xin(48) = xin(48) + dxij*xin(47)
                                      yin(48) = yin(48) + dyij*yin(47)
                                      zin(48) = zin(48) + dzij*zin(47)

                                      ! i3 = i4 =   47
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   46

                                      xin(47) = xin(47) + dxij*xin(46)
                                      yin(47) = yin(47) + dyij*yin(46)
                                      zin(47) = zin(47) + dzij*zin(46)

                                      ! i3 = i4 =   46
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   47

                                      xin(48) = xin(48) + dxij*xin(47)
                                      yin(48) = yin(48) + dyij*yin(47)
                                      zin(48) = zin(48) + dzij*zin(47)

                                      ! i3 = i4 =   47
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   34

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   34

                                      ! do ni = 1,    3

                                      xin(34) = xin(37) + dxij*xin(33)
                                      yin(34) = yin(37) + dyij*yin(33)
                                      zin(34) = zin(37) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   38

                                      ! ni =    2

                                      xin(38) = xin(41) + dxij*xin(37)
                                      yin(38) = yin(41) + dyij*yin(37)
                                      zin(38) = zin(41) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    3

                                      xin(42) = xin(45) + dxij*xin(41)
                                      yin(42) = yin(45) + dyij*yin(41)
                                      zin(42) = zin(45) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   35

                                      ! nj =    2

                                      ! i4 = i3 =   35

                                      ! do ni = 1,    3

                                      xin(35) = xin(38) + dxij*xin(34)
                                      yin(35) = yin(38) + dyij*yin(34)
                                      zin(35) = zin(38) + dzij*zin(34)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    2

                                      xin(39) = xin(42) + dxij*xin(38)
                                      yin(39) = yin(42) + dyij*yin(38)
                                      zin(39) = zin(42) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    3

                                      xin(43) = xin(46) + dxij*xin(42)
                                      yin(43) = yin(46) + dyij*yin(42)
                                      zin(43) = zin(46) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   36

                                      ! nj =    3

                                      ! i4 = i3 =   36

                                      ! do ni = 1,    3

                                      xin(36) = xin(39) + dxij*xin(35)
                                      yin(36) = yin(39) + dyij*yin(35)
                                      zin(36) = zin(39) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    2

                                      xin(40) = xin(43) + dxij*xin(39)
                                      yin(40) = yin(43) + dyij*yin(39)
                                      zin(40) = zin(43) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    3

                                      xin(44) = xin(47) + dxij*xin(43)
                                      yin(44) = yin(47) + dyij*yin(43)
                                      zin(44) = zin(47) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   37

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   48

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

                                      ! i1 = in(1) =   49

                                      xin(49) = 1.0_dp
                                      yin(49) = 1.0_dp
                                      zin(49) = f00

                                      ! i2 = in(2) =   53
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(53) = xc00
                                      yin(53) = yc00
                                      zin(53) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   53

                                      ! do n = 2,   6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   57
                                      ! i3 =   49
                                      ! i4 =   53

                                      xin(57) = c10*xin(49) + xc00*xin(53)
                                      yin(57) = c10*yin(49) + yc00*yin(53)
                                      zin(57) = c10*zin(49) + zc00*zin(53)

                                      ! i3 = i4 =   53
                                      ! i4 = i5 =   57

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   53
                                      ! i4 =   57

                                      xin(61) = c10*xin(53) + xc00*xin(57)
                                      yin(61) = c10*yin(53) + yc00*yin(57)
                                      zin(61) = c10*zin(53) + zc00*zin(57)

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   62
                                      ! i3 =   57
                                      ! i4 =   61

                                      xin(62) = c10*xin(57) + xc00*xin(61)
                                      yin(62) = c10*yin(57) + yc00*yin(61)
                                      zin(62) = c10*zin(57) + zc00*zin(61)

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   62

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   63
                                      ! i3 =   61
                                      ! i4 =   62

                                      xin(63) = c10*xin(61) + xc00*xin(62)
                                      yin(63) = c10*yin(61) + yc00*yin(62)
                                      zin(63) = c10*zin(61) + zc00*zin(62)

                                      ! i3 = i4 =   62
                                      ! i4 = i5 =   63

                                      ! n =    6

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   64
                                      ! i3 =   62
                                      ! i4 =   63

                                      xin(64) = c10*xin(62) + xc00*xin(63)
                                      yin(64) = c10*yin(62) + yc00*yin(63)
                                      zin(64) = c10*zin(62) + zc00*zin(63)

                                      ! i3 = i4 =   63
                                      ! i4 = i5 =   64

                                      ! n =    7

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   64

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(64) = xin(64) + dxij*xin(63)
                                      yin(64) = yin(64) + dyij*yin(63)
                                      zin(64) = zin(64) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   62

                                      xin(63) = xin(63) + dxij*xin(62)
                                      yin(63) = yin(63) + dyij*yin(62)
                                      zin(63) = zin(63) + dzij*zin(62)

                                      ! i3 = i4 =   62
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   61

                                      xin(62) = xin(62) + dxij*xin(61)
                                      yin(62) = yin(62) + dyij*yin(61)
                                      zin(62) = zin(62) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(64) = xin(64) + dxij*xin(63)
                                      yin(64) = yin(64) + dyij*yin(63)
                                      zin(64) = zin(64) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    5

                                      ! i4 = in(nn)+km =   62

                                      xin(63) = xin(63) + dxij*xin(62)
                                      yin(63) = yin(63) + dyij*yin(62)
                                      zin(63) = zin(63) + dzij*zin(62)

                                      ! i3 = i4 =   62
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    6

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(64) = xin(64) + dxij*xin(63)
                                      yin(64) = yin(64) + dyij*yin(63)
                                      zin(64) = zin(64) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    5

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   50

                                      ! do nj = 1,    3

                                      ! i4 = i3 =   50

                                      ! do ni = 1,    3

                                      xin(50) = xin(53) + dxij*xin(49)
                                      yin(50) = yin(53) + dyij*yin(49)
                                      zin(50) = zin(53) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    2

                                      xin(54) = xin(57) + dxij*xin(53)
                                      yin(54) = yin(57) + dyij*yin(53)
                                      zin(54) = zin(57) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    3

                                      xin(58) = xin(61) + dxij*xin(57)
                                      yin(58) = yin(61) + dyij*yin(57)
                                      zin(58) = zin(61) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   62

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   51

                                      ! nj =    2

                                      ! i4 = i3 =   51

                                      ! do ni = 1,    3

                                      xin(51) = xin(54) + dxij*xin(50)
                                      yin(51) = yin(54) + dyij*yin(50)
                                      zin(51) = zin(54) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    2

                                      xin(55) = xin(58) + dxij*xin(54)
                                      yin(55) = yin(58) + dyij*yin(54)
                                      zin(55) = zin(58) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    3

                                      xin(59) = xin(62) + dxij*xin(58)
                                      yin(59) = yin(62) + dyij*yin(58)
                                      zin(59) = zin(62) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   52

                                      ! nj =    3

                                      ! i4 = i3 =   52

                                      ! do ni = 1,    3

                                      xin(52) = xin(55) + dxij*xin(51)
                                      yin(52) = yin(55) + dyij*yin(51)
                                      zin(52) = zin(55) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   56

                                      ! ni =    2

                                      xin(56) = xin(59) + dxij*xin(55)
                                      yin(56) = yin(59) + dyij*yin(55)
                                      zin(56) = zin(59) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    3

                                      xin(60) = xin(63) + dxij*xin(59)
                                      yin(60) = yin(63) + dyij*yin(59)
                                      zin(60) = zin(63) + dzij*zin(59)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   53

                                      ! nj =    4

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   64

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d33bra(  1)*d00ket(  1)*(xin(  16)*yin(   1)*zin(   1)+xin(  32)*yin(  17)*zin(  17)+xin(  48)*yin(  33)*zin(  33)+xin(  64)*yin(  49)*zin(  49))
          eri_value(    2)=eri_value(    2)+d33bra(  2)*d00ket(  1)*(xin(  13)*yin(   4)*zin(   1)+xin(  29)*yin(  20)*zin(  17)+xin(  45)*yin(  36)*zin(  33)+xin(  61)*yin(  52)*zin(  49))
          eri_value(    3)=eri_value(    3)+d33bra(  3)*d00ket(  1)*(xin(  13)*yin(   1)*zin(   4)+xin(  29)*yin(  17)*zin(  20)+xin(  45)*yin(  33)*zin(  36)+xin(  61)*yin(  49)*zin(  52))
          eri_value(    4)=eri_value(    4)+d33bra(  4)*d00ket(  1)*(xin(  15)*yin(   2)*zin(   1)+xin(  31)*yin(  18)*zin(  17)+xin(  47)*yin(  34)*zin(  33)+xin(  63)*yin(  50)*zin(  49))
          eri_value(    5)=eri_value(    5)+d33bra(  5)*d00ket(  1)*(xin(  15)*yin(   1)*zin(   2)+xin(  31)*yin(  17)*zin(  18)+xin(  47)*yin(  33)*zin(  34)+xin(  63)*yin(  49)*zin(  50))
          eri_value(    6)=eri_value(    6)+d33bra(  6)*d00ket(  1)*(xin(  14)*yin(   3)*zin(   1)+xin(  30)*yin(  19)*zin(  17)+xin(  46)*yin(  35)*zin(  33)+xin(  62)*yin(  51)*zin(  49))
          eri_value(    7)=eri_value(    7)+d33bra(  7)*d00ket(  1)*(xin(  13)*yin(   3)*zin(   2)+xin(  29)*yin(  19)*zin(  18)+xin(  45)*yin(  35)*zin(  34)+xin(  61)*yin(  51)*zin(  50))
          eri_value(    8)=eri_value(    8)+d33bra(  8)*d00ket(  1)*(xin(  14)*yin(   1)*zin(   3)+xin(  30)*yin(  17)*zin(  19)+xin(  46)*yin(  33)*zin(  35)+xin(  62)*yin(  49)*zin(  51))
          eri_value(    9)=eri_value(    9)+d33bra(  9)*d00ket(  1)*(xin(  13)*yin(   2)*zin(   3)+xin(  29)*yin(  18)*zin(  19)+xin(  45)*yin(  34)*zin(  35)+xin(  61)*yin(  50)*zin(  51))
          eri_value(   10)=eri_value(   10)+d33bra( 10)*d00ket(  1)*(xin(  14)*yin(   2)*zin(   2)+xin(  30)*yin(  18)*zin(  18)+xin(  46)*yin(  34)*zin(  34)+xin(  62)*yin(  50)*zin(  50))
          eri_value(   11)=eri_value(   11)+d33bra( 11)*d00ket(  1)*(xin(   4)*yin(  13)*zin(   1)+xin(  20)*yin(  29)*zin(  17)+xin(  36)*yin(  45)*zin(  33)+xin(  52)*yin(  61)*zin(  49))
          eri_value(   12)=eri_value(   12)+d33bra( 12)*d00ket(  1)*(xin(   1)*yin(  16)*zin(   1)+xin(  17)*yin(  32)*zin(  17)+xin(  33)*yin(  48)*zin(  33)+xin(  49)*yin(  64)*zin(  49))
          eri_value(   13)=eri_value(   13)+d33bra( 13)*d00ket(  1)*(xin(   1)*yin(  13)*zin(   4)+xin(  17)*yin(  29)*zin(  20)+xin(  33)*yin(  45)*zin(  36)+xin(  49)*yin(  61)*zin(  52))
          eri_value(   14)=eri_value(   14)+d33bra( 14)*d00ket(  1)*(xin(   3)*yin(  14)*zin(   1)+xin(  19)*yin(  30)*zin(  17)+xin(  35)*yin(  46)*zin(  33)+xin(  51)*yin(  62)*zin(  49))
          eri_value(   15)=eri_value(   15)+d33bra( 15)*d00ket(  1)*(xin(   3)*yin(  13)*zin(   2)+xin(  19)*yin(  29)*zin(  18)+xin(  35)*yin(  45)*zin(  34)+xin(  51)*yin(  61)*zin(  50))
          eri_value(   16)=eri_value(   16)+d33bra( 16)*d00ket(  1)*(xin(   2)*yin(  15)*zin(   1)+xin(  18)*yin(  31)*zin(  17)+xin(  34)*yin(  47)*zin(  33)+xin(  50)*yin(  63)*zin(  49))
          eri_value(   17)=eri_value(   17)+d33bra( 17)*d00ket(  1)*(xin(   1)*yin(  15)*zin(   2)+xin(  17)*yin(  31)*zin(  18)+xin(  33)*yin(  47)*zin(  34)+xin(  49)*yin(  63)*zin(  50))
          eri_value(   18)=eri_value(   18)+d33bra( 18)*d00ket(  1)*(xin(   2)*yin(  13)*zin(   3)+xin(  18)*yin(  29)*zin(  19)+xin(  34)*yin(  45)*zin(  35)+xin(  50)*yin(  61)*zin(  51))
          eri_value(   19)=eri_value(   19)+d33bra( 19)*d00ket(  1)*(xin(   1)*yin(  14)*zin(   3)+xin(  17)*yin(  30)*zin(  19)+xin(  33)*yin(  46)*zin(  35)+xin(  49)*yin(  62)*zin(  51))
          eri_value(   20)=eri_value(   20)+d33bra( 20)*d00ket(  1)*(xin(   2)*yin(  14)*zin(   2)+xin(  18)*yin(  30)*zin(  18)+xin(  34)*yin(  46)*zin(  34)+xin(  50)*yin(  62)*zin(  50))
          eri_value(   21)=eri_value(   21)+d33bra( 21)*d00ket(  1)*(xin(   4)*yin(   1)*zin(  13)+xin(  20)*yin(  17)*zin(  29)+xin(  36)*yin(  33)*zin(  45)+xin(  52)*yin(  49)*zin(  61))
          eri_value(   22)=eri_value(   22)+d33bra( 22)*d00ket(  1)*(xin(   1)*yin(   4)*zin(  13)+xin(  17)*yin(  20)*zin(  29)+xin(  33)*yin(  36)*zin(  45)+xin(  49)*yin(  52)*zin(  61))
          eri_value(   23)=eri_value(   23)+d33bra( 23)*d00ket(  1)*(xin(   1)*yin(   1)*zin(  16)+xin(  17)*yin(  17)*zin(  32)+xin(  33)*yin(  33)*zin(  48)+xin(  49)*yin(  49)*zin(  64))
          eri_value(   24)=eri_value(   24)+d33bra( 24)*d00ket(  1)*(xin(   3)*yin(   2)*zin(  13)+xin(  19)*yin(  18)*zin(  29)+xin(  35)*yin(  34)*zin(  45)+xin(  51)*yin(  50)*zin(  61))
          eri_value(   25)=eri_value(   25)+d33bra( 25)*d00ket(  1)*(xin(   3)*yin(   1)*zin(  14)+xin(  19)*yin(  17)*zin(  30)+xin(  35)*yin(  33)*zin(  46)+xin(  51)*yin(  49)*zin(  62))
          eri_value(   26)=eri_value(   26)+d33bra( 26)*d00ket(  1)*(xin(   2)*yin(   3)*zin(  13)+xin(  18)*yin(  19)*zin(  29)+xin(  34)*yin(  35)*zin(  45)+xin(  50)*yin(  51)*zin(  61))
          eri_value(   27)=eri_value(   27)+d33bra( 27)*d00ket(  1)*(xin(   1)*yin(   3)*zin(  14)+xin(  17)*yin(  19)*zin(  30)+xin(  33)*yin(  35)*zin(  46)+xin(  49)*yin(  51)*zin(  62))
          eri_value(   28)=eri_value(   28)+d33bra( 28)*d00ket(  1)*(xin(   2)*yin(   1)*zin(  15)+xin(  18)*yin(  17)*zin(  31)+xin(  34)*yin(  33)*zin(  47)+xin(  50)*yin(  49)*zin(  63))
          eri_value(   29)=eri_value(   29)+d33bra( 29)*d00ket(  1)*(xin(   1)*yin(   2)*zin(  15)+xin(  17)*yin(  18)*zin(  31)+xin(  33)*yin(  34)*zin(  47)+xin(  49)*yin(  50)*zin(  63))
          eri_value(   30)=eri_value(   30)+d33bra( 30)*d00ket(  1)*(xin(   2)*yin(   2)*zin(  14)+xin(  18)*yin(  18)*zin(  30)+xin(  34)*yin(  34)*zin(  46)+xin(  50)*yin(  50)*zin(  62))
          eri_value(   31)=eri_value(   31)+d33bra( 31)*d00ket(  1)*(xin(  12)*yin(   5)*zin(   1)+xin(  28)*yin(  21)*zin(  17)+xin(  44)*yin(  37)*zin(  33)+xin(  60)*yin(  53)*zin(  49))
          eri_value(   32)=eri_value(   32)+d33bra( 32)*d00ket(  1)*(xin(   9)*yin(   8)*zin(   1)+xin(  25)*yin(  24)*zin(  17)+xin(  41)*yin(  40)*zin(  33)+xin(  57)*yin(  56)*zin(  49))
          eri_value(   33)=eri_value(   33)+d33bra( 33)*d00ket(  1)*(xin(   9)*yin(   5)*zin(   4)+xin(  25)*yin(  21)*zin(  20)+xin(  41)*yin(  37)*zin(  36)+xin(  57)*yin(  53)*zin(  52))
          eri_value(   34)=eri_value(   34)+d33bra( 34)*d00ket(  1)*(xin(  11)*yin(   6)*zin(   1)+xin(  27)*yin(  22)*zin(  17)+xin(  43)*yin(  38)*zin(  33)+xin(  59)*yin(  54)*zin(  49))
          eri_value(   35)=eri_value(   35)+d33bra( 35)*d00ket(  1)*(xin(  11)*yin(   5)*zin(   2)+xin(  27)*yin(  21)*zin(  18)+xin(  43)*yin(  37)*zin(  34)+xin(  59)*yin(  53)*zin(  50))
          eri_value(   36)=eri_value(   36)+d33bra( 36)*d00ket(  1)*(xin(  10)*yin(   7)*zin(   1)+xin(  26)*yin(  23)*zin(  17)+xin(  42)*yin(  39)*zin(  33)+xin(  58)*yin(  55)*zin(  49))
          eri_value(   37)=eri_value(   37)+d33bra( 37)*d00ket(  1)*(xin(   9)*yin(   7)*zin(   2)+xin(  25)*yin(  23)*zin(  18)+xin(  41)*yin(  39)*zin(  34)+xin(  57)*yin(  55)*zin(  50))
          eri_value(   38)=eri_value(   38)+d33bra( 38)*d00ket(  1)*(xin(  10)*yin(   5)*zin(   3)+xin(  26)*yin(  21)*zin(  19)+xin(  42)*yin(  37)*zin(  35)+xin(  58)*yin(  53)*zin(  51))
          eri_value(   39)=eri_value(   39)+d33bra( 39)*d00ket(  1)*(xin(   9)*yin(   6)*zin(   3)+xin(  25)*yin(  22)*zin(  19)+xin(  41)*yin(  38)*zin(  35)+xin(  57)*yin(  54)*zin(  51))
          eri_value(   40)=eri_value(   40)+d33bra( 40)*d00ket(  1)*(xin(  10)*yin(   6)*zin(   2)+xin(  26)*yin(  22)*zin(  18)+xin(  42)*yin(  38)*zin(  34)+xin(  58)*yin(  54)*zin(  50))
          eri_value(   41)=eri_value(   41)+d33bra( 41)*d00ket(  1)*(xin(  12)*yin(   1)*zin(   5)+xin(  28)*yin(  17)*zin(  21)+xin(  44)*yin(  33)*zin(  37)+xin(  60)*yin(  49)*zin(  53))
          eri_value(   42)=eri_value(   42)+d33bra( 42)*d00ket(  1)*(xin(   9)*yin(   4)*zin(   5)+xin(  25)*yin(  20)*zin(  21)+xin(  41)*yin(  36)*zin(  37)+xin(  57)*yin(  52)*zin(  53))
          eri_value(   43)=eri_value(   43)+d33bra( 43)*d00ket(  1)*(xin(   9)*yin(   1)*zin(   8)+xin(  25)*yin(  17)*zin(  24)+xin(  41)*yin(  33)*zin(  40)+xin(  57)*yin(  49)*zin(  56))
          eri_value(   44)=eri_value(   44)+d33bra( 44)*d00ket(  1)*(xin(  11)*yin(   2)*zin(   5)+xin(  27)*yin(  18)*zin(  21)+xin(  43)*yin(  34)*zin(  37)+xin(  59)*yin(  50)*zin(  53))
          eri_value(   45)=eri_value(   45)+d33bra( 45)*d00ket(  1)*(xin(  11)*yin(   1)*zin(   6)+xin(  27)*yin(  17)*zin(  22)+xin(  43)*yin(  33)*zin(  38)+xin(  59)*yin(  49)*zin(  54))
          eri_value(   46)=eri_value(   46)+d33bra( 46)*d00ket(  1)*(xin(  10)*yin(   3)*zin(   5)+xin(  26)*yin(  19)*zin(  21)+xin(  42)*yin(  35)*zin(  37)+xin(  58)*yin(  51)*zin(  53))
          eri_value(   47)=eri_value(   47)+d33bra( 47)*d00ket(  1)*(xin(   9)*yin(   3)*zin(   6)+xin(  25)*yin(  19)*zin(  22)+xin(  41)*yin(  35)*zin(  38)+xin(  57)*yin(  51)*zin(  54))
          eri_value(   48)=eri_value(   48)+d33bra( 48)*d00ket(  1)*(xin(  10)*yin(   1)*zin(   7)+xin(  26)*yin(  17)*zin(  23)+xin(  42)*yin(  33)*zin(  39)+xin(  58)*yin(  49)*zin(  55))
          eri_value(   49)=eri_value(   49)+d33bra( 49)*d00ket(  1)*(xin(   9)*yin(   2)*zin(   7)+xin(  25)*yin(  18)*zin(  23)+xin(  41)*yin(  34)*zin(  39)+xin(  57)*yin(  50)*zin(  55))
          eri_value(   50)=eri_value(   50)+d33bra( 50)*d00ket(  1)*(xin(  10)*yin(   2)*zin(   6)+xin(  26)*yin(  18)*zin(  22)+xin(  42)*yin(  34)*zin(  38)+xin(  58)*yin(  50)*zin(  54))
          eri_value(   51)=eri_value(   51)+d33bra( 51)*d00ket(  1)*(xin(   8)*yin(   9)*zin(   1)+xin(  24)*yin(  25)*zin(  17)+xin(  40)*yin(  41)*zin(  33)+xin(  56)*yin(  57)*zin(  49))
          eri_value(   52)=eri_value(   52)+d33bra( 52)*d00ket(  1)*(xin(   5)*yin(  12)*zin(   1)+xin(  21)*yin(  28)*zin(  17)+xin(  37)*yin(  44)*zin(  33)+xin(  53)*yin(  60)*zin(  49))
          eri_value(   53)=eri_value(   53)+d33bra( 53)*d00ket(  1)*(xin(   5)*yin(   9)*zin(   4)+xin(  21)*yin(  25)*zin(  20)+xin(  37)*yin(  41)*zin(  36)+xin(  53)*yin(  57)*zin(  52))
          eri_value(   54)=eri_value(   54)+d33bra( 54)*d00ket(  1)*(xin(   7)*yin(  10)*zin(   1)+xin(  23)*yin(  26)*zin(  17)+xin(  39)*yin(  42)*zin(  33)+xin(  55)*yin(  58)*zin(  49))
          eri_value(   55)=eri_value(   55)+d33bra( 55)*d00ket(  1)*(xin(   7)*yin(   9)*zin(   2)+xin(  23)*yin(  25)*zin(  18)+xin(  39)*yin(  41)*zin(  34)+xin(  55)*yin(  57)*zin(  50))
          eri_value(   56)=eri_value(   56)+d33bra( 56)*d00ket(  1)*(xin(   6)*yin(  11)*zin(   1)+xin(  22)*yin(  27)*zin(  17)+xin(  38)*yin(  43)*zin(  33)+xin(  54)*yin(  59)*zin(  49))
          eri_value(   57)=eri_value(   57)+d33bra( 57)*d00ket(  1)*(xin(   5)*yin(  11)*zin(   2)+xin(  21)*yin(  27)*zin(  18)+xin(  37)*yin(  43)*zin(  34)+xin(  53)*yin(  59)*zin(  50))
          eri_value(   58)=eri_value(   58)+d33bra( 58)*d00ket(  1)*(xin(   6)*yin(   9)*zin(   3)+xin(  22)*yin(  25)*zin(  19)+xin(  38)*yin(  41)*zin(  35)+xin(  54)*yin(  57)*zin(  51))
          eri_value(   59)=eri_value(   59)+d33bra( 59)*d00ket(  1)*(xin(   5)*yin(  10)*zin(   3)+xin(  21)*yin(  26)*zin(  19)+xin(  37)*yin(  42)*zin(  35)+xin(  53)*yin(  58)*zin(  51))
          eri_value(   60)=eri_value(   60)+d33bra( 60)*d00ket(  1)*(xin(   6)*yin(  10)*zin(   2)+xin(  22)*yin(  26)*zin(  18)+xin(  38)*yin(  42)*zin(  34)+xin(  54)*yin(  58)*zin(  50))
          eri_value(   61)=eri_value(   61)+d33bra( 61)*d00ket(  1)*(xin(   4)*yin(   9)*zin(   5)+xin(  20)*yin(  25)*zin(  21)+xin(  36)*yin(  41)*zin(  37)+xin(  52)*yin(  57)*zin(  53))
          eri_value(   62)=eri_value(   62)+d33bra( 62)*d00ket(  1)*(xin(   1)*yin(  12)*zin(   5)+xin(  17)*yin(  28)*zin(  21)+xin(  33)*yin(  44)*zin(  37)+xin(  49)*yin(  60)*zin(  53))
          eri_value(   63)=eri_value(   63)+d33bra( 63)*d00ket(  1)*(xin(   1)*yin(   9)*zin(   8)+xin(  17)*yin(  25)*zin(  24)+xin(  33)*yin(  41)*zin(  40)+xin(  49)*yin(  57)*zin(  56))
          eri_value(   64)=eri_value(   64)+d33bra( 64)*d00ket(  1)*(xin(   3)*yin(  10)*zin(   5)+xin(  19)*yin(  26)*zin(  21)+xin(  35)*yin(  42)*zin(  37)+xin(  51)*yin(  58)*zin(  53))
          eri_value(   65)=eri_value(   65)+d33bra( 65)*d00ket(  1)*(xin(   3)*yin(   9)*zin(   6)+xin(  19)*yin(  25)*zin(  22)+xin(  35)*yin(  41)*zin(  38)+xin(  51)*yin(  57)*zin(  54))
          eri_value(   66)=eri_value(   66)+d33bra( 66)*d00ket(  1)*(xin(   2)*yin(  11)*zin(   5)+xin(  18)*yin(  27)*zin(  21)+xin(  34)*yin(  43)*zin(  37)+xin(  50)*yin(  59)*zin(  53))
          eri_value(   67)=eri_value(   67)+d33bra( 67)*d00ket(  1)*(xin(   1)*yin(  11)*zin(   6)+xin(  17)*yin(  27)*zin(  22)+xin(  33)*yin(  43)*zin(  38)+xin(  49)*yin(  59)*zin(  54))
          eri_value(   68)=eri_value(   68)+d33bra( 68)*d00ket(  1)*(xin(   2)*yin(   9)*zin(   7)+xin(  18)*yin(  25)*zin(  23)+xin(  34)*yin(  41)*zin(  39)+xin(  50)*yin(  57)*zin(  55))
          eri_value(   69)=eri_value(   69)+d33bra( 69)*d00ket(  1)*(xin(   1)*yin(  10)*zin(   7)+xin(  17)*yin(  26)*zin(  23)+xin(  33)*yin(  42)*zin(  39)+xin(  49)*yin(  58)*zin(  55))
          eri_value(   70)=eri_value(   70)+d33bra( 70)*d00ket(  1)*(xin(   2)*yin(  10)*zin(   6)+xin(  18)*yin(  26)*zin(  22)+xin(  34)*yin(  42)*zin(  38)+xin(  50)*yin(  58)*zin(  54))
          eri_value(   71)=eri_value(   71)+d33bra( 71)*d00ket(  1)*(xin(   8)*yin(   1)*zin(   9)+xin(  24)*yin(  17)*zin(  25)+xin(  40)*yin(  33)*zin(  41)+xin(  56)*yin(  49)*zin(  57))
          eri_value(   72)=eri_value(   72)+d33bra( 72)*d00ket(  1)*(xin(   5)*yin(   4)*zin(   9)+xin(  21)*yin(  20)*zin(  25)+xin(  37)*yin(  36)*zin(  41)+xin(  53)*yin(  52)*zin(  57))
          eri_value(   73)=eri_value(   73)+d33bra( 73)*d00ket(  1)*(xin(   5)*yin(   1)*zin(  12)+xin(  21)*yin(  17)*zin(  28)+xin(  37)*yin(  33)*zin(  44)+xin(  53)*yin(  49)*zin(  60))
          eri_value(   74)=eri_value(   74)+d33bra( 74)*d00ket(  1)*(xin(   7)*yin(   2)*zin(   9)+xin(  23)*yin(  18)*zin(  25)+xin(  39)*yin(  34)*zin(  41)+xin(  55)*yin(  50)*zin(  57))
          eri_value(   75)=eri_value(   75)+d33bra( 75)*d00ket(  1)*(xin(   7)*yin(   1)*zin(  10)+xin(  23)*yin(  17)*zin(  26)+xin(  39)*yin(  33)*zin(  42)+xin(  55)*yin(  49)*zin(  58))
          eri_value(   76)=eri_value(   76)+d33bra( 76)*d00ket(  1)*(xin(   6)*yin(   3)*zin(   9)+xin(  22)*yin(  19)*zin(  25)+xin(  38)*yin(  35)*zin(  41)+xin(  54)*yin(  51)*zin(  57))
          eri_value(   77)=eri_value(   77)+d33bra( 77)*d00ket(  1)*(xin(   5)*yin(   3)*zin(  10)+xin(  21)*yin(  19)*zin(  26)+xin(  37)*yin(  35)*zin(  42)+xin(  53)*yin(  51)*zin(  58))
          eri_value(   78)=eri_value(   78)+d33bra( 78)*d00ket(  1)*(xin(   6)*yin(   1)*zin(  11)+xin(  22)*yin(  17)*zin(  27)+xin(  38)*yin(  33)*zin(  43)+xin(  54)*yin(  49)*zin(  59))
          eri_value(   79)=eri_value(   79)+d33bra( 79)*d00ket(  1)*(xin(   5)*yin(   2)*zin(  11)+xin(  21)*yin(  18)*zin(  27)+xin(  37)*yin(  34)*zin(  43)+xin(  53)*yin(  50)*zin(  59))
          eri_value(   80)=eri_value(   80)+d33bra( 80)*d00ket(  1)*(xin(   6)*yin(   2)*zin(  10)+xin(  22)*yin(  18)*zin(  26)+xin(  38)*yin(  34)*zin(  42)+xin(  54)*yin(  50)*zin(  58))
          eri_value(   81)=eri_value(   81)+d33bra( 81)*d00ket(  1)*(xin(   4)*yin(   5)*zin(   9)+xin(  20)*yin(  21)*zin(  25)+xin(  36)*yin(  37)*zin(  41)+xin(  52)*yin(  53)*zin(  57))
          eri_value(   82)=eri_value(   82)+d33bra( 82)*d00ket(  1)*(xin(   1)*yin(   8)*zin(   9)+xin(  17)*yin(  24)*zin(  25)+xin(  33)*yin(  40)*zin(  41)+xin(  49)*yin(  56)*zin(  57))
          eri_value(   83)=eri_value(   83)+d33bra( 83)*d00ket(  1)*(xin(   1)*yin(   5)*zin(  12)+xin(  17)*yin(  21)*zin(  28)+xin(  33)*yin(  37)*zin(  44)+xin(  49)*yin(  53)*zin(  60))
          eri_value(   84)=eri_value(   84)+d33bra( 84)*d00ket(  1)*(xin(   3)*yin(   6)*zin(   9)+xin(  19)*yin(  22)*zin(  25)+xin(  35)*yin(  38)*zin(  41)+xin(  51)*yin(  54)*zin(  57))
          eri_value(   85)=eri_value(   85)+d33bra( 85)*d00ket(  1)*(xin(   3)*yin(   5)*zin(  10)+xin(  19)*yin(  21)*zin(  26)+xin(  35)*yin(  37)*zin(  42)+xin(  51)*yin(  53)*zin(  58))
          eri_value(   86)=eri_value(   86)+d33bra( 86)*d00ket(  1)*(xin(   2)*yin(   7)*zin(   9)+xin(  18)*yin(  23)*zin(  25)+xin(  34)*yin(  39)*zin(  41)+xin(  50)*yin(  55)*zin(  57))
          eri_value(   87)=eri_value(   87)+d33bra( 87)*d00ket(  1)*(xin(   1)*yin(   7)*zin(  10)+xin(  17)*yin(  23)*zin(  26)+xin(  33)*yin(  39)*zin(  42)+xin(  49)*yin(  55)*zin(  58))
          eri_value(   88)=eri_value(   88)+d33bra( 88)*d00ket(  1)*(xin(   2)*yin(   5)*zin(  11)+xin(  18)*yin(  21)*zin(  27)+xin(  34)*yin(  37)*zin(  43)+xin(  50)*yin(  53)*zin(  59))
          eri_value(   89)=eri_value(   89)+d33bra( 89)*d00ket(  1)*(xin(   1)*yin(   6)*zin(  11)+xin(  17)*yin(  22)*zin(  27)+xin(  33)*yin(  38)*zin(  43)+xin(  49)*yin(  54)*zin(  59))
          eri_value(   90)=eri_value(   90)+d33bra( 90)*d00ket(  1)*(xin(   2)*yin(   6)*zin(  10)+xin(  18)*yin(  22)*zin(  26)+xin(  34)*yin(  38)*zin(  42)+xin(  50)*yin(  54)*zin(  58))
          eri_value(   91)=eri_value(   91)+d33bra( 91)*d00ket(  1)*(xin(   8)*yin(   5)*zin(   5)+xin(  24)*yin(  21)*zin(  21)+xin(  40)*yin(  37)*zin(  37)+xin(  56)*yin(  53)*zin(  53))
          eri_value(   92)=eri_value(   92)+d33bra( 92)*d00ket(  1)*(xin(   5)*yin(   8)*zin(   5)+xin(  21)*yin(  24)*zin(  21)+xin(  37)*yin(  40)*zin(  37)+xin(  53)*yin(  56)*zin(  53))
          eri_value(   93)=eri_value(   93)+d33bra( 93)*d00ket(  1)*(xin(   5)*yin(   5)*zin(   8)+xin(  21)*yin(  21)*zin(  24)+xin(  37)*yin(  37)*zin(  40)+xin(  53)*yin(  53)*zin(  56))
          eri_value(   94)=eri_value(   94)+d33bra( 94)*d00ket(  1)*(xin(   7)*yin(   6)*zin(   5)+xin(  23)*yin(  22)*zin(  21)+xin(  39)*yin(  38)*zin(  37)+xin(  55)*yin(  54)*zin(  53))
          eri_value(   95)=eri_value(   95)+d33bra( 95)*d00ket(  1)*(xin(   7)*yin(   5)*zin(   6)+xin(  23)*yin(  21)*zin(  22)+xin(  39)*yin(  37)*zin(  38)+xin(  55)*yin(  53)*zin(  54))
          eri_value(   96)=eri_value(   96)+d33bra( 96)*d00ket(  1)*(xin(   6)*yin(   7)*zin(   5)+xin(  22)*yin(  23)*zin(  21)+xin(  38)*yin(  39)*zin(  37)+xin(  54)*yin(  55)*zin(  53))
          eri_value(   97)=eri_value(   97)+d33bra( 97)*d00ket(  1)*(xin(   5)*yin(   7)*zin(   6)+xin(  21)*yin(  23)*zin(  22)+xin(  37)*yin(  39)*zin(  38)+xin(  53)*yin(  55)*zin(  54))
          eri_value(   98)=eri_value(   98)+d33bra( 98)*d00ket(  1)*(xin(   6)*yin(   5)*zin(   7)+xin(  22)*yin(  21)*zin(  23)+xin(  38)*yin(  37)*zin(  39)+xin(  54)*yin(  53)*zin(  55))
          eri_value(   99)=eri_value(   99)+d33bra( 99)*d00ket(  1)*(xin(   5)*yin(   6)*zin(   7)+xin(  21)*yin(  22)*zin(  23)+xin(  37)*yin(  38)*zin(  39)+xin(  53)*yin(  54)*zin(  55))
          eri_value(  100)=eri_value(  100)+d33bra(100)*d00ket(  1)*(xin(   6)*yin(   6)*zin(   6)+xin(  22)*yin(  22)*zin(  22)+xin(  38)*yin(  38)*zin(  38)+xin(  54)*yin(  54)*zin(  54))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxj2 = 10
                                  iandj = ish .eq. jsh
                                  maxl = 1
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    if (iandj) maxj2 = i

                                    ii1 = i + loci
                                    ip = (i - 1)*10 ! Stride between functions in i

                                    do j = 1, maxj2

                                      maxl2 = maxl

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*1 + ip ! Add stride between functions in j

                                      do k = 1, 1 ! # of cartesians in k

                                        if (kandl) maxl2 = k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*1 + ijp ! Add stride between functions in k

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
                              deallocate (n00ket)
                              deallocate (xint00ket)

                              end subroutine int3300
                              end submodule
