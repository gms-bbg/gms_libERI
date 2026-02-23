! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3210_impl
contains
  module subroutine int3210(df_pair, sp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, sp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n01ket(:)
    real(dp), allocatable :: xint23bra(:), xint01ket(:)
    integer(kind=int64) :: ndfbra, nspket
    real(dp) :: scutdfbra, scutspket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2
    integer(kind=int64) :: n, i1, i3, i4, i5, nm, nn, km, nj, ni, nl, nk
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
    real(dp) :: xin(96), yin(96), zin(96)
    real(dp) :: eri_value(180)
    real(dp) :: d23bra(60), d01ket(3)
    integer(kind=int64) :: ix(10), jx(6), kx(3), lx(1)
    integer(kind=int64) :: iy(10), jy(6), ky(3), ly(1)
    integer(kind=int64) :: iz(10), jz(6), kz(3), lz(1)
    integer(kind=int64) :: in(6), in1(6), kn(2)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(3), kly(3), klz(3)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 7
    in1(3) = 13
    in1(4) = 19
    in1(5) = 21
    in1(6) = 23

    kn(1) = 0
    kn(2) = 1

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 1
    kx(2) = 0
    kx(3) = 0

    jx(1) = 4
    jx(2) = 0
    jx(3) = 0
    jx(4) = 2
    jx(5) = 2
    jx(6) = 0

    ix(1) = 19
    ix(2) = 1
    ix(3) = 1
    ix(4) = 13
    ix(5) = 13
    ix(6) = 7
    ix(7) = 1
    ix(8) = 7
    ix(9) = 1
    ix(10) = 7

    ! y-arrays

    ly(1) = 0

    ky(1) = 0
    ky(2) = 1
    ky(3) = 0

    jy(1) = 0
    jy(2) = 4
    jy(3) = 0
    jy(4) = 2
    jy(5) = 0
    jy(6) = 2

    iy(1) = 1
    iy(2) = 19
    iy(3) = 1
    iy(4) = 7
    iy(5) = 1
    iy(6) = 13
    iy(7) = 13
    iy(8) = 1
    iy(9) = 7
    iy(10) = 7

    ! z-arrays

    lz(1) = 0

    kz(1) = 0
    kz(2) = 0
    kz(3) = 1

    jz(1) = 0
    jz(2) = 0
    jz(3) = 4
    jz(4) = 0
    jz(5) = 2
    jz(6) = 2

    iz(1) = 1
    iz(2) = 1
    iz(3) = 19
    iz(4) = 1
    iz(5) = 7
    iz(6) = 1
    iz(7) = 7
    iz(8) = 13
    iz(9) = 13
    iz(10) = 7

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 23
    ijx(2) = 19
    ijx(3) = 19
    ijx(4) = 21
    ijx(5) = 21
    ijx(6) = 19
    ijx(7) = 5
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 3
    ijx(11) = 3
    ijx(12) = 1
    ijx(13) = 5
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 3
    ijx(17) = 3
    ijx(18) = 1
    ijx(19) = 17
    ijx(20) = 13
    ijx(21) = 13
    ijx(22) = 15
    ijx(23) = 15
    ijx(24) = 13
    ijx(25) = 17
    ijx(26) = 13
    ijx(27) = 13
    ijx(28) = 15
    ijx(29) = 15
    ijx(30) = 13
    ijx(31) = 11
    ijx(32) = 7
    ijx(33) = 7
    ijx(34) = 9
    ijx(35) = 9
    ijx(36) = 7
    ijx(37) = 5
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 3
    ijx(41) = 3
    ijx(42) = 1
    ijx(43) = 11
    ijx(44) = 7
    ijx(45) = 7
    ijx(46) = 9
    ijx(47) = 9
    ijx(48) = 7
    ijx(49) = 5
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 3
    ijx(53) = 3
    ijx(54) = 1
    ijx(55) = 11
    ijx(56) = 7
    ijx(57) = 7
    ijx(58) = 9
    ijx(59) = 9
    ijx(60) = 7

    ijy(1) = 1
    ijy(2) = 5
    ijy(3) = 1
    ijy(4) = 3
    ijy(5) = 1
    ijy(6) = 3
    ijy(7) = 19
    ijy(8) = 23
    ijy(9) = 19
    ijy(10) = 21
    ijy(11) = 19
    ijy(12) = 21
    ijy(13) = 1
    ijy(14) = 5
    ijy(15) = 1
    ijy(16) = 3
    ijy(17) = 1
    ijy(18) = 3
    ijy(19) = 7
    ijy(20) = 11
    ijy(21) = 7
    ijy(22) = 9
    ijy(23) = 7
    ijy(24) = 9
    ijy(25) = 1
    ijy(26) = 5
    ijy(27) = 1
    ijy(28) = 3
    ijy(29) = 1
    ijy(30) = 3
    ijy(31) = 13
    ijy(32) = 17
    ijy(33) = 13
    ijy(34) = 15
    ijy(35) = 13
    ijy(36) = 15
    ijy(37) = 13
    ijy(38) = 17
    ijy(39) = 13
    ijy(40) = 15
    ijy(41) = 13
    ijy(42) = 15
    ijy(43) = 1
    ijy(44) = 5
    ijy(45) = 1
    ijy(46) = 3
    ijy(47) = 1
    ijy(48) = 3
    ijy(49) = 7
    ijy(50) = 11
    ijy(51) = 7
    ijy(52) = 9
    ijy(53) = 7
    ijy(54) = 9
    ijy(55) = 7
    ijy(56) = 11
    ijy(57) = 7
    ijy(58) = 9
    ijy(59) = 7
    ijy(60) = 9

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 5
    ijz(4) = 1
    ijz(5) = 3
    ijz(6) = 3
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 5
    ijz(10) = 1
    ijz(11) = 3
    ijz(12) = 3
    ijz(13) = 19
    ijz(14) = 19
    ijz(15) = 23
    ijz(16) = 19
    ijz(17) = 21
    ijz(18) = 21
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 5
    ijz(22) = 1
    ijz(23) = 3
    ijz(24) = 3
    ijz(25) = 7
    ijz(26) = 7
    ijz(27) = 11
    ijz(28) = 7
    ijz(29) = 9
    ijz(30) = 9
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 5
    ijz(34) = 1
    ijz(35) = 3
    ijz(36) = 3
    ijz(37) = 7
    ijz(38) = 7
    ijz(39) = 11
    ijz(40) = 7
    ijz(41) = 9
    ijz(42) = 9
    ijz(43) = 13
    ijz(44) = 13
    ijz(45) = 17
    ijz(46) = 13
    ijz(47) = 15
    ijz(48) = 15
    ijz(49) = 13
    ijz(50) = 13
    ijz(51) = 17
    ijz(52) = 13
    ijz(53) = 15
    ijz(54) = 15
    ijz(55) = 7
    ijz(56) = 7
    ijz(57) = 11
    ijz(58) = 7
    ijz(59) = 9
    ijz(60) = 9

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 1
    klx(2) = 0
    klx(3) = 0

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n01ket(res%n_s_shl*res%n_p_shl))
    allocate (xint01ket(res%n_s_shl*res%n_p_shl))

    ! Start screening

    scutdfbra = cutoff_schwarz/maxval(df_pair%xints)
    ndfbra = 0
    do ij = 1, res%n_d_shl*res%n_f_shl
      if (df_pair%xints(ij) .ge. scutdfbra) then
        ndfbra = ndfbra + 1
        xint23bra(ndfbra) = df_pair%xints(ij)
        n23bra(ndfbra) = ij
      end if
    end do

    scutspket = cutoff_schwarz/maxval(sp_pair%xints)
    nspket = 0
    do ij = 1, res%n_s_shl*res%n_p_shl
      if (sp_pair%xints(ij) .ge. scutspket) then
        nspket = nspket + 1
        xint01ket(nspket) = sp_pair%xints(ij)
        n01ket(nspket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((ndfbra*nspket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nspket
    ntile = int(ndfbra*nspket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nspket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint01ket, n01ket, df_pair, sp_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d01ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,nm,nn,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint01ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n01ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_p_shl) + 1
                lsh_tmp = (kl - 1)/res%n_p_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
                ksh = res%i_p_shl(ksh_tmp)
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

                  t_expon_cd = sp_pair%t_expon_ab(sp_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = sp_pair%expon_b(sp_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = sp_pair%expon_a(sp_pair%pair_loc(kl) + ket_loop)
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

                  d01ket(1) = sp_pair%d_coeff_alt(sp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d01ket(2) = sp_pair%d_coeff_alt(sp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d01ket(3) = sp_pair%d_coeff_alt(sp_pair%pair_loc(kl) + ket_loop)*twopi_5_2

                  bra_loop = 0

                  do i = 1, ijtop

                    bra_loop = bra_loop + 1

                    t_expon_ab = df_pair%t_expon_ab(df_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = df_pair%expon_b(df_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = df_pair%expon_a(df_pair%pair_loc(ij) + bra_loop)
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

                    d23bra(1) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(2) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(3) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(4) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(5) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(6) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(7) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(8) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(9) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(10) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(11) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(12) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(13) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(14) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(15) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(16) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(17) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(18) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(19) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(20) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(21) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(22) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(23) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(24) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(25) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(26) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(27) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(28) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(29) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(30) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(31) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(32) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(33) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(34) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(35) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(36) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(37) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(38) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(39) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(40) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(41) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(42) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(43) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(44) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(45) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(46) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(47) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(48) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(49) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(50) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(51) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(52) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(53) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(54) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(55) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(56) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(57) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(58) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3
                    d23bra(59) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3
                    d23bra(60) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3

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

                                      ! i2 = in(2) =    7
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(7) = xc00
                                      yin(7) = yc00
                                      zin(7) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =    8
                                      ! i2 =    7

                                      xin(8) = xcp00*xin(7) + cp10
                                      yin(8) = ycp00*yin(7) + cp10
                                      zin(8) = zcp00*zin(7) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    7

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   13
                                      ! i3 =    1
                                      ! i4 =    7

                                      xin(13) = c10*xin(1) + xc00*xin(7)
                                      yin(13) = c10*yin(1) + yc00*yin(7)
                                      zin(13) = c10*zin(1) + zc00*zin(7)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   14
                                      ! i5 =   13
                                      ! i4 =    7

                                      xin(14) = xcp00*xin(13) + cp10*xin(7)
                                      yin(14) = ycp00*yin(13) + cp10*yin(7)
                                      zin(14) = zcp00*zin(13) + cp10*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =   13

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   19
                                      ! i3 =    7
                                      ! i4 =   13

                                      xin(19) = c10*xin(7) + xc00*xin(13)
                                      yin(19) = c10*yin(7) + yc00*yin(13)
                                      zin(19) = c10*zin(7) + zc00*zin(13)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   20
                                      ! i5 =   19
                                      ! i4 =   13

                                      xin(20) = xcp00*xin(19) + cp10*xin(13)
                                      yin(20) = ycp00*yin(19) + cp10*yin(13)
                                      zin(20) = zcp00*zin(19) + cp10*zin(13)

                                      ! ------------------

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   19

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   21
                                      ! i3 =   13
                                      ! i4 =   19

                                      xin(21) = c10*xin(13) + xc00*xin(19)
                                      yin(21) = c10*yin(13) + yc00*yin(19)
                                      zin(21) = c10*zin(13) + zc00*zin(19)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   22
                                      ! i5 =   21
                                      ! i4 =   19

                                      xin(22) = xcp00*xin(21) + cp10*xin(19)
                                      yin(22) = ycp00*yin(21) + cp10*yin(19)
                                      zin(22) = zcp00*zin(21) + cp10*zin(19)

                                      ! ------------------

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   21

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   23
                                      ! i3 =   19
                                      ! i4 =   21

                                      xin(23) = c10*xin(19) + xc00*xin(21)
                                      yin(23) = c10*yin(19) + yc00*yin(21)
                                      zin(23) = c10*zin(19) + zc00*zin(21)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   24
                                      ! i5 =   23
                                      ! i4 =   21

                                      xin(24) = xcp00*xin(23) + cp10*xin(21)
                                      yin(24) = ycp00*yin(23) + cp10*yin(21)
                                      zin(24) = zcp00*zin(23) + cp10*zin(21)

                                      ! ------------------

                                      ! i3 = i4 =   21
                                      ! i4 = i5 =   23

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   23

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   23

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   21

                                      xin(23) = xin(23) + dxij*xin(21)
                                      yin(23) = yin(23) + dyij*yin(21)
                                      zin(23) = zin(23) + dzij*zin(21)

                                      ! i3 = i4 =   21
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   19

                                      xin(21) = xin(21) + dxij*xin(19)
                                      yin(21) = yin(21) + dyij*yin(19)
                                      zin(21) = zin(21) + dzij*zin(19)

                                      ! i3 = i4 =   19
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   23

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   21

                                      xin(23) = xin(23) + dxij*xin(21)
                                      yin(23) = yin(23) + dyij*yin(21)
                                      zin(23) = zin(23) + dzij*zin(21)

                                      ! i3 = i4 =   21
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    3

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    3

                                      ! do ni = 1,    3

                                      xin(3) = xin(7) + dxij*xin(1)
                                      yin(3) = yin(7) + dyij*yin(1)
                                      zin(3) = zin(7) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    9

                                      ! ni =    2

                                      xin(9) = xin(13) + dxij*xin(7)
                                      yin(9) = yin(13) + dyij*yin(7)
                                      zin(9) = zin(13) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   15

                                      ! ni =    3

                                      xin(15) = xin(19) + dxij*xin(13)
                                      yin(15) = yin(19) + dyij*yin(13)
                                      zin(15) = zin(19) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    5

                                      ! nj =    2

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(9) + dxij*xin(3)
                                      yin(5) = yin(9) + dyij*yin(3)
                                      zin(5) = zin(9) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   11

                                      ! ni =    2

                                      xin(11) = xin(15) + dxij*xin(9)
                                      yin(11) = yin(15) + dyij*yin(9)
                                      zin(11) = zin(15) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    3

                                      xin(17) = xin(21) + dxij*xin(15)
                                      yin(17) = yin(21) + dyij*yin(15)
                                      zin(17) = zin(21) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    7

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   22

                                      xin(24) = xin(24) + dxij*xin(22)
                                      yin(24) = yin(24) + dyij*yin(22)
                                      zin(24) = zin(24) + dzij*zin(22)

                                      ! i3 = i4 =   22
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   20

                                      xin(22) = xin(22) + dxij*xin(20)
                                      yin(22) = yin(22) + dyij*yin(20)
                                      zin(22) = zin(22) + dzij*zin(20)

                                      ! i3 = i4 =   20
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   22

                                      xin(24) = xin(24) + dxij*xin(22)
                                      yin(24) = yin(24) + dyij*yin(22)
                                      zin(24) = zin(24) + dzij*zin(22)

                                      ! i3 = i4 =   22
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    4

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    4

                                      ! do ni = 1,    3

                                      xin(4) = xin(8) + dxij*xin(2)
                                      yin(4) = yin(8) + dyij*yin(2)
                                      zin(4) = zin(8) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   10

                                      ! ni =    2

                                      xin(10) = xin(14) + dxij*xin(8)
                                      yin(10) = yin(14) + dyij*yin(8)
                                      zin(10) = zin(14) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    3

                                      xin(16) = xin(20) + dxij*xin(14)
                                      yin(16) = yin(20) + dyij*yin(14)
                                      zin(16) = zin(20) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    6

                                      ! nj =    2

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    3

                                      xin(6) = xin(10) + dxij*xin(4)
                                      yin(6) = yin(10) + dyij*yin(4)
                                      zin(6) = zin(10) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   12

                                      ! ni =    2

                                      xin(12) = xin(16) + dxij*xin(10)
                                      yin(12) = yin(16) + dyij*yin(10)
                                      zin(12) = zin(16) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    3

                                      xin(18) = xin(22) + dxij*xin(16)
                                      yin(18) = yin(22) + dyij*yin(16)
                                      zin(18) = zin(22) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    8

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   24

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

                                      ! i1 = in(1) =   25

                                      xin(25) = 1.0_dp
                                      yin(25) = 1.0_dp
                                      zin(25) = f00

                                      ! i2 = in(2) =   31
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(31) = xc00
                                      yin(31) = yc00
                                      zin(31) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   26

                                      xin(26) = xcp00
                                      yin(26) = ycp00
                                      zin(26) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   32
                                      ! i2 =   31

                                      xin(32) = xcp00*xin(31) + cp10
                                      yin(32) = ycp00*yin(31) + cp10
                                      zin(32) = zcp00*zin(31) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   31

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   37
                                      ! i3 =   25
                                      ! i4 =   31

                                      xin(37) = c10*xin(25) + xc00*xin(31)
                                      yin(37) = c10*yin(25) + yc00*yin(31)
                                      zin(37) = c10*zin(25) + zc00*zin(31)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   38
                                      ! i5 =   37
                                      ! i4 =   31

                                      xin(38) = xcp00*xin(37) + cp10*xin(31)
                                      yin(38) = ycp00*yin(37) + cp10*yin(31)
                                      zin(38) = zcp00*zin(37) + cp10*zin(31)

                                      ! ------------------

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   37

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   43
                                      ! i3 =   31
                                      ! i4 =   37

                                      xin(43) = c10*xin(31) + xc00*xin(37)
                                      yin(43) = c10*yin(31) + yc00*yin(37)
                                      zin(43) = c10*zin(31) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   44
                                      ! i5 =   43
                                      ! i4 =   37

                                      xin(44) = xcp00*xin(43) + cp10*xin(37)
                                      yin(44) = ycp00*yin(43) + cp10*yin(37)
                                      zin(44) = zcp00*zin(43) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   45
                                      ! i3 =   37
                                      ! i4 =   43

                                      xin(45) = c10*xin(37) + xc00*xin(43)
                                      yin(45) = c10*yin(37) + yc00*yin(43)
                                      zin(45) = c10*zin(37) + zc00*zin(43)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   46
                                      ! i5 =   45
                                      ! i4 =   43

                                      xin(46) = xcp00*xin(45) + cp10*xin(43)
                                      yin(46) = ycp00*yin(45) + cp10*yin(43)
                                      zin(46) = zcp00*zin(45) + cp10*zin(43)

                                      ! ------------------

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   45

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   47
                                      ! i3 =   43
                                      ! i4 =   45

                                      xin(47) = c10*xin(43) + xc00*xin(45)
                                      yin(47) = c10*yin(43) + yc00*yin(45)
                                      zin(47) = c10*zin(43) + zc00*zin(45)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   48
                                      ! i5 =   47
                                      ! i4 =   45

                                      xin(48) = xcp00*xin(47) + cp10*xin(45)
                                      yin(48) = ycp00*yin(47) + cp10*yin(45)
                                      zin(48) = zcp00*zin(47) + cp10*zin(45)

                                      ! ------------------

                                      ! i3 = i4 =   45
                                      ! i4 = i5 =   47

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   47

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(47) = xin(47) + dxij*xin(45)
                                      yin(47) = yin(47) + dyij*yin(45)
                                      zin(47) = zin(47) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   43

                                      xin(45) = xin(45) + dxij*xin(43)
                                      yin(45) = yin(45) + dyij*yin(43)
                                      zin(45) = zin(45) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(47) = xin(47) + dxij*xin(45)
                                      yin(47) = yin(47) + dyij*yin(45)
                                      zin(47) = zin(47) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   27

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   27

                                      ! do ni = 1,    3

                                      xin(27) = xin(31) + dxij*xin(25)
                                      yin(27) = yin(31) + dyij*yin(25)
                                      zin(27) = zin(31) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    2

                                      xin(33) = xin(37) + dxij*xin(31)
                                      yin(33) = yin(37) + dyij*yin(31)
                                      zin(33) = zin(37) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    3

                                      xin(39) = xin(43) + dxij*xin(37)
                                      yin(39) = yin(43) + dyij*yin(37)
                                      zin(39) = zin(43) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   29

                                      ! nj =    2

                                      ! i4 = i3 =   29

                                      ! do ni = 1,    3

                                      xin(29) = xin(33) + dxij*xin(27)
                                      yin(29) = yin(33) + dyij*yin(27)
                                      zin(29) = zin(33) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    2

                                      xin(35) = xin(39) + dxij*xin(33)
                                      yin(35) = yin(39) + dyij*yin(33)
                                      zin(35) = zin(39) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    3

                                      xin(41) = xin(45) + dxij*xin(39)
                                      yin(41) = yin(45) + dyij*yin(39)
                                      zin(41) = zin(45) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   31

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   46

                                      xin(48) = xin(48) + dxij*xin(46)
                                      yin(48) = yin(48) + dyij*yin(46)
                                      zin(48) = zin(48) + dzij*zin(46)

                                      ! i3 = i4 =   46
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   44

                                      xin(46) = xin(46) + dxij*xin(44)
                                      yin(46) = yin(46) + dyij*yin(44)
                                      zin(46) = zin(46) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   46

                                      xin(48) = xin(48) + dxij*xin(46)
                                      yin(48) = yin(48) + dyij*yin(46)
                                      zin(48) = zin(48) + dzij*zin(46)

                                      ! i3 = i4 =   46
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   28

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   28

                                      ! do ni = 1,    3

                                      xin(28) = xin(32) + dxij*xin(26)
                                      yin(28) = yin(32) + dyij*yin(26)
                                      zin(28) = zin(32) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    2

                                      xin(34) = xin(38) + dxij*xin(32)
                                      yin(34) = yin(38) + dyij*yin(32)
                                      zin(34) = zin(38) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    3

                                      xin(40) = xin(44) + dxij*xin(38)
                                      yin(40) = yin(44) + dyij*yin(38)
                                      zin(40) = zin(44) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   30

                                      ! nj =    2

                                      ! i4 = i3 =   30

                                      ! do ni = 1,    3

                                      xin(30) = xin(34) + dxij*xin(28)
                                      yin(30) = yin(34) + dyij*yin(28)
                                      zin(30) = zin(34) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    2

                                      xin(36) = xin(40) + dxij*xin(34)
                                      yin(36) = yin(40) + dyij*yin(34)
                                      zin(36) = zin(40) + dzij*zin(34)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    3

                                      xin(42) = xin(46) + dxij*xin(40)
                                      yin(42) = yin(46) + dyij*yin(40)
                                      zin(42) = zin(46) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   32

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   48

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

                                      ! i1 = in(1) =   49

                                      xin(49) = 1.0_dp
                                      yin(49) = 1.0_dp
                                      zin(49) = f00

                                      ! i2 = in(2) =   55
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(55) = xc00
                                      yin(55) = yc00
                                      zin(55) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   50

                                      xin(50) = xcp00
                                      yin(50) = ycp00
                                      zin(50) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   56
                                      ! i2 =   55

                                      xin(56) = xcp00*xin(55) + cp10
                                      yin(56) = ycp00*yin(55) + cp10
                                      zin(56) = zcp00*zin(55) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   55

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   49
                                      ! i4 =   55

                                      xin(61) = c10*xin(49) + xc00*xin(55)
                                      yin(61) = c10*yin(49) + yc00*yin(55)
                                      zin(61) = c10*zin(49) + zc00*zin(55)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   62
                                      ! i5 =   61
                                      ! i4 =   55

                                      xin(62) = xcp00*xin(61) + cp10*xin(55)
                                      yin(62) = ycp00*yin(61) + cp10*yin(55)
                                      zin(62) = zcp00*zin(61) + cp10*zin(55)

                                      ! ------------------

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   61

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   67
                                      ! i3 =   55
                                      ! i4 =   61

                                      xin(67) = c10*xin(55) + xc00*xin(61)
                                      yin(67) = c10*yin(55) + yc00*yin(61)
                                      zin(67) = c10*zin(55) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   68
                                      ! i5 =   67
                                      ! i4 =   61

                                      xin(68) = xcp00*xin(67) + cp10*xin(61)
                                      yin(68) = ycp00*yin(67) + cp10*yin(61)
                                      zin(68) = zcp00*zin(67) + cp10*zin(61)

                                      ! ------------------

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   67

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   69
                                      ! i3 =   61
                                      ! i4 =   67

                                      xin(69) = c10*xin(61) + xc00*xin(67)
                                      yin(69) = c10*yin(61) + yc00*yin(67)
                                      zin(69) = c10*zin(61) + zc00*zin(67)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   70
                                      ! i5 =   69
                                      ! i4 =   67

                                      xin(70) = xcp00*xin(69) + cp10*xin(67)
                                      yin(70) = ycp00*yin(69) + cp10*yin(67)
                                      zin(70) = zcp00*zin(69) + cp10*zin(67)

                                      ! ------------------

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   69

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   71
                                      ! i3 =   67
                                      ! i4 =   69

                                      xin(71) = c10*xin(67) + xc00*xin(69)
                                      yin(71) = c10*yin(67) + yc00*yin(69)
                                      zin(71) = c10*zin(67) + zc00*zin(69)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   72
                                      ! i5 =   71
                                      ! i4 =   69

                                      xin(72) = xcp00*xin(71) + cp10*xin(69)
                                      yin(72) = ycp00*yin(71) + cp10*yin(69)
                                      zin(72) = zcp00*zin(71) + cp10*zin(69)

                                      ! ------------------

                                      ! i3 = i4 =   69
                                      ! i4 = i5 =   71

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   71

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   69

                                      xin(71) = xin(71) + dxij*xin(69)
                                      yin(71) = yin(71) + dyij*yin(69)
                                      zin(71) = zin(71) + dzij*zin(69)

                                      ! i3 = i4 =   69
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   67

                                      xin(69) = xin(69) + dxij*xin(67)
                                      yin(69) = yin(69) + dyij*yin(67)
                                      zin(69) = zin(69) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   69

                                      xin(71) = xin(71) + dxij*xin(69)
                                      yin(71) = yin(71) + dyij*yin(69)
                                      zin(71) = zin(71) + dzij*zin(69)

                                      ! i3 = i4 =   69
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   51

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   51

                                      ! do ni = 1,    3

                                      xin(51) = xin(55) + dxij*xin(49)
                                      yin(51) = yin(55) + dyij*yin(49)
                                      zin(51) = zin(55) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    2

                                      xin(57) = xin(61) + dxij*xin(55)
                                      yin(57) = yin(61) + dyij*yin(55)
                                      zin(57) = zin(61) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    3

                                      xin(63) = xin(67) + dxij*xin(61)
                                      yin(63) = yin(67) + dyij*yin(61)
                                      zin(63) = zin(67) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   53

                                      ! nj =    2

                                      ! i4 = i3 =   53

                                      ! do ni = 1,    3

                                      xin(53) = xin(57) + dxij*xin(51)
                                      yin(53) = yin(57) + dyij*yin(51)
                                      zin(53) = zin(57) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    2

                                      xin(59) = xin(63) + dxij*xin(57)
                                      yin(59) = yin(63) + dyij*yin(57)
                                      zin(59) = zin(63) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    3

                                      xin(65) = xin(69) + dxij*xin(63)
                                      yin(65) = yin(69) + dyij*yin(63)
                                      zin(65) = zin(69) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   55

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   70

                                      xin(72) = xin(72) + dxij*xin(70)
                                      yin(72) = yin(72) + dyij*yin(70)
                                      zin(72) = zin(72) + dzij*zin(70)

                                      ! i3 = i4 =   70
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   68

                                      xin(70) = xin(70) + dxij*xin(68)
                                      yin(70) = yin(70) + dyij*yin(68)
                                      zin(70) = zin(70) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   70

                                      xin(72) = xin(72) + dxij*xin(70)
                                      yin(72) = yin(72) + dyij*yin(70)
                                      zin(72) = zin(72) + dzij*zin(70)

                                      ! i3 = i4 =   70
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   52

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   52

                                      ! do ni = 1,    3

                                      xin(52) = xin(56) + dxij*xin(50)
                                      yin(52) = yin(56) + dyij*yin(50)
                                      zin(52) = zin(56) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    2

                                      xin(58) = xin(62) + dxij*xin(56)
                                      yin(58) = yin(62) + dyij*yin(56)
                                      zin(58) = zin(62) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    3

                                      xin(64) = xin(68) + dxij*xin(62)
                                      yin(64) = yin(68) + dyij*yin(62)
                                      zin(64) = zin(68) + dzij*zin(62)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   54

                                      ! nj =    2

                                      ! i4 = i3 =   54

                                      ! do ni = 1,    3

                                      xin(54) = xin(58) + dxij*xin(52)
                                      yin(54) = yin(58) + dyij*yin(52)
                                      zin(54) = zin(58) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    2

                                      xin(60) = xin(64) + dxij*xin(58)
                                      yin(60) = yin(64) + dyij*yin(58)
                                      zin(60) = zin(64) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    3

                                      xin(66) = xin(70) + dxij*xin(64)
                                      yin(66) = yin(70) + dyij*yin(64)
                                      zin(66) = zin(70) + dzij*zin(64)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   56

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   72

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

                                      ! i1 = in(1) =   73

                                      xin(73) = 1.0_dp
                                      yin(73) = 1.0_dp
                                      zin(73) = f00

                                      ! i2 = in(2) =   79
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(79) = xc00
                                      yin(79) = yc00
                                      zin(79) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   74

                                      xin(74) = xcp00
                                      yin(74) = ycp00
                                      zin(74) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   80
                                      ! i2 =   79

                                      xin(80) = xcp00*xin(79) + cp10
                                      yin(80) = ycp00*yin(79) + cp10
                                      zin(80) = zcp00*zin(79) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   79

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   85
                                      ! i3 =   73
                                      ! i4 =   79

                                      xin(85) = c10*xin(73) + xc00*xin(79)
                                      yin(85) = c10*yin(73) + yc00*yin(79)
                                      zin(85) = c10*zin(73) + zc00*zin(79)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   86
                                      ! i5 =   85
                                      ! i4 =   79

                                      xin(86) = xcp00*xin(85) + cp10*xin(79)
                                      yin(86) = ycp00*yin(85) + cp10*yin(79)
                                      zin(86) = zcp00*zin(85) + cp10*zin(79)

                                      ! ------------------

                                      ! i3 = i4 =   79
                                      ! i4 = i5 =   85

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   79
                                      ! i4 =   85

                                      xin(91) = c10*xin(79) + xc00*xin(85)
                                      yin(91) = c10*yin(79) + yc00*yin(85)
                                      zin(91) = c10*zin(79) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   92
                                      ! i5 =   91
                                      ! i4 =   85

                                      xin(92) = xcp00*xin(91) + cp10*xin(85)
                                      yin(92) = ycp00*yin(91) + cp10*yin(85)
                                      zin(92) = zcp00*zin(91) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   93
                                      ! i3 =   85
                                      ! i4 =   91

                                      xin(93) = c10*xin(85) + xc00*xin(91)
                                      yin(93) = c10*yin(85) + yc00*yin(91)
                                      zin(93) = c10*zin(85) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   94
                                      ! i5 =   93
                                      ! i4 =   91

                                      xin(94) = xcp00*xin(93) + cp10*xin(91)
                                      yin(94) = ycp00*yin(93) + cp10*yin(91)
                                      zin(94) = zcp00*zin(93) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   93

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   95
                                      ! i3 =   91
                                      ! i4 =   93

                                      xin(95) = c10*xin(91) + xc00*xin(93)
                                      yin(95) = c10*yin(91) + yc00*yin(93)
                                      zin(95) = c10*zin(91) + zc00*zin(93)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   96
                                      ! i5 =   95
                                      ! i4 =   93

                                      xin(96) = xcp00*xin(95) + cp10*xin(93)
                                      yin(96) = ycp00*yin(95) + cp10*yin(93)
                                      zin(96) = zcp00*zin(95) + cp10*zin(93)

                                      ! ------------------

                                      ! i3 = i4 =   93
                                      ! i4 = i5 =   95

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   95

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   93

                                      xin(95) = xin(95) + dxij*xin(93)
                                      yin(95) = yin(95) + dyij*yin(93)
                                      zin(95) = zin(95) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   91

                                      xin(93) = xin(93) + dxij*xin(91)
                                      yin(93) = yin(93) + dyij*yin(91)
                                      zin(93) = zin(93) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   93

                                      xin(95) = xin(95) + dxij*xin(93)
                                      yin(95) = yin(95) + dyij*yin(93)
                                      zin(95) = zin(95) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   75

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   75

                                      ! do ni = 1,    3

                                      xin(75) = xin(79) + dxij*xin(73)
                                      yin(75) = yin(79) + dyij*yin(73)
                                      zin(75) = zin(79) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    2

                                      xin(81) = xin(85) + dxij*xin(79)
                                      yin(81) = yin(85) + dyij*yin(79)
                                      zin(81) = zin(85) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    3

                                      xin(87) = xin(91) + dxij*xin(85)
                                      yin(87) = yin(91) + dyij*yin(85)
                                      zin(87) = zin(91) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   77

                                      ! nj =    2

                                      ! i4 = i3 =   77

                                      ! do ni = 1,    3

                                      xin(77) = xin(81) + dxij*xin(75)
                                      yin(77) = yin(81) + dyij*yin(75)
                                      zin(77) = zin(81) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    2

                                      xin(83) = xin(87) + dxij*xin(81)
                                      yin(83) = yin(87) + dyij*yin(81)
                                      zin(83) = zin(87) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    3

                                      xin(89) = xin(93) + dxij*xin(87)
                                      yin(89) = yin(93) + dyij*yin(87)
                                      zin(89) = zin(93) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   79

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   94

                                      xin(96) = xin(96) + dxij*xin(94)
                                      yin(96) = yin(96) + dyij*yin(94)
                                      zin(96) = zin(96) + dzij*zin(94)

                                      ! i3 = i4 =   94
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   92

                                      xin(94) = xin(94) + dxij*xin(92)
                                      yin(94) = yin(94) + dyij*yin(92)
                                      zin(94) = zin(94) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   94

                                      xin(96) = xin(96) + dxij*xin(94)
                                      yin(96) = yin(96) + dyij*yin(94)
                                      zin(96) = zin(96) + dzij*zin(94)

                                      ! i3 = i4 =   94
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   76

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   76

                                      ! do ni = 1,    3

                                      xin(76) = xin(80) + dxij*xin(74)
                                      yin(76) = yin(80) + dyij*yin(74)
                                      zin(76) = zin(80) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   82

                                      ! ni =    2

                                      xin(82) = xin(86) + dxij*xin(80)
                                      yin(82) = yin(86) + dyij*yin(80)
                                      zin(82) = zin(86) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    3

                                      xin(88) = xin(92) + dxij*xin(86)
                                      yin(88) = yin(92) + dyij*yin(86)
                                      zin(88) = zin(92) + dzij*zin(86)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   78

                                      ! nj =    2

                                      ! i4 = i3 =   78

                                      ! do ni = 1,    3

                                      xin(78) = xin(82) + dxij*xin(76)
                                      yin(78) = yin(82) + dyij*yin(76)
                                      zin(78) = zin(82) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    2

                                      xin(84) = xin(88) + dxij*xin(82)
                                      yin(84) = yin(88) + dyij*yin(82)
                                      zin(84) = zin(88) + dzij*zin(82)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    3

                                      xin(90) = xin(94) + dxij*xin(88)
                                      yin(90) = yin(94) + dyij*yin(88)
                                      zin(90) = zin(94) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   80

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d23bra(  1)*d01ket(  1)*(xin(  24)*yin(   1)*zin(   1)+xin(  48)*yin(  25)*zin(  25)+xin(  72)*yin(  49)*zin(  49)+xin(  96)*yin(  73)*zin(  73))
          eri_value(    2)=eri_value(    2)+d23bra(  1)*d01ket(  2)*(xin(  23)*yin(   2)*zin(   1)+xin(  47)*yin(  26)*zin(  25)+xin(  71)*yin(  50)*zin(  49)+xin(  95)*yin(  74)*zin(  73))
          eri_value(    3)=eri_value(    3)+d23bra(  1)*d01ket(  3)*(xin(  23)*yin(   1)*zin(   2)+xin(  47)*yin(  25)*zin(  26)+xin(  71)*yin(  49)*zin(  50)+xin(  95)*yin(  73)*zin(  74))
          eri_value(    4)=eri_value(    4)+d23bra(  2)*d01ket(  1)*(xin(  20)*yin(   5)*zin(   1)+xin(  44)*yin(  29)*zin(  25)+xin(  68)*yin(  53)*zin(  49)+xin(  92)*yin(  77)*zin(  73))
          eri_value(    5)=eri_value(    5)+d23bra(  2)*d01ket(  2)*(xin(  19)*yin(   6)*zin(   1)+xin(  43)*yin(  30)*zin(  25)+xin(  67)*yin(  54)*zin(  49)+xin(  91)*yin(  78)*zin(  73))
          eri_value(    6)=eri_value(    6)+d23bra(  2)*d01ket(  3)*(xin(  19)*yin(   5)*zin(   2)+xin(  43)*yin(  29)*zin(  26)+xin(  67)*yin(  53)*zin(  50)+xin(  91)*yin(  77)*zin(  74))
          eri_value(    7)=eri_value(    7)+d23bra(  3)*d01ket(  1)*(xin(  20)*yin(   1)*zin(   5)+xin(  44)*yin(  25)*zin(  29)+xin(  68)*yin(  49)*zin(  53)+xin(  92)*yin(  73)*zin(  77))
          eri_value(    8)=eri_value(    8)+d23bra(  3)*d01ket(  2)*(xin(  19)*yin(   2)*zin(   5)+xin(  43)*yin(  26)*zin(  29)+xin(  67)*yin(  50)*zin(  53)+xin(  91)*yin(  74)*zin(  77))
          eri_value(    9)=eri_value(    9)+d23bra(  3)*d01ket(  3)*(xin(  19)*yin(   1)*zin(   6)+xin(  43)*yin(  25)*zin(  30)+xin(  67)*yin(  49)*zin(  54)+xin(  91)*yin(  73)*zin(  78))
          eri_value(   10)=eri_value(   10)+d23bra(  4)*d01ket(  1)*(xin(  22)*yin(   3)*zin(   1)+xin(  46)*yin(  27)*zin(  25)+xin(  70)*yin(  51)*zin(  49)+xin(  94)*yin(  75)*zin(  73))
          eri_value(   11)=eri_value(   11)+d23bra(  4)*d01ket(  2)*(xin(  21)*yin(   4)*zin(   1)+xin(  45)*yin(  28)*zin(  25)+xin(  69)*yin(  52)*zin(  49)+xin(  93)*yin(  76)*zin(  73))
          eri_value(   12)=eri_value(   12)+d23bra(  4)*d01ket(  3)*(xin(  21)*yin(   3)*zin(   2)+xin(  45)*yin(  27)*zin(  26)+xin(  69)*yin(  51)*zin(  50)+xin(  93)*yin(  75)*zin(  74))
          eri_value(   13)=eri_value(   13)+d23bra(  5)*d01ket(  1)*(xin(  22)*yin(   1)*zin(   3)+xin(  46)*yin(  25)*zin(  27)+xin(  70)*yin(  49)*zin(  51)+xin(  94)*yin(  73)*zin(  75))
          eri_value(   14)=eri_value(   14)+d23bra(  5)*d01ket(  2)*(xin(  21)*yin(   2)*zin(   3)+xin(  45)*yin(  26)*zin(  27)+xin(  69)*yin(  50)*zin(  51)+xin(  93)*yin(  74)*zin(  75))
          eri_value(   15)=eri_value(   15)+d23bra(  5)*d01ket(  3)*(xin(  21)*yin(   1)*zin(   4)+xin(  45)*yin(  25)*zin(  28)+xin(  69)*yin(  49)*zin(  52)+xin(  93)*yin(  73)*zin(  76))
          eri_value(   16)=eri_value(   16)+d23bra(  6)*d01ket(  1)*(xin(  20)*yin(   3)*zin(   3)+xin(  44)*yin(  27)*zin(  27)+xin(  68)*yin(  51)*zin(  51)+xin(  92)*yin(  75)*zin(  75))
          eri_value(   17)=eri_value(   17)+d23bra(  6)*d01ket(  2)*(xin(  19)*yin(   4)*zin(   3)+xin(  43)*yin(  28)*zin(  27)+xin(  67)*yin(  52)*zin(  51)+xin(  91)*yin(  76)*zin(  75))
          eri_value(   18)=eri_value(   18)+d23bra(  6)*d01ket(  3)*(xin(  19)*yin(   3)*zin(   4)+xin(  43)*yin(  27)*zin(  28)+xin(  67)*yin(  51)*zin(  52)+xin(  91)*yin(  75)*zin(  76))
          eri_value(   19)=eri_value(   19)+d23bra(  7)*d01ket(  1)*(xin(   6)*yin(  19)*zin(   1)+xin(  30)*yin(  43)*zin(  25)+xin(  54)*yin(  67)*zin(  49)+xin(  78)*yin(  91)*zin(  73))
          eri_value(   20)=eri_value(   20)+d23bra(  7)*d01ket(  2)*(xin(   5)*yin(  20)*zin(   1)+xin(  29)*yin(  44)*zin(  25)+xin(  53)*yin(  68)*zin(  49)+xin(  77)*yin(  92)*zin(  73))
          eri_value(   21)=eri_value(   21)+d23bra(  7)*d01ket(  3)*(xin(   5)*yin(  19)*zin(   2)+xin(  29)*yin(  43)*zin(  26)+xin(  53)*yin(  67)*zin(  50)+xin(  77)*yin(  91)*zin(  74))
          eri_value(   22)=eri_value(   22)+d23bra(  8)*d01ket(  1)*(xin(   2)*yin(  23)*zin(   1)+xin(  26)*yin(  47)*zin(  25)+xin(  50)*yin(  71)*zin(  49)+xin(  74)*yin(  95)*zin(  73))
          eri_value(   23)=eri_value(   23)+d23bra(  8)*d01ket(  2)*(xin(   1)*yin(  24)*zin(   1)+xin(  25)*yin(  48)*zin(  25)+xin(  49)*yin(  72)*zin(  49)+xin(  73)*yin(  96)*zin(  73))
          eri_value(   24)=eri_value(   24)+d23bra(  8)*d01ket(  3)*(xin(   1)*yin(  23)*zin(   2)+xin(  25)*yin(  47)*zin(  26)+xin(  49)*yin(  71)*zin(  50)+xin(  73)*yin(  95)*zin(  74))
          eri_value(   25)=eri_value(   25)+d23bra(  9)*d01ket(  1)*(xin(   2)*yin(  19)*zin(   5)+xin(  26)*yin(  43)*zin(  29)+xin(  50)*yin(  67)*zin(  53)+xin(  74)*yin(  91)*zin(  77))
          eri_value(   26)=eri_value(   26)+d23bra(  9)*d01ket(  2)*(xin(   1)*yin(  20)*zin(   5)+xin(  25)*yin(  44)*zin(  29)+xin(  49)*yin(  68)*zin(  53)+xin(  73)*yin(  92)*zin(  77))
          eri_value(   27)=eri_value(   27)+d23bra(  9)*d01ket(  3)*(xin(   1)*yin(  19)*zin(   6)+xin(  25)*yin(  43)*zin(  30)+xin(  49)*yin(  67)*zin(  54)+xin(  73)*yin(  91)*zin(  78))
          eri_value(   28)=eri_value(   28)+d23bra( 10)*d01ket(  1)*(xin(   4)*yin(  21)*zin(   1)+xin(  28)*yin(  45)*zin(  25)+xin(  52)*yin(  69)*zin(  49)+xin(  76)*yin(  93)*zin(  73))
          eri_value(   29)=eri_value(   29)+d23bra( 10)*d01ket(  2)*(xin(   3)*yin(  22)*zin(   1)+xin(  27)*yin(  46)*zin(  25)+xin(  51)*yin(  70)*zin(  49)+xin(  75)*yin(  94)*zin(  73))
          eri_value(   30)=eri_value(   30)+d23bra( 10)*d01ket(  3)*(xin(   3)*yin(  21)*zin(   2)+xin(  27)*yin(  45)*zin(  26)+xin(  51)*yin(  69)*zin(  50)+xin(  75)*yin(  93)*zin(  74))
          eri_value(   31)=eri_value(   31)+d23bra( 11)*d01ket(  1)*(xin(   4)*yin(  19)*zin(   3)+xin(  28)*yin(  43)*zin(  27)+xin(  52)*yin(  67)*zin(  51)+xin(  76)*yin(  91)*zin(  75))
          eri_value(   32)=eri_value(   32)+d23bra( 11)*d01ket(  2)*(xin(   3)*yin(  20)*zin(   3)+xin(  27)*yin(  44)*zin(  27)+xin(  51)*yin(  68)*zin(  51)+xin(  75)*yin(  92)*zin(  75))
          eri_value(   33)=eri_value(   33)+d23bra( 11)*d01ket(  3)*(xin(   3)*yin(  19)*zin(   4)+xin(  27)*yin(  43)*zin(  28)+xin(  51)*yin(  67)*zin(  52)+xin(  75)*yin(  91)*zin(  76))
          eri_value(   34)=eri_value(   34)+d23bra( 12)*d01ket(  1)*(xin(   2)*yin(  21)*zin(   3)+xin(  26)*yin(  45)*zin(  27)+xin(  50)*yin(  69)*zin(  51)+xin(  74)*yin(  93)*zin(  75))
          eri_value(   35)=eri_value(   35)+d23bra( 12)*d01ket(  2)*(xin(   1)*yin(  22)*zin(   3)+xin(  25)*yin(  46)*zin(  27)+xin(  49)*yin(  70)*zin(  51)+xin(  73)*yin(  94)*zin(  75))
          eri_value(   36)=eri_value(   36)+d23bra( 12)*d01ket(  3)*(xin(   1)*yin(  21)*zin(   4)+xin(  25)*yin(  45)*zin(  28)+xin(  49)*yin(  69)*zin(  52)+xin(  73)*yin(  93)*zin(  76))
          eri_value(   37)=eri_value(   37)+d23bra( 13)*d01ket(  1)*(xin(   6)*yin(   1)*zin(  19)+xin(  30)*yin(  25)*zin(  43)+xin(  54)*yin(  49)*zin(  67)+xin(  78)*yin(  73)*zin(  91))
          eri_value(   38)=eri_value(   38)+d23bra( 13)*d01ket(  2)*(xin(   5)*yin(   2)*zin(  19)+xin(  29)*yin(  26)*zin(  43)+xin(  53)*yin(  50)*zin(  67)+xin(  77)*yin(  74)*zin(  91))
          eri_value(   39)=eri_value(   39)+d23bra( 13)*d01ket(  3)*(xin(   5)*yin(   1)*zin(  20)+xin(  29)*yin(  25)*zin(  44)+xin(  53)*yin(  49)*zin(  68)+xin(  77)*yin(  73)*zin(  92))
          eri_value(   40)=eri_value(   40)+d23bra( 14)*d01ket(  1)*(xin(   2)*yin(   5)*zin(  19)+xin(  26)*yin(  29)*zin(  43)+xin(  50)*yin(  53)*zin(  67)+xin(  74)*yin(  77)*zin(  91))
          eri_value(   41)=eri_value(   41)+d23bra( 14)*d01ket(  2)*(xin(   1)*yin(   6)*zin(  19)+xin(  25)*yin(  30)*zin(  43)+xin(  49)*yin(  54)*zin(  67)+xin(  73)*yin(  78)*zin(  91))
          eri_value(   42)=eri_value(   42)+d23bra( 14)*d01ket(  3)*(xin(   1)*yin(   5)*zin(  20)+xin(  25)*yin(  29)*zin(  44)+xin(  49)*yin(  53)*zin(  68)+xin(  73)*yin(  77)*zin(  92))
          eri_value(   43)=eri_value(   43)+d23bra( 15)*d01ket(  1)*(xin(   2)*yin(   1)*zin(  23)+xin(  26)*yin(  25)*zin(  47)+xin(  50)*yin(  49)*zin(  71)+xin(  74)*yin(  73)*zin(  95))
          eri_value(   44)=eri_value(   44)+d23bra( 15)*d01ket(  2)*(xin(   1)*yin(   2)*zin(  23)+xin(  25)*yin(  26)*zin(  47)+xin(  49)*yin(  50)*zin(  71)+xin(  73)*yin(  74)*zin(  95))
          eri_value(   45)=eri_value(   45)+d23bra( 15)*d01ket(  3)*(xin(   1)*yin(   1)*zin(  24)+xin(  25)*yin(  25)*zin(  48)+xin(  49)*yin(  49)*zin(  72)+xin(  73)*yin(  73)*zin(  96))
          eri_value(   46)=eri_value(   46)+d23bra( 16)*d01ket(  1)*(xin(   4)*yin(   3)*zin(  19)+xin(  28)*yin(  27)*zin(  43)+xin(  52)*yin(  51)*zin(  67)+xin(  76)*yin(  75)*zin(  91))
          eri_value(   47)=eri_value(   47)+d23bra( 16)*d01ket(  2)*(xin(   3)*yin(   4)*zin(  19)+xin(  27)*yin(  28)*zin(  43)+xin(  51)*yin(  52)*zin(  67)+xin(  75)*yin(  76)*zin(  91))
          eri_value(   48)=eri_value(   48)+d23bra( 16)*d01ket(  3)*(xin(   3)*yin(   3)*zin(  20)+xin(  27)*yin(  27)*zin(  44)+xin(  51)*yin(  51)*zin(  68)+xin(  75)*yin(  75)*zin(  92))
          eri_value(   49)=eri_value(   49)+d23bra( 17)*d01ket(  1)*(xin(   4)*yin(   1)*zin(  21)+xin(  28)*yin(  25)*zin(  45)+xin(  52)*yin(  49)*zin(  69)+xin(  76)*yin(  73)*zin(  93))
          eri_value(   50)=eri_value(   50)+d23bra( 17)*d01ket(  2)*(xin(   3)*yin(   2)*zin(  21)+xin(  27)*yin(  26)*zin(  45)+xin(  51)*yin(  50)*zin(  69)+xin(  75)*yin(  74)*zin(  93))
          eri_value(   51)=eri_value(   51)+d23bra( 17)*d01ket(  3)*(xin(   3)*yin(   1)*zin(  22)+xin(  27)*yin(  25)*zin(  46)+xin(  51)*yin(  49)*zin(  70)+xin(  75)*yin(  73)*zin(  94))
          eri_value(   52)=eri_value(   52)+d23bra( 18)*d01ket(  1)*(xin(   2)*yin(   3)*zin(  21)+xin(  26)*yin(  27)*zin(  45)+xin(  50)*yin(  51)*zin(  69)+xin(  74)*yin(  75)*zin(  93))
          eri_value(   53)=eri_value(   53)+d23bra( 18)*d01ket(  2)*(xin(   1)*yin(   4)*zin(  21)+xin(  25)*yin(  28)*zin(  45)+xin(  49)*yin(  52)*zin(  69)+xin(  73)*yin(  76)*zin(  93))
          eri_value(   54)=eri_value(   54)+d23bra( 18)*d01ket(  3)*(xin(   1)*yin(   3)*zin(  22)+xin(  25)*yin(  27)*zin(  46)+xin(  49)*yin(  51)*zin(  70)+xin(  73)*yin(  75)*zin(  94))
          eri_value(   55)=eri_value(   55)+d23bra( 19)*d01ket(  1)*(xin(  18)*yin(   7)*zin(   1)+xin(  42)*yin(  31)*zin(  25)+xin(  66)*yin(  55)*zin(  49)+xin(  90)*yin(  79)*zin(  73))
          eri_value(   56)=eri_value(   56)+d23bra( 19)*d01ket(  2)*(xin(  17)*yin(   8)*zin(   1)+xin(  41)*yin(  32)*zin(  25)+xin(  65)*yin(  56)*zin(  49)+xin(  89)*yin(  80)*zin(  73))
          eri_value(   57)=eri_value(   57)+d23bra( 19)*d01ket(  3)*(xin(  17)*yin(   7)*zin(   2)+xin(  41)*yin(  31)*zin(  26)+xin(  65)*yin(  55)*zin(  50)+xin(  89)*yin(  79)*zin(  74))
          eri_value(   58)=eri_value(   58)+d23bra( 20)*d01ket(  1)*(xin(  14)*yin(  11)*zin(   1)+xin(  38)*yin(  35)*zin(  25)+xin(  62)*yin(  59)*zin(  49)+xin(  86)*yin(  83)*zin(  73))
          eri_value(   59)=eri_value(   59)+d23bra( 20)*d01ket(  2)*(xin(  13)*yin(  12)*zin(   1)+xin(  37)*yin(  36)*zin(  25)+xin(  61)*yin(  60)*zin(  49)+xin(  85)*yin(  84)*zin(  73))
          eri_value(   60)=eri_value(   60)+d23bra( 20)*d01ket(  3)*(xin(  13)*yin(  11)*zin(   2)+xin(  37)*yin(  35)*zin(  26)+xin(  61)*yin(  59)*zin(  50)+xin(  85)*yin(  83)*zin(  74))
          eri_value(   61)=eri_value(   61)+d23bra( 21)*d01ket(  1)*(xin(  14)*yin(   7)*zin(   5)+xin(  38)*yin(  31)*zin(  29)+xin(  62)*yin(  55)*zin(  53)+xin(  86)*yin(  79)*zin(  77))
          eri_value(   62)=eri_value(   62)+d23bra( 21)*d01ket(  2)*(xin(  13)*yin(   8)*zin(   5)+xin(  37)*yin(  32)*zin(  29)+xin(  61)*yin(  56)*zin(  53)+xin(  85)*yin(  80)*zin(  77))
          eri_value(   63)=eri_value(   63)+d23bra( 21)*d01ket(  3)*(xin(  13)*yin(   7)*zin(   6)+xin(  37)*yin(  31)*zin(  30)+xin(  61)*yin(  55)*zin(  54)+xin(  85)*yin(  79)*zin(  78))
          eri_value(   64)=eri_value(   64)+d23bra( 22)*d01ket(  1)*(xin(  16)*yin(   9)*zin(   1)+xin(  40)*yin(  33)*zin(  25)+xin(  64)*yin(  57)*zin(  49)+xin(  88)*yin(  81)*zin(  73))
          eri_value(   65)=eri_value(   65)+d23bra( 22)*d01ket(  2)*(xin(  15)*yin(  10)*zin(   1)+xin(  39)*yin(  34)*zin(  25)+xin(  63)*yin(  58)*zin(  49)+xin(  87)*yin(  82)*zin(  73))
          eri_value(   66)=eri_value(   66)+d23bra( 22)*d01ket(  3)*(xin(  15)*yin(   9)*zin(   2)+xin(  39)*yin(  33)*zin(  26)+xin(  63)*yin(  57)*zin(  50)+xin(  87)*yin(  81)*zin(  74))
          eri_value(   67)=eri_value(   67)+d23bra( 23)*d01ket(  1)*(xin(  16)*yin(   7)*zin(   3)+xin(  40)*yin(  31)*zin(  27)+xin(  64)*yin(  55)*zin(  51)+xin(  88)*yin(  79)*zin(  75))
          eri_value(   68)=eri_value(   68)+d23bra( 23)*d01ket(  2)*(xin(  15)*yin(   8)*zin(   3)+xin(  39)*yin(  32)*zin(  27)+xin(  63)*yin(  56)*zin(  51)+xin(  87)*yin(  80)*zin(  75))
          eri_value(   69)=eri_value(   69)+d23bra( 23)*d01ket(  3)*(xin(  15)*yin(   7)*zin(   4)+xin(  39)*yin(  31)*zin(  28)+xin(  63)*yin(  55)*zin(  52)+xin(  87)*yin(  79)*zin(  76))
          eri_value(   70)=eri_value(   70)+d23bra( 24)*d01ket(  1)*(xin(  14)*yin(   9)*zin(   3)+xin(  38)*yin(  33)*zin(  27)+xin(  62)*yin(  57)*zin(  51)+xin(  86)*yin(  81)*zin(  75))
          eri_value(   71)=eri_value(   71)+d23bra( 24)*d01ket(  2)*(xin(  13)*yin(  10)*zin(   3)+xin(  37)*yin(  34)*zin(  27)+xin(  61)*yin(  58)*zin(  51)+xin(  85)*yin(  82)*zin(  75))
          eri_value(   72)=eri_value(   72)+d23bra( 24)*d01ket(  3)*(xin(  13)*yin(   9)*zin(   4)+xin(  37)*yin(  33)*zin(  28)+xin(  61)*yin(  57)*zin(  52)+xin(  85)*yin(  81)*zin(  76))
          eri_value(   73)=eri_value(   73)+d23bra( 25)*d01ket(  1)*(xin(  18)*yin(   1)*zin(   7)+xin(  42)*yin(  25)*zin(  31)+xin(  66)*yin(  49)*zin(  55)+xin(  90)*yin(  73)*zin(  79))
          eri_value(   74)=eri_value(   74)+d23bra( 25)*d01ket(  2)*(xin(  17)*yin(   2)*zin(   7)+xin(  41)*yin(  26)*zin(  31)+xin(  65)*yin(  50)*zin(  55)+xin(  89)*yin(  74)*zin(  79))
          eri_value(   75)=eri_value(   75)+d23bra( 25)*d01ket(  3)*(xin(  17)*yin(   1)*zin(   8)+xin(  41)*yin(  25)*zin(  32)+xin(  65)*yin(  49)*zin(  56)+xin(  89)*yin(  73)*zin(  80))
          eri_value(   76)=eri_value(   76)+d23bra( 26)*d01ket(  1)*(xin(  14)*yin(   5)*zin(   7)+xin(  38)*yin(  29)*zin(  31)+xin(  62)*yin(  53)*zin(  55)+xin(  86)*yin(  77)*zin(  79))
          eri_value(   77)=eri_value(   77)+d23bra( 26)*d01ket(  2)*(xin(  13)*yin(   6)*zin(   7)+xin(  37)*yin(  30)*zin(  31)+xin(  61)*yin(  54)*zin(  55)+xin(  85)*yin(  78)*zin(  79))
          eri_value(   78)=eri_value(   78)+d23bra( 26)*d01ket(  3)*(xin(  13)*yin(   5)*zin(   8)+xin(  37)*yin(  29)*zin(  32)+xin(  61)*yin(  53)*zin(  56)+xin(  85)*yin(  77)*zin(  80))
          eri_value(   79)=eri_value(   79)+d23bra( 27)*d01ket(  1)*(xin(  14)*yin(   1)*zin(  11)+xin(  38)*yin(  25)*zin(  35)+xin(  62)*yin(  49)*zin(  59)+xin(  86)*yin(  73)*zin(  83))
          eri_value(   80)=eri_value(   80)+d23bra( 27)*d01ket(  2)*(xin(  13)*yin(   2)*zin(  11)+xin(  37)*yin(  26)*zin(  35)+xin(  61)*yin(  50)*zin(  59)+xin(  85)*yin(  74)*zin(  83))
          eri_value(   81)=eri_value(   81)+d23bra( 27)*d01ket(  3)*(xin(  13)*yin(   1)*zin(  12)+xin(  37)*yin(  25)*zin(  36)+xin(  61)*yin(  49)*zin(  60)+xin(  85)*yin(  73)*zin(  84))
          eri_value(   82)=eri_value(   82)+d23bra( 28)*d01ket(  1)*(xin(  16)*yin(   3)*zin(   7)+xin(  40)*yin(  27)*zin(  31)+xin(  64)*yin(  51)*zin(  55)+xin(  88)*yin(  75)*zin(  79))
          eri_value(   83)=eri_value(   83)+d23bra( 28)*d01ket(  2)*(xin(  15)*yin(   4)*zin(   7)+xin(  39)*yin(  28)*zin(  31)+xin(  63)*yin(  52)*zin(  55)+xin(  87)*yin(  76)*zin(  79))
          eri_value(   84)=eri_value(   84)+d23bra( 28)*d01ket(  3)*(xin(  15)*yin(   3)*zin(   8)+xin(  39)*yin(  27)*zin(  32)+xin(  63)*yin(  51)*zin(  56)+xin(  87)*yin(  75)*zin(  80))
          eri_value(   85)=eri_value(   85)+d23bra( 29)*d01ket(  1)*(xin(  16)*yin(   1)*zin(   9)+xin(  40)*yin(  25)*zin(  33)+xin(  64)*yin(  49)*zin(  57)+xin(  88)*yin(  73)*zin(  81))
          eri_value(   86)=eri_value(   86)+d23bra( 29)*d01ket(  2)*(xin(  15)*yin(   2)*zin(   9)+xin(  39)*yin(  26)*zin(  33)+xin(  63)*yin(  50)*zin(  57)+xin(  87)*yin(  74)*zin(  81))
          eri_value(   87)=eri_value(   87)+d23bra( 29)*d01ket(  3)*(xin(  15)*yin(   1)*zin(  10)+xin(  39)*yin(  25)*zin(  34)+xin(  63)*yin(  49)*zin(  58)+xin(  87)*yin(  73)*zin(  82))
          eri_value(   88)=eri_value(   88)+d23bra( 30)*d01ket(  1)*(xin(  14)*yin(   3)*zin(   9)+xin(  38)*yin(  27)*zin(  33)+xin(  62)*yin(  51)*zin(  57)+xin(  86)*yin(  75)*zin(  81))
          eri_value(   89)=eri_value(   89)+d23bra( 30)*d01ket(  2)*(xin(  13)*yin(   4)*zin(   9)+xin(  37)*yin(  28)*zin(  33)+xin(  61)*yin(  52)*zin(  57)+xin(  85)*yin(  76)*zin(  81))
          eri_value(   90)=eri_value(   90)+d23bra( 30)*d01ket(  3)*(xin(  13)*yin(   3)*zin(  10)+xin(  37)*yin(  27)*zin(  34)+xin(  61)*yin(  51)*zin(  58)+xin(  85)*yin(  75)*zin(  82))
          eri_value(   91)=eri_value(   91)+d23bra( 31)*d01ket(  1)*(xin(  12)*yin(  13)*zin(   1)+xin(  36)*yin(  37)*zin(  25)+xin(  60)*yin(  61)*zin(  49)+xin(  84)*yin(  85)*zin(  73))
          eri_value(   92)=eri_value(   92)+d23bra( 31)*d01ket(  2)*(xin(  11)*yin(  14)*zin(   1)+xin(  35)*yin(  38)*zin(  25)+xin(  59)*yin(  62)*zin(  49)+xin(  83)*yin(  86)*zin(  73))
          eri_value(   93)=eri_value(   93)+d23bra( 31)*d01ket(  3)*(xin(  11)*yin(  13)*zin(   2)+xin(  35)*yin(  37)*zin(  26)+xin(  59)*yin(  61)*zin(  50)+xin(  83)*yin(  85)*zin(  74))
          eri_value(   94)=eri_value(   94)+d23bra( 32)*d01ket(  1)*(xin(   8)*yin(  17)*zin(   1)+xin(  32)*yin(  41)*zin(  25)+xin(  56)*yin(  65)*zin(  49)+xin(  80)*yin(  89)*zin(  73))
          eri_value(   95)=eri_value(   95)+d23bra( 32)*d01ket(  2)*(xin(   7)*yin(  18)*zin(   1)+xin(  31)*yin(  42)*zin(  25)+xin(  55)*yin(  66)*zin(  49)+xin(  79)*yin(  90)*zin(  73))
          eri_value(   96)=eri_value(   96)+d23bra( 32)*d01ket(  3)*(xin(   7)*yin(  17)*zin(   2)+xin(  31)*yin(  41)*zin(  26)+xin(  55)*yin(  65)*zin(  50)+xin(  79)*yin(  89)*zin(  74))
          eri_value(   97)=eri_value(   97)+d23bra( 33)*d01ket(  1)*(xin(   8)*yin(  13)*zin(   5)+xin(  32)*yin(  37)*zin(  29)+xin(  56)*yin(  61)*zin(  53)+xin(  80)*yin(  85)*zin(  77))
          eri_value(   98)=eri_value(   98)+d23bra( 33)*d01ket(  2)*(xin(   7)*yin(  14)*zin(   5)+xin(  31)*yin(  38)*zin(  29)+xin(  55)*yin(  62)*zin(  53)+xin(  79)*yin(  86)*zin(  77))
          eri_value(   99)=eri_value(   99)+d23bra( 33)*d01ket(  3)*(xin(   7)*yin(  13)*zin(   6)+xin(  31)*yin(  37)*zin(  30)+xin(  55)*yin(  61)*zin(  54)+xin(  79)*yin(  85)*zin(  78))
          eri_value(  100)=eri_value(  100)+d23bra( 34)*d01ket(  1)*(xin(  10)*yin(  15)*zin(   1)+xin(  34)*yin(  39)*zin(  25)+xin(  58)*yin(  63)*zin(  49)+xin(  82)*yin(  87)*zin(  73))
          eri_value(  101)=eri_value(  101)+d23bra( 34)*d01ket(  2)*(xin(   9)*yin(  16)*zin(   1)+xin(  33)*yin(  40)*zin(  25)+xin(  57)*yin(  64)*zin(  49)+xin(  81)*yin(  88)*zin(  73))
          eri_value(  102)=eri_value(  102)+d23bra( 34)*d01ket(  3)*(xin(   9)*yin(  15)*zin(   2)+xin(  33)*yin(  39)*zin(  26)+xin(  57)*yin(  63)*zin(  50)+xin(  81)*yin(  87)*zin(  74))
          eri_value(  103)=eri_value(  103)+d23bra( 35)*d01ket(  1)*(xin(  10)*yin(  13)*zin(   3)+xin(  34)*yin(  37)*zin(  27)+xin(  58)*yin(  61)*zin(  51)+xin(  82)*yin(  85)*zin(  75))
          eri_value(  104)=eri_value(  104)+d23bra( 35)*d01ket(  2)*(xin(   9)*yin(  14)*zin(   3)+xin(  33)*yin(  38)*zin(  27)+xin(  57)*yin(  62)*zin(  51)+xin(  81)*yin(  86)*zin(  75))
          eri_value(  105)=eri_value(  105)+d23bra( 35)*d01ket(  3)*(xin(   9)*yin(  13)*zin(   4)+xin(  33)*yin(  37)*zin(  28)+xin(  57)*yin(  61)*zin(  52)+xin(  81)*yin(  85)*zin(  76))
          eri_value(  106)=eri_value(  106)+d23bra( 36)*d01ket(  1)*(xin(   8)*yin(  15)*zin(   3)+xin(  32)*yin(  39)*zin(  27)+xin(  56)*yin(  63)*zin(  51)+xin(  80)*yin(  87)*zin(  75))
          eri_value(  107)=eri_value(  107)+d23bra( 36)*d01ket(  2)*(xin(   7)*yin(  16)*zin(   3)+xin(  31)*yin(  40)*zin(  27)+xin(  55)*yin(  64)*zin(  51)+xin(  79)*yin(  88)*zin(  75))
          eri_value(  108)=eri_value(  108)+d23bra( 36)*d01ket(  3)*(xin(   7)*yin(  15)*zin(   4)+xin(  31)*yin(  39)*zin(  28)+xin(  55)*yin(  63)*zin(  52)+xin(  79)*yin(  87)*zin(  76))
          eri_value(  109)=eri_value(  109)+d23bra( 37)*d01ket(  1)*(xin(   6)*yin(  13)*zin(   7)+xin(  30)*yin(  37)*zin(  31)+xin(  54)*yin(  61)*zin(  55)+xin(  78)*yin(  85)*zin(  79))
          eri_value(  110)=eri_value(  110)+d23bra( 37)*d01ket(  2)*(xin(   5)*yin(  14)*zin(   7)+xin(  29)*yin(  38)*zin(  31)+xin(  53)*yin(  62)*zin(  55)+xin(  77)*yin(  86)*zin(  79))
          eri_value(  111)=eri_value(  111)+d23bra( 37)*d01ket(  3)*(xin(   5)*yin(  13)*zin(   8)+xin(  29)*yin(  37)*zin(  32)+xin(  53)*yin(  61)*zin(  56)+xin(  77)*yin(  85)*zin(  80))
          eri_value(  112)=eri_value(  112)+d23bra( 38)*d01ket(  1)*(xin(   2)*yin(  17)*zin(   7)+xin(  26)*yin(  41)*zin(  31)+xin(  50)*yin(  65)*zin(  55)+xin(  74)*yin(  89)*zin(  79))
          eri_value(  113)=eri_value(  113)+d23bra( 38)*d01ket(  2)*(xin(   1)*yin(  18)*zin(   7)+xin(  25)*yin(  42)*zin(  31)+xin(  49)*yin(  66)*zin(  55)+xin(  73)*yin(  90)*zin(  79))
          eri_value(  114)=eri_value(  114)+d23bra( 38)*d01ket(  3)*(xin(   1)*yin(  17)*zin(   8)+xin(  25)*yin(  41)*zin(  32)+xin(  49)*yin(  65)*zin(  56)+xin(  73)*yin(  89)*zin(  80))
          eri_value(  115)=eri_value(  115)+d23bra( 39)*d01ket(  1)*(xin(   2)*yin(  13)*zin(  11)+xin(  26)*yin(  37)*zin(  35)+xin(  50)*yin(  61)*zin(  59)+xin(  74)*yin(  85)*zin(  83))
          eri_value(  116)=eri_value(  116)+d23bra( 39)*d01ket(  2)*(xin(   1)*yin(  14)*zin(  11)+xin(  25)*yin(  38)*zin(  35)+xin(  49)*yin(  62)*zin(  59)+xin(  73)*yin(  86)*zin(  83))
          eri_value(  117)=eri_value(  117)+d23bra( 39)*d01ket(  3)*(xin(   1)*yin(  13)*zin(  12)+xin(  25)*yin(  37)*zin(  36)+xin(  49)*yin(  61)*zin(  60)+xin(  73)*yin(  85)*zin(  84))
          eri_value(  118)=eri_value(  118)+d23bra( 40)*d01ket(  1)*(xin(   4)*yin(  15)*zin(   7)+xin(  28)*yin(  39)*zin(  31)+xin(  52)*yin(  63)*zin(  55)+xin(  76)*yin(  87)*zin(  79))
          eri_value(  119)=eri_value(  119)+d23bra( 40)*d01ket(  2)*(xin(   3)*yin(  16)*zin(   7)+xin(  27)*yin(  40)*zin(  31)+xin(  51)*yin(  64)*zin(  55)+xin(  75)*yin(  88)*zin(  79))
          eri_value(  120)=eri_value(  120)+d23bra( 40)*d01ket(  3)*(xin(   3)*yin(  15)*zin(   8)+xin(  27)*yin(  39)*zin(  32)+xin(  51)*yin(  63)*zin(  56)+xin(  75)*yin(  87)*zin(  80))
          eri_value(  121)=eri_value(  121)+d23bra( 41)*d01ket(  1)*(xin(   4)*yin(  13)*zin(   9)+xin(  28)*yin(  37)*zin(  33)+xin(  52)*yin(  61)*zin(  57)+xin(  76)*yin(  85)*zin(  81))
          eri_value(  122)=eri_value(  122)+d23bra( 41)*d01ket(  2)*(xin(   3)*yin(  14)*zin(   9)+xin(  27)*yin(  38)*zin(  33)+xin(  51)*yin(  62)*zin(  57)+xin(  75)*yin(  86)*zin(  81))
          eri_value(  123)=eri_value(  123)+d23bra( 41)*d01ket(  3)*(xin(   3)*yin(  13)*zin(  10)+xin(  27)*yin(  37)*zin(  34)+xin(  51)*yin(  61)*zin(  58)+xin(  75)*yin(  85)*zin(  82))
          eri_value(  124)=eri_value(  124)+d23bra( 42)*d01ket(  1)*(xin(   2)*yin(  15)*zin(   9)+xin(  26)*yin(  39)*zin(  33)+xin(  50)*yin(  63)*zin(  57)+xin(  74)*yin(  87)*zin(  81))
          eri_value(  125)=eri_value(  125)+d23bra( 42)*d01ket(  2)*(xin(   1)*yin(  16)*zin(   9)+xin(  25)*yin(  40)*zin(  33)+xin(  49)*yin(  64)*zin(  57)+xin(  73)*yin(  88)*zin(  81))
          eri_value(  126)=eri_value(  126)+d23bra( 42)*d01ket(  3)*(xin(   1)*yin(  15)*zin(  10)+xin(  25)*yin(  39)*zin(  34)+xin(  49)*yin(  63)*zin(  58)+xin(  73)*yin(  87)*zin(  82))
          eri_value(  127)=eri_value(  127)+d23bra( 43)*d01ket(  1)*(xin(  12)*yin(   1)*zin(  13)+xin(  36)*yin(  25)*zin(  37)+xin(  60)*yin(  49)*zin(  61)+xin(  84)*yin(  73)*zin(  85))
          eri_value(  128)=eri_value(  128)+d23bra( 43)*d01ket(  2)*(xin(  11)*yin(   2)*zin(  13)+xin(  35)*yin(  26)*zin(  37)+xin(  59)*yin(  50)*zin(  61)+xin(  83)*yin(  74)*zin(  85))
          eri_value(  129)=eri_value(  129)+d23bra( 43)*d01ket(  3)*(xin(  11)*yin(   1)*zin(  14)+xin(  35)*yin(  25)*zin(  38)+xin(  59)*yin(  49)*zin(  62)+xin(  83)*yin(  73)*zin(  86))
          eri_value(  130)=eri_value(  130)+d23bra( 44)*d01ket(  1)*(xin(   8)*yin(   5)*zin(  13)+xin(  32)*yin(  29)*zin(  37)+xin(  56)*yin(  53)*zin(  61)+xin(  80)*yin(  77)*zin(  85))
          eri_value(  131)=eri_value(  131)+d23bra( 44)*d01ket(  2)*(xin(   7)*yin(   6)*zin(  13)+xin(  31)*yin(  30)*zin(  37)+xin(  55)*yin(  54)*zin(  61)+xin(  79)*yin(  78)*zin(  85))
          eri_value(  132)=eri_value(  132)+d23bra( 44)*d01ket(  3)*(xin(   7)*yin(   5)*zin(  14)+xin(  31)*yin(  29)*zin(  38)+xin(  55)*yin(  53)*zin(  62)+xin(  79)*yin(  77)*zin(  86))
          eri_value(  133)=eri_value(  133)+d23bra( 45)*d01ket(  1)*(xin(   8)*yin(   1)*zin(  17)+xin(  32)*yin(  25)*zin(  41)+xin(  56)*yin(  49)*zin(  65)+xin(  80)*yin(  73)*zin(  89))
          eri_value(  134)=eri_value(  134)+d23bra( 45)*d01ket(  2)*(xin(   7)*yin(   2)*zin(  17)+xin(  31)*yin(  26)*zin(  41)+xin(  55)*yin(  50)*zin(  65)+xin(  79)*yin(  74)*zin(  89))
          eri_value(  135)=eri_value(  135)+d23bra( 45)*d01ket(  3)*(xin(   7)*yin(   1)*zin(  18)+xin(  31)*yin(  25)*zin(  42)+xin(  55)*yin(  49)*zin(  66)+xin(  79)*yin(  73)*zin(  90))
          eri_value(  136)=eri_value(  136)+d23bra( 46)*d01ket(  1)*(xin(  10)*yin(   3)*zin(  13)+xin(  34)*yin(  27)*zin(  37)+xin(  58)*yin(  51)*zin(  61)+xin(  82)*yin(  75)*zin(  85))
          eri_value(  137)=eri_value(  137)+d23bra( 46)*d01ket(  2)*(xin(   9)*yin(   4)*zin(  13)+xin(  33)*yin(  28)*zin(  37)+xin(  57)*yin(  52)*zin(  61)+xin(  81)*yin(  76)*zin(  85))
          eri_value(  138)=eri_value(  138)+d23bra( 46)*d01ket(  3)*(xin(   9)*yin(   3)*zin(  14)+xin(  33)*yin(  27)*zin(  38)+xin(  57)*yin(  51)*zin(  62)+xin(  81)*yin(  75)*zin(  86))
          eri_value(  139)=eri_value(  139)+d23bra( 47)*d01ket(  1)*(xin(  10)*yin(   1)*zin(  15)+xin(  34)*yin(  25)*zin(  39)+xin(  58)*yin(  49)*zin(  63)+xin(  82)*yin(  73)*zin(  87))
          eri_value(  140)=eri_value(  140)+d23bra( 47)*d01ket(  2)*(xin(   9)*yin(   2)*zin(  15)+xin(  33)*yin(  26)*zin(  39)+xin(  57)*yin(  50)*zin(  63)+xin(  81)*yin(  74)*zin(  87))
          eri_value(  141)=eri_value(  141)+d23bra( 47)*d01ket(  3)*(xin(   9)*yin(   1)*zin(  16)+xin(  33)*yin(  25)*zin(  40)+xin(  57)*yin(  49)*zin(  64)+xin(  81)*yin(  73)*zin(  88))
          eri_value(  142)=eri_value(  142)+d23bra( 48)*d01ket(  1)*(xin(   8)*yin(   3)*zin(  15)+xin(  32)*yin(  27)*zin(  39)+xin(  56)*yin(  51)*zin(  63)+xin(  80)*yin(  75)*zin(  87))
          eri_value(  143)=eri_value(  143)+d23bra( 48)*d01ket(  2)*(xin(   7)*yin(   4)*zin(  15)+xin(  31)*yin(  28)*zin(  39)+xin(  55)*yin(  52)*zin(  63)+xin(  79)*yin(  76)*zin(  87))
          eri_value(  144)=eri_value(  144)+d23bra( 48)*d01ket(  3)*(xin(   7)*yin(   3)*zin(  16)+xin(  31)*yin(  27)*zin(  40)+xin(  55)*yin(  51)*zin(  64)+xin(  79)*yin(  75)*zin(  88))
          eri_value(  145)=eri_value(  145)+d23bra( 49)*d01ket(  1)*(xin(   6)*yin(   7)*zin(  13)+xin(  30)*yin(  31)*zin(  37)+xin(  54)*yin(  55)*zin(  61)+xin(  78)*yin(  79)*zin(  85))
          eri_value(  146)=eri_value(  146)+d23bra( 49)*d01ket(  2)*(xin(   5)*yin(   8)*zin(  13)+xin(  29)*yin(  32)*zin(  37)+xin(  53)*yin(  56)*zin(  61)+xin(  77)*yin(  80)*zin(  85))
          eri_value(  147)=eri_value(  147)+d23bra( 49)*d01ket(  3)*(xin(   5)*yin(   7)*zin(  14)+xin(  29)*yin(  31)*zin(  38)+xin(  53)*yin(  55)*zin(  62)+xin(  77)*yin(  79)*zin(  86))
          eri_value(  148)=eri_value(  148)+d23bra( 50)*d01ket(  1)*(xin(   2)*yin(  11)*zin(  13)+xin(  26)*yin(  35)*zin(  37)+xin(  50)*yin(  59)*zin(  61)+xin(  74)*yin(  83)*zin(  85))
          eri_value(  149)=eri_value(  149)+d23bra( 50)*d01ket(  2)*(xin(   1)*yin(  12)*zin(  13)+xin(  25)*yin(  36)*zin(  37)+xin(  49)*yin(  60)*zin(  61)+xin(  73)*yin(  84)*zin(  85))
          eri_value(  150)=eri_value(  150)+d23bra( 50)*d01ket(  3)*(xin(   1)*yin(  11)*zin(  14)+xin(  25)*yin(  35)*zin(  38)+xin(  49)*yin(  59)*zin(  62)+xin(  73)*yin(  83)*zin(  86))
          eri_value(  151)=eri_value(  151)+d23bra( 51)*d01ket(  1)*(xin(   2)*yin(   7)*zin(  17)+xin(  26)*yin(  31)*zin(  41)+xin(  50)*yin(  55)*zin(  65)+xin(  74)*yin(  79)*zin(  89))
          eri_value(  152)=eri_value(  152)+d23bra( 51)*d01ket(  2)*(xin(   1)*yin(   8)*zin(  17)+xin(  25)*yin(  32)*zin(  41)+xin(  49)*yin(  56)*zin(  65)+xin(  73)*yin(  80)*zin(  89))
          eri_value(  153)=eri_value(  153)+d23bra( 51)*d01ket(  3)*(xin(   1)*yin(   7)*zin(  18)+xin(  25)*yin(  31)*zin(  42)+xin(  49)*yin(  55)*zin(  66)+xin(  73)*yin(  79)*zin(  90))
          eri_value(  154)=eri_value(  154)+d23bra( 52)*d01ket(  1)*(xin(   4)*yin(   9)*zin(  13)+xin(  28)*yin(  33)*zin(  37)+xin(  52)*yin(  57)*zin(  61)+xin(  76)*yin(  81)*zin(  85))
          eri_value(  155)=eri_value(  155)+d23bra( 52)*d01ket(  2)*(xin(   3)*yin(  10)*zin(  13)+xin(  27)*yin(  34)*zin(  37)+xin(  51)*yin(  58)*zin(  61)+xin(  75)*yin(  82)*zin(  85))
          eri_value(  156)=eri_value(  156)+d23bra( 52)*d01ket(  3)*(xin(   3)*yin(   9)*zin(  14)+xin(  27)*yin(  33)*zin(  38)+xin(  51)*yin(  57)*zin(  62)+xin(  75)*yin(  81)*zin(  86))
          eri_value(  157)=eri_value(  157)+d23bra( 53)*d01ket(  1)*(xin(   4)*yin(   7)*zin(  15)+xin(  28)*yin(  31)*zin(  39)+xin(  52)*yin(  55)*zin(  63)+xin(  76)*yin(  79)*zin(  87))
          eri_value(  158)=eri_value(  158)+d23bra( 53)*d01ket(  2)*(xin(   3)*yin(   8)*zin(  15)+xin(  27)*yin(  32)*zin(  39)+xin(  51)*yin(  56)*zin(  63)+xin(  75)*yin(  80)*zin(  87))
          eri_value(  159)=eri_value(  159)+d23bra( 53)*d01ket(  3)*(xin(   3)*yin(   7)*zin(  16)+xin(  27)*yin(  31)*zin(  40)+xin(  51)*yin(  55)*zin(  64)+xin(  75)*yin(  79)*zin(  88))
          eri_value(  160)=eri_value(  160)+d23bra( 54)*d01ket(  1)*(xin(   2)*yin(   9)*zin(  15)+xin(  26)*yin(  33)*zin(  39)+xin(  50)*yin(  57)*zin(  63)+xin(  74)*yin(  81)*zin(  87))
          eri_value(  161)=eri_value(  161)+d23bra( 54)*d01ket(  2)*(xin(   1)*yin(  10)*zin(  15)+xin(  25)*yin(  34)*zin(  39)+xin(  49)*yin(  58)*zin(  63)+xin(  73)*yin(  82)*zin(  87))
          eri_value(  162)=eri_value(  162)+d23bra( 54)*d01ket(  3)*(xin(   1)*yin(   9)*zin(  16)+xin(  25)*yin(  33)*zin(  40)+xin(  49)*yin(  57)*zin(  64)+xin(  73)*yin(  81)*zin(  88))
          eri_value(  163)=eri_value(  163)+d23bra( 55)*d01ket(  1)*(xin(  12)*yin(   7)*zin(   7)+xin(  36)*yin(  31)*zin(  31)+xin(  60)*yin(  55)*zin(  55)+xin(  84)*yin(  79)*zin(  79))
          eri_value(  164)=eri_value(  164)+d23bra( 55)*d01ket(  2)*(xin(  11)*yin(   8)*zin(   7)+xin(  35)*yin(  32)*zin(  31)+xin(  59)*yin(  56)*zin(  55)+xin(  83)*yin(  80)*zin(  79))
          eri_value(  165)=eri_value(  165)+d23bra( 55)*d01ket(  3)*(xin(  11)*yin(   7)*zin(   8)+xin(  35)*yin(  31)*zin(  32)+xin(  59)*yin(  55)*zin(  56)+xin(  83)*yin(  79)*zin(  80))
          eri_value(  166)=eri_value(  166)+d23bra( 56)*d01ket(  1)*(xin(   8)*yin(  11)*zin(   7)+xin(  32)*yin(  35)*zin(  31)+xin(  56)*yin(  59)*zin(  55)+xin(  80)*yin(  83)*zin(  79))
          eri_value(  167)=eri_value(  167)+d23bra( 56)*d01ket(  2)*(xin(   7)*yin(  12)*zin(   7)+xin(  31)*yin(  36)*zin(  31)+xin(  55)*yin(  60)*zin(  55)+xin(  79)*yin(  84)*zin(  79))
          eri_value(  168)=eri_value(  168)+d23bra( 56)*d01ket(  3)*(xin(   7)*yin(  11)*zin(   8)+xin(  31)*yin(  35)*zin(  32)+xin(  55)*yin(  59)*zin(  56)+xin(  79)*yin(  83)*zin(  80))
          eri_value(  169)=eri_value(  169)+d23bra( 57)*d01ket(  1)*(xin(   8)*yin(   7)*zin(  11)+xin(  32)*yin(  31)*zin(  35)+xin(  56)*yin(  55)*zin(  59)+xin(  80)*yin(  79)*zin(  83))
          eri_value(  170)=eri_value(  170)+d23bra( 57)*d01ket(  2)*(xin(   7)*yin(   8)*zin(  11)+xin(  31)*yin(  32)*zin(  35)+xin(  55)*yin(  56)*zin(  59)+xin(  79)*yin(  80)*zin(  83))
          eri_value(  171)=eri_value(  171)+d23bra( 57)*d01ket(  3)*(xin(   7)*yin(   7)*zin(  12)+xin(  31)*yin(  31)*zin(  36)+xin(  55)*yin(  55)*zin(  60)+xin(  79)*yin(  79)*zin(  84))
          eri_value(  172)=eri_value(  172)+d23bra( 58)*d01ket(  1)*(xin(  10)*yin(   9)*zin(   7)+xin(  34)*yin(  33)*zin(  31)+xin(  58)*yin(  57)*zin(  55)+xin(  82)*yin(  81)*zin(  79))
          eri_value(  173)=eri_value(  173)+d23bra( 58)*d01ket(  2)*(xin(   9)*yin(  10)*zin(   7)+xin(  33)*yin(  34)*zin(  31)+xin(  57)*yin(  58)*zin(  55)+xin(  81)*yin(  82)*zin(  79))
          eri_value(  174)=eri_value(  174)+d23bra( 58)*d01ket(  3)*(xin(   9)*yin(   9)*zin(   8)+xin(  33)*yin(  33)*zin(  32)+xin(  57)*yin(  57)*zin(  56)+xin(  81)*yin(  81)*zin(  80))
          eri_value(  175)=eri_value(  175)+d23bra( 59)*d01ket(  1)*(xin(  10)*yin(   7)*zin(   9)+xin(  34)*yin(  31)*zin(  33)+xin(  58)*yin(  55)*zin(  57)+xin(  82)*yin(  79)*zin(  81))
          eri_value(  176)=eri_value(  176)+d23bra( 59)*d01ket(  2)*(xin(   9)*yin(   8)*zin(   9)+xin(  33)*yin(  32)*zin(  33)+xin(  57)*yin(  56)*zin(  57)+xin(  81)*yin(  80)*zin(  81))
          eri_value(  177)=eri_value(  177)+d23bra( 59)*d01ket(  3)*(xin(   9)*yin(   7)*zin(  10)+xin(  33)*yin(  31)*zin(  34)+xin(  57)*yin(  55)*zin(  58)+xin(  81)*yin(  79)*zin(  82))
          eri_value(  178)=eri_value(  178)+d23bra( 60)*d01ket(  1)*(xin(   8)*yin(   9)*zin(   9)+xin(  32)*yin(  33)*zin(  33)+xin(  56)*yin(  57)*zin(  57)+xin(  80)*yin(  81)*zin(  81))
          eri_value(  179)=eri_value(  179)+d23bra( 60)*d01ket(  2)*(xin(   7)*yin(  10)*zin(   9)+xin(  31)*yin(  34)*zin(  33)+xin(  55)*yin(  58)*zin(  57)+xin(  79)*yin(  82)*zin(  81))
          eri_value(  180)=eri_value(  180)+d23bra( 60)*d01ket(  3)*(xin(   7)*yin(   9)*zin(  10)+xin(  31)*yin(  33)*zin(  34)+xin(  55)*yin(  57)*zin(  58)+xin(  79)*yin(  81)*zin(  82))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*18 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*3 + ip ! Add stride between functions in j

                                      do k = 1, 3 ! # of cartesians in k

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

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n01ket)
                              deallocate (xint01ket)

                              end subroutine int3210
                              end submodule
