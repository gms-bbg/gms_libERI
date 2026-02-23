! The total angular momentum of this class is:           7
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3220_impl
contains
  module subroutine int3220(df_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n02ket(:)
    real(dp), allocatable :: xint23bra(:), xint02ket(:)
    integer(kind=int64) :: ndfbra, nsdket
    real(dp) :: scutdfbra, scutsdket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(144), yin(144), zin(144)
    real(dp) :: eri_value(360)
    real(dp) :: d23bra(60), d02ket(6)
    integer(kind=int64) :: ix(10), jx(6), kx(6), lx(1)
    integer(kind=int64) :: iy(10), jy(6), ky(6), ly(1)
    integer(kind=int64) :: iz(10), jz(6), kz(6), lz(1)
    integer(kind=int64) :: in(6), in1(6), kn(3)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(6), kly(6), klz(6)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 10
    in1(3) = 19
    in1(4) = 28
    in1(5) = 31
    in1(6) = 34

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

    jx(1) = 6
    jx(2) = 0
    jx(3) = 0
    jx(4) = 3
    jx(5) = 3
    jx(6) = 0

    ix(1) = 28
    ix(2) = 1
    ix(3) = 1
    ix(4) = 19
    ix(5) = 19
    ix(6) = 10
    ix(7) = 1
    ix(8) = 10
    ix(9) = 1
    ix(10) = 10

    ! y-arrays

    ly(1) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 1

    jy(1) = 0
    jy(2) = 6
    jy(3) = 0
    jy(4) = 3
    jy(5) = 0
    jy(6) = 3

    iy(1) = 1
    iy(2) = 28
    iy(3) = 1
    iy(4) = 10
    iy(5) = 1
    iy(6) = 19
    iy(7) = 19
    iy(8) = 1
    iy(9) = 10
    iy(10) = 10

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
    jz(3) = 6
    jz(4) = 0
    jz(5) = 3
    jz(6) = 3

    iz(1) = 1
    iz(2) = 1
    iz(3) = 28
    iz(4) = 1
    iz(5) = 10
    iz(6) = 1
    iz(7) = 10
    iz(8) = 19
    iz(9) = 19
    iz(10) = 10

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 34
    ijx(2) = 28
    ijx(3) = 28
    ijx(4) = 31
    ijx(5) = 31
    ijx(6) = 28
    ijx(7) = 7
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 4
    ijx(11) = 4
    ijx(12) = 1
    ijx(13) = 7
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 4
    ijx(17) = 4
    ijx(18) = 1
    ijx(19) = 25
    ijx(20) = 19
    ijx(21) = 19
    ijx(22) = 22
    ijx(23) = 22
    ijx(24) = 19
    ijx(25) = 25
    ijx(26) = 19
    ijx(27) = 19
    ijx(28) = 22
    ijx(29) = 22
    ijx(30) = 19
    ijx(31) = 16
    ijx(32) = 10
    ijx(33) = 10
    ijx(34) = 13
    ijx(35) = 13
    ijx(36) = 10
    ijx(37) = 7
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 4
    ijx(41) = 4
    ijx(42) = 1
    ijx(43) = 16
    ijx(44) = 10
    ijx(45) = 10
    ijx(46) = 13
    ijx(47) = 13
    ijx(48) = 10
    ijx(49) = 7
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 4
    ijx(53) = 4
    ijx(54) = 1
    ijx(55) = 16
    ijx(56) = 10
    ijx(57) = 10
    ijx(58) = 13
    ijx(59) = 13
    ijx(60) = 10

    ijy(1) = 1
    ijy(2) = 7
    ijy(3) = 1
    ijy(4) = 4
    ijy(5) = 1
    ijy(6) = 4
    ijy(7) = 28
    ijy(8) = 34
    ijy(9) = 28
    ijy(10) = 31
    ijy(11) = 28
    ijy(12) = 31
    ijy(13) = 1
    ijy(14) = 7
    ijy(15) = 1
    ijy(16) = 4
    ijy(17) = 1
    ijy(18) = 4
    ijy(19) = 10
    ijy(20) = 16
    ijy(21) = 10
    ijy(22) = 13
    ijy(23) = 10
    ijy(24) = 13
    ijy(25) = 1
    ijy(26) = 7
    ijy(27) = 1
    ijy(28) = 4
    ijy(29) = 1
    ijy(30) = 4
    ijy(31) = 19
    ijy(32) = 25
    ijy(33) = 19
    ijy(34) = 22
    ijy(35) = 19
    ijy(36) = 22
    ijy(37) = 19
    ijy(38) = 25
    ijy(39) = 19
    ijy(40) = 22
    ijy(41) = 19
    ijy(42) = 22
    ijy(43) = 1
    ijy(44) = 7
    ijy(45) = 1
    ijy(46) = 4
    ijy(47) = 1
    ijy(48) = 4
    ijy(49) = 10
    ijy(50) = 16
    ijy(51) = 10
    ijy(52) = 13
    ijy(53) = 10
    ijy(54) = 13
    ijy(55) = 10
    ijy(56) = 16
    ijy(57) = 10
    ijy(58) = 13
    ijy(59) = 10
    ijy(60) = 13

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 7
    ijz(4) = 1
    ijz(5) = 4
    ijz(6) = 4
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 7
    ijz(10) = 1
    ijz(11) = 4
    ijz(12) = 4
    ijz(13) = 28
    ijz(14) = 28
    ijz(15) = 34
    ijz(16) = 28
    ijz(17) = 31
    ijz(18) = 31
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 7
    ijz(22) = 1
    ijz(23) = 4
    ijz(24) = 4
    ijz(25) = 10
    ijz(26) = 10
    ijz(27) = 16
    ijz(28) = 10
    ijz(29) = 13
    ijz(30) = 13
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 7
    ijz(34) = 1
    ijz(35) = 4
    ijz(36) = 4
    ijz(37) = 10
    ijz(38) = 10
    ijz(39) = 16
    ijz(40) = 10
    ijz(41) = 13
    ijz(42) = 13
    ijz(43) = 19
    ijz(44) = 19
    ijz(45) = 25
    ijz(46) = 19
    ijz(47) = 22
    ijz(48) = 22
    ijz(49) = 19
    ijz(50) = 19
    ijz(51) = 25
    ijz(52) = 19
    ijz(53) = 22
    ijz(54) = 22
    ijz(55) = 10
    ijz(56) = 10
    ijz(57) = 16
    ijz(58) = 10
    ijz(59) = 13
    ijz(60) = 13

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

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))

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

    if ((ndfbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nsdket
    ntile = int(ndfbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nsdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint02ket, n02ket, df_pair, sd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d02ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint02ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n02ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
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

                                      ! i2 = in(2) =   10
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(10) = xc00
                                      yin(10) = yc00
                                      zin(10) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   11
                                      ! i2 =   10

                                      xin(11) = xcp00*xin(10) + cp10
                                      yin(11) = ycp00*yin(10) + cp10
                                      zin(11) = zcp00*zin(10) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   10

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   19
                                      ! i3 =    1
                                      ! i4 =   10

                                      xin(19) = c10*xin(1) + xc00*xin(10)
                                      yin(19) = c10*yin(1) + yc00*yin(10)
                                      zin(19) = c10*zin(1) + zc00*zin(10)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   20
                                      ! i5 =   19
                                      ! i4 =   10

                                      xin(20) = xcp00*xin(19) + cp10*xin(10)
                                      yin(20) = ycp00*yin(19) + cp10*yin(10)
                                      zin(20) = zcp00*zin(19) + cp10*zin(10)

                                      ! ------------------

                                      ! i3 = i4 =   10
                                      ! i4 = i5 =   19

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   28
                                      ! i3 =   10
                                      ! i4 =   19

                                      xin(28) = c10*xin(10) + xc00*xin(19)
                                      yin(28) = c10*yin(10) + yc00*yin(19)
                                      zin(28) = c10*zin(10) + zc00*zin(19)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   29
                                      ! i5 =   28
                                      ! i4 =   19

                                      xin(29) = xcp00*xin(28) + cp10*xin(19)
                                      yin(29) = ycp00*yin(28) + cp10*yin(19)
                                      zin(29) = zcp00*zin(28) + cp10*zin(19)

                                      ! ------------------

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   28

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   31
                                      ! i3 =   19
                                      ! i4 =   28

                                      xin(31) = c10*xin(19) + xc00*xin(28)
                                      yin(31) = c10*yin(19) + yc00*yin(28)
                                      zin(31) = c10*zin(19) + zc00*zin(28)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   32
                                      ! i5 =   31
                                      ! i4 =   28

                                      xin(32) = xcp00*xin(31) + cp10*xin(28)
                                      yin(32) = ycp00*yin(31) + cp10*yin(28)
                                      zin(32) = zcp00*zin(31) + cp10*zin(28)

                                      ! ------------------

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   31

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   34
                                      ! i3 =   28
                                      ! i4 =   31

                                      xin(34) = c10*xin(28) + xc00*xin(31)
                                      yin(34) = c10*yin(28) + yc00*yin(31)
                                      zin(34) = c10*zin(28) + zc00*zin(31)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   35
                                      ! i5 =   34
                                      ! i4 =   31

                                      xin(35) = xcp00*xin(34) + cp10*xin(31)
                                      yin(35) = ycp00*yin(34) + cp10*yin(31)
                                      zin(35) = zcp00*zin(34) + cp10*zin(31)

                                      ! ------------------

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   34

                                      ! n =    6

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

                                      ! i3 = i2+kn(n+1) =   12

                                      xin(12) = xc00*xin(3) + c01*xin(2)
                                      yin(12) = yc00*yin(3) + c01*yin(2)
                                      zin(12) = zc00*zin(3) + c01*zin(2)

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
                                      ! i4 = i2 =   10

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   19

                                      xin(21) = c10*xin(3) + xc00*xin(12) + c01*xin(11)
                                      yin(21) = c10*yin(3) + yc00*yin(12) + c01*yin(11)
                                      zin(21) = c10*zin(3) + zc00*zin(12) + c01*zin(11)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   10
                                      ! i4 = i5 =   19

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   28

                                      xin(30) = c10*xin(12) + xc00*xin(21) + c01*xin(20)
                                      yin(30) = c10*yin(12) + yc00*yin(21) + c01*yin(20)
                                      zin(30) = c10*zin(12) + zc00*zin(21) + c01*zin(20)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   28

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   31

                                      xin(33) = c10*xin(21) + xc00*xin(30) + c01*xin(29)
                                      yin(33) = c10*yin(21) + yc00*yin(30) + c01*yin(29)
                                      zin(33) = c10*zin(21) + zc00*zin(30) + c01*zin(29)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   31

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   34

                                      xin(36) = c10*xin(30) + xc00*xin(33) + c01*xin(32)
                                      yin(36) = c10*yin(30) + yc00*yin(33) + c01*yin(32)
                                      zin(36) = c10*zin(30) + zc00*zin(33) + c01*zin(32)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   34

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   34

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   34

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(34) = xin(34) + dxij*xin(31)
                                      yin(34) = yin(34) + dyij*yin(31)
                                      zin(34) = zin(34) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   28

                                      xin(31) = xin(31) + dxij*xin(28)
                                      yin(31) = yin(31) + dyij*yin(28)
                                      zin(31) = zin(31) + dzij*zin(28)

                                      ! i3 = i4 =   28
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   34

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(34) = xin(34) + dxij*xin(31)
                                      yin(34) = yin(34) + dyij*yin(31)
                                      zin(34) = zin(34) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    4

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    4

                                      ! do ni = 1,    3

                                      xin(4) = xin(10) + dxij*xin(1)
                                      yin(4) = yin(10) + dyij*yin(1)
                                      zin(4) = zin(10) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   13

                                      ! ni =    2

                                      xin(13) = xin(19) + dxij*xin(10)
                                      yin(13) = yin(19) + dyij*yin(10)
                                      zin(13) = zin(19) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    3

                                      xin(22) = xin(28) + dxij*xin(19)
                                      yin(22) = yin(28) + dyij*yin(19)
                                      zin(22) = zin(28) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    7

                                      ! nj =    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(13) + dxij*xin(4)
                                      yin(7) = yin(13) + dyij*yin(4)
                                      zin(7) = zin(13) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    2

                                      xin(16) = xin(22) + dxij*xin(13)
                                      yin(16) = yin(22) + dyij*yin(13)
                                      zin(16) = zin(22) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! ni =    3

                                      xin(25) = xin(31) + dxij*xin(22)
                                      yin(25) = yin(31) + dyij*yin(22)
                                      zin(25) = zin(31) + dzij*zin(22)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   35

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   32

                                      xin(35) = xin(35) + dxij*xin(32)
                                      yin(35) = yin(35) + dyij*yin(32)
                                      zin(35) = zin(35) + dzij*zin(32)

                                      ! i3 = i4 =   32
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   29

                                      xin(32) = xin(32) + dxij*xin(29)
                                      yin(32) = yin(32) + dyij*yin(29)
                                      zin(32) = zin(32) + dzij*zin(29)

                                      ! i3 = i4 =   29
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   35

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   32

                                      xin(35) = xin(35) + dxij*xin(32)
                                      yin(35) = yin(35) + dyij*yin(32)
                                      zin(35) = zin(35) + dzij*zin(32)

                                      ! i3 = i4 =   32
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(11) + dxij*xin(2)
                                      yin(5) = yin(11) + dyij*yin(2)
                                      zin(5) = zin(11) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   14

                                      ! ni =    2

                                      xin(14) = xin(20) + dxij*xin(11)
                                      yin(14) = yin(20) + dyij*yin(11)
                                      zin(14) = zin(20) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    3

                                      xin(23) = xin(29) + dxij*xin(20)
                                      yin(23) = yin(29) + dyij*yin(20)
                                      zin(23) = zin(29) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    8

                                      ! nj =    2

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

                                      xin(8) = xin(14) + dxij*xin(5)
                                      yin(8) = yin(14) + dyij*yin(5)
                                      zin(8) = zin(14) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    2

                                      xin(17) = xin(23) + dxij*xin(14)
                                      yin(17) = yin(23) + dyij*yin(14)
                                      zin(17) = zin(23) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   26

                                      ! ni =    3

                                      xin(26) = xin(32) + dxij*xin(23)
                                      yin(26) = yin(32) + dyij*yin(23)
                                      zin(26) = zin(32) + dzij*zin(23)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   33

                                      xin(36) = xin(36) + dxij*xin(33)
                                      yin(36) = yin(36) + dyij*yin(33)
                                      zin(36) = zin(36) + dzij*zin(33)

                                      ! i3 = i4 =   33
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   30

                                      xin(33) = xin(33) + dxij*xin(30)
                                      yin(33) = yin(33) + dyij*yin(30)
                                      zin(33) = zin(33) + dzij*zin(30)

                                      ! i3 = i4 =   30
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   33

                                      xin(36) = xin(36) + dxij*xin(33)
                                      yin(36) = yin(36) + dyij*yin(33)
                                      zin(36) = zin(36) + dzij*zin(33)

                                      ! i3 = i4 =   33
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    3

                                      xin(6) = xin(12) + dxij*xin(3)
                                      yin(6) = yin(12) + dyij*yin(3)
                                      zin(6) = zin(12) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   15

                                      ! ni =    2

                                      xin(15) = xin(21) + dxij*xin(12)
                                      yin(15) = yin(21) + dyij*yin(12)
                                      zin(15) = zin(21) + dzij*zin(12)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    3

                                      xin(24) = xin(30) + dxij*xin(21)
                                      yin(24) = yin(30) + dyij*yin(21)
                                      zin(24) = zin(30) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(15) + dxij*xin(6)
                                      yin(9) = yin(15) + dyij*yin(6)
                                      zin(9) = zin(15) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    2

                                      xin(18) = xin(24) + dxij*xin(15)
                                      yin(18) = yin(24) + dyij*yin(15)
                                      zin(18) = zin(24) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   27

                                      ! ni =    3

                                      xin(27) = xin(33) + dxij*xin(24)
                                      yin(27) = yin(33) + dyij*yin(24)
                                      zin(27) = zin(33) + dzij*zin(24)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   36

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

                                      ! i1 = in(1) =   37

                                      xin(37) = 1.0_dp
                                      yin(37) = 1.0_dp
                                      zin(37) = f00

                                      ! i2 = in(2) =   46
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(46) = xc00
                                      yin(46) = yc00
                                      zin(46) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   38

                                      xin(38) = xcp00
                                      yin(38) = ycp00
                                      zin(38) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   47
                                      ! i2 =   46

                                      xin(47) = xcp00*xin(46) + cp10
                                      yin(47) = ycp00*yin(46) + cp10
                                      zin(47) = zcp00*zin(46) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   46

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   55
                                      ! i3 =   37
                                      ! i4 =   46

                                      xin(55) = c10*xin(37) + xc00*xin(46)
                                      yin(55) = c10*yin(37) + yc00*yin(46)
                                      zin(55) = c10*zin(37) + zc00*zin(46)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   56
                                      ! i5 =   55
                                      ! i4 =   46

                                      xin(56) = xcp00*xin(55) + cp10*xin(46)
                                      yin(56) = ycp00*yin(55) + cp10*yin(46)
                                      zin(56) = zcp00*zin(55) + cp10*zin(46)

                                      ! ------------------

                                      ! i3 = i4 =   46
                                      ! i4 = i5 =   55

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   64
                                      ! i3 =   46
                                      ! i4 =   55

                                      xin(64) = c10*xin(46) + xc00*xin(55)
                                      yin(64) = c10*yin(46) + yc00*yin(55)
                                      zin(64) = c10*zin(46) + zc00*zin(55)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   65
                                      ! i5 =   64
                                      ! i4 =   55

                                      xin(65) = xcp00*xin(64) + cp10*xin(55)
                                      yin(65) = ycp00*yin(64) + cp10*yin(55)
                                      zin(65) = zcp00*zin(64) + cp10*zin(55)

                                      ! ------------------

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   64

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   67
                                      ! i3 =   55
                                      ! i4 =   64

                                      xin(67) = c10*xin(55) + xc00*xin(64)
                                      yin(67) = c10*yin(55) + yc00*yin(64)
                                      zin(67) = c10*zin(55) + zc00*zin(64)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   68
                                      ! i5 =   67
                                      ! i4 =   64

                                      xin(68) = xcp00*xin(67) + cp10*xin(64)
                                      yin(68) = ycp00*yin(67) + cp10*yin(64)
                                      zin(68) = zcp00*zin(67) + cp10*zin(64)

                                      ! ------------------

                                      ! i3 = i4 =   64
                                      ! i4 = i5 =   67

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   70
                                      ! i3 =   64
                                      ! i4 =   67

                                      xin(70) = c10*xin(64) + xc00*xin(67)
                                      yin(70) = c10*yin(64) + yc00*yin(67)
                                      zin(70) = c10*zin(64) + zc00*zin(67)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   71
                                      ! i5 =   70
                                      ! i4 =   67

                                      xin(71) = xcp00*xin(70) + cp10*xin(67)
                                      yin(71) = ycp00*yin(70) + cp10*yin(67)
                                      zin(71) = zcp00*zin(70) + cp10*zin(67)

                                      ! ------------------

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   70

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   37
                                      ! i4 = i1+k2 =   38

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   39
                                      ! i3 =   37
                                      ! i4 =   38

                                      xin(39) = cp01*xin(37) + xcp00*xin(38)
                                      yin(39) = cp01*yin(37) + ycp00*yin(38)
                                      zin(39) = cp01*zin(37) + zcp00*zin(38)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   48

                                      xin(48) = xc00*xin(39) + c01*xin(38)
                                      yin(48) = yc00*yin(39) + c01*yin(38)
                                      zin(48) = zc00*zin(39) + c01*zin(38)

                                      ! ------------------

                                      ! i3 = i4 =   38
                                      ! i4 = i5 =   39

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   46

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   55

                                      xin(57) = c10*xin(39) + xc00*xin(48) + c01*xin(47)
                                      yin(57) = c10*yin(39) + yc00*yin(48) + c01*yin(47)
                                      zin(57) = c10*zin(39) + zc00*zin(48) + c01*zin(47)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   46
                                      ! i4 = i5 =   55

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   64

                                      xin(66) = c10*xin(48) + xc00*xin(57) + c01*xin(56)
                                      yin(66) = c10*yin(48) + yc00*yin(57) + c01*yin(56)
                                      zin(66) = c10*zin(48) + zc00*zin(57) + c01*zin(56)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   64

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   67

                                      xin(69) = c10*xin(57) + xc00*xin(66) + c01*xin(65)
                                      yin(69) = c10*yin(57) + yc00*yin(66) + c01*yin(65)
                                      zin(69) = c10*zin(57) + zc00*zin(66) + c01*zin(65)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   64
                                      ! i4 = i5 =   67

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   70

                                      xin(72) = c10*xin(66) + xc00*xin(69) + c01*xin(68)
                                      yin(72) = c10*yin(66) + yc00*yin(69) + c01*yin(68)
                                      zin(72) = c10*zin(66) + zc00*zin(69) + c01*zin(68)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   70

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   70

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   70

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   67

                                      xin(70) = xin(70) + dxij*xin(67)
                                      yin(70) = yin(70) + dyij*yin(67)
                                      zin(70) = zin(70) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   64

                                      xin(67) = xin(67) + dxij*xin(64)
                                      yin(67) = yin(67) + dyij*yin(64)
                                      zin(67) = zin(67) + dzij*zin(64)

                                      ! i3 = i4 =   64
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   70

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   67

                                      xin(70) = xin(70) + dxij*xin(67)
                                      yin(70) = yin(70) + dyij*yin(67)
                                      zin(70) = zin(70) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   40

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   40

                                      ! do ni = 1,    3

                                      xin(40) = xin(46) + dxij*xin(37)
                                      yin(40) = yin(46) + dyij*yin(37)
                                      zin(40) = zin(46) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! ni =    2

                                      xin(49) = xin(55) + dxij*xin(46)
                                      yin(49) = yin(55) + dyij*yin(46)
                                      zin(49) = zin(55) + dzij*zin(46)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    3

                                      xin(58) = xin(64) + dxij*xin(55)
                                      yin(58) = yin(64) + dyij*yin(55)
                                      zin(58) = zin(64) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   43

                                      ! nj =    2

                                      ! i4 = i3 =   43

                                      ! do ni = 1,    3

                                      xin(43) = xin(49) + dxij*xin(40)
                                      yin(43) = yin(49) + dyij*yin(40)
                                      zin(43) = zin(49) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   52

                                      ! ni =    2

                                      xin(52) = xin(58) + dxij*xin(49)
                                      yin(52) = yin(58) + dyij*yin(49)
                                      zin(52) = zin(58) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    3

                                      xin(61) = xin(67) + dxij*xin(58)
                                      yin(61) = yin(67) + dyij*yin(58)
                                      zin(61) = zin(67) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   46

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   68

                                      xin(71) = xin(71) + dxij*xin(68)
                                      yin(71) = yin(71) + dyij*yin(68)
                                      zin(71) = zin(71) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   65

                                      xin(68) = xin(68) + dxij*xin(65)
                                      yin(68) = yin(68) + dyij*yin(65)
                                      zin(68) = zin(68) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   68

                                      xin(71) = xin(71) + dxij*xin(68)
                                      yin(71) = yin(71) + dyij*yin(68)
                                      zin(71) = zin(71) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   41

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   41

                                      ! do ni = 1,    3

                                      xin(41) = xin(47) + dxij*xin(38)
                                      yin(41) = yin(47) + dyij*yin(38)
                                      zin(41) = zin(47) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   50

                                      ! ni =    2

                                      xin(50) = xin(56) + dxij*xin(47)
                                      yin(50) = yin(56) + dyij*yin(47)
                                      zin(50) = zin(56) + dzij*zin(47)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    3

                                      xin(59) = xin(65) + dxij*xin(56)
                                      yin(59) = yin(65) + dyij*yin(56)
                                      zin(59) = zin(65) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   68

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   44

                                      ! nj =    2

                                      ! i4 = i3 =   44

                                      ! do ni = 1,    3

                                      xin(44) = xin(50) + dxij*xin(41)
                                      yin(44) = yin(50) + dyij*yin(41)
                                      zin(44) = zin(50) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    2

                                      xin(53) = xin(59) + dxij*xin(50)
                                      yin(53) = yin(59) + dyij*yin(50)
                                      zin(53) = zin(59) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   62

                                      ! ni =    3

                                      xin(62) = xin(68) + dxij*xin(59)
                                      yin(62) = yin(68) + dyij*yin(59)
                                      zin(62) = zin(68) + dzij*zin(59)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   47

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   69

                                      xin(72) = xin(72) + dxij*xin(69)
                                      yin(72) = yin(72) + dyij*yin(69)
                                      zin(72) = zin(72) + dzij*zin(69)

                                      ! i3 = i4 =   69
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   66

                                      xin(69) = xin(69) + dxij*xin(66)
                                      yin(69) = yin(69) + dyij*yin(66)
                                      zin(69) = zin(69) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   69

                                      xin(72) = xin(72) + dxij*xin(69)
                                      yin(72) = yin(72) + dyij*yin(69)
                                      zin(72) = zin(72) + dzij*zin(69)

                                      ! i3 = i4 =   69
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   42

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   42

                                      ! do ni = 1,    3

                                      xin(42) = xin(48) + dxij*xin(39)
                                      yin(42) = yin(48) + dyij*yin(39)
                                      zin(42) = zin(48) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   51

                                      ! ni =    2

                                      xin(51) = xin(57) + dxij*xin(48)
                                      yin(51) = yin(57) + dyij*yin(48)
                                      zin(51) = zin(57) + dzij*zin(48)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    3

                                      xin(60) = xin(66) + dxij*xin(57)
                                      yin(60) = yin(66) + dyij*yin(57)
                                      zin(60) = zin(66) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   45

                                      ! nj =    2

                                      ! i4 = i3 =   45

                                      ! do ni = 1,    3

                                      xin(45) = xin(51) + dxij*xin(42)
                                      yin(45) = yin(51) + dyij*yin(42)
                                      zin(45) = zin(51) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    2

                                      xin(54) = xin(60) + dxij*xin(51)
                                      yin(54) = yin(60) + dyij*yin(51)
                                      zin(54) = zin(60) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    3

                                      xin(63) = xin(69) + dxij*xin(60)
                                      yin(63) = yin(69) + dyij*yin(60)
                                      zin(63) = zin(69) + dzij*zin(60)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   48

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   72

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

                                      ! i1 = in(1) =   73

                                      xin(73) = 1.0_dp
                                      yin(73) = 1.0_dp
                                      zin(73) = f00

                                      ! i2 = in(2) =   82
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(82) = xc00
                                      yin(82) = yc00
                                      zin(82) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   74

                                      xin(74) = xcp00
                                      yin(74) = ycp00
                                      zin(74) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   83
                                      ! i2 =   82

                                      xin(83) = xcp00*xin(82) + cp10
                                      yin(83) = ycp00*yin(82) + cp10
                                      zin(83) = zcp00*zin(82) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   82

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   73
                                      ! i4 =   82

                                      xin(91) = c10*xin(73) + xc00*xin(82)
                                      yin(91) = c10*yin(73) + yc00*yin(82)
                                      zin(91) = c10*zin(73) + zc00*zin(82)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   92
                                      ! i5 =   91
                                      ! i4 =   82

                                      xin(92) = xcp00*xin(91) + cp10*xin(82)
                                      yin(92) = ycp00*yin(91) + cp10*yin(82)
                                      zin(92) = zcp00*zin(91) + cp10*zin(82)

                                      ! ------------------

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  100
                                      ! i3 =   82
                                      ! i4 =   91

                                      xin(100) = c10*xin(82) + xc00*xin(91)
                                      yin(100) = c10*yin(82) + yc00*yin(91)
                                      zin(100) = c10*zin(82) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  101
                                      ! i5 =  100
                                      ! i4 =   91

                                      xin(101) = xcp00*xin(100) + cp10*xin(91)
                                      yin(101) = ycp00*yin(100) + cp10*yin(91)
                                      zin(101) = zcp00*zin(100) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  103
                                      ! i3 =   91
                                      ! i4 =  100

                                      xin(103) = c10*xin(91) + xc00*xin(100)
                                      yin(103) = c10*yin(91) + yc00*yin(100)
                                      zin(103) = c10*zin(91) + zc00*zin(100)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  104
                                      ! i5 =  103
                                      ! i4 =  100

                                      xin(104) = xcp00*xin(103) + cp10*xin(100)
                                      yin(104) = ycp00*yin(103) + cp10*yin(100)
                                      zin(104) = zcp00*zin(103) + cp10*zin(100)

                                      ! ------------------

                                      ! i3 = i4 =  100
                                      ! i4 = i5 =  103

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  106
                                      ! i3 =  100
                                      ! i4 =  103

                                      xin(106) = c10*xin(100) + xc00*xin(103)
                                      yin(106) = c10*yin(100) + yc00*yin(103)
                                      zin(106) = c10*zin(100) + zc00*zin(103)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  107
                                      ! i5 =  106
                                      ! i4 =  103

                                      xin(107) = xcp00*xin(106) + cp10*xin(103)
                                      yin(107) = ycp00*yin(106) + cp10*yin(103)
                                      zin(107) = zcp00*zin(106) + cp10*zin(103)

                                      ! ------------------

                                      ! i3 = i4 =  103
                                      ! i4 = i5 =  106

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   73
                                      ! i4 = i1+k2 =   74

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   75
                                      ! i3 =   73
                                      ! i4 =   74

                                      xin(75) = cp01*xin(73) + xcp00*xin(74)
                                      yin(75) = cp01*yin(73) + ycp00*yin(74)
                                      zin(75) = cp01*zin(73) + zcp00*zin(74)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   84

                                      xin(84) = xc00*xin(75) + c01*xin(74)
                                      yin(84) = yc00*yin(75) + c01*yin(74)
                                      zin(84) = zc00*zin(75) + c01*zin(74)

                                      ! ------------------

                                      ! i3 = i4 =   74
                                      ! i4 = i5 =   75

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   82

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   91

                                      xin(93) = c10*xin(75) + xc00*xin(84) + c01*xin(83)
                                      yin(93) = c10*yin(75) + yc00*yin(84) + c01*yin(83)
                                      zin(93) = c10*zin(75) + zc00*zin(84) + c01*zin(83)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  100

                                      xin(102) = c10*xin(84) + xc00*xin(93) + c01*xin(92)
                                      yin(102) = c10*yin(84) + yc00*yin(93) + c01*yin(92)
                                      zin(102) = c10*zin(84) + zc00*zin(93) + c01*zin(92)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  103

                                      xin(105) = c10*xin(93) + xc00*xin(102) + c01*xin(101)
                                      yin(105) = c10*yin(93) + yc00*yin(102) + c01*yin(101)
                                      zin(105) = c10*zin(93) + zc00*zin(102) + c01*zin(101)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  100
                                      ! i4 = i5 =  103

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  106

                                      xin(108) = c10*xin(102) + xc00*xin(105) + c01*xin(104)
                                      yin(108) = c10*yin(102) + yc00*yin(105) + c01*yin(104)
                                      zin(108) = c10*zin(102) + zc00*zin(105) + c01*zin(104)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  103
                                      ! i4 = i5 =  106

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  106

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  103

                                      xin(106) = xin(106) + dxij*xin(103)
                                      yin(106) = yin(106) + dyij*yin(103)
                                      zin(106) = zin(106) + dzij*zin(103)

                                      ! i3 = i4 =  103
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  100

                                      xin(103) = xin(103) + dxij*xin(100)
                                      yin(103) = yin(103) + dyij*yin(100)
                                      zin(103) = zin(103) + dzij*zin(100)

                                      ! i3 = i4 =  100
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  103

                                      xin(106) = xin(106) + dxij*xin(103)
                                      yin(106) = yin(106) + dyij*yin(103)
                                      zin(106) = zin(106) + dzij*zin(103)

                                      ! i3 = i4 =  103
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   76

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   76

                                      ! do ni = 1,    3

                                      xin(76) = xin(82) + dxij*xin(73)
                                      yin(76) = yin(82) + dyij*yin(73)
                                      zin(76) = zin(82) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    2

                                      xin(85) = xin(91) + dxij*xin(82)
                                      yin(85) = yin(91) + dyij*yin(82)
                                      zin(85) = zin(91) + dzij*zin(82)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    3

                                      xin(94) = xin(100) + dxij*xin(91)
                                      yin(94) = yin(100) + dyij*yin(91)
                                      zin(94) = zin(100) + dzij*zin(91)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   79

                                      ! nj =    2

                                      ! i4 = i3 =   79

                                      ! do ni = 1,    3

                                      xin(79) = xin(85) + dxij*xin(76)
                                      yin(79) = yin(85) + dyij*yin(76)
                                      zin(79) = zin(85) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    2

                                      xin(88) = xin(94) + dxij*xin(85)
                                      yin(88) = yin(94) + dyij*yin(85)
                                      zin(88) = zin(94) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! ni =    3

                                      xin(97) = xin(103) + dxij*xin(94)
                                      yin(97) = yin(103) + dyij*yin(94)
                                      zin(97) = zin(103) + dzij*zin(94)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   82

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  104

                                      xin(107) = xin(107) + dxij*xin(104)
                                      yin(107) = yin(107) + dyij*yin(104)
                                      zin(107) = zin(107) + dzij*zin(104)

                                      ! i3 = i4 =  104
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  101

                                      xin(104) = xin(104) + dxij*xin(101)
                                      yin(104) = yin(104) + dyij*yin(101)
                                      zin(104) = zin(104) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  104

                                      xin(107) = xin(107) + dxij*xin(104)
                                      yin(107) = yin(107) + dyij*yin(104)
                                      zin(107) = zin(107) + dzij*zin(104)

                                      ! i3 = i4 =  104
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   77

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   77

                                      ! do ni = 1,    3

                                      xin(77) = xin(83) + dxij*xin(74)
                                      yin(77) = yin(83) + dyij*yin(74)
                                      zin(77) = zin(83) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   86

                                      ! ni =    2

                                      xin(86) = xin(92) + dxij*xin(83)
                                      yin(86) = yin(92) + dyij*yin(83)
                                      zin(86) = zin(92) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    3

                                      xin(95) = xin(101) + dxij*xin(92)
                                      yin(95) = yin(101) + dyij*yin(92)
                                      zin(95) = zin(101) + dzij*zin(92)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  104

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   80

                                      ! nj =    2

                                      ! i4 = i3 =   80

                                      ! do ni = 1,    3

                                      xin(80) = xin(86) + dxij*xin(77)
                                      yin(80) = yin(86) + dyij*yin(77)
                                      zin(80) = zin(86) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    2

                                      xin(89) = xin(95) + dxij*xin(86)
                                      yin(89) = yin(95) + dyij*yin(86)
                                      zin(89) = zin(95) + dzij*zin(86)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   98

                                      ! ni =    3

                                      xin(98) = xin(104) + dxij*xin(95)
                                      yin(98) = yin(104) + dyij*yin(95)
                                      zin(98) = zin(104) + dzij*zin(95)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   83

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  105

                                      xin(108) = xin(108) + dxij*xin(105)
                                      yin(108) = yin(108) + dyij*yin(105)
                                      zin(108) = zin(108) + dzij*zin(105)

                                      ! i3 = i4 =  105
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  102

                                      xin(105) = xin(105) + dxij*xin(102)
                                      yin(105) = yin(105) + dyij*yin(102)
                                      zin(105) = zin(105) + dzij*zin(102)

                                      ! i3 = i4 =  102
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  105

                                      xin(108) = xin(108) + dxij*xin(105)
                                      yin(108) = yin(108) + dyij*yin(105)
                                      zin(108) = zin(108) + dzij*zin(105)

                                      ! i3 = i4 =  105
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   78

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   78

                                      ! do ni = 1,    3

                                      xin(78) = xin(84) + dxij*xin(75)
                                      yin(78) = yin(84) + dyij*yin(75)
                                      zin(78) = zin(84) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    2

                                      xin(87) = xin(93) + dxij*xin(84)
                                      yin(87) = yin(93) + dyij*yin(84)
                                      zin(87) = zin(93) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    3

                                      xin(96) = xin(102) + dxij*xin(93)
                                      yin(96) = yin(102) + dyij*yin(93)
                                      zin(96) = zin(102) + dzij*zin(93)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   81

                                      ! nj =    2

                                      ! i4 = i3 =   81

                                      ! do ni = 1,    3

                                      xin(81) = xin(87) + dxij*xin(78)
                                      yin(81) = yin(87) + dyij*yin(78)
                                      zin(81) = zin(87) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    2

                                      xin(90) = xin(96) + dxij*xin(87)
                                      yin(90) = yin(96) + dyij*yin(87)
                                      zin(90) = zin(96) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   99

                                      ! ni =    3

                                      xin(99) = xin(105) + dxij*xin(96)
                                      yin(99) = yin(105) + dyij*yin(96)
                                      zin(99) = zin(105) + dzij*zin(96)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   84

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  108

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

                                      ! i1 = in(1) =  109

                                      xin(109) = 1.0_dp
                                      yin(109) = 1.0_dp
                                      zin(109) = f00

                                      ! i2 = in(2) =  118
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(118) = xc00
                                      yin(118) = yc00
                                      zin(118) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  110

                                      xin(110) = xcp00
                                      yin(110) = ycp00
                                      zin(110) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  119
                                      ! i2 =  118

                                      xin(119) = xcp00*xin(118) + cp10
                                      yin(119) = ycp00*yin(118) + cp10
                                      zin(119) = zcp00*zin(118) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  118

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  127
                                      ! i3 =  109
                                      ! i4 =  118

                                      xin(127) = c10*xin(109) + xc00*xin(118)
                                      yin(127) = c10*yin(109) + yc00*yin(118)
                                      zin(127) = c10*zin(109) + zc00*zin(118)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  128
                                      ! i5 =  127
                                      ! i4 =  118

                                      xin(128) = xcp00*xin(127) + cp10*xin(118)
                                      yin(128) = ycp00*yin(127) + cp10*yin(118)
                                      zin(128) = zcp00*zin(127) + cp10*zin(118)

                                      ! ------------------

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  136
                                      ! i3 =  118
                                      ! i4 =  127

                                      xin(136) = c10*xin(118) + xc00*xin(127)
                                      yin(136) = c10*yin(118) + yc00*yin(127)
                                      zin(136) = c10*zin(118) + zc00*zin(127)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  137
                                      ! i5 =  136
                                      ! i4 =  127

                                      xin(137) = xcp00*xin(136) + cp10*xin(127)
                                      yin(137) = ycp00*yin(136) + cp10*yin(127)
                                      zin(137) = zcp00*zin(136) + cp10*zin(127)

                                      ! ------------------

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  139
                                      ! i3 =  127
                                      ! i4 =  136

                                      xin(139) = c10*xin(127) + xc00*xin(136)
                                      yin(139) = c10*yin(127) + yc00*yin(136)
                                      zin(139) = c10*zin(127) + zc00*zin(136)

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

                                      ! n =    5

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

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  109
                                      ! i4 = i1+k2 =  110

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  111
                                      ! i3 =  109
                                      ! i4 =  110

                                      xin(111) = cp01*xin(109) + xcp00*xin(110)
                                      yin(111) = cp01*yin(109) + ycp00*yin(110)
                                      zin(111) = cp01*zin(109) + zcp00*zin(110)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  120

                                      xin(120) = xc00*xin(111) + c01*xin(110)
                                      yin(120) = yc00*yin(111) + c01*yin(110)
                                      zin(120) = zc00*zin(111) + c01*zin(110)

                                      ! ------------------

                                      ! i3 = i4 =  110
                                      ! i4 = i5 =  111

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  118

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  127

                                      xin(129) = c10*xin(111) + xc00*xin(120) + c01*xin(119)
                                      yin(129) = c10*yin(111) + yc00*yin(120) + c01*yin(119)
                                      zin(129) = c10*zin(111) + zc00*zin(120) + c01*zin(119)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  118
                                      ! i4 = i5 =  127

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  136

                                      xin(138) = c10*xin(120) + xc00*xin(129) + c01*xin(128)
                                      yin(138) = c10*yin(120) + yc00*yin(129) + c01*yin(128)
                                      zin(138) = c10*zin(120) + zc00*zin(129) + c01*zin(128)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  136

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  139

                                      xin(141) = c10*xin(129) + xc00*xin(138) + c01*xin(137)
                                      yin(141) = c10*yin(129) + yc00*yin(138) + c01*yin(137)
                                      zin(141) = c10*zin(129) + zc00*zin(138) + c01*zin(137)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  139

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  142

                                      xin(144) = c10*xin(138) + xc00*xin(141) + c01*xin(140)
                                      yin(144) = c10*yin(138) + yc00*yin(141) + c01*yin(140)
                                      zin(144) = c10*zin(138) + zc00*zin(141) + c01*zin(140)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  139
                                      ! i4 = i5 =  142

                                      ! nn =    6

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

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(142) = xin(142) + dxij*xin(139)
                                      yin(142) = yin(142) + dyij*yin(139)
                                      zin(142) = zin(142) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  136

                                      xin(139) = xin(139) + dxij*xin(136)
                                      yin(139) = yin(139) + dyij*yin(136)
                                      zin(139) = zin(139) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(142) = xin(142) + dxij*xin(139)
                                      yin(142) = yin(142) + dyij*yin(139)
                                      zin(142) = zin(142) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  112

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  112

                                      ! do ni = 1,    3

                                      xin(112) = xin(118) + dxij*xin(109)
                                      yin(112) = yin(118) + dyij*yin(109)
                                      zin(112) = zin(118) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! ni =    2

                                      xin(121) = xin(127) + dxij*xin(118)
                                      yin(121) = yin(127) + dyij*yin(118)
                                      zin(121) = zin(127) + dzij*zin(118)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                      ! ni =    3

                                      xin(130) = xin(136) + dxij*xin(127)
                                      yin(130) = yin(136) + dyij*yin(127)
                                      zin(130) = zin(136) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  115

                                      ! nj =    2

                                      ! i4 = i3 =  115

                                      ! do ni = 1,    3

                                      xin(115) = xin(121) + dxij*xin(112)
                                      yin(115) = yin(121) + dyij*yin(112)
                                      zin(115) = zin(121) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  124

                                      ! ni =    2

                                      xin(124) = xin(130) + dxij*xin(121)
                                      yin(124) = yin(130) + dyij*yin(121)
                                      zin(124) = zin(130) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    3

                                      xin(133) = xin(139) + dxij*xin(130)
                                      yin(133) = yin(139) + dyij*yin(130)
                                      zin(133) = zin(139) + dzij*zin(130)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  118

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(143) = xin(143) + dxij*xin(140)
                                      yin(143) = yin(143) + dyij*yin(140)
                                      zin(143) = zin(143) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  137

                                      xin(140) = xin(140) + dxij*xin(137)
                                      yin(140) = yin(140) + dyij*yin(137)
                                      zin(140) = zin(140) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(143) = xin(143) + dxij*xin(140)
                                      yin(143) = yin(143) + dyij*yin(140)
                                      zin(143) = zin(143) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  113

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  113

                                      ! do ni = 1,    3

                                      xin(113) = xin(119) + dxij*xin(110)
                                      yin(113) = yin(119) + dyij*yin(110)
                                      zin(113) = zin(119) + dzij*zin(110)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  122

                                      ! ni =    2

                                      xin(122) = xin(128) + dxij*xin(119)
                                      yin(122) = yin(128) + dyij*yin(119)
                                      zin(122) = zin(128) + dzij*zin(119)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      xin(131) = xin(137) + dxij*xin(128)
                                      yin(131) = yin(137) + dyij*yin(128)
                                      zin(131) = zin(137) + dzij*zin(128)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  140

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  116

                                      ! nj =    2

                                      ! i4 = i3 =  116

                                      ! do ni = 1,    3

                                      xin(116) = xin(122) + dxij*xin(113)
                                      yin(116) = yin(122) + dyij*yin(113)
                                      zin(116) = zin(122) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    2

                                      xin(125) = xin(131) + dxij*xin(122)
                                      yin(125) = yin(131) + dyij*yin(122)
                                      zin(125) = zin(131) + dzij*zin(122)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  134

                                      ! ni =    3

                                      xin(134) = xin(140) + dxij*xin(131)
                                      yin(134) = yin(140) + dyij*yin(131)
                                      zin(134) = zin(140) + dzij*zin(131)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  119

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  141

                                      xin(144) = xin(144) + dxij*xin(141)
                                      yin(144) = yin(144) + dyij*yin(141)
                                      zin(144) = zin(144) + dzij*zin(141)

                                      ! i3 = i4 =  141
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  138

                                      xin(141) = xin(141) + dxij*xin(138)
                                      yin(141) = yin(141) + dyij*yin(138)
                                      zin(141) = zin(141) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  141

                                      xin(144) = xin(144) + dxij*xin(141)
                                      yin(144) = yin(144) + dyij*yin(141)
                                      zin(144) = zin(144) + dzij*zin(141)

                                      ! i3 = i4 =  141
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  114

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  114

                                      ! do ni = 1,    3

                                      xin(114) = xin(120) + dxij*xin(111)
                                      yin(114) = yin(120) + dyij*yin(111)
                                      zin(114) = zin(120) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  123

                                      ! ni =    2

                                      xin(123) = xin(129) + dxij*xin(120)
                                      yin(123) = yin(129) + dyij*yin(120)
                                      zin(123) = zin(129) + dzij*zin(120)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    3

                                      xin(132) = xin(138) + dxij*xin(129)
                                      yin(132) = yin(138) + dyij*yin(129)
                                      zin(132) = zin(138) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  117

                                      ! nj =    2

                                      ! i4 = i3 =  117

                                      ! do ni = 1,    3

                                      xin(117) = xin(123) + dxij*xin(114)
                                      yin(117) = yin(123) + dyij*yin(114)
                                      zin(117) = zin(123) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    2

                                      xin(126) = xin(132) + dxij*xin(123)
                                      yin(126) = yin(132) + dyij*yin(123)
                                      zin(126) = zin(132) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    3

                                      xin(135) = xin(141) + dxij*xin(132)
                                      yin(135) = yin(141) + dyij*yin(132)
                                      zin(135) = zin(141) + dzij*zin(132)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  120

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 360! loop over all integrals

                                        l = n - 6*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d23bra(j)*d02ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 36)*yin(my + 36)*zin(mz + 36) & ! root  2
                                                        + xin(mx + 72)*yin(my + 72)*zin(mz + 72) & ! root  3
                                                        + xin(mx + 108)*yin(my + 108)*zin(mz + 108)) ! root  4

                                        j = int(n/6) + 1 ! index for the next bra cartesian pair

                                      end do

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
                                    ip = (i - 1)*36 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

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

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n02ket)
                              deallocate (xint02ket)

                              end subroutine int3220
                              end submodule
