! The total angular momentum of this class is:           5
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3200_impl
contains
  module subroutine int3200(df_pair, ss_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, ss_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n00ket(:)
    real(dp), allocatable :: xint23bra(:), xint00ket(:)
    integer(kind=int64) :: ndfbra, nssket
    real(dp) :: scutdfbra, scutssket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
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
    real(dp) :: roots(3), wghts(3)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(3), wts(3), alpha(3), beta(3), wrk(3)
    real(dp) :: xin(36), yin(36), zin(36)
    real(dp) :: eri_value(60)
    real(dp) :: d23bra(60), d00ket(1)
    integer(kind=int64) :: ix(10), jx(6), kx(1), lx(1)
    integer(kind=int64) :: iy(10), jy(6), ky(1), ly(1)
    integer(kind=int64) :: iz(10), jz(6), kz(1), lz(1)
    integer(kind=int64) :: in(6), in1(6), kn(1)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(1), kly(1), klz(1)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 4
    in1(3) = 7
    in1(4) = 10
    in1(5) = 11
    in1(6) = 12

    kn(1) = 0

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 0

    jx(1) = 2
    jx(2) = 0
    jx(3) = 0
    jx(4) = 1
    jx(5) = 1
    jx(6) = 0

    ix(1) = 10
    ix(2) = 1
    ix(3) = 1
    ix(4) = 7
    ix(5) = 7
    ix(6) = 4
    ix(7) = 1
    ix(8) = 4
    ix(9) = 1
    ix(10) = 4

    ! y-arrays

    ly(1) = 0

    ky(1) = 0

    jy(1) = 0
    jy(2) = 2
    jy(3) = 0
    jy(4) = 1
    jy(5) = 0
    jy(6) = 1

    iy(1) = 1
    iy(2) = 10
    iy(3) = 1
    iy(4) = 4
    iy(5) = 1
    iy(6) = 7
    iy(7) = 7
    iy(8) = 1
    iy(9) = 4
    iy(10) = 4

    ! z-arrays

    lz(1) = 0

    kz(1) = 0

    jz(1) = 0
    jz(2) = 0
    jz(3) = 2
    jz(4) = 0
    jz(5) = 1
    jz(6) = 1

    iz(1) = 1
    iz(2) = 1
    iz(3) = 10
    iz(4) = 1
    iz(5) = 4
    iz(6) = 1
    iz(7) = 4
    iz(8) = 7
    iz(9) = 7
    iz(10) = 4

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 12
    ijx(2) = 10
    ijx(3) = 10
    ijx(4) = 11
    ijx(5) = 11
    ijx(6) = 10
    ijx(7) = 3
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 2
    ijx(11) = 2
    ijx(12) = 1
    ijx(13) = 3
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 2
    ijx(17) = 2
    ijx(18) = 1
    ijx(19) = 9
    ijx(20) = 7
    ijx(21) = 7
    ijx(22) = 8
    ijx(23) = 8
    ijx(24) = 7
    ijx(25) = 9
    ijx(26) = 7
    ijx(27) = 7
    ijx(28) = 8
    ijx(29) = 8
    ijx(30) = 7
    ijx(31) = 6
    ijx(32) = 4
    ijx(33) = 4
    ijx(34) = 5
    ijx(35) = 5
    ijx(36) = 4
    ijx(37) = 3
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 2
    ijx(41) = 2
    ijx(42) = 1
    ijx(43) = 6
    ijx(44) = 4
    ijx(45) = 4
    ijx(46) = 5
    ijx(47) = 5
    ijx(48) = 4
    ijx(49) = 3
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 2
    ijx(53) = 2
    ijx(54) = 1
    ijx(55) = 6
    ijx(56) = 4
    ijx(57) = 4
    ijx(58) = 5
    ijx(59) = 5
    ijx(60) = 4

    ijy(1) = 1
    ijy(2) = 3
    ijy(3) = 1
    ijy(4) = 2
    ijy(5) = 1
    ijy(6) = 2
    ijy(7) = 10
    ijy(8) = 12
    ijy(9) = 10
    ijy(10) = 11
    ijy(11) = 10
    ijy(12) = 11
    ijy(13) = 1
    ijy(14) = 3
    ijy(15) = 1
    ijy(16) = 2
    ijy(17) = 1
    ijy(18) = 2
    ijy(19) = 4
    ijy(20) = 6
    ijy(21) = 4
    ijy(22) = 5
    ijy(23) = 4
    ijy(24) = 5
    ijy(25) = 1
    ijy(26) = 3
    ijy(27) = 1
    ijy(28) = 2
    ijy(29) = 1
    ijy(30) = 2
    ijy(31) = 7
    ijy(32) = 9
    ijy(33) = 7
    ijy(34) = 8
    ijy(35) = 7
    ijy(36) = 8
    ijy(37) = 7
    ijy(38) = 9
    ijy(39) = 7
    ijy(40) = 8
    ijy(41) = 7
    ijy(42) = 8
    ijy(43) = 1
    ijy(44) = 3
    ijy(45) = 1
    ijy(46) = 2
    ijy(47) = 1
    ijy(48) = 2
    ijy(49) = 4
    ijy(50) = 6
    ijy(51) = 4
    ijy(52) = 5
    ijy(53) = 4
    ijy(54) = 5
    ijy(55) = 4
    ijy(56) = 6
    ijy(57) = 4
    ijy(58) = 5
    ijy(59) = 4
    ijy(60) = 5

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 3
    ijz(4) = 1
    ijz(5) = 2
    ijz(6) = 2
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 3
    ijz(10) = 1
    ijz(11) = 2
    ijz(12) = 2
    ijz(13) = 10
    ijz(14) = 10
    ijz(15) = 12
    ijz(16) = 10
    ijz(17) = 11
    ijz(18) = 11
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 3
    ijz(22) = 1
    ijz(23) = 2
    ijz(24) = 2
    ijz(25) = 4
    ijz(26) = 4
    ijz(27) = 6
    ijz(28) = 4
    ijz(29) = 5
    ijz(30) = 5
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 3
    ijz(34) = 1
    ijz(35) = 2
    ijz(36) = 2
    ijz(37) = 4
    ijz(38) = 4
    ijz(39) = 6
    ijz(40) = 4
    ijz(41) = 5
    ijz(42) = 5
    ijz(43) = 7
    ijz(44) = 7
    ijz(45) = 9
    ijz(46) = 7
    ijz(47) = 8
    ijz(48) = 8
    ijz(49) = 7
    ijz(50) = 7
    ijz(51) = 9
    ijz(52) = 7
    ijz(53) = 8
    ijz(54) = 8
    ijz(55) = 4
    ijz(56) = 4
    ijz(57) = 6
    ijz(58) = 4
    ijz(59) = 5
    ijz(60) = 5

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 0

    kly(1) = 0

    klz(1) = 0

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n00ket(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00ket(res%n_s_shl*(res%n_s_shl + 1)/2))

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

    if ((ndfbra*nssket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nssket
    ntile = int(ndfbra*nssket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nssket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint00ket, n00ket, df_pair, ss_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d00ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,nm,nn,km,nj,ni) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint00ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n00ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
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

                    if (xx .ge. 43.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.1901635091934877D+00
                      rts(2) = factr*0.1784492748543251D+01
                      rts(3) = factr*0.5525343742263263D+01

                      wts(1) = factw*0.7246295952243917D+00
                      wts(2) = factw*0.1570673203228566D+00
                      wts(3) = factw*0.4530009905508823D-02

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

                      do kk = 1, 2

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
                        wrk(3) = 0.0D+00
                        do 100 kk = 2, 3

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 3

                            jj = 0

105                         do 110 m = l, 3
                              if (m .eq. 3) go to 120
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

                                do 300 ii = 2, 3

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 3
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

                                    do 310 kk = 1, 3
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 3
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

                                      ! i2 = in(2) =    4
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(4) = xc00
                                      yin(4) = yc00
                                      zin(4) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    4

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    7
                                      ! i3 =    1
                                      ! i4 =    4

                                      xin(7) = c10*xin(1) + xc00*xin(4)
                                      yin(7) = c10*yin(1) + yc00*yin(4)
                                      zin(7) = c10*zin(1) + zc00*zin(4)

                                      ! i3 = i4 =    4
                                      ! i4 = i5 =    7

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   10
                                      ! i3 =    4
                                      ! i4 =    7

                                      xin(10) = c10*xin(4) + xc00*xin(7)
                                      yin(10) = c10*yin(4) + yc00*yin(7)
                                      zin(10) = c10*zin(4) + zc00*zin(7)

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =   10

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   11
                                      ! i3 =    7
                                      ! i4 =   10

                                      xin(11) = c10*xin(7) + xc00*xin(10)
                                      yin(11) = c10*yin(7) + yc00*yin(10)
                                      zin(11) = c10*zin(7) + zc00*zin(10)

                                      ! i3 = i4 =   10
                                      ! i4 = i5 =   11

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   12
                                      ! i3 =   10
                                      ! i4 =   11

                                      xin(12) = c10*xin(10) + xc00*xin(11)
                                      yin(12) = c10*yin(10) + yc00*yin(11)
                                      zin(12) = c10*zin(10) + zc00*zin(11)

                                      ! i3 = i4 =   11
                                      ! i4 = i5 =   12

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   12

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   12

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   11

                                      xin(12) = xin(12) + dxij*xin(11)
                                      yin(12) = yin(12) + dyij*yin(11)
                                      zin(12) = zin(12) + dzij*zin(11)

                                      ! i3 = i4 =   11
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   10

                                      xin(11) = xin(11) + dxij*xin(10)
                                      yin(11) = yin(11) + dyij*yin(10)
                                      zin(11) = zin(11) + dzij*zin(10)

                                      ! i3 = i4 =   10
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   12

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   11

                                      xin(12) = xin(12) + dxij*xin(11)
                                      yin(12) = yin(12) + dyij*yin(11)
                                      zin(12) = zin(12) + dzij*zin(11)

                                      ! i3 = i4 =   11
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    2

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    2

                                      ! do ni = 1,    3

                                      xin(2) = xin(4) + dxij*xin(1)
                                      yin(2) = yin(4) + dyij*yin(1)
                                      zin(2) = zin(4) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    5

                                      ! ni =    2

                                      xin(5) = xin(7) + dxij*xin(4)
                                      yin(5) = yin(7) + dyij*yin(4)
                                      zin(5) = zin(7) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    8

                                      ! ni =    3

                                      xin(8) = xin(10) + dxij*xin(7)
                                      yin(8) = yin(10) + dyij*yin(7)
                                      zin(8) = zin(10) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   11

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    3

                                      ! nj =    2

                                      ! i4 = i3 =    3

                                      ! do ni = 1,    3

                                      xin(3) = xin(5) + dxij*xin(2)
                                      yin(3) = yin(5) + dyij*yin(2)
                                      zin(3) = zin(5) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    6

                                      ! ni =    2

                                      xin(6) = xin(8) + dxij*xin(5)
                                      yin(6) = yin(8) + dyij*yin(5)
                                      zin(6) = zin(8) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    9

                                      ! ni =    3

                                      xin(9) = xin(11) + dxij*xin(8)
                                      yin(9) = yin(11) + dyij*yin(8)
                                      zin(9) = zin(11) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   12

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    4

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   12

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

                                      ! i1 = in(1) =   13

                                      xin(13) = 1.0_dp
                                      yin(13) = 1.0_dp
                                      zin(13) = f00

                                      ! i2 = in(2) =   16
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(16) = xc00
                                      yin(16) = yc00
                                      zin(16) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   13
                                      ! i4 = i2 =   16

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   19
                                      ! i3 =   13
                                      ! i4 =   16

                                      xin(19) = c10*xin(13) + xc00*xin(16)
                                      yin(19) = c10*yin(13) + yc00*yin(16)
                                      zin(19) = c10*zin(13) + zc00*zin(16)

                                      ! i3 = i4 =   16
                                      ! i4 = i5 =   19

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   22
                                      ! i3 =   16
                                      ! i4 =   19

                                      xin(22) = c10*xin(16) + xc00*xin(19)
                                      yin(22) = c10*yin(16) + yc00*yin(19)
                                      zin(22) = c10*zin(16) + zc00*zin(19)

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   22

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   23
                                      ! i3 =   19
                                      ! i4 =   22

                                      xin(23) = c10*xin(19) + xc00*xin(22)
                                      yin(23) = c10*yin(19) + yc00*yin(22)
                                      zin(23) = c10*zin(19) + zc00*zin(22)

                                      ! i3 = i4 =   22
                                      ! i4 = i5 =   23

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   24
                                      ! i3 =   22
                                      ! i4 =   23

                                      xin(24) = c10*xin(22) + xc00*xin(23)
                                      yin(24) = c10*yin(22) + yc00*yin(23)
                                      zin(24) = c10*zin(22) + zc00*zin(23)

                                      ! i3 = i4 =   23
                                      ! i4 = i5 =   24

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   24

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   23

                                      xin(24) = xin(24) + dxij*xin(23)
                                      yin(24) = yin(24) + dyij*yin(23)
                                      zin(24) = zin(24) + dzij*zin(23)

                                      ! i3 = i4 =   23
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   22

                                      xin(23) = xin(23) + dxij*xin(22)
                                      yin(23) = yin(23) + dyij*yin(22)
                                      zin(23) = zin(23) + dzij*zin(22)

                                      ! i3 = i4 =   22
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   23

                                      xin(24) = xin(24) + dxij*xin(23)
                                      yin(24) = yin(24) + dyij*yin(23)
                                      zin(24) = zin(24) + dzij*zin(23)

                                      ! i3 = i4 =   23
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   14

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   14

                                      ! do ni = 1,    3

                                      xin(14) = xin(16) + dxij*xin(13)
                                      yin(14) = yin(16) + dyij*yin(13)
                                      zin(14) = zin(16) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    2

                                      xin(17) = xin(19) + dxij*xin(16)
                                      yin(17) = yin(19) + dyij*yin(16)
                                      zin(17) = zin(19) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   20

                                      ! ni =    3

                                      xin(20) = xin(22) + dxij*xin(19)
                                      yin(20) = yin(22) + dyij*yin(19)
                                      zin(20) = zin(22) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    2

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    3

                                      xin(15) = xin(17) + dxij*xin(14)
                                      yin(15) = yin(17) + dyij*yin(14)
                                      zin(15) = zin(17) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    2

                                      xin(18) = xin(20) + dxij*xin(17)
                                      yin(18) = yin(20) + dyij*yin(17)
                                      zin(18) = zin(20) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    3

                                      xin(21) = xin(23) + dxij*xin(20)
                                      yin(21) = yin(23) + dyij*yin(20)
                                      zin(21) = zin(23) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   16

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   24

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

                                      ! i1 = in(1) =   25

                                      xin(25) = 1.0_dp
                                      yin(25) = 1.0_dp
                                      zin(25) = f00

                                      ! i2 = in(2) =   28
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(28) = xc00
                                      yin(28) = yc00
                                      zin(28) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   28

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   31
                                      ! i3 =   25
                                      ! i4 =   28

                                      xin(31) = c10*xin(25) + xc00*xin(28)
                                      yin(31) = c10*yin(25) + yc00*yin(28)
                                      zin(31) = c10*zin(25) + zc00*zin(28)

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   31

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   34
                                      ! i3 =   28
                                      ! i4 =   31

                                      xin(34) = c10*xin(28) + xc00*xin(31)
                                      yin(34) = c10*yin(28) + yc00*yin(31)
                                      zin(34) = c10*zin(28) + zc00*zin(31)

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   34

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   35
                                      ! i3 =   31
                                      ! i4 =   34

                                      xin(35) = c10*xin(31) + xc00*xin(34)
                                      yin(35) = c10*yin(31) + yc00*yin(34)
                                      zin(35) = c10*zin(31) + zc00*zin(34)

                                      ! i3 = i4 =   34
                                      ! i4 = i5 =   35

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   36
                                      ! i3 =   34
                                      ! i4 =   35

                                      xin(36) = c10*xin(34) + xc00*xin(35)
                                      yin(36) = c10*yin(34) + yc00*yin(35)
                                      zin(36) = c10*zin(34) + zc00*zin(35)

                                      ! i3 = i4 =   35
                                      ! i4 = i5 =   36

                                      ! n =    6

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   36

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   35

                                      xin(36) = xin(36) + dxij*xin(35)
                                      yin(36) = yin(36) + dyij*yin(35)
                                      zin(36) = zin(36) + dzij*zin(35)

                                      ! i3 = i4 =   35
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   34

                                      xin(35) = xin(35) + dxij*xin(34)
                                      yin(35) = yin(35) + dyij*yin(34)
                                      zin(35) = zin(35) + dzij*zin(34)

                                      ! i3 = i4 =   34
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   35

                                      xin(36) = xin(36) + dxij*xin(35)
                                      yin(36) = yin(36) + dyij*yin(35)
                                      zin(36) = zin(36) + dzij*zin(35)

                                      ! i3 = i4 =   35
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   26

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   26

                                      ! do ni = 1,    3

                                      xin(26) = xin(28) + dxij*xin(25)
                                      yin(26) = yin(28) + dyij*yin(25)
                                      zin(26) = zin(28) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    2

                                      xin(29) = xin(31) + dxij*xin(28)
                                      yin(29) = yin(31) + dyij*yin(28)
                                      zin(29) = zin(31) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    3

                                      xin(32) = xin(34) + dxij*xin(31)
                                      yin(32) = yin(34) + dyij*yin(31)
                                      zin(32) = zin(34) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    2

                                      ! i4 = i3 =   27

                                      ! do ni = 1,    3

                                      xin(27) = xin(29) + dxij*xin(26)
                                      yin(27) = yin(29) + dyij*yin(26)
                                      zin(27) = zin(29) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    2

                                      xin(30) = xin(32) + dxij*xin(29)
                                      yin(30) = yin(32) + dyij*yin(29)
                                      zin(30) = zin(32) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    3

                                      xin(33) = xin(35) + dxij*xin(32)
                                      yin(33) = yin(35) + dyij*yin(32)
                                      zin(33) = zin(35) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   28

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   36

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

       eri_value(1) = eri_value(1) + d23bra(1)*d00ket(1)*(xin(12)*yin(1)*zin(1) + xin(24)*yin(13)*zin(13) + xin(36)*yin(25)*zin(25))
       eri_value(2) = eri_value(2) + d23bra(2)*d00ket(1)*(xin(10)*yin(3)*zin(1) + xin(22)*yin(15)*zin(13) + xin(34)*yin(27)*zin(25))
       eri_value(3) = eri_value(3) + d23bra(3)*d00ket(1)*(xin(10)*yin(1)*zin(3) + xin(22)*yin(13)*zin(15) + xin(34)*yin(25)*zin(27))
       eri_value(4) = eri_value(4) + d23bra(4)*d00ket(1)*(xin(11)*yin(2)*zin(1) + xin(23)*yin(14)*zin(13) + xin(35)*yin(26)*zin(25))
       eri_value(5) = eri_value(5) + d23bra(5)*d00ket(1)*(xin(11)*yin(1)*zin(2) + xin(23)*yin(13)*zin(14) + xin(35)*yin(25)*zin(26))
       eri_value(6) = eri_value(6) + d23bra(6)*d00ket(1)*(xin(10)*yin(2)*zin(2) + xin(22)*yin(14)*zin(14) + xin(34)*yin(26)*zin(26))
       eri_value(7) = eri_value(7) + d23bra(7)*d00ket(1)*(xin(3)*yin(10)*zin(1) + xin(15)*yin(22)*zin(13) + xin(27)*yin(34)*zin(25))
       eri_value(8) = eri_value(8) + d23bra(8)*d00ket(1)*(xin(1)*yin(12)*zin(1) + xin(13)*yin(24)*zin(13) + xin(25)*yin(36)*zin(25))
       eri_value(9) = eri_value(9) + d23bra(9)*d00ket(1)*(xin(1)*yin(10)*zin(3) + xin(13)*yin(22)*zin(15) + xin(25)*yin(34)*zin(27))
    eri_value(10) = eri_value(10) + d23bra(10)*d00ket(1)*(xin(2)*yin(11)*zin(1) + xin(14)*yin(23)*zin(13) + xin(26)*yin(35)*zin(25))
    eri_value(11) = eri_value(11) + d23bra(11)*d00ket(1)*(xin(2)*yin(10)*zin(2) + xin(14)*yin(22)*zin(14) + xin(26)*yin(34)*zin(26))
    eri_value(12) = eri_value(12) + d23bra(12)*d00ket(1)*(xin(1)*yin(11)*zin(2) + xin(13)*yin(23)*zin(14) + xin(25)*yin(35)*zin(26))
    eri_value(13) = eri_value(13) + d23bra(13)*d00ket(1)*(xin(3)*yin(1)*zin(10) + xin(15)*yin(13)*zin(22) + xin(27)*yin(25)*zin(34))
    eri_value(14) = eri_value(14) + d23bra(14)*d00ket(1)*(xin(1)*yin(3)*zin(10) + xin(13)*yin(15)*zin(22) + xin(25)*yin(27)*zin(34))
    eri_value(15) = eri_value(15) + d23bra(15)*d00ket(1)*(xin(1)*yin(1)*zin(12) + xin(13)*yin(13)*zin(24) + xin(25)*yin(25)*zin(36))
    eri_value(16) = eri_value(16) + d23bra(16)*d00ket(1)*(xin(2)*yin(2)*zin(10) + xin(14)*yin(14)*zin(22) + xin(26)*yin(26)*zin(34))
    eri_value(17) = eri_value(17) + d23bra(17)*d00ket(1)*(xin(2)*yin(1)*zin(11) + xin(14)*yin(13)*zin(23) + xin(26)*yin(25)*zin(35))
    eri_value(18) = eri_value(18) + d23bra(18)*d00ket(1)*(xin(1)*yin(2)*zin(11) + xin(13)*yin(14)*zin(23) + xin(25)*yin(26)*zin(35))
     eri_value(19) = eri_value(19) + d23bra(19)*d00ket(1)*(xin(9)*yin(4)*zin(1) + xin(21)*yin(16)*zin(13) + xin(33)*yin(28)*zin(25))
     eri_value(20) = eri_value(20) + d23bra(20)*d00ket(1)*(xin(7)*yin(6)*zin(1) + xin(19)*yin(18)*zin(13) + xin(31)*yin(30)*zin(25))
     eri_value(21) = eri_value(21) + d23bra(21)*d00ket(1)*(xin(7)*yin(4)*zin(3) + xin(19)*yin(16)*zin(15) + xin(31)*yin(28)*zin(27))
     eri_value(22) = eri_value(22) + d23bra(22)*d00ket(1)*(xin(8)*yin(5)*zin(1) + xin(20)*yin(17)*zin(13) + xin(32)*yin(29)*zin(25))
     eri_value(23) = eri_value(23) + d23bra(23)*d00ket(1)*(xin(8)*yin(4)*zin(2) + xin(20)*yin(16)*zin(14) + xin(32)*yin(28)*zin(26))
     eri_value(24) = eri_value(24) + d23bra(24)*d00ket(1)*(xin(7)*yin(5)*zin(2) + xin(19)*yin(17)*zin(14) + xin(31)*yin(29)*zin(26))
     eri_value(25) = eri_value(25) + d23bra(25)*d00ket(1)*(xin(9)*yin(1)*zin(4) + xin(21)*yin(13)*zin(16) + xin(33)*yin(25)*zin(28))
     eri_value(26) = eri_value(26) + d23bra(26)*d00ket(1)*(xin(7)*yin(3)*zin(4) + xin(19)*yin(15)*zin(16) + xin(31)*yin(27)*zin(28))
     eri_value(27) = eri_value(27) + d23bra(27)*d00ket(1)*(xin(7)*yin(1)*zin(6) + xin(19)*yin(13)*zin(18) + xin(31)*yin(25)*zin(30))
     eri_value(28) = eri_value(28) + d23bra(28)*d00ket(1)*(xin(8)*yin(2)*zin(4) + xin(20)*yin(14)*zin(16) + xin(32)*yin(26)*zin(28))
     eri_value(29) = eri_value(29) + d23bra(29)*d00ket(1)*(xin(8)*yin(1)*zin(5) + xin(20)*yin(13)*zin(17) + xin(32)*yin(25)*zin(29))
     eri_value(30) = eri_value(30) + d23bra(30)*d00ket(1)*(xin(7)*yin(2)*zin(5) + xin(19)*yin(14)*zin(17) + xin(31)*yin(26)*zin(29))
     eri_value(31) = eri_value(31) + d23bra(31)*d00ket(1)*(xin(6)*yin(7)*zin(1) + xin(18)*yin(19)*zin(13) + xin(30)*yin(31)*zin(25))
     eri_value(32) = eri_value(32) + d23bra(32)*d00ket(1)*(xin(4)*yin(9)*zin(1) + xin(16)*yin(21)*zin(13) + xin(28)*yin(33)*zin(25))
     eri_value(33) = eri_value(33) + d23bra(33)*d00ket(1)*(xin(4)*yin(7)*zin(3) + xin(16)*yin(19)*zin(15) + xin(28)*yin(31)*zin(27))
     eri_value(34) = eri_value(34) + d23bra(34)*d00ket(1)*(xin(5)*yin(8)*zin(1) + xin(17)*yin(20)*zin(13) + xin(29)*yin(32)*zin(25))
     eri_value(35) = eri_value(35) + d23bra(35)*d00ket(1)*(xin(5)*yin(7)*zin(2) + xin(17)*yin(19)*zin(14) + xin(29)*yin(31)*zin(26))
     eri_value(36) = eri_value(36) + d23bra(36)*d00ket(1)*(xin(4)*yin(8)*zin(2) + xin(16)*yin(20)*zin(14) + xin(28)*yin(32)*zin(26))
     eri_value(37) = eri_value(37) + d23bra(37)*d00ket(1)*(xin(3)*yin(7)*zin(4) + xin(15)*yin(19)*zin(16) + xin(27)*yin(31)*zin(28))
     eri_value(38) = eri_value(38) + d23bra(38)*d00ket(1)*(xin(1)*yin(9)*zin(4) + xin(13)*yin(21)*zin(16) + xin(25)*yin(33)*zin(28))
     eri_value(39) = eri_value(39) + d23bra(39)*d00ket(1)*(xin(1)*yin(7)*zin(6) + xin(13)*yin(19)*zin(18) + xin(25)*yin(31)*zin(30))
     eri_value(40) = eri_value(40) + d23bra(40)*d00ket(1)*(xin(2)*yin(8)*zin(4) + xin(14)*yin(20)*zin(16) + xin(26)*yin(32)*zin(28))
     eri_value(41) = eri_value(41) + d23bra(41)*d00ket(1)*(xin(2)*yin(7)*zin(5) + xin(14)*yin(19)*zin(17) + xin(26)*yin(31)*zin(29))
     eri_value(42) = eri_value(42) + d23bra(42)*d00ket(1)*(xin(1)*yin(8)*zin(5) + xin(13)*yin(20)*zin(17) + xin(25)*yin(32)*zin(29))
     eri_value(43) = eri_value(43) + d23bra(43)*d00ket(1)*(xin(6)*yin(1)*zin(7) + xin(18)*yin(13)*zin(19) + xin(30)*yin(25)*zin(31))
     eri_value(44) = eri_value(44) + d23bra(44)*d00ket(1)*(xin(4)*yin(3)*zin(7) + xin(16)*yin(15)*zin(19) + xin(28)*yin(27)*zin(31))
     eri_value(45) = eri_value(45) + d23bra(45)*d00ket(1)*(xin(4)*yin(1)*zin(9) + xin(16)*yin(13)*zin(21) + xin(28)*yin(25)*zin(33))
     eri_value(46) = eri_value(46) + d23bra(46)*d00ket(1)*(xin(5)*yin(2)*zin(7) + xin(17)*yin(14)*zin(19) + xin(29)*yin(26)*zin(31))
     eri_value(47) = eri_value(47) + d23bra(47)*d00ket(1)*(xin(5)*yin(1)*zin(8) + xin(17)*yin(13)*zin(20) + xin(29)*yin(25)*zin(32))
     eri_value(48) = eri_value(48) + d23bra(48)*d00ket(1)*(xin(4)*yin(2)*zin(8) + xin(16)*yin(14)*zin(20) + xin(28)*yin(26)*zin(32))
     eri_value(49) = eri_value(49) + d23bra(49)*d00ket(1)*(xin(3)*yin(4)*zin(7) + xin(15)*yin(16)*zin(19) + xin(27)*yin(28)*zin(31))
     eri_value(50) = eri_value(50) + d23bra(50)*d00ket(1)*(xin(1)*yin(6)*zin(7) + xin(13)*yin(18)*zin(19) + xin(25)*yin(30)*zin(31))
     eri_value(51) = eri_value(51) + d23bra(51)*d00ket(1)*(xin(1)*yin(4)*zin(9) + xin(13)*yin(16)*zin(21) + xin(25)*yin(28)*zin(33))
     eri_value(52) = eri_value(52) + d23bra(52)*d00ket(1)*(xin(2)*yin(5)*zin(7) + xin(14)*yin(17)*zin(19) + xin(26)*yin(29)*zin(31))
     eri_value(53) = eri_value(53) + d23bra(53)*d00ket(1)*(xin(2)*yin(4)*zin(8) + xin(14)*yin(16)*zin(20) + xin(26)*yin(28)*zin(32))
     eri_value(54) = eri_value(54) + d23bra(54)*d00ket(1)*(xin(1)*yin(5)*zin(8) + xin(13)*yin(17)*zin(20) + xin(25)*yin(29)*zin(32))
     eri_value(55) = eri_value(55) + d23bra(55)*d00ket(1)*(xin(6)*yin(4)*zin(4) + xin(18)*yin(16)*zin(16) + xin(30)*yin(28)*zin(28))
     eri_value(56) = eri_value(56) + d23bra(56)*d00ket(1)*(xin(4)*yin(6)*zin(4) + xin(16)*yin(18)*zin(16) + xin(28)*yin(30)*zin(28))
     eri_value(57) = eri_value(57) + d23bra(57)*d00ket(1)*(xin(4)*yin(4)*zin(6) + xin(16)*yin(16)*zin(18) + xin(28)*yin(28)*zin(30))
     eri_value(58) = eri_value(58) + d23bra(58)*d00ket(1)*(xin(5)*yin(5)*zin(4) + xin(17)*yin(17)*zin(16) + xin(29)*yin(29)*zin(28))
     eri_value(59) = eri_value(59) + d23bra(59)*d00ket(1)*(xin(5)*yin(4)*zin(5) + xin(17)*yin(16)*zin(17) + xin(29)*yin(28)*zin(29))
     eri_value(60) = eri_value(60) + d23bra(60)*d00ket(1)*(xin(4)*yin(5)*zin(5) + xin(16)*yin(17)*zin(17) + xin(28)*yin(29)*zin(29))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxl = 1
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*6 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

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

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n00ket)
                              deallocate (xint00ket)

                              end subroutine int3200
                              end submodule
