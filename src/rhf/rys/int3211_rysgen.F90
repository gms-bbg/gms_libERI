! The total angular momentum of this class is:           7
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3211_impl
contains
  module subroutine int3211(df_pair, pp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, pp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n11ket(:)
    real(dp), allocatable :: xint23bra(:), xint11ket(:)
    integer(kind=int64) :: ndfbra, nppket
    real(dp) :: scutdfbra, scutppket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(192), yin(192), zin(192)
    real(dp) :: eri_value(540)
    real(dp) :: d23bra(60), d11ket(9)
    integer(kind=int64) :: ix(10), jx(6), kx(3), lx(3)
    integer(kind=int64) :: iy(10), jy(6), ky(3), ly(3)
    integer(kind=int64) :: iz(10), jz(6), kz(3), lz(3)
    integer(kind=int64) :: in(6), in1(6), kn(3)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(9), kly(9), klz(9)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 37
    in1(5) = 41
    in1(6) = 45

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

    jx(1) = 8
    jx(2) = 0
    jx(3) = 0
    jx(4) = 4
    jx(5) = 4
    jx(6) = 0

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
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0

    jy(1) = 0
    jy(2) = 8
    jy(3) = 0
    jy(4) = 4
    jy(5) = 0
    jy(6) = 4

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
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 8
    jz(4) = 0
    jz(5) = 4
    jz(6) = 4

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

    ijx(1) = 45
    ijx(2) = 37
    ijx(3) = 37
    ijx(4) = 41
    ijx(5) = 41
    ijx(6) = 37
    ijx(7) = 9
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 5
    ijx(11) = 5
    ijx(12) = 1
    ijx(13) = 9
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 5
    ijx(17) = 5
    ijx(18) = 1
    ijx(19) = 33
    ijx(20) = 25
    ijx(21) = 25
    ijx(22) = 29
    ijx(23) = 29
    ijx(24) = 25
    ijx(25) = 33
    ijx(26) = 25
    ijx(27) = 25
    ijx(28) = 29
    ijx(29) = 29
    ijx(30) = 25
    ijx(31) = 21
    ijx(32) = 13
    ijx(33) = 13
    ijx(34) = 17
    ijx(35) = 17
    ijx(36) = 13
    ijx(37) = 9
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 5
    ijx(41) = 5
    ijx(42) = 1
    ijx(43) = 21
    ijx(44) = 13
    ijx(45) = 13
    ijx(46) = 17
    ijx(47) = 17
    ijx(48) = 13
    ijx(49) = 9
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 5
    ijx(53) = 5
    ijx(54) = 1
    ijx(55) = 21
    ijx(56) = 13
    ijx(57) = 13
    ijx(58) = 17
    ijx(59) = 17
    ijx(60) = 13

    ijy(1) = 1
    ijy(2) = 9
    ijy(3) = 1
    ijy(4) = 5
    ijy(5) = 1
    ijy(6) = 5
    ijy(7) = 37
    ijy(8) = 45
    ijy(9) = 37
    ijy(10) = 41
    ijy(11) = 37
    ijy(12) = 41
    ijy(13) = 1
    ijy(14) = 9
    ijy(15) = 1
    ijy(16) = 5
    ijy(17) = 1
    ijy(18) = 5
    ijy(19) = 13
    ijy(20) = 21
    ijy(21) = 13
    ijy(22) = 17
    ijy(23) = 13
    ijy(24) = 17
    ijy(25) = 1
    ijy(26) = 9
    ijy(27) = 1
    ijy(28) = 5
    ijy(29) = 1
    ijy(30) = 5
    ijy(31) = 25
    ijy(32) = 33
    ijy(33) = 25
    ijy(34) = 29
    ijy(35) = 25
    ijy(36) = 29
    ijy(37) = 25
    ijy(38) = 33
    ijy(39) = 25
    ijy(40) = 29
    ijy(41) = 25
    ijy(42) = 29
    ijy(43) = 1
    ijy(44) = 9
    ijy(45) = 1
    ijy(46) = 5
    ijy(47) = 1
    ijy(48) = 5
    ijy(49) = 13
    ijy(50) = 21
    ijy(51) = 13
    ijy(52) = 17
    ijy(53) = 13
    ijy(54) = 17
    ijy(55) = 13
    ijy(56) = 21
    ijy(57) = 13
    ijy(58) = 17
    ijy(59) = 13
    ijy(60) = 17

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 9
    ijz(4) = 1
    ijz(5) = 5
    ijz(6) = 5
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 9
    ijz(10) = 1
    ijz(11) = 5
    ijz(12) = 5
    ijz(13) = 37
    ijz(14) = 37
    ijz(15) = 45
    ijz(16) = 37
    ijz(17) = 41
    ijz(18) = 41
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 9
    ijz(22) = 1
    ijz(23) = 5
    ijz(24) = 5
    ijz(25) = 13
    ijz(26) = 13
    ijz(27) = 21
    ijz(28) = 13
    ijz(29) = 17
    ijz(30) = 17
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 9
    ijz(34) = 1
    ijz(35) = 5
    ijz(36) = 5
    ijz(37) = 13
    ijz(38) = 13
    ijz(39) = 21
    ijz(40) = 13
    ijz(41) = 17
    ijz(42) = 17
    ijz(43) = 25
    ijz(44) = 25
    ijz(45) = 33
    ijz(46) = 25
    ijz(47) = 29
    ijz(48) = 29
    ijz(49) = 25
    ijz(50) = 25
    ijz(51) = 33
    ijz(52) = 25
    ijz(53) = 29
    ijz(54) = 29
    ijz(55) = 13
    ijz(56) = 13
    ijz(57) = 21
    ijz(58) = 13
    ijz(59) = 17
    ijz(60) = 17

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

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11ket(res%n_p_shl*(res%n_p_shl + 1)/2))

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

    if ((ndfbra*nppket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nppket
    ntile = int(ndfbra*nppket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nppket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint11ket, n11ket, df_pair, pp_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d11ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint11ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n11ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
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

                                      ! i2 = in(2) =   13
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(13) = xc00
                                      yin(13) = yc00
                                      zin(13) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   15
                                      ! i2 =   13

                                      xin(15) = xcp00*xin(13) + cp10
                                      yin(15) = ycp00*yin(13) + cp10
                                      zin(15) = zcp00*zin(13) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   25
                                      ! i3 =    1
                                      ! i4 =   13

                                      xin(25) = c10*xin(1) + xc00*xin(13)
                                      yin(25) = c10*yin(1) + yc00*yin(13)
                                      zin(25) = c10*zin(1) + zc00*zin(13)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   27
                                      ! i5 =   25
                                      ! i4 =   13

                                      xin(27) = xcp00*xin(25) + cp10*xin(13)
                                      yin(27) = ycp00*yin(25) + cp10*yin(13)
                                      zin(27) = zcp00*zin(25) + cp10*zin(13)

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

                                      ! i3 = i5 + k2 =   39
                                      ! i5 =   37
                                      ! i4 =   25

                                      xin(39) = xcp00*xin(37) + cp10*xin(25)
                                      yin(39) = ycp00*yin(37) + cp10*yin(25)
                                      zin(39) = zcp00*zin(37) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   41
                                      ! i3 =   25
                                      ! i4 =   37

                                      xin(41) = c10*xin(25) + xc00*xin(37)
                                      yin(41) = c10*yin(25) + yc00*yin(37)
                                      zin(41) = c10*zin(25) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   43
                                      ! i5 =   41
                                      ! i4 =   37

                                      xin(43) = xcp00*xin(41) + cp10*xin(37)
                                      yin(43) = ycp00*yin(41) + cp10*yin(37)
                                      zin(43) = zcp00*zin(41) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   45
                                      ! i3 =   37
                                      ! i4 =   41

                                      xin(45) = c10*xin(37) + xc00*xin(41)
                                      yin(45) = c10*yin(37) + yc00*yin(41)
                                      zin(45) = c10*zin(37) + zc00*zin(41)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   47
                                      ! i5 =   45
                                      ! i4 =   41

                                      xin(47) = xcp00*xin(45) + cp10*xin(41)
                                      yin(47) = ycp00*yin(45) + cp10*yin(41)
                                      zin(47) = zcp00*zin(45) + cp10*zin(41)

                                      ! ------------------

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! n =    6

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

                                      ! i3 = i2+kn(n+1) =   16

                                      xin(16) = xc00*xin(4) + c01*xin(3)
                                      yin(16) = yc00*yin(4) + c01*yin(3)
                                      zin(16) = zc00*zin(4) + c01*zin(3)

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
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   25

                                      xin(28) = c10*xin(4) + xc00*xin(16) + c01*xin(15)
                                      yin(28) = c10*yin(4) + yc00*yin(16) + c01*yin(15)
                                      zin(28) = c10*zin(4) + zc00*zin(16) + c01*zin(15)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(40) = c10*xin(16) + xc00*xin(28) + c01*xin(27)
                                      yin(40) = c10*yin(16) + yc00*yin(28) + c01*yin(27)
                                      zin(40) = c10*zin(16) + zc00*zin(28) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   41

                                      xin(44) = c10*xin(28) + xc00*xin(40) + c01*xin(39)
                                      yin(44) = c10*yin(28) + yc00*yin(40) + c01*yin(39)
                                      zin(44) = c10*zin(28) + zc00*zin(40) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   45

                                      xin(48) = c10*xin(40) + xc00*xin(44) + c01*xin(43)
                                      yin(48) = c10*yin(40) + yc00*yin(44) + c01*yin(43)
                                      zin(48) = c10*zin(40) + zc00*zin(44) + c01*zin(43)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   45

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(45) = xin(45) + dxij*xin(41)
                                      yin(45) = yin(45) + dyij*yin(41)
                                      zin(45) = zin(45) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   37

                                      xin(41) = xin(41) + dxij*xin(37)
                                      yin(41) = yin(41) + dyij*yin(37)
                                      zin(41) = zin(41) + dzij*zin(37)

                                      ! i3 = i4 =   37
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(45) = xin(45) + dxij*xin(41)
                                      yin(45) = yin(45) + dyij*yin(41)
                                      zin(45) = zin(45) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(13) + dxij*xin(1)
                                      yin(5) = yin(13) + dyij*yin(1)
                                      zin(5) = zin(13) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    2

                                      xin(17) = xin(25) + dxij*xin(13)
                                      yin(17) = yin(25) + dyij*yin(13)
                                      zin(17) = zin(25) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    3

                                      xin(29) = xin(37) + dxij*xin(25)
                                      yin(29) = yin(37) + dyij*yin(25)
                                      zin(29) = zin(37) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(17) + dxij*xin(5)
                                      yin(9) = yin(17) + dyij*yin(5)
                                      zin(9) = zin(17) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    2

                                      xin(21) = xin(29) + dxij*xin(17)
                                      yin(21) = yin(29) + dyij*yin(17)
                                      zin(21) = zin(29) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    3

                                      xin(33) = xin(41) + dxij*xin(29)
                                      yin(33) = yin(41) + dyij*yin(29)
                                      zin(33) = zin(41) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(47) = xin(47) + dxij*xin(43)
                                      yin(47) = yin(47) + dyij*yin(43)
                                      zin(47) = zin(47) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   39

                                      xin(43) = xin(43) + dxij*xin(39)
                                      yin(43) = yin(43) + dyij*yin(39)
                                      zin(43) = zin(43) + dzij*zin(39)

                                      ! i3 = i4 =   39
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(47) = xin(47) + dxij*xin(43)
                                      yin(47) = yin(47) + dyij*yin(43)
                                      zin(47) = zin(47) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(15) + dxij*xin(3)
                                      yin(7) = yin(15) + dyij*yin(3)
                                      zin(7) = zin(15) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   19

                                      ! ni =    2

                                      xin(19) = xin(27) + dxij*xin(15)
                                      yin(19) = yin(27) + dyij*yin(15)
                                      zin(19) = zin(27) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    3

                                      xin(31) = xin(39) + dxij*xin(27)
                                      yin(31) = yin(39) + dyij*yin(27)
                                      zin(31) = zin(39) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(19) + dxij*xin(7)
                                      yin(11) = yin(19) + dyij*yin(7)
                                      zin(11) = zin(19) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    2

                                      xin(23) = xin(31) + dxij*xin(19)
                                      yin(23) = yin(31) + dyij*yin(19)
                                      zin(23) = zin(31) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    3

                                      xin(35) = xin(43) + dxij*xin(31)
                                      yin(35) = yin(43) + dyij*yin(31)
                                      zin(35) = zin(43) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(48) = xin(48) + dxij*xin(44)
                                      yin(48) = yin(48) + dyij*yin(44)
                                      zin(48) = zin(48) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   40

                                      xin(44) = xin(44) + dxij*xin(40)
                                      yin(44) = yin(44) + dyij*yin(40)
                                      zin(44) = zin(44) + dzij*zin(40)

                                      ! i3 = i4 =   40
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(48) = xin(48) + dxij*xin(44)
                                      yin(48) = yin(48) + dyij*yin(44)
                                      zin(48) = zin(48) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

                                      xin(8) = xin(16) + dxij*xin(4)
                                      yin(8) = yin(16) + dyij*yin(4)
                                      zin(8) = zin(16) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   20

                                      ! ni =    2

                                      xin(20) = xin(28) + dxij*xin(16)
                                      yin(20) = yin(28) + dyij*yin(16)
                                      zin(20) = zin(28) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    3

                                      xin(32) = xin(40) + dxij*xin(28)
                                      yin(32) = yin(40) + dyij*yin(28)
                                      zin(32) = zin(40) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(20) + dxij*xin(8)
                                      yin(12) = yin(20) + dyij*yin(8)
                                      zin(12) = zin(20) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    2

                                      xin(24) = xin(32) + dxij*xin(20)
                                      yin(24) = yin(32) + dyij*yin(20)
                                      zin(24) = zin(32) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    3

                                      xin(36) = xin(44) + dxij*xin(32)
                                      yin(36) = yin(44) + dyij*yin(32)
                                      zin(36) = zin(44) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   16

                                      ! nj =    3

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   13

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

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
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(61) = xc00
                                      yin(61) = yc00
                                      zin(61) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   51

                                      xin(51) = xcp00
                                      yin(51) = ycp00
                                      zin(51) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   63
                                      ! i2 =   61

                                      xin(63) = xcp00*xin(61) + cp10
                                      yin(63) = ycp00*yin(61) + cp10
                                      zin(63) = zcp00*zin(61) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   73
                                      ! i3 =   49
                                      ! i4 =   61

                                      xin(73) = c10*xin(49) + xc00*xin(61)
                                      yin(73) = c10*yin(49) + yc00*yin(61)
                                      zin(73) = c10*zin(49) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   75
                                      ! i5 =   73
                                      ! i4 =   61

                                      xin(75) = xcp00*xin(73) + cp10*xin(61)
                                      yin(75) = ycp00*yin(73) + cp10*yin(61)
                                      zin(75) = zcp00*zin(73) + cp10*zin(61)

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

                                      ! i3 = i5 + k2 =   87
                                      ! i5 =   85
                                      ! i4 =   73

                                      xin(87) = xcp00*xin(85) + cp10*xin(73)
                                      yin(87) = ycp00*yin(85) + cp10*yin(73)
                                      zin(87) = zcp00*zin(85) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   89
                                      ! i3 =   73
                                      ! i4 =   85

                                      xin(89) = c10*xin(73) + xc00*xin(85)
                                      yin(89) = c10*yin(73) + yc00*yin(85)
                                      zin(89) = c10*zin(73) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   91
                                      ! i5 =   89
                                      ! i4 =   85

                                      xin(91) = xcp00*xin(89) + cp10*xin(85)
                                      yin(91) = ycp00*yin(89) + cp10*yin(85)
                                      zin(91) = zcp00*zin(89) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   89

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   93
                                      ! i3 =   85
                                      ! i4 =   89

                                      xin(93) = c10*xin(85) + xc00*xin(89)
                                      yin(93) = c10*yin(85) + yc00*yin(89)
                                      zin(93) = c10*zin(85) + zc00*zin(89)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   95
                                      ! i5 =   93
                                      ! i4 =   89

                                      xin(95) = xcp00*xin(93) + cp10*xin(89)
                                      yin(95) = ycp00*yin(93) + cp10*yin(89)
                                      zin(95) = zcp00*zin(93) + cp10*zin(89)

                                      ! ------------------

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   49
                                      ! i4 = i1+k2 =   51

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   52
                                      ! i3 =   49
                                      ! i4 =   51

                                      xin(52) = cp01*xin(49) + xcp00*xin(51)
                                      yin(52) = cp01*yin(49) + ycp00*yin(51)
                                      zin(52) = cp01*zin(49) + zcp00*zin(51)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   64

                                      xin(64) = xc00*xin(52) + c01*xin(51)
                                      yin(64) = yc00*yin(52) + c01*yin(51)
                                      zin(64) = zc00*zin(52) + c01*zin(51)

                                      ! ------------------

                                      ! i3 = i4 =   51
                                      ! i4 = i5 =   52

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   73

                                      xin(76) = c10*xin(52) + xc00*xin(64) + c01*xin(63)
                                      yin(76) = c10*yin(52) + yc00*yin(64) + c01*yin(63)
                                      zin(76) = c10*zin(52) + zc00*zin(64) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(88) = c10*xin(64) + xc00*xin(76) + c01*xin(75)
                                      yin(88) = c10*yin(64) + yc00*yin(76) + c01*yin(75)
                                      zin(88) = c10*zin(64) + zc00*zin(76) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   89

                                      xin(92) = c10*xin(76) + xc00*xin(88) + c01*xin(87)
                                      yin(92) = c10*yin(76) + yc00*yin(88) + c01*yin(87)
                                      zin(92) = c10*zin(76) + zc00*zin(88) + c01*zin(87)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   89

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   93

                                      xin(96) = c10*xin(88) + xc00*xin(92) + c01*xin(91)
                                      yin(96) = c10*yin(88) + yc00*yin(92) + c01*yin(91)
                                      zin(96) = c10*zin(88) + zc00*zin(92) + c01*zin(91)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   93

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   85

                                      xin(89) = xin(89) + dxij*xin(85)
                                      yin(89) = yin(89) + dyij*yin(85)
                                      zin(89) = zin(89) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   53

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   53

                                      ! do ni = 1,    3

                                      xin(53) = xin(61) + dxij*xin(49)
                                      yin(53) = yin(61) + dyij*yin(49)
                                      zin(53) = zin(61) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    2

                                      xin(65) = xin(73) + dxij*xin(61)
                                      yin(65) = yin(73) + dyij*yin(61)
                                      zin(65) = zin(73) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   77

                                      ! ni =    3

                                      xin(77) = xin(85) + dxij*xin(73)
                                      yin(77) = yin(85) + dyij*yin(73)
                                      zin(77) = zin(85) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   57

                                      ! nj =    2

                                      ! i4 = i3 =   57

                                      ! do ni = 1,    3

                                      xin(57) = xin(65) + dxij*xin(53)
                                      yin(57) = yin(65) + dyij*yin(53)
                                      zin(57) = zin(65) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    2

                                      xin(69) = xin(77) + dxij*xin(65)
                                      yin(69) = yin(77) + dyij*yin(65)
                                      zin(69) = zin(77) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(89) + dxij*xin(77)
                                      yin(81) = yin(89) + dyij*yin(77)
                                      zin(81) = zin(89) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   61

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   87

                                      xin(91) = xin(91) + dxij*xin(87)
                                      yin(91) = yin(91) + dyij*yin(87)
                                      zin(91) = zin(91) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   55

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   55

                                      ! do ni = 1,    3

                                      xin(55) = xin(63) + dxij*xin(51)
                                      yin(55) = yin(63) + dyij*yin(51)
                                      zin(55) = zin(63) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    2

                                      xin(67) = xin(75) + dxij*xin(63)
                                      yin(67) = yin(75) + dyij*yin(63)
                                      zin(67) = zin(75) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(87) + dxij*xin(75)
                                      yin(79) = yin(87) + dyij*yin(75)
                                      zin(79) = zin(87) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   59

                                      ! nj =    2

                                      ! i4 = i3 =   59

                                      ! do ni = 1,    3

                                      xin(59) = xin(67) + dxij*xin(55)
                                      yin(59) = yin(67) + dyij*yin(55)
                                      zin(59) = zin(67) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    2

                                      xin(71) = xin(79) + dxij*xin(67)
                                      yin(71) = yin(79) + dyij*yin(67)
                                      zin(71) = zin(79) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    3

                                      xin(83) = xin(91) + dxij*xin(79)
                                      yin(83) = yin(91) + dyij*yin(79)
                                      zin(83) = zin(91) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   63

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   88

                                      xin(92) = xin(92) + dxij*xin(88)
                                      yin(92) = yin(92) + dyij*yin(88)
                                      zin(92) = zin(92) + dzij*zin(88)

                                      ! i3 = i4 =   88
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   56

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   56

                                      ! do ni = 1,    3

                                      xin(56) = xin(64) + dxij*xin(52)
                                      yin(56) = yin(64) + dyij*yin(52)
                                      zin(56) = zin(64) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   68

                                      ! ni =    2

                                      xin(68) = xin(76) + dxij*xin(64)
                                      yin(68) = yin(76) + dyij*yin(64)
                                      zin(68) = zin(76) + dzij*zin(64)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    3

                                      xin(80) = xin(88) + dxij*xin(76)
                                      yin(80) = yin(88) + dyij*yin(76)
                                      zin(80) = zin(88) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   92

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   60

                                      ! nj =    2

                                      ! i4 = i3 =   60

                                      ! do ni = 1,    3

                                      xin(60) = xin(68) + dxij*xin(56)
                                      yin(60) = yin(68) + dyij*yin(56)
                                      zin(60) = zin(68) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    2

                                      xin(72) = xin(80) + dxij*xin(68)
                                      yin(72) = yin(80) + dyij*yin(68)
                                      zin(72) = zin(80) + dzij*zin(68)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    3

                                      xin(84) = xin(92) + dxij*xin(80)
                                      yin(84) = yin(92) + dyij*yin(80)
                                      zin(84) = zin(92) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   64

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   49

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

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
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(109) = xc00
                                      yin(109) = yc00
                                      zin(109) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   99

                                      xin(99) = xcp00
                                      yin(99) = ycp00
                                      zin(99) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  111
                                      ! i2 =  109

                                      xin(111) = xcp00*xin(109) + cp10
                                      yin(111) = ycp00*yin(109) + cp10
                                      zin(111) = zcp00*zin(109) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =   97
                                      ! i4 =  109

                                      xin(121) = c10*xin(97) + xc00*xin(109)
                                      yin(121) = c10*yin(97) + yc00*yin(109)
                                      zin(121) = c10*zin(97) + zc00*zin(109)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  123
                                      ! i5 =  121
                                      ! i4 =  109

                                      xin(123) = xcp00*xin(121) + cp10*xin(109)
                                      yin(123) = ycp00*yin(121) + cp10*yin(109)
                                      zin(123) = zcp00*zin(121) + cp10*zin(109)

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

                                      ! i3 = i5 + k2 =  135
                                      ! i5 =  133
                                      ! i4 =  121

                                      xin(135) = xcp00*xin(133) + cp10*xin(121)
                                      yin(135) = ycp00*yin(133) + cp10*yin(121)
                                      zin(135) = zcp00*zin(133) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  137
                                      ! i3 =  121
                                      ! i4 =  133

                                      xin(137) = c10*xin(121) + xc00*xin(133)
                                      yin(137) = c10*yin(121) + yc00*yin(133)
                                      zin(137) = c10*zin(121) + zc00*zin(133)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  139
                                      ! i5 =  137
                                      ! i4 =  133

                                      xin(139) = xcp00*xin(137) + cp10*xin(133)
                                      yin(139) = ycp00*yin(137) + cp10*yin(133)
                                      zin(139) = zcp00*zin(137) + cp10*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  141
                                      ! i3 =  133
                                      ! i4 =  137

                                      xin(141) = c10*xin(133) + xc00*xin(137)
                                      yin(141) = c10*yin(133) + yc00*yin(137)
                                      zin(141) = c10*zin(133) + zc00*zin(137)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  143
                                      ! i5 =  141
                                      ! i4 =  137

                                      xin(143) = xcp00*xin(141) + cp10*xin(137)
                                      yin(143) = ycp00*yin(141) + cp10*yin(137)
                                      zin(143) = zcp00*zin(141) + cp10*zin(137)

                                      ! ------------------

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   99

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  100
                                      ! i3 =   97
                                      ! i4 =   99

                                      xin(100) = cp01*xin(97) + xcp00*xin(99)
                                      yin(100) = cp01*yin(97) + ycp00*yin(99)
                                      zin(100) = cp01*zin(97) + zcp00*zin(99)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  112

                                      xin(112) = xc00*xin(100) + c01*xin(99)
                                      yin(112) = yc00*yin(100) + c01*yin(99)
                                      zin(112) = zc00*zin(100) + c01*zin(99)

                                      ! ------------------

                                      ! i3 = i4 =   99
                                      ! i4 = i5 =  100

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  121

                                      xin(124) = c10*xin(100) + xc00*xin(112) + c01*xin(111)
                                      yin(124) = c10*yin(100) + yc00*yin(112) + c01*yin(111)
                                      zin(124) = c10*zin(100) + zc00*zin(112) + c01*zin(111)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(136) = c10*xin(112) + xc00*xin(124) + c01*xin(123)
                                      yin(136) = c10*yin(112) + yc00*yin(124) + c01*yin(123)
                                      zin(136) = c10*zin(112) + zc00*zin(124) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(140) = c10*xin(124) + xc00*xin(136) + c01*xin(135)
                                      yin(140) = c10*yin(124) + yc00*yin(136) + c01*yin(135)
                                      zin(140) = c10*zin(124) + zc00*zin(136) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  141

                                      xin(144) = c10*xin(136) + xc00*xin(140) + c01*xin(139)
                                      yin(144) = c10*yin(136) + yc00*yin(140) + c01*yin(139)
                                      zin(144) = c10*zin(136) + zc00*zin(140) + c01*zin(139)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  141

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  133

                                      xin(137) = xin(137) + dxij*xin(133)
                                      yin(137) = yin(137) + dyij*yin(133)
                                      zin(137) = zin(137) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  101

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  101

                                      ! do ni = 1,    3

                                      xin(101) = xin(109) + dxij*xin(97)
                                      yin(101) = yin(109) + dyij*yin(97)
                                      zin(101) = zin(109) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! ni =    2

                                      xin(113) = xin(121) + dxij*xin(109)
                                      yin(113) = yin(121) + dyij*yin(109)
                                      zin(113) = zin(121) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    3

                                      xin(125) = xin(133) + dxij*xin(121)
                                      yin(125) = yin(133) + dyij*yin(121)
                                      zin(125) = zin(133) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  105

                                      ! nj =    2

                                      ! i4 = i3 =  105

                                      ! do ni = 1,    3

                                      xin(105) = xin(113) + dxij*xin(101)
                                      yin(105) = yin(113) + dyij*yin(101)
                                      zin(105) = zin(113) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    2

                                      xin(117) = xin(125) + dxij*xin(113)
                                      yin(117) = yin(125) + dyij*yin(113)
                                      zin(117) = zin(125) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    3

                                      xin(129) = xin(137) + dxij*xin(125)
                                      yin(129) = yin(137) + dyij*yin(125)
                                      zin(129) = zin(137) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  109

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  135

                                      xin(139) = xin(139) + dxij*xin(135)
                                      yin(139) = yin(139) + dyij*yin(135)
                                      zin(139) = zin(139) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  103

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(111) + dxij*xin(99)
                                      yin(103) = yin(111) + dyij*yin(99)
                                      zin(103) = zin(111) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    2

                                      xin(115) = xin(123) + dxij*xin(111)
                                      yin(115) = yin(123) + dyij*yin(111)
                                      zin(115) = zin(123) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    3

                                      xin(127) = xin(135) + dxij*xin(123)
                                      yin(127) = yin(135) + dyij*yin(123)
                                      zin(127) = zin(135) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  107

                                      ! nj =    2

                                      ! i4 = i3 =  107

                                      ! do ni = 1,    3

                                      xin(107) = xin(115) + dxij*xin(103)
                                      yin(107) = yin(115) + dyij*yin(103)
                                      zin(107) = zin(115) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    2

                                      xin(119) = xin(127) + dxij*xin(115)
                                      yin(119) = yin(127) + dyij*yin(115)
                                      zin(119) = zin(127) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      xin(131) = xin(139) + dxij*xin(127)
                                      yin(131) = yin(139) + dyij*yin(127)
                                      zin(131) = zin(139) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  111

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  136

                                      xin(140) = xin(140) + dxij*xin(136)
                                      yin(140) = yin(140) + dyij*yin(136)
                                      zin(140) = zin(140) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  104

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  104

                                      ! do ni = 1,    3

                                      xin(104) = xin(112) + dxij*xin(100)
                                      yin(104) = yin(112) + dyij*yin(100)
                                      zin(104) = zin(112) + dzij*zin(100)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  116

                                      ! ni =    2

                                      xin(116) = xin(124) + dxij*xin(112)
                                      yin(116) = yin(124) + dyij*yin(112)
                                      zin(116) = zin(124) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    3

                                      xin(128) = xin(136) + dxij*xin(124)
                                      yin(128) = yin(136) + dyij*yin(124)
                                      zin(128) = zin(136) + dzij*zin(124)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  140

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  108

                                      ! nj =    2

                                      ! i4 = i3 =  108

                                      ! do ni = 1,    3

                                      xin(108) = xin(116) + dxij*xin(104)
                                      yin(108) = yin(116) + dyij*yin(104)
                                      zin(108) = zin(116) + dzij*zin(104)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    2

                                      xin(120) = xin(128) + dxij*xin(116)
                                      yin(120) = yin(128) + dyij*yin(116)
                                      zin(120) = zin(128) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    3

                                      xin(132) = xin(140) + dxij*xin(128)
                                      yin(132) = yin(140) + dyij*yin(128)
                                      zin(132) = zin(140) + dzij*zin(128)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  112

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   97

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

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
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(157) = xc00
                                      yin(157) = yc00
                                      zin(157) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  147

                                      xin(147) = xcp00
                                      yin(147) = ycp00
                                      zin(147) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  159
                                      ! i2 =  157

                                      xin(159) = xcp00*xin(157) + cp10
                                      yin(159) = ycp00*yin(157) + cp10
                                      zin(159) = zcp00*zin(157) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  169
                                      ! i3 =  145
                                      ! i4 =  157

                                      xin(169) = c10*xin(145) + xc00*xin(157)
                                      yin(169) = c10*yin(145) + yc00*yin(157)
                                      zin(169) = c10*zin(145) + zc00*zin(157)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  171
                                      ! i5 =  169
                                      ! i4 =  157

                                      xin(171) = xcp00*xin(169) + cp10*xin(157)
                                      yin(171) = ycp00*yin(169) + cp10*yin(157)
                                      zin(171) = zcp00*zin(169) + cp10*zin(157)

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

                                      ! i3 = i5 + k2 =  183
                                      ! i5 =  181
                                      ! i4 =  169

                                      xin(183) = xcp00*xin(181) + cp10*xin(169)
                                      yin(183) = ycp00*yin(181) + cp10*yin(169)
                                      zin(183) = zcp00*zin(181) + cp10*zin(169)

                                      ! ------------------

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  185
                                      ! i3 =  169
                                      ! i4 =  181

                                      xin(185) = c10*xin(169) + xc00*xin(181)
                                      yin(185) = c10*yin(169) + yc00*yin(181)
                                      zin(185) = c10*zin(169) + zc00*zin(181)

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

                                      ! n =    5

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

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  147

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  148
                                      ! i3 =  145
                                      ! i4 =  147

                                      xin(148) = cp01*xin(145) + xcp00*xin(147)
                                      yin(148) = cp01*yin(145) + ycp00*yin(147)
                                      zin(148) = cp01*zin(145) + zcp00*zin(147)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  160

                                      xin(160) = xc00*xin(148) + c01*xin(147)
                                      yin(160) = yc00*yin(148) + c01*yin(147)
                                      zin(160) = zc00*zin(148) + c01*zin(147)

                                      ! ------------------

                                      ! i3 = i4 =  147
                                      ! i4 = i5 =  148

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  169

                                      xin(172) = c10*xin(148) + xc00*xin(160) + c01*xin(159)
                                      yin(172) = c10*yin(148) + yc00*yin(160) + c01*yin(159)
                                      zin(172) = c10*zin(148) + zc00*zin(160) + c01*zin(159)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(184) = c10*xin(160) + xc00*xin(172) + c01*xin(171)
                                      yin(184) = c10*yin(160) + yc00*yin(172) + c01*yin(171)
                                      zin(184) = c10*zin(160) + zc00*zin(172) + c01*zin(171)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  185

                                      xin(188) = c10*xin(172) + xc00*xin(184) + c01*xin(183)
                                      yin(188) = c10*yin(172) + yc00*yin(184) + c01*yin(183)
                                      zin(188) = c10*zin(172) + zc00*zin(184) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  189

                                      xin(192) = c10*xin(184) + xc00*xin(188) + c01*xin(187)
                                      yin(192) = c10*yin(184) + yc00*yin(188) + c01*yin(187)
                                      zin(192) = c10*zin(184) + zc00*zin(188) + c01*zin(187)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! nn =    6

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

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  181

                                      xin(185) = xin(185) + dxij*xin(181)
                                      yin(185) = yin(185) + dyij*yin(181)
                                      zin(185) = zin(185) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  149

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  149

                                      ! do ni = 1,    3

                                      xin(149) = xin(157) + dxij*xin(145)
                                      yin(149) = yin(157) + dyij*yin(145)
                                      zin(149) = zin(157) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    2

                                      xin(161) = xin(169) + dxij*xin(157)
                                      yin(161) = yin(169) + dyij*yin(157)
                                      zin(161) = zin(169) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    3

                                      xin(173) = xin(181) + dxij*xin(169)
                                      yin(173) = yin(181) + dyij*yin(169)
                                      zin(173) = zin(181) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  153

                                      ! nj =    2

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(161) + dxij*xin(149)
                                      yin(153) = yin(161) + dyij*yin(149)
                                      zin(153) = zin(161) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    2

                                      xin(165) = xin(173) + dxij*xin(161)
                                      yin(165) = yin(173) + dyij*yin(161)
                                      zin(165) = zin(173) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    3

                                      xin(177) = xin(185) + dxij*xin(173)
                                      yin(177) = yin(185) + dyij*yin(173)
                                      zin(177) = zin(185) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  183

                                      xin(187) = xin(187) + dxij*xin(183)
                                      yin(187) = yin(187) + dyij*yin(183)
                                      zin(187) = zin(187) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  151

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(159) + dxij*xin(147)
                                      yin(151) = yin(159) + dyij*yin(147)
                                      zin(151) = zin(159) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! ni =    2

                                      xin(163) = xin(171) + dxij*xin(159)
                                      yin(163) = yin(171) + dyij*yin(159)
                                      zin(163) = zin(171) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(183) + dxij*xin(171)
                                      yin(175) = yin(183) + dyij*yin(171)
                                      zin(175) = zin(183) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  155

                                      ! nj =    2

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(163) + dxij*xin(151)
                                      yin(155) = yin(163) + dyij*yin(151)
                                      zin(155) = zin(163) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    2

                                      xin(167) = xin(175) + dxij*xin(163)
                                      yin(167) = yin(175) + dyij*yin(163)
                                      zin(167) = zin(175) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    3

                                      xin(179) = xin(187) + dxij*xin(175)
                                      yin(179) = yin(187) + dyij*yin(175)
                                      zin(179) = zin(187) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  184

                                      xin(188) = xin(188) + dxij*xin(184)
                                      yin(188) = yin(188) + dyij*yin(184)
                                      zin(188) = zin(188) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  152

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  152

                                      ! do ni = 1,    3

                                      xin(152) = xin(160) + dxij*xin(148)
                                      yin(152) = yin(160) + dyij*yin(148)
                                      zin(152) = zin(160) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  164

                                      ! ni =    2

                                      xin(164) = xin(172) + dxij*xin(160)
                                      yin(164) = yin(172) + dyij*yin(160)
                                      zin(164) = zin(172) + dzij*zin(160)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  176

                                      ! ni =    3

                                      xin(176) = xin(184) + dxij*xin(172)
                                      yin(176) = yin(184) + dyij*yin(172)
                                      zin(176) = zin(184) + dzij*zin(172)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  188

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  156

                                      ! nj =    2

                                      ! i4 = i3 =  156

                                      ! do ni = 1,    3

                                      xin(156) = xin(164) + dxij*xin(152)
                                      yin(156) = yin(164) + dyij*yin(152)
                                      zin(156) = zin(164) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    2

                                      xin(168) = xin(176) + dxij*xin(164)
                                      yin(168) = yin(176) + dyij*yin(164)
                                      zin(168) = zin(176) + dzij*zin(164)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    3

                                      xin(180) = xin(188) + dxij*xin(176)
                                      yin(180) = yin(188) + dyij*yin(176)
                                      zin(180) = zin(188) + dzij*zin(176)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  160

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =  145

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d23bra(  1)*d11ket(  1)*(xin(  48)*yin(   1)*zin(   1)+xin(  96)*yin(  49)*zin(  49)+xin( 144)*yin(  97)*zin(  97)+xin( 192)*yin( 145)*zin( 145))
          eri_value(    2)=eri_value(    2)+d23bra(  1)*d11ket(  2)*(xin(  47)*yin(   2)*zin(   1)+xin(  95)*yin(  50)*zin(  49)+xin( 143)*yin(  98)*zin(  97)+xin( 191)*yin( 146)*zin( 145))
          eri_value(    3)=eri_value(    3)+d23bra(  1)*d11ket(  3)*(xin(  47)*yin(   1)*zin(   2)+xin(  95)*yin(  49)*zin(  50)+xin( 143)*yin(  97)*zin(  98)+xin( 191)*yin( 145)*zin( 146))
          eri_value(    4)=eri_value(    4)+d23bra(  1)*d11ket(  4)*(xin(  46)*yin(   3)*zin(   1)+xin(  94)*yin(  51)*zin(  49)+xin( 142)*yin(  99)*zin(  97)+xin( 190)*yin( 147)*zin( 145))
          eri_value(    5)=eri_value(    5)+d23bra(  1)*d11ket(  5)*(xin(  45)*yin(   4)*zin(   1)+xin(  93)*yin(  52)*zin(  49)+xin( 141)*yin( 100)*zin(  97)+xin( 189)*yin( 148)*zin( 145))
          eri_value(    6)=eri_value(    6)+d23bra(  1)*d11ket(  6)*(xin(  45)*yin(   3)*zin(   2)+xin(  93)*yin(  51)*zin(  50)+xin( 141)*yin(  99)*zin(  98)+xin( 189)*yin( 147)*zin( 146))
          eri_value(    7)=eri_value(    7)+d23bra(  1)*d11ket(  7)*(xin(  46)*yin(   1)*zin(   3)+xin(  94)*yin(  49)*zin(  51)+xin( 142)*yin(  97)*zin(  99)+xin( 190)*yin( 145)*zin( 147))
          eri_value(    8)=eri_value(    8)+d23bra(  1)*d11ket(  8)*(xin(  45)*yin(   2)*zin(   3)+xin(  93)*yin(  50)*zin(  51)+xin( 141)*yin(  98)*zin(  99)+xin( 189)*yin( 146)*zin( 147))
          eri_value(    9)=eri_value(    9)+d23bra(  1)*d11ket(  9)*(xin(  45)*yin(   1)*zin(   4)+xin(  93)*yin(  49)*zin(  52)+xin( 141)*yin(  97)*zin( 100)+xin( 189)*yin( 145)*zin( 148))
          eri_value(   10)=eri_value(   10)+d23bra(  2)*d11ket(  1)*(xin(  40)*yin(   9)*zin(   1)+xin(  88)*yin(  57)*zin(  49)+xin( 136)*yin( 105)*zin(  97)+xin( 184)*yin( 153)*zin( 145))
          eri_value(   11)=eri_value(   11)+d23bra(  2)*d11ket(  2)*(xin(  39)*yin(  10)*zin(   1)+xin(  87)*yin(  58)*zin(  49)+xin( 135)*yin( 106)*zin(  97)+xin( 183)*yin( 154)*zin( 145))
          eri_value(   12)=eri_value(   12)+d23bra(  2)*d11ket(  3)*(xin(  39)*yin(   9)*zin(   2)+xin(  87)*yin(  57)*zin(  50)+xin( 135)*yin( 105)*zin(  98)+xin( 183)*yin( 153)*zin( 146))
          eri_value(   13)=eri_value(   13)+d23bra(  2)*d11ket(  4)*(xin(  38)*yin(  11)*zin(   1)+xin(  86)*yin(  59)*zin(  49)+xin( 134)*yin( 107)*zin(  97)+xin( 182)*yin( 155)*zin( 145))
          eri_value(   14)=eri_value(   14)+d23bra(  2)*d11ket(  5)*(xin(  37)*yin(  12)*zin(   1)+xin(  85)*yin(  60)*zin(  49)+xin( 133)*yin( 108)*zin(  97)+xin( 181)*yin( 156)*zin( 145))
          eri_value(   15)=eri_value(   15)+d23bra(  2)*d11ket(  6)*(xin(  37)*yin(  11)*zin(   2)+xin(  85)*yin(  59)*zin(  50)+xin( 133)*yin( 107)*zin(  98)+xin( 181)*yin( 155)*zin( 146))
          eri_value(   16)=eri_value(   16)+d23bra(  2)*d11ket(  7)*(xin(  38)*yin(   9)*zin(   3)+xin(  86)*yin(  57)*zin(  51)+xin( 134)*yin( 105)*zin(  99)+xin( 182)*yin( 153)*zin( 147))
          eri_value(   17)=eri_value(   17)+d23bra(  2)*d11ket(  8)*(xin(  37)*yin(  10)*zin(   3)+xin(  85)*yin(  58)*zin(  51)+xin( 133)*yin( 106)*zin(  99)+xin( 181)*yin( 154)*zin( 147))
          eri_value(   18)=eri_value(   18)+d23bra(  2)*d11ket(  9)*(xin(  37)*yin(   9)*zin(   4)+xin(  85)*yin(  57)*zin(  52)+xin( 133)*yin( 105)*zin( 100)+xin( 181)*yin( 153)*zin( 148))
          eri_value(   19)=eri_value(   19)+d23bra(  3)*d11ket(  1)*(xin(  40)*yin(   1)*zin(   9)+xin(  88)*yin(  49)*zin(  57)+xin( 136)*yin(  97)*zin( 105)+xin( 184)*yin( 145)*zin( 153))
          eri_value(   20)=eri_value(   20)+d23bra(  3)*d11ket(  2)*(xin(  39)*yin(   2)*zin(   9)+xin(  87)*yin(  50)*zin(  57)+xin( 135)*yin(  98)*zin( 105)+xin( 183)*yin( 146)*zin( 153))
          eri_value(   21)=eri_value(   21)+d23bra(  3)*d11ket(  3)*(xin(  39)*yin(   1)*zin(  10)+xin(  87)*yin(  49)*zin(  58)+xin( 135)*yin(  97)*zin( 106)+xin( 183)*yin( 145)*zin( 154))
          eri_value(   22)=eri_value(   22)+d23bra(  3)*d11ket(  4)*(xin(  38)*yin(   3)*zin(   9)+xin(  86)*yin(  51)*zin(  57)+xin( 134)*yin(  99)*zin( 105)+xin( 182)*yin( 147)*zin( 153))
          eri_value(   23)=eri_value(   23)+d23bra(  3)*d11ket(  5)*(xin(  37)*yin(   4)*zin(   9)+xin(  85)*yin(  52)*zin(  57)+xin( 133)*yin( 100)*zin( 105)+xin( 181)*yin( 148)*zin( 153))
          eri_value(   24)=eri_value(   24)+d23bra(  3)*d11ket(  6)*(xin(  37)*yin(   3)*zin(  10)+xin(  85)*yin(  51)*zin(  58)+xin( 133)*yin(  99)*zin( 106)+xin( 181)*yin( 147)*zin( 154))
          eri_value(   25)=eri_value(   25)+d23bra(  3)*d11ket(  7)*(xin(  38)*yin(   1)*zin(  11)+xin(  86)*yin(  49)*zin(  59)+xin( 134)*yin(  97)*zin( 107)+xin( 182)*yin( 145)*zin( 155))
          eri_value(   26)=eri_value(   26)+d23bra(  3)*d11ket(  8)*(xin(  37)*yin(   2)*zin(  11)+xin(  85)*yin(  50)*zin(  59)+xin( 133)*yin(  98)*zin( 107)+xin( 181)*yin( 146)*zin( 155))
          eri_value(   27)=eri_value(   27)+d23bra(  3)*d11ket(  9)*(xin(  37)*yin(   1)*zin(  12)+xin(  85)*yin(  49)*zin(  60)+xin( 133)*yin(  97)*zin( 108)+xin( 181)*yin( 145)*zin( 156))
          eri_value(   28)=eri_value(   28)+d23bra(  4)*d11ket(  1)*(xin(  44)*yin(   5)*zin(   1)+xin(  92)*yin(  53)*zin(  49)+xin( 140)*yin( 101)*zin(  97)+xin( 188)*yin( 149)*zin( 145))
          eri_value(   29)=eri_value(   29)+d23bra(  4)*d11ket(  2)*(xin(  43)*yin(   6)*zin(   1)+xin(  91)*yin(  54)*zin(  49)+xin( 139)*yin( 102)*zin(  97)+xin( 187)*yin( 150)*zin( 145))
          eri_value(   30)=eri_value(   30)+d23bra(  4)*d11ket(  3)*(xin(  43)*yin(   5)*zin(   2)+xin(  91)*yin(  53)*zin(  50)+xin( 139)*yin( 101)*zin(  98)+xin( 187)*yin( 149)*zin( 146))
          eri_value(   31)=eri_value(   31)+d23bra(  4)*d11ket(  4)*(xin(  42)*yin(   7)*zin(   1)+xin(  90)*yin(  55)*zin(  49)+xin( 138)*yin( 103)*zin(  97)+xin( 186)*yin( 151)*zin( 145))
          eri_value(   32)=eri_value(   32)+d23bra(  4)*d11ket(  5)*(xin(  41)*yin(   8)*zin(   1)+xin(  89)*yin(  56)*zin(  49)+xin( 137)*yin( 104)*zin(  97)+xin( 185)*yin( 152)*zin( 145))
          eri_value(   33)=eri_value(   33)+d23bra(  4)*d11ket(  6)*(xin(  41)*yin(   7)*zin(   2)+xin(  89)*yin(  55)*zin(  50)+xin( 137)*yin( 103)*zin(  98)+xin( 185)*yin( 151)*zin( 146))
          eri_value(   34)=eri_value(   34)+d23bra(  4)*d11ket(  7)*(xin(  42)*yin(   5)*zin(   3)+xin(  90)*yin(  53)*zin(  51)+xin( 138)*yin( 101)*zin(  99)+xin( 186)*yin( 149)*zin( 147))
          eri_value(   35)=eri_value(   35)+d23bra(  4)*d11ket(  8)*(xin(  41)*yin(   6)*zin(   3)+xin(  89)*yin(  54)*zin(  51)+xin( 137)*yin( 102)*zin(  99)+xin( 185)*yin( 150)*zin( 147))
          eri_value(   36)=eri_value(   36)+d23bra(  4)*d11ket(  9)*(xin(  41)*yin(   5)*zin(   4)+xin(  89)*yin(  53)*zin(  52)+xin( 137)*yin( 101)*zin( 100)+xin( 185)*yin( 149)*zin( 148))
          eri_value(   37)=eri_value(   37)+d23bra(  5)*d11ket(  1)*(xin(  44)*yin(   1)*zin(   5)+xin(  92)*yin(  49)*zin(  53)+xin( 140)*yin(  97)*zin( 101)+xin( 188)*yin( 145)*zin( 149))
          eri_value(   38)=eri_value(   38)+d23bra(  5)*d11ket(  2)*(xin(  43)*yin(   2)*zin(   5)+xin(  91)*yin(  50)*zin(  53)+xin( 139)*yin(  98)*zin( 101)+xin( 187)*yin( 146)*zin( 149))
          eri_value(   39)=eri_value(   39)+d23bra(  5)*d11ket(  3)*(xin(  43)*yin(   1)*zin(   6)+xin(  91)*yin(  49)*zin(  54)+xin( 139)*yin(  97)*zin( 102)+xin( 187)*yin( 145)*zin( 150))
          eri_value(   40)=eri_value(   40)+d23bra(  5)*d11ket(  4)*(xin(  42)*yin(   3)*zin(   5)+xin(  90)*yin(  51)*zin(  53)+xin( 138)*yin(  99)*zin( 101)+xin( 186)*yin( 147)*zin( 149))
          eri_value(   41)=eri_value(   41)+d23bra(  5)*d11ket(  5)*(xin(  41)*yin(   4)*zin(   5)+xin(  89)*yin(  52)*zin(  53)+xin( 137)*yin( 100)*zin( 101)+xin( 185)*yin( 148)*zin( 149))
          eri_value(   42)=eri_value(   42)+d23bra(  5)*d11ket(  6)*(xin(  41)*yin(   3)*zin(   6)+xin(  89)*yin(  51)*zin(  54)+xin( 137)*yin(  99)*zin( 102)+xin( 185)*yin( 147)*zin( 150))
          eri_value(   43)=eri_value(   43)+d23bra(  5)*d11ket(  7)*(xin(  42)*yin(   1)*zin(   7)+xin(  90)*yin(  49)*zin(  55)+xin( 138)*yin(  97)*zin( 103)+xin( 186)*yin( 145)*zin( 151))
          eri_value(   44)=eri_value(   44)+d23bra(  5)*d11ket(  8)*(xin(  41)*yin(   2)*zin(   7)+xin(  89)*yin(  50)*zin(  55)+xin( 137)*yin(  98)*zin( 103)+xin( 185)*yin( 146)*zin( 151))
          eri_value(   45)=eri_value(   45)+d23bra(  5)*d11ket(  9)*(xin(  41)*yin(   1)*zin(   8)+xin(  89)*yin(  49)*zin(  56)+xin( 137)*yin(  97)*zin( 104)+xin( 185)*yin( 145)*zin( 152))
          eri_value(   46)=eri_value(   46)+d23bra(  6)*d11ket(  1)*(xin(  40)*yin(   5)*zin(   5)+xin(  88)*yin(  53)*zin(  53)+xin( 136)*yin( 101)*zin( 101)+xin( 184)*yin( 149)*zin( 149))
          eri_value(   47)=eri_value(   47)+d23bra(  6)*d11ket(  2)*(xin(  39)*yin(   6)*zin(   5)+xin(  87)*yin(  54)*zin(  53)+xin( 135)*yin( 102)*zin( 101)+xin( 183)*yin( 150)*zin( 149))
          eri_value(   48)=eri_value(   48)+d23bra(  6)*d11ket(  3)*(xin(  39)*yin(   5)*zin(   6)+xin(  87)*yin(  53)*zin(  54)+xin( 135)*yin( 101)*zin( 102)+xin( 183)*yin( 149)*zin( 150))
          eri_value(   49)=eri_value(   49)+d23bra(  6)*d11ket(  4)*(xin(  38)*yin(   7)*zin(   5)+xin(  86)*yin(  55)*zin(  53)+xin( 134)*yin( 103)*zin( 101)+xin( 182)*yin( 151)*zin( 149))
          eri_value(   50)=eri_value(   50)+d23bra(  6)*d11ket(  5)*(xin(  37)*yin(   8)*zin(   5)+xin(  85)*yin(  56)*zin(  53)+xin( 133)*yin( 104)*zin( 101)+xin( 181)*yin( 152)*zin( 149))
          eri_value(   51)=eri_value(   51)+d23bra(  6)*d11ket(  6)*(xin(  37)*yin(   7)*zin(   6)+xin(  85)*yin(  55)*zin(  54)+xin( 133)*yin( 103)*zin( 102)+xin( 181)*yin( 151)*zin( 150))
          eri_value(   52)=eri_value(   52)+d23bra(  6)*d11ket(  7)*(xin(  38)*yin(   5)*zin(   7)+xin(  86)*yin(  53)*zin(  55)+xin( 134)*yin( 101)*zin( 103)+xin( 182)*yin( 149)*zin( 151))
          eri_value(   53)=eri_value(   53)+d23bra(  6)*d11ket(  8)*(xin(  37)*yin(   6)*zin(   7)+xin(  85)*yin(  54)*zin(  55)+xin( 133)*yin( 102)*zin( 103)+xin( 181)*yin( 150)*zin( 151))
          eri_value(   54)=eri_value(   54)+d23bra(  6)*d11ket(  9)*(xin(  37)*yin(   5)*zin(   8)+xin(  85)*yin(  53)*zin(  56)+xin( 133)*yin( 101)*zin( 104)+xin( 181)*yin( 149)*zin( 152))
          eri_value(   55)=eri_value(   55)+d23bra(  7)*d11ket(  1)*(xin(  12)*yin(  37)*zin(   1)+xin(  60)*yin(  85)*zin(  49)+xin( 108)*yin( 133)*zin(  97)+xin( 156)*yin( 181)*zin( 145))
          eri_value(   56)=eri_value(   56)+d23bra(  7)*d11ket(  2)*(xin(  11)*yin(  38)*zin(   1)+xin(  59)*yin(  86)*zin(  49)+xin( 107)*yin( 134)*zin(  97)+xin( 155)*yin( 182)*zin( 145))
          eri_value(   57)=eri_value(   57)+d23bra(  7)*d11ket(  3)*(xin(  11)*yin(  37)*zin(   2)+xin(  59)*yin(  85)*zin(  50)+xin( 107)*yin( 133)*zin(  98)+xin( 155)*yin( 181)*zin( 146))
          eri_value(   58)=eri_value(   58)+d23bra(  7)*d11ket(  4)*(xin(  10)*yin(  39)*zin(   1)+xin(  58)*yin(  87)*zin(  49)+xin( 106)*yin( 135)*zin(  97)+xin( 154)*yin( 183)*zin( 145))
          eri_value(   59)=eri_value(   59)+d23bra(  7)*d11ket(  5)*(xin(   9)*yin(  40)*zin(   1)+xin(  57)*yin(  88)*zin(  49)+xin( 105)*yin( 136)*zin(  97)+xin( 153)*yin( 184)*zin( 145))
          eri_value(   60)=eri_value(   60)+d23bra(  7)*d11ket(  6)*(xin(   9)*yin(  39)*zin(   2)+xin(  57)*yin(  87)*zin(  50)+xin( 105)*yin( 135)*zin(  98)+xin( 153)*yin( 183)*zin( 146))
          eri_value(   61)=eri_value(   61)+d23bra(  7)*d11ket(  7)*(xin(  10)*yin(  37)*zin(   3)+xin(  58)*yin(  85)*zin(  51)+xin( 106)*yin( 133)*zin(  99)+xin( 154)*yin( 181)*zin( 147))
          eri_value(   62)=eri_value(   62)+d23bra(  7)*d11ket(  8)*(xin(   9)*yin(  38)*zin(   3)+xin(  57)*yin(  86)*zin(  51)+xin( 105)*yin( 134)*zin(  99)+xin( 153)*yin( 182)*zin( 147))
          eri_value(   63)=eri_value(   63)+d23bra(  7)*d11ket(  9)*(xin(   9)*yin(  37)*zin(   4)+xin(  57)*yin(  85)*zin(  52)+xin( 105)*yin( 133)*zin( 100)+xin( 153)*yin( 181)*zin( 148))
          eri_value(   64)=eri_value(   64)+d23bra(  8)*d11ket(  1)*(xin(   4)*yin(  45)*zin(   1)+xin(  52)*yin(  93)*zin(  49)+xin( 100)*yin( 141)*zin(  97)+xin( 148)*yin( 189)*zin( 145))
          eri_value(   65)=eri_value(   65)+d23bra(  8)*d11ket(  2)*(xin(   3)*yin(  46)*zin(   1)+xin(  51)*yin(  94)*zin(  49)+xin(  99)*yin( 142)*zin(  97)+xin( 147)*yin( 190)*zin( 145))
          eri_value(   66)=eri_value(   66)+d23bra(  8)*d11ket(  3)*(xin(   3)*yin(  45)*zin(   2)+xin(  51)*yin(  93)*zin(  50)+xin(  99)*yin( 141)*zin(  98)+xin( 147)*yin( 189)*zin( 146))
          eri_value(   67)=eri_value(   67)+d23bra(  8)*d11ket(  4)*(xin(   2)*yin(  47)*zin(   1)+xin(  50)*yin(  95)*zin(  49)+xin(  98)*yin( 143)*zin(  97)+xin( 146)*yin( 191)*zin( 145))
          eri_value(   68)=eri_value(   68)+d23bra(  8)*d11ket(  5)*(xin(   1)*yin(  48)*zin(   1)+xin(  49)*yin(  96)*zin(  49)+xin(  97)*yin( 144)*zin(  97)+xin( 145)*yin( 192)*zin( 145))
          eri_value(   69)=eri_value(   69)+d23bra(  8)*d11ket(  6)*(xin(   1)*yin(  47)*zin(   2)+xin(  49)*yin(  95)*zin(  50)+xin(  97)*yin( 143)*zin(  98)+xin( 145)*yin( 191)*zin( 146))
          eri_value(   70)=eri_value(   70)+d23bra(  8)*d11ket(  7)*(xin(   2)*yin(  45)*zin(   3)+xin(  50)*yin(  93)*zin(  51)+xin(  98)*yin( 141)*zin(  99)+xin( 146)*yin( 189)*zin( 147))
          eri_value(   71)=eri_value(   71)+d23bra(  8)*d11ket(  8)*(xin(   1)*yin(  46)*zin(   3)+xin(  49)*yin(  94)*zin(  51)+xin(  97)*yin( 142)*zin(  99)+xin( 145)*yin( 190)*zin( 147))
          eri_value(   72)=eri_value(   72)+d23bra(  8)*d11ket(  9)*(xin(   1)*yin(  45)*zin(   4)+xin(  49)*yin(  93)*zin(  52)+xin(  97)*yin( 141)*zin( 100)+xin( 145)*yin( 189)*zin( 148))
          eri_value(   73)=eri_value(   73)+d23bra(  9)*d11ket(  1)*(xin(   4)*yin(  37)*zin(   9)+xin(  52)*yin(  85)*zin(  57)+xin( 100)*yin( 133)*zin( 105)+xin( 148)*yin( 181)*zin( 153))
          eri_value(   74)=eri_value(   74)+d23bra(  9)*d11ket(  2)*(xin(   3)*yin(  38)*zin(   9)+xin(  51)*yin(  86)*zin(  57)+xin(  99)*yin( 134)*zin( 105)+xin( 147)*yin( 182)*zin( 153))
          eri_value(   75)=eri_value(   75)+d23bra(  9)*d11ket(  3)*(xin(   3)*yin(  37)*zin(  10)+xin(  51)*yin(  85)*zin(  58)+xin(  99)*yin( 133)*zin( 106)+xin( 147)*yin( 181)*zin( 154))
          eri_value(   76)=eri_value(   76)+d23bra(  9)*d11ket(  4)*(xin(   2)*yin(  39)*zin(   9)+xin(  50)*yin(  87)*zin(  57)+xin(  98)*yin( 135)*zin( 105)+xin( 146)*yin( 183)*zin( 153))
          eri_value(   77)=eri_value(   77)+d23bra(  9)*d11ket(  5)*(xin(   1)*yin(  40)*zin(   9)+xin(  49)*yin(  88)*zin(  57)+xin(  97)*yin( 136)*zin( 105)+xin( 145)*yin( 184)*zin( 153))
          eri_value(   78)=eri_value(   78)+d23bra(  9)*d11ket(  6)*(xin(   1)*yin(  39)*zin(  10)+xin(  49)*yin(  87)*zin(  58)+xin(  97)*yin( 135)*zin( 106)+xin( 145)*yin( 183)*zin( 154))
          eri_value(   79)=eri_value(   79)+d23bra(  9)*d11ket(  7)*(xin(   2)*yin(  37)*zin(  11)+xin(  50)*yin(  85)*zin(  59)+xin(  98)*yin( 133)*zin( 107)+xin( 146)*yin( 181)*zin( 155))
          eri_value(   80)=eri_value(   80)+d23bra(  9)*d11ket(  8)*(xin(   1)*yin(  38)*zin(  11)+xin(  49)*yin(  86)*zin(  59)+xin(  97)*yin( 134)*zin( 107)+xin( 145)*yin( 182)*zin( 155))
          eri_value(   81)=eri_value(   81)+d23bra(  9)*d11ket(  9)*(xin(   1)*yin(  37)*zin(  12)+xin(  49)*yin(  85)*zin(  60)+xin(  97)*yin( 133)*zin( 108)+xin( 145)*yin( 181)*zin( 156))
          eri_value(   82)=eri_value(   82)+d23bra( 10)*d11ket(  1)*(xin(   8)*yin(  41)*zin(   1)+xin(  56)*yin(  89)*zin(  49)+xin( 104)*yin( 137)*zin(  97)+xin( 152)*yin( 185)*zin( 145))
          eri_value(   83)=eri_value(   83)+d23bra( 10)*d11ket(  2)*(xin(   7)*yin(  42)*zin(   1)+xin(  55)*yin(  90)*zin(  49)+xin( 103)*yin( 138)*zin(  97)+xin( 151)*yin( 186)*zin( 145))
          eri_value(   84)=eri_value(   84)+d23bra( 10)*d11ket(  3)*(xin(   7)*yin(  41)*zin(   2)+xin(  55)*yin(  89)*zin(  50)+xin( 103)*yin( 137)*zin(  98)+xin( 151)*yin( 185)*zin( 146))
          eri_value(   85)=eri_value(   85)+d23bra( 10)*d11ket(  4)*(xin(   6)*yin(  43)*zin(   1)+xin(  54)*yin(  91)*zin(  49)+xin( 102)*yin( 139)*zin(  97)+xin( 150)*yin( 187)*zin( 145))
          eri_value(   86)=eri_value(   86)+d23bra( 10)*d11ket(  5)*(xin(   5)*yin(  44)*zin(   1)+xin(  53)*yin(  92)*zin(  49)+xin( 101)*yin( 140)*zin(  97)+xin( 149)*yin( 188)*zin( 145))
          eri_value(   87)=eri_value(   87)+d23bra( 10)*d11ket(  6)*(xin(   5)*yin(  43)*zin(   2)+xin(  53)*yin(  91)*zin(  50)+xin( 101)*yin( 139)*zin(  98)+xin( 149)*yin( 187)*zin( 146))
          eri_value(   88)=eri_value(   88)+d23bra( 10)*d11ket(  7)*(xin(   6)*yin(  41)*zin(   3)+xin(  54)*yin(  89)*zin(  51)+xin( 102)*yin( 137)*zin(  99)+xin( 150)*yin( 185)*zin( 147))
          eri_value(   89)=eri_value(   89)+d23bra( 10)*d11ket(  8)*(xin(   5)*yin(  42)*zin(   3)+xin(  53)*yin(  90)*zin(  51)+xin( 101)*yin( 138)*zin(  99)+xin( 149)*yin( 186)*zin( 147))
          eri_value(   90)=eri_value(   90)+d23bra( 10)*d11ket(  9)*(xin(   5)*yin(  41)*zin(   4)+xin(  53)*yin(  89)*zin(  52)+xin( 101)*yin( 137)*zin( 100)+xin( 149)*yin( 185)*zin( 148))
          eri_value(   91)=eri_value(   91)+d23bra( 11)*d11ket(  1)*(xin(   8)*yin(  37)*zin(   5)+xin(  56)*yin(  85)*zin(  53)+xin( 104)*yin( 133)*zin( 101)+xin( 152)*yin( 181)*zin( 149))
          eri_value(   92)=eri_value(   92)+d23bra( 11)*d11ket(  2)*(xin(   7)*yin(  38)*zin(   5)+xin(  55)*yin(  86)*zin(  53)+xin( 103)*yin( 134)*zin( 101)+xin( 151)*yin( 182)*zin( 149))
          eri_value(   93)=eri_value(   93)+d23bra( 11)*d11ket(  3)*(xin(   7)*yin(  37)*zin(   6)+xin(  55)*yin(  85)*zin(  54)+xin( 103)*yin( 133)*zin( 102)+xin( 151)*yin( 181)*zin( 150))
          eri_value(   94)=eri_value(   94)+d23bra( 11)*d11ket(  4)*(xin(   6)*yin(  39)*zin(   5)+xin(  54)*yin(  87)*zin(  53)+xin( 102)*yin( 135)*zin( 101)+xin( 150)*yin( 183)*zin( 149))
          eri_value(   95)=eri_value(   95)+d23bra( 11)*d11ket(  5)*(xin(   5)*yin(  40)*zin(   5)+xin(  53)*yin(  88)*zin(  53)+xin( 101)*yin( 136)*zin( 101)+xin( 149)*yin( 184)*zin( 149))
          eri_value(   96)=eri_value(   96)+d23bra( 11)*d11ket(  6)*(xin(   5)*yin(  39)*zin(   6)+xin(  53)*yin(  87)*zin(  54)+xin( 101)*yin( 135)*zin( 102)+xin( 149)*yin( 183)*zin( 150))
          eri_value(   97)=eri_value(   97)+d23bra( 11)*d11ket(  7)*(xin(   6)*yin(  37)*zin(   7)+xin(  54)*yin(  85)*zin(  55)+xin( 102)*yin( 133)*zin( 103)+xin( 150)*yin( 181)*zin( 151))
          eri_value(   98)=eri_value(   98)+d23bra( 11)*d11ket(  8)*(xin(   5)*yin(  38)*zin(   7)+xin(  53)*yin(  86)*zin(  55)+xin( 101)*yin( 134)*zin( 103)+xin( 149)*yin( 182)*zin( 151))
          eri_value(   99)=eri_value(   99)+d23bra( 11)*d11ket(  9)*(xin(   5)*yin(  37)*zin(   8)+xin(  53)*yin(  85)*zin(  56)+xin( 101)*yin( 133)*zin( 104)+xin( 149)*yin( 181)*zin( 152))
          eri_value(  100)=eri_value(  100)+d23bra( 12)*d11ket(  1)*(xin(   4)*yin(  41)*zin(   5)+xin(  52)*yin(  89)*zin(  53)+xin( 100)*yin( 137)*zin( 101)+xin( 148)*yin( 185)*zin( 149))
          eri_value(  101)=eri_value(  101)+d23bra( 12)*d11ket(  2)*(xin(   3)*yin(  42)*zin(   5)+xin(  51)*yin(  90)*zin(  53)+xin(  99)*yin( 138)*zin( 101)+xin( 147)*yin( 186)*zin( 149))
          eri_value(  102)=eri_value(  102)+d23bra( 12)*d11ket(  3)*(xin(   3)*yin(  41)*zin(   6)+xin(  51)*yin(  89)*zin(  54)+xin(  99)*yin( 137)*zin( 102)+xin( 147)*yin( 185)*zin( 150))
          eri_value(  103)=eri_value(  103)+d23bra( 12)*d11ket(  4)*(xin(   2)*yin(  43)*zin(   5)+xin(  50)*yin(  91)*zin(  53)+xin(  98)*yin( 139)*zin( 101)+xin( 146)*yin( 187)*zin( 149))
          eri_value(  104)=eri_value(  104)+d23bra( 12)*d11ket(  5)*(xin(   1)*yin(  44)*zin(   5)+xin(  49)*yin(  92)*zin(  53)+xin(  97)*yin( 140)*zin( 101)+xin( 145)*yin( 188)*zin( 149))
          eri_value(  105)=eri_value(  105)+d23bra( 12)*d11ket(  6)*(xin(   1)*yin(  43)*zin(   6)+xin(  49)*yin(  91)*zin(  54)+xin(  97)*yin( 139)*zin( 102)+xin( 145)*yin( 187)*zin( 150))
          eri_value(  106)=eri_value(  106)+d23bra( 12)*d11ket(  7)*(xin(   2)*yin(  41)*zin(   7)+xin(  50)*yin(  89)*zin(  55)+xin(  98)*yin( 137)*zin( 103)+xin( 146)*yin( 185)*zin( 151))
          eri_value(  107)=eri_value(  107)+d23bra( 12)*d11ket(  8)*(xin(   1)*yin(  42)*zin(   7)+xin(  49)*yin(  90)*zin(  55)+xin(  97)*yin( 138)*zin( 103)+xin( 145)*yin( 186)*zin( 151))
          eri_value(  108)=eri_value(  108)+d23bra( 12)*d11ket(  9)*(xin(   1)*yin(  41)*zin(   8)+xin(  49)*yin(  89)*zin(  56)+xin(  97)*yin( 137)*zin( 104)+xin( 145)*yin( 185)*zin( 152))
          eri_value(  109)=eri_value(  109)+d23bra( 13)*d11ket(  1)*(xin(  12)*yin(   1)*zin(  37)+xin(  60)*yin(  49)*zin(  85)+xin( 108)*yin(  97)*zin( 133)+xin( 156)*yin( 145)*zin( 181))
          eri_value(  110)=eri_value(  110)+d23bra( 13)*d11ket(  2)*(xin(  11)*yin(   2)*zin(  37)+xin(  59)*yin(  50)*zin(  85)+xin( 107)*yin(  98)*zin( 133)+xin( 155)*yin( 146)*zin( 181))
          eri_value(  111)=eri_value(  111)+d23bra( 13)*d11ket(  3)*(xin(  11)*yin(   1)*zin(  38)+xin(  59)*yin(  49)*zin(  86)+xin( 107)*yin(  97)*zin( 134)+xin( 155)*yin( 145)*zin( 182))
          eri_value(  112)=eri_value(  112)+d23bra( 13)*d11ket(  4)*(xin(  10)*yin(   3)*zin(  37)+xin(  58)*yin(  51)*zin(  85)+xin( 106)*yin(  99)*zin( 133)+xin( 154)*yin( 147)*zin( 181))
          eri_value(  113)=eri_value(  113)+d23bra( 13)*d11ket(  5)*(xin(   9)*yin(   4)*zin(  37)+xin(  57)*yin(  52)*zin(  85)+xin( 105)*yin( 100)*zin( 133)+xin( 153)*yin( 148)*zin( 181))
          eri_value(  114)=eri_value(  114)+d23bra( 13)*d11ket(  6)*(xin(   9)*yin(   3)*zin(  38)+xin(  57)*yin(  51)*zin(  86)+xin( 105)*yin(  99)*zin( 134)+xin( 153)*yin( 147)*zin( 182))
          eri_value(  115)=eri_value(  115)+d23bra( 13)*d11ket(  7)*(xin(  10)*yin(   1)*zin(  39)+xin(  58)*yin(  49)*zin(  87)+xin( 106)*yin(  97)*zin( 135)+xin( 154)*yin( 145)*zin( 183))
          eri_value(  116)=eri_value(  116)+d23bra( 13)*d11ket(  8)*(xin(   9)*yin(   2)*zin(  39)+xin(  57)*yin(  50)*zin(  87)+xin( 105)*yin(  98)*zin( 135)+xin( 153)*yin( 146)*zin( 183))
          eri_value(  117)=eri_value(  117)+d23bra( 13)*d11ket(  9)*(xin(   9)*yin(   1)*zin(  40)+xin(  57)*yin(  49)*zin(  88)+xin( 105)*yin(  97)*zin( 136)+xin( 153)*yin( 145)*zin( 184))
          eri_value(  118)=eri_value(  118)+d23bra( 14)*d11ket(  1)*(xin(   4)*yin(   9)*zin(  37)+xin(  52)*yin(  57)*zin(  85)+xin( 100)*yin( 105)*zin( 133)+xin( 148)*yin( 153)*zin( 181))
          eri_value(  119)=eri_value(  119)+d23bra( 14)*d11ket(  2)*(xin(   3)*yin(  10)*zin(  37)+xin(  51)*yin(  58)*zin(  85)+xin(  99)*yin( 106)*zin( 133)+xin( 147)*yin( 154)*zin( 181))
          eri_value(  120)=eri_value(  120)+d23bra( 14)*d11ket(  3)*(xin(   3)*yin(   9)*zin(  38)+xin(  51)*yin(  57)*zin(  86)+xin(  99)*yin( 105)*zin( 134)+xin( 147)*yin( 153)*zin( 182))
          eri_value(  121)=eri_value(  121)+d23bra( 14)*d11ket(  4)*(xin(   2)*yin(  11)*zin(  37)+xin(  50)*yin(  59)*zin(  85)+xin(  98)*yin( 107)*zin( 133)+xin( 146)*yin( 155)*zin( 181))
          eri_value(  122)=eri_value(  122)+d23bra( 14)*d11ket(  5)*(xin(   1)*yin(  12)*zin(  37)+xin(  49)*yin(  60)*zin(  85)+xin(  97)*yin( 108)*zin( 133)+xin( 145)*yin( 156)*zin( 181))
          eri_value(  123)=eri_value(  123)+d23bra( 14)*d11ket(  6)*(xin(   1)*yin(  11)*zin(  38)+xin(  49)*yin(  59)*zin(  86)+xin(  97)*yin( 107)*zin( 134)+xin( 145)*yin( 155)*zin( 182))
          eri_value(  124)=eri_value(  124)+d23bra( 14)*d11ket(  7)*(xin(   2)*yin(   9)*zin(  39)+xin(  50)*yin(  57)*zin(  87)+xin(  98)*yin( 105)*zin( 135)+xin( 146)*yin( 153)*zin( 183))
          eri_value(  125)=eri_value(  125)+d23bra( 14)*d11ket(  8)*(xin(   1)*yin(  10)*zin(  39)+xin(  49)*yin(  58)*zin(  87)+xin(  97)*yin( 106)*zin( 135)+xin( 145)*yin( 154)*zin( 183))
          eri_value(  126)=eri_value(  126)+d23bra( 14)*d11ket(  9)*(xin(   1)*yin(   9)*zin(  40)+xin(  49)*yin(  57)*zin(  88)+xin(  97)*yin( 105)*zin( 136)+xin( 145)*yin( 153)*zin( 184))
          eri_value(  127)=eri_value(  127)+d23bra( 15)*d11ket(  1)*(xin(   4)*yin(   1)*zin(  45)+xin(  52)*yin(  49)*zin(  93)+xin( 100)*yin(  97)*zin( 141)+xin( 148)*yin( 145)*zin( 189))
          eri_value(  128)=eri_value(  128)+d23bra( 15)*d11ket(  2)*(xin(   3)*yin(   2)*zin(  45)+xin(  51)*yin(  50)*zin(  93)+xin(  99)*yin(  98)*zin( 141)+xin( 147)*yin( 146)*zin( 189))
          eri_value(  129)=eri_value(  129)+d23bra( 15)*d11ket(  3)*(xin(   3)*yin(   1)*zin(  46)+xin(  51)*yin(  49)*zin(  94)+xin(  99)*yin(  97)*zin( 142)+xin( 147)*yin( 145)*zin( 190))
          eri_value(  130)=eri_value(  130)+d23bra( 15)*d11ket(  4)*(xin(   2)*yin(   3)*zin(  45)+xin(  50)*yin(  51)*zin(  93)+xin(  98)*yin(  99)*zin( 141)+xin( 146)*yin( 147)*zin( 189))
          eri_value(  131)=eri_value(  131)+d23bra( 15)*d11ket(  5)*(xin(   1)*yin(   4)*zin(  45)+xin(  49)*yin(  52)*zin(  93)+xin(  97)*yin( 100)*zin( 141)+xin( 145)*yin( 148)*zin( 189))
          eri_value(  132)=eri_value(  132)+d23bra( 15)*d11ket(  6)*(xin(   1)*yin(   3)*zin(  46)+xin(  49)*yin(  51)*zin(  94)+xin(  97)*yin(  99)*zin( 142)+xin( 145)*yin( 147)*zin( 190))
          eri_value(  133)=eri_value(  133)+d23bra( 15)*d11ket(  7)*(xin(   2)*yin(   1)*zin(  47)+xin(  50)*yin(  49)*zin(  95)+xin(  98)*yin(  97)*zin( 143)+xin( 146)*yin( 145)*zin( 191))
          eri_value(  134)=eri_value(  134)+d23bra( 15)*d11ket(  8)*(xin(   1)*yin(   2)*zin(  47)+xin(  49)*yin(  50)*zin(  95)+xin(  97)*yin(  98)*zin( 143)+xin( 145)*yin( 146)*zin( 191))
          eri_value(  135)=eri_value(  135)+d23bra( 15)*d11ket(  9)*(xin(   1)*yin(   1)*zin(  48)+xin(  49)*yin(  49)*zin(  96)+xin(  97)*yin(  97)*zin( 144)+xin( 145)*yin( 145)*zin( 192))
          eri_value(  136)=eri_value(  136)+d23bra( 16)*d11ket(  1)*(xin(   8)*yin(   5)*zin(  37)+xin(  56)*yin(  53)*zin(  85)+xin( 104)*yin( 101)*zin( 133)+xin( 152)*yin( 149)*zin( 181))
          eri_value(  137)=eri_value(  137)+d23bra( 16)*d11ket(  2)*(xin(   7)*yin(   6)*zin(  37)+xin(  55)*yin(  54)*zin(  85)+xin( 103)*yin( 102)*zin( 133)+xin( 151)*yin( 150)*zin( 181))
          eri_value(  138)=eri_value(  138)+d23bra( 16)*d11ket(  3)*(xin(   7)*yin(   5)*zin(  38)+xin(  55)*yin(  53)*zin(  86)+xin( 103)*yin( 101)*zin( 134)+xin( 151)*yin( 149)*zin( 182))
          eri_value(  139)=eri_value(  139)+d23bra( 16)*d11ket(  4)*(xin(   6)*yin(   7)*zin(  37)+xin(  54)*yin(  55)*zin(  85)+xin( 102)*yin( 103)*zin( 133)+xin( 150)*yin( 151)*zin( 181))
          eri_value(  140)=eri_value(  140)+d23bra( 16)*d11ket(  5)*(xin(   5)*yin(   8)*zin(  37)+xin(  53)*yin(  56)*zin(  85)+xin( 101)*yin( 104)*zin( 133)+xin( 149)*yin( 152)*zin( 181))
          eri_value(  141)=eri_value(  141)+d23bra( 16)*d11ket(  6)*(xin(   5)*yin(   7)*zin(  38)+xin(  53)*yin(  55)*zin(  86)+xin( 101)*yin( 103)*zin( 134)+xin( 149)*yin( 151)*zin( 182))
          eri_value(  142)=eri_value(  142)+d23bra( 16)*d11ket(  7)*(xin(   6)*yin(   5)*zin(  39)+xin(  54)*yin(  53)*zin(  87)+xin( 102)*yin( 101)*zin( 135)+xin( 150)*yin( 149)*zin( 183))
          eri_value(  143)=eri_value(  143)+d23bra( 16)*d11ket(  8)*(xin(   5)*yin(   6)*zin(  39)+xin(  53)*yin(  54)*zin(  87)+xin( 101)*yin( 102)*zin( 135)+xin( 149)*yin( 150)*zin( 183))
          eri_value(  144)=eri_value(  144)+d23bra( 16)*d11ket(  9)*(xin(   5)*yin(   5)*zin(  40)+xin(  53)*yin(  53)*zin(  88)+xin( 101)*yin( 101)*zin( 136)+xin( 149)*yin( 149)*zin( 184))
          eri_value(  145)=eri_value(  145)+d23bra( 17)*d11ket(  1)*(xin(   8)*yin(   1)*zin(  41)+xin(  56)*yin(  49)*zin(  89)+xin( 104)*yin(  97)*zin( 137)+xin( 152)*yin( 145)*zin( 185))
          eri_value(  146)=eri_value(  146)+d23bra( 17)*d11ket(  2)*(xin(   7)*yin(   2)*zin(  41)+xin(  55)*yin(  50)*zin(  89)+xin( 103)*yin(  98)*zin( 137)+xin( 151)*yin( 146)*zin( 185))
          eri_value(  147)=eri_value(  147)+d23bra( 17)*d11ket(  3)*(xin(   7)*yin(   1)*zin(  42)+xin(  55)*yin(  49)*zin(  90)+xin( 103)*yin(  97)*zin( 138)+xin( 151)*yin( 145)*zin( 186))
          eri_value(  148)=eri_value(  148)+d23bra( 17)*d11ket(  4)*(xin(   6)*yin(   3)*zin(  41)+xin(  54)*yin(  51)*zin(  89)+xin( 102)*yin(  99)*zin( 137)+xin( 150)*yin( 147)*zin( 185))
          eri_value(  149)=eri_value(  149)+d23bra( 17)*d11ket(  5)*(xin(   5)*yin(   4)*zin(  41)+xin(  53)*yin(  52)*zin(  89)+xin( 101)*yin( 100)*zin( 137)+xin( 149)*yin( 148)*zin( 185))
          eri_value(  150)=eri_value(  150)+d23bra( 17)*d11ket(  6)*(xin(   5)*yin(   3)*zin(  42)+xin(  53)*yin(  51)*zin(  90)+xin( 101)*yin(  99)*zin( 138)+xin( 149)*yin( 147)*zin( 186))
          eri_value(  151)=eri_value(  151)+d23bra( 17)*d11ket(  7)*(xin(   6)*yin(   1)*zin(  43)+xin(  54)*yin(  49)*zin(  91)+xin( 102)*yin(  97)*zin( 139)+xin( 150)*yin( 145)*zin( 187))
          eri_value(  152)=eri_value(  152)+d23bra( 17)*d11ket(  8)*(xin(   5)*yin(   2)*zin(  43)+xin(  53)*yin(  50)*zin(  91)+xin( 101)*yin(  98)*zin( 139)+xin( 149)*yin( 146)*zin( 187))
          eri_value(  153)=eri_value(  153)+d23bra( 17)*d11ket(  9)*(xin(   5)*yin(   1)*zin(  44)+xin(  53)*yin(  49)*zin(  92)+xin( 101)*yin(  97)*zin( 140)+xin( 149)*yin( 145)*zin( 188))
          eri_value(  154)=eri_value(  154)+d23bra( 18)*d11ket(  1)*(xin(   4)*yin(   5)*zin(  41)+xin(  52)*yin(  53)*zin(  89)+xin( 100)*yin( 101)*zin( 137)+xin( 148)*yin( 149)*zin( 185))
          eri_value(  155)=eri_value(  155)+d23bra( 18)*d11ket(  2)*(xin(   3)*yin(   6)*zin(  41)+xin(  51)*yin(  54)*zin(  89)+xin(  99)*yin( 102)*zin( 137)+xin( 147)*yin( 150)*zin( 185))
          eri_value(  156)=eri_value(  156)+d23bra( 18)*d11ket(  3)*(xin(   3)*yin(   5)*zin(  42)+xin(  51)*yin(  53)*zin(  90)+xin(  99)*yin( 101)*zin( 138)+xin( 147)*yin( 149)*zin( 186))
          eri_value(  157)=eri_value(  157)+d23bra( 18)*d11ket(  4)*(xin(   2)*yin(   7)*zin(  41)+xin(  50)*yin(  55)*zin(  89)+xin(  98)*yin( 103)*zin( 137)+xin( 146)*yin( 151)*zin( 185))
          eri_value(  158)=eri_value(  158)+d23bra( 18)*d11ket(  5)*(xin(   1)*yin(   8)*zin(  41)+xin(  49)*yin(  56)*zin(  89)+xin(  97)*yin( 104)*zin( 137)+xin( 145)*yin( 152)*zin( 185))
          eri_value(  159)=eri_value(  159)+d23bra( 18)*d11ket(  6)*(xin(   1)*yin(   7)*zin(  42)+xin(  49)*yin(  55)*zin(  90)+xin(  97)*yin( 103)*zin( 138)+xin( 145)*yin( 151)*zin( 186))
          eri_value(  160)=eri_value(  160)+d23bra( 18)*d11ket(  7)*(xin(   2)*yin(   5)*zin(  43)+xin(  50)*yin(  53)*zin(  91)+xin(  98)*yin( 101)*zin( 139)+xin( 146)*yin( 149)*zin( 187))
          eri_value(  161)=eri_value(  161)+d23bra( 18)*d11ket(  8)*(xin(   1)*yin(   6)*zin(  43)+xin(  49)*yin(  54)*zin(  91)+xin(  97)*yin( 102)*zin( 139)+xin( 145)*yin( 150)*zin( 187))
          eri_value(  162)=eri_value(  162)+d23bra( 18)*d11ket(  9)*(xin(   1)*yin(   5)*zin(  44)+xin(  49)*yin(  53)*zin(  92)+xin(  97)*yin( 101)*zin( 140)+xin( 145)*yin( 149)*zin( 188))
          eri_value(  163)=eri_value(  163)+d23bra( 19)*d11ket(  1)*(xin(  36)*yin(  13)*zin(   1)+xin(  84)*yin(  61)*zin(  49)+xin( 132)*yin( 109)*zin(  97)+xin( 180)*yin( 157)*zin( 145))
          eri_value(  164)=eri_value(  164)+d23bra( 19)*d11ket(  2)*(xin(  35)*yin(  14)*zin(   1)+xin(  83)*yin(  62)*zin(  49)+xin( 131)*yin( 110)*zin(  97)+xin( 179)*yin( 158)*zin( 145))
          eri_value(  165)=eri_value(  165)+d23bra( 19)*d11ket(  3)*(xin(  35)*yin(  13)*zin(   2)+xin(  83)*yin(  61)*zin(  50)+xin( 131)*yin( 109)*zin(  98)+xin( 179)*yin( 157)*zin( 146))
          eri_value(  166)=eri_value(  166)+d23bra( 19)*d11ket(  4)*(xin(  34)*yin(  15)*zin(   1)+xin(  82)*yin(  63)*zin(  49)+xin( 130)*yin( 111)*zin(  97)+xin( 178)*yin( 159)*zin( 145))
          eri_value(  167)=eri_value(  167)+d23bra( 19)*d11ket(  5)*(xin(  33)*yin(  16)*zin(   1)+xin(  81)*yin(  64)*zin(  49)+xin( 129)*yin( 112)*zin(  97)+xin( 177)*yin( 160)*zin( 145))
          eri_value(  168)=eri_value(  168)+d23bra( 19)*d11ket(  6)*(xin(  33)*yin(  15)*zin(   2)+xin(  81)*yin(  63)*zin(  50)+xin( 129)*yin( 111)*zin(  98)+xin( 177)*yin( 159)*zin( 146))
          eri_value(  169)=eri_value(  169)+d23bra( 19)*d11ket(  7)*(xin(  34)*yin(  13)*zin(   3)+xin(  82)*yin(  61)*zin(  51)+xin( 130)*yin( 109)*zin(  99)+xin( 178)*yin( 157)*zin( 147))
          eri_value(  170)=eri_value(  170)+d23bra( 19)*d11ket(  8)*(xin(  33)*yin(  14)*zin(   3)+xin(  81)*yin(  62)*zin(  51)+xin( 129)*yin( 110)*zin(  99)+xin( 177)*yin( 158)*zin( 147))
          eri_value(  171)=eri_value(  171)+d23bra( 19)*d11ket(  9)*(xin(  33)*yin(  13)*zin(   4)+xin(  81)*yin(  61)*zin(  52)+xin( 129)*yin( 109)*zin( 100)+xin( 177)*yin( 157)*zin( 148))
          eri_value(  172)=eri_value(  172)+d23bra( 20)*d11ket(  1)*(xin(  28)*yin(  21)*zin(   1)+xin(  76)*yin(  69)*zin(  49)+xin( 124)*yin( 117)*zin(  97)+xin( 172)*yin( 165)*zin( 145))
          eri_value(  173)=eri_value(  173)+d23bra( 20)*d11ket(  2)*(xin(  27)*yin(  22)*zin(   1)+xin(  75)*yin(  70)*zin(  49)+xin( 123)*yin( 118)*zin(  97)+xin( 171)*yin( 166)*zin( 145))
          eri_value(  174)=eri_value(  174)+d23bra( 20)*d11ket(  3)*(xin(  27)*yin(  21)*zin(   2)+xin(  75)*yin(  69)*zin(  50)+xin( 123)*yin( 117)*zin(  98)+xin( 171)*yin( 165)*zin( 146))
          eri_value(  175)=eri_value(  175)+d23bra( 20)*d11ket(  4)*(xin(  26)*yin(  23)*zin(   1)+xin(  74)*yin(  71)*zin(  49)+xin( 122)*yin( 119)*zin(  97)+xin( 170)*yin( 167)*zin( 145))
          eri_value(  176)=eri_value(  176)+d23bra( 20)*d11ket(  5)*(xin(  25)*yin(  24)*zin(   1)+xin(  73)*yin(  72)*zin(  49)+xin( 121)*yin( 120)*zin(  97)+xin( 169)*yin( 168)*zin( 145))
          eri_value(  177)=eri_value(  177)+d23bra( 20)*d11ket(  6)*(xin(  25)*yin(  23)*zin(   2)+xin(  73)*yin(  71)*zin(  50)+xin( 121)*yin( 119)*zin(  98)+xin( 169)*yin( 167)*zin( 146))
          eri_value(  178)=eri_value(  178)+d23bra( 20)*d11ket(  7)*(xin(  26)*yin(  21)*zin(   3)+xin(  74)*yin(  69)*zin(  51)+xin( 122)*yin( 117)*zin(  99)+xin( 170)*yin( 165)*zin( 147))
          eri_value(  179)=eri_value(  179)+d23bra( 20)*d11ket(  8)*(xin(  25)*yin(  22)*zin(   3)+xin(  73)*yin(  70)*zin(  51)+xin( 121)*yin( 118)*zin(  99)+xin( 169)*yin( 166)*zin( 147))
          eri_value(  180)=eri_value(  180)+d23bra( 20)*d11ket(  9)*(xin(  25)*yin(  21)*zin(   4)+xin(  73)*yin(  69)*zin(  52)+xin( 121)*yin( 117)*zin( 100)+xin( 169)*yin( 165)*zin( 148))
          eri_value(  181)=eri_value(  181)+d23bra( 21)*d11ket(  1)*(xin(  28)*yin(  13)*zin(   9)+xin(  76)*yin(  61)*zin(  57)+xin( 124)*yin( 109)*zin( 105)+xin( 172)*yin( 157)*zin( 153))
          eri_value(  182)=eri_value(  182)+d23bra( 21)*d11ket(  2)*(xin(  27)*yin(  14)*zin(   9)+xin(  75)*yin(  62)*zin(  57)+xin( 123)*yin( 110)*zin( 105)+xin( 171)*yin( 158)*zin( 153))
          eri_value(  183)=eri_value(  183)+d23bra( 21)*d11ket(  3)*(xin(  27)*yin(  13)*zin(  10)+xin(  75)*yin(  61)*zin(  58)+xin( 123)*yin( 109)*zin( 106)+xin( 171)*yin( 157)*zin( 154))
          eri_value(  184)=eri_value(  184)+d23bra( 21)*d11ket(  4)*(xin(  26)*yin(  15)*zin(   9)+xin(  74)*yin(  63)*zin(  57)+xin( 122)*yin( 111)*zin( 105)+xin( 170)*yin( 159)*zin( 153))
          eri_value(  185)=eri_value(  185)+d23bra( 21)*d11ket(  5)*(xin(  25)*yin(  16)*zin(   9)+xin(  73)*yin(  64)*zin(  57)+xin( 121)*yin( 112)*zin( 105)+xin( 169)*yin( 160)*zin( 153))
          eri_value(  186)=eri_value(  186)+d23bra( 21)*d11ket(  6)*(xin(  25)*yin(  15)*zin(  10)+xin(  73)*yin(  63)*zin(  58)+xin( 121)*yin( 111)*zin( 106)+xin( 169)*yin( 159)*zin( 154))
          eri_value(  187)=eri_value(  187)+d23bra( 21)*d11ket(  7)*(xin(  26)*yin(  13)*zin(  11)+xin(  74)*yin(  61)*zin(  59)+xin( 122)*yin( 109)*zin( 107)+xin( 170)*yin( 157)*zin( 155))
          eri_value(  188)=eri_value(  188)+d23bra( 21)*d11ket(  8)*(xin(  25)*yin(  14)*zin(  11)+xin(  73)*yin(  62)*zin(  59)+xin( 121)*yin( 110)*zin( 107)+xin( 169)*yin( 158)*zin( 155))
          eri_value(  189)=eri_value(  189)+d23bra( 21)*d11ket(  9)*(xin(  25)*yin(  13)*zin(  12)+xin(  73)*yin(  61)*zin(  60)+xin( 121)*yin( 109)*zin( 108)+xin( 169)*yin( 157)*zin( 156))
          eri_value(  190)=eri_value(  190)+d23bra( 22)*d11ket(  1)*(xin(  32)*yin(  17)*zin(   1)+xin(  80)*yin(  65)*zin(  49)+xin( 128)*yin( 113)*zin(  97)+xin( 176)*yin( 161)*zin( 145))
          eri_value(  191)=eri_value(  191)+d23bra( 22)*d11ket(  2)*(xin(  31)*yin(  18)*zin(   1)+xin(  79)*yin(  66)*zin(  49)+xin( 127)*yin( 114)*zin(  97)+xin( 175)*yin( 162)*zin( 145))
          eri_value(  192)=eri_value(  192)+d23bra( 22)*d11ket(  3)*(xin(  31)*yin(  17)*zin(   2)+xin(  79)*yin(  65)*zin(  50)+xin( 127)*yin( 113)*zin(  98)+xin( 175)*yin( 161)*zin( 146))
          eri_value(  193)=eri_value(  193)+d23bra( 22)*d11ket(  4)*(xin(  30)*yin(  19)*zin(   1)+xin(  78)*yin(  67)*zin(  49)+xin( 126)*yin( 115)*zin(  97)+xin( 174)*yin( 163)*zin( 145))
          eri_value(  194)=eri_value(  194)+d23bra( 22)*d11ket(  5)*(xin(  29)*yin(  20)*zin(   1)+xin(  77)*yin(  68)*zin(  49)+xin( 125)*yin( 116)*zin(  97)+xin( 173)*yin( 164)*zin( 145))
          eri_value(  195)=eri_value(  195)+d23bra( 22)*d11ket(  6)*(xin(  29)*yin(  19)*zin(   2)+xin(  77)*yin(  67)*zin(  50)+xin( 125)*yin( 115)*zin(  98)+xin( 173)*yin( 163)*zin( 146))
          eri_value(  196)=eri_value(  196)+d23bra( 22)*d11ket(  7)*(xin(  30)*yin(  17)*zin(   3)+xin(  78)*yin(  65)*zin(  51)+xin( 126)*yin( 113)*zin(  99)+xin( 174)*yin( 161)*zin( 147))
          eri_value(  197)=eri_value(  197)+d23bra( 22)*d11ket(  8)*(xin(  29)*yin(  18)*zin(   3)+xin(  77)*yin(  66)*zin(  51)+xin( 125)*yin( 114)*zin(  99)+xin( 173)*yin( 162)*zin( 147))
          eri_value(  198)=eri_value(  198)+d23bra( 22)*d11ket(  9)*(xin(  29)*yin(  17)*zin(   4)+xin(  77)*yin(  65)*zin(  52)+xin( 125)*yin( 113)*zin( 100)+xin( 173)*yin( 161)*zin( 148))
          eri_value(  199)=eri_value(  199)+d23bra( 23)*d11ket(  1)*(xin(  32)*yin(  13)*zin(   5)+xin(  80)*yin(  61)*zin(  53)+xin( 128)*yin( 109)*zin( 101)+xin( 176)*yin( 157)*zin( 149))
          eri_value(  200)=eri_value(  200)+d23bra( 23)*d11ket(  2)*(xin(  31)*yin(  14)*zin(   5)+xin(  79)*yin(  62)*zin(  53)+xin( 127)*yin( 110)*zin( 101)+xin( 175)*yin( 158)*zin( 149))
          eri_value(  201)=eri_value(  201)+d23bra( 23)*d11ket(  3)*(xin(  31)*yin(  13)*zin(   6)+xin(  79)*yin(  61)*zin(  54)+xin( 127)*yin( 109)*zin( 102)+xin( 175)*yin( 157)*zin( 150))
          eri_value(  202)=eri_value(  202)+d23bra( 23)*d11ket(  4)*(xin(  30)*yin(  15)*zin(   5)+xin(  78)*yin(  63)*zin(  53)+xin( 126)*yin( 111)*zin( 101)+xin( 174)*yin( 159)*zin( 149))
          eri_value(  203)=eri_value(  203)+d23bra( 23)*d11ket(  5)*(xin(  29)*yin(  16)*zin(   5)+xin(  77)*yin(  64)*zin(  53)+xin( 125)*yin( 112)*zin( 101)+xin( 173)*yin( 160)*zin( 149))
          eri_value(  204)=eri_value(  204)+d23bra( 23)*d11ket(  6)*(xin(  29)*yin(  15)*zin(   6)+xin(  77)*yin(  63)*zin(  54)+xin( 125)*yin( 111)*zin( 102)+xin( 173)*yin( 159)*zin( 150))
          eri_value(  205)=eri_value(  205)+d23bra( 23)*d11ket(  7)*(xin(  30)*yin(  13)*zin(   7)+xin(  78)*yin(  61)*zin(  55)+xin( 126)*yin( 109)*zin( 103)+xin( 174)*yin( 157)*zin( 151))
          eri_value(  206)=eri_value(  206)+d23bra( 23)*d11ket(  8)*(xin(  29)*yin(  14)*zin(   7)+xin(  77)*yin(  62)*zin(  55)+xin( 125)*yin( 110)*zin( 103)+xin( 173)*yin( 158)*zin( 151))
          eri_value(  207)=eri_value(  207)+d23bra( 23)*d11ket(  9)*(xin(  29)*yin(  13)*zin(   8)+xin(  77)*yin(  61)*zin(  56)+xin( 125)*yin( 109)*zin( 104)+xin( 173)*yin( 157)*zin( 152))
          eri_value(  208)=eri_value(  208)+d23bra( 24)*d11ket(  1)*(xin(  28)*yin(  17)*zin(   5)+xin(  76)*yin(  65)*zin(  53)+xin( 124)*yin( 113)*zin( 101)+xin( 172)*yin( 161)*zin( 149))
          eri_value(  209)=eri_value(  209)+d23bra( 24)*d11ket(  2)*(xin(  27)*yin(  18)*zin(   5)+xin(  75)*yin(  66)*zin(  53)+xin( 123)*yin( 114)*zin( 101)+xin( 171)*yin( 162)*zin( 149))
          eri_value(  210)=eri_value(  210)+d23bra( 24)*d11ket(  3)*(xin(  27)*yin(  17)*zin(   6)+xin(  75)*yin(  65)*zin(  54)+xin( 123)*yin( 113)*zin( 102)+xin( 171)*yin( 161)*zin( 150))
          eri_value(  211)=eri_value(  211)+d23bra( 24)*d11ket(  4)*(xin(  26)*yin(  19)*zin(   5)+xin(  74)*yin(  67)*zin(  53)+xin( 122)*yin( 115)*zin( 101)+xin( 170)*yin( 163)*zin( 149))
          eri_value(  212)=eri_value(  212)+d23bra( 24)*d11ket(  5)*(xin(  25)*yin(  20)*zin(   5)+xin(  73)*yin(  68)*zin(  53)+xin( 121)*yin( 116)*zin( 101)+xin( 169)*yin( 164)*zin( 149))
          eri_value(  213)=eri_value(  213)+d23bra( 24)*d11ket(  6)*(xin(  25)*yin(  19)*zin(   6)+xin(  73)*yin(  67)*zin(  54)+xin( 121)*yin( 115)*zin( 102)+xin( 169)*yin( 163)*zin( 150))
          eri_value(  214)=eri_value(  214)+d23bra( 24)*d11ket(  7)*(xin(  26)*yin(  17)*zin(   7)+xin(  74)*yin(  65)*zin(  55)+xin( 122)*yin( 113)*zin( 103)+xin( 170)*yin( 161)*zin( 151))
          eri_value(  215)=eri_value(  215)+d23bra( 24)*d11ket(  8)*(xin(  25)*yin(  18)*zin(   7)+xin(  73)*yin(  66)*zin(  55)+xin( 121)*yin( 114)*zin( 103)+xin( 169)*yin( 162)*zin( 151))
          eri_value(  216)=eri_value(  216)+d23bra( 24)*d11ket(  9)*(xin(  25)*yin(  17)*zin(   8)+xin(  73)*yin(  65)*zin(  56)+xin( 121)*yin( 113)*zin( 104)+xin( 169)*yin( 161)*zin( 152))
          eri_value(  217)=eri_value(  217)+d23bra( 25)*d11ket(  1)*(xin(  36)*yin(   1)*zin(  13)+xin(  84)*yin(  49)*zin(  61)+xin( 132)*yin(  97)*zin( 109)+xin( 180)*yin( 145)*zin( 157))
          eri_value(  218)=eri_value(  218)+d23bra( 25)*d11ket(  2)*(xin(  35)*yin(   2)*zin(  13)+xin(  83)*yin(  50)*zin(  61)+xin( 131)*yin(  98)*zin( 109)+xin( 179)*yin( 146)*zin( 157))
          eri_value(  219)=eri_value(  219)+d23bra( 25)*d11ket(  3)*(xin(  35)*yin(   1)*zin(  14)+xin(  83)*yin(  49)*zin(  62)+xin( 131)*yin(  97)*zin( 110)+xin( 179)*yin( 145)*zin( 158))
          eri_value(  220)=eri_value(  220)+d23bra( 25)*d11ket(  4)*(xin(  34)*yin(   3)*zin(  13)+xin(  82)*yin(  51)*zin(  61)+xin( 130)*yin(  99)*zin( 109)+xin( 178)*yin( 147)*zin( 157))
          eri_value(  221)=eri_value(  221)+d23bra( 25)*d11ket(  5)*(xin(  33)*yin(   4)*zin(  13)+xin(  81)*yin(  52)*zin(  61)+xin( 129)*yin( 100)*zin( 109)+xin( 177)*yin( 148)*zin( 157))
          eri_value(  222)=eri_value(  222)+d23bra( 25)*d11ket(  6)*(xin(  33)*yin(   3)*zin(  14)+xin(  81)*yin(  51)*zin(  62)+xin( 129)*yin(  99)*zin( 110)+xin( 177)*yin( 147)*zin( 158))
          eri_value(  223)=eri_value(  223)+d23bra( 25)*d11ket(  7)*(xin(  34)*yin(   1)*zin(  15)+xin(  82)*yin(  49)*zin(  63)+xin( 130)*yin(  97)*zin( 111)+xin( 178)*yin( 145)*zin( 159))
          eri_value(  224)=eri_value(  224)+d23bra( 25)*d11ket(  8)*(xin(  33)*yin(   2)*zin(  15)+xin(  81)*yin(  50)*zin(  63)+xin( 129)*yin(  98)*zin( 111)+xin( 177)*yin( 146)*zin( 159))
          eri_value(  225)=eri_value(  225)+d23bra( 25)*d11ket(  9)*(xin(  33)*yin(   1)*zin(  16)+xin(  81)*yin(  49)*zin(  64)+xin( 129)*yin(  97)*zin( 112)+xin( 177)*yin( 145)*zin( 160))
          eri_value(  226)=eri_value(  226)+d23bra( 26)*d11ket(  1)*(xin(  28)*yin(   9)*zin(  13)+xin(  76)*yin(  57)*zin(  61)+xin( 124)*yin( 105)*zin( 109)+xin( 172)*yin( 153)*zin( 157))
          eri_value(  227)=eri_value(  227)+d23bra( 26)*d11ket(  2)*(xin(  27)*yin(  10)*zin(  13)+xin(  75)*yin(  58)*zin(  61)+xin( 123)*yin( 106)*zin( 109)+xin( 171)*yin( 154)*zin( 157))
          eri_value(  228)=eri_value(  228)+d23bra( 26)*d11ket(  3)*(xin(  27)*yin(   9)*zin(  14)+xin(  75)*yin(  57)*zin(  62)+xin( 123)*yin( 105)*zin( 110)+xin( 171)*yin( 153)*zin( 158))
          eri_value(  229)=eri_value(  229)+d23bra( 26)*d11ket(  4)*(xin(  26)*yin(  11)*zin(  13)+xin(  74)*yin(  59)*zin(  61)+xin( 122)*yin( 107)*zin( 109)+xin( 170)*yin( 155)*zin( 157))
          eri_value(  230)=eri_value(  230)+d23bra( 26)*d11ket(  5)*(xin(  25)*yin(  12)*zin(  13)+xin(  73)*yin(  60)*zin(  61)+xin( 121)*yin( 108)*zin( 109)+xin( 169)*yin( 156)*zin( 157))
          eri_value(  231)=eri_value(  231)+d23bra( 26)*d11ket(  6)*(xin(  25)*yin(  11)*zin(  14)+xin(  73)*yin(  59)*zin(  62)+xin( 121)*yin( 107)*zin( 110)+xin( 169)*yin( 155)*zin( 158))
          eri_value(  232)=eri_value(  232)+d23bra( 26)*d11ket(  7)*(xin(  26)*yin(   9)*zin(  15)+xin(  74)*yin(  57)*zin(  63)+xin( 122)*yin( 105)*zin( 111)+xin( 170)*yin( 153)*zin( 159))
          eri_value(  233)=eri_value(  233)+d23bra( 26)*d11ket(  8)*(xin(  25)*yin(  10)*zin(  15)+xin(  73)*yin(  58)*zin(  63)+xin( 121)*yin( 106)*zin( 111)+xin( 169)*yin( 154)*zin( 159))
          eri_value(  234)=eri_value(  234)+d23bra( 26)*d11ket(  9)*(xin(  25)*yin(   9)*zin(  16)+xin(  73)*yin(  57)*zin(  64)+xin( 121)*yin( 105)*zin( 112)+xin( 169)*yin( 153)*zin( 160))
          eri_value(  235)=eri_value(  235)+d23bra( 27)*d11ket(  1)*(xin(  28)*yin(   1)*zin(  21)+xin(  76)*yin(  49)*zin(  69)+xin( 124)*yin(  97)*zin( 117)+xin( 172)*yin( 145)*zin( 165))
          eri_value(  236)=eri_value(  236)+d23bra( 27)*d11ket(  2)*(xin(  27)*yin(   2)*zin(  21)+xin(  75)*yin(  50)*zin(  69)+xin( 123)*yin(  98)*zin( 117)+xin( 171)*yin( 146)*zin( 165))
          eri_value(  237)=eri_value(  237)+d23bra( 27)*d11ket(  3)*(xin(  27)*yin(   1)*zin(  22)+xin(  75)*yin(  49)*zin(  70)+xin( 123)*yin(  97)*zin( 118)+xin( 171)*yin( 145)*zin( 166))
          eri_value(  238)=eri_value(  238)+d23bra( 27)*d11ket(  4)*(xin(  26)*yin(   3)*zin(  21)+xin(  74)*yin(  51)*zin(  69)+xin( 122)*yin(  99)*zin( 117)+xin( 170)*yin( 147)*zin( 165))
          eri_value(  239)=eri_value(  239)+d23bra( 27)*d11ket(  5)*(xin(  25)*yin(   4)*zin(  21)+xin(  73)*yin(  52)*zin(  69)+xin( 121)*yin( 100)*zin( 117)+xin( 169)*yin( 148)*zin( 165))
          eri_value(  240)=eri_value(  240)+d23bra( 27)*d11ket(  6)*(xin(  25)*yin(   3)*zin(  22)+xin(  73)*yin(  51)*zin(  70)+xin( 121)*yin(  99)*zin( 118)+xin( 169)*yin( 147)*zin( 166))
          eri_value(  241)=eri_value(  241)+d23bra( 27)*d11ket(  7)*(xin(  26)*yin(   1)*zin(  23)+xin(  74)*yin(  49)*zin(  71)+xin( 122)*yin(  97)*zin( 119)+xin( 170)*yin( 145)*zin( 167))
          eri_value(  242)=eri_value(  242)+d23bra( 27)*d11ket(  8)*(xin(  25)*yin(   2)*zin(  23)+xin(  73)*yin(  50)*zin(  71)+xin( 121)*yin(  98)*zin( 119)+xin( 169)*yin( 146)*zin( 167))
          eri_value(  243)=eri_value(  243)+d23bra( 27)*d11ket(  9)*(xin(  25)*yin(   1)*zin(  24)+xin(  73)*yin(  49)*zin(  72)+xin( 121)*yin(  97)*zin( 120)+xin( 169)*yin( 145)*zin( 168))
          eri_value(  244)=eri_value(  244)+d23bra( 28)*d11ket(  1)*(xin(  32)*yin(   5)*zin(  13)+xin(  80)*yin(  53)*zin(  61)+xin( 128)*yin( 101)*zin( 109)+xin( 176)*yin( 149)*zin( 157))
          eri_value(  245)=eri_value(  245)+d23bra( 28)*d11ket(  2)*(xin(  31)*yin(   6)*zin(  13)+xin(  79)*yin(  54)*zin(  61)+xin( 127)*yin( 102)*zin( 109)+xin( 175)*yin( 150)*zin( 157))
          eri_value(  246)=eri_value(  246)+d23bra( 28)*d11ket(  3)*(xin(  31)*yin(   5)*zin(  14)+xin(  79)*yin(  53)*zin(  62)+xin( 127)*yin( 101)*zin( 110)+xin( 175)*yin( 149)*zin( 158))
          eri_value(  247)=eri_value(  247)+d23bra( 28)*d11ket(  4)*(xin(  30)*yin(   7)*zin(  13)+xin(  78)*yin(  55)*zin(  61)+xin( 126)*yin( 103)*zin( 109)+xin( 174)*yin( 151)*zin( 157))
          eri_value(  248)=eri_value(  248)+d23bra( 28)*d11ket(  5)*(xin(  29)*yin(   8)*zin(  13)+xin(  77)*yin(  56)*zin(  61)+xin( 125)*yin( 104)*zin( 109)+xin( 173)*yin( 152)*zin( 157))
          eri_value(  249)=eri_value(  249)+d23bra( 28)*d11ket(  6)*(xin(  29)*yin(   7)*zin(  14)+xin(  77)*yin(  55)*zin(  62)+xin( 125)*yin( 103)*zin( 110)+xin( 173)*yin( 151)*zin( 158))
          eri_value(  250)=eri_value(  250)+d23bra( 28)*d11ket(  7)*(xin(  30)*yin(   5)*zin(  15)+xin(  78)*yin(  53)*zin(  63)+xin( 126)*yin( 101)*zin( 111)+xin( 174)*yin( 149)*zin( 159))
          eri_value(  251)=eri_value(  251)+d23bra( 28)*d11ket(  8)*(xin(  29)*yin(   6)*zin(  15)+xin(  77)*yin(  54)*zin(  63)+xin( 125)*yin( 102)*zin( 111)+xin( 173)*yin( 150)*zin( 159))
          eri_value(  252)=eri_value(  252)+d23bra( 28)*d11ket(  9)*(xin(  29)*yin(   5)*zin(  16)+xin(  77)*yin(  53)*zin(  64)+xin( 125)*yin( 101)*zin( 112)+xin( 173)*yin( 149)*zin( 160))
          eri_value(  253)=eri_value(  253)+d23bra( 29)*d11ket(  1)*(xin(  32)*yin(   1)*zin(  17)+xin(  80)*yin(  49)*zin(  65)+xin( 128)*yin(  97)*zin( 113)+xin( 176)*yin( 145)*zin( 161))
          eri_value(  254)=eri_value(  254)+d23bra( 29)*d11ket(  2)*(xin(  31)*yin(   2)*zin(  17)+xin(  79)*yin(  50)*zin(  65)+xin( 127)*yin(  98)*zin( 113)+xin( 175)*yin( 146)*zin( 161))
          eri_value(  255)=eri_value(  255)+d23bra( 29)*d11ket(  3)*(xin(  31)*yin(   1)*zin(  18)+xin(  79)*yin(  49)*zin(  66)+xin( 127)*yin(  97)*zin( 114)+xin( 175)*yin( 145)*zin( 162))
          eri_value(  256)=eri_value(  256)+d23bra( 29)*d11ket(  4)*(xin(  30)*yin(   3)*zin(  17)+xin(  78)*yin(  51)*zin(  65)+xin( 126)*yin(  99)*zin( 113)+xin( 174)*yin( 147)*zin( 161))
          eri_value(  257)=eri_value(  257)+d23bra( 29)*d11ket(  5)*(xin(  29)*yin(   4)*zin(  17)+xin(  77)*yin(  52)*zin(  65)+xin( 125)*yin( 100)*zin( 113)+xin( 173)*yin( 148)*zin( 161))
          eri_value(  258)=eri_value(  258)+d23bra( 29)*d11ket(  6)*(xin(  29)*yin(   3)*zin(  18)+xin(  77)*yin(  51)*zin(  66)+xin( 125)*yin(  99)*zin( 114)+xin( 173)*yin( 147)*zin( 162))
          eri_value(  259)=eri_value(  259)+d23bra( 29)*d11ket(  7)*(xin(  30)*yin(   1)*zin(  19)+xin(  78)*yin(  49)*zin(  67)+xin( 126)*yin(  97)*zin( 115)+xin( 174)*yin( 145)*zin( 163))
          eri_value(  260)=eri_value(  260)+d23bra( 29)*d11ket(  8)*(xin(  29)*yin(   2)*zin(  19)+xin(  77)*yin(  50)*zin(  67)+xin( 125)*yin(  98)*zin( 115)+xin( 173)*yin( 146)*zin( 163))
          eri_value(  261)=eri_value(  261)+d23bra( 29)*d11ket(  9)*(xin(  29)*yin(   1)*zin(  20)+xin(  77)*yin(  49)*zin(  68)+xin( 125)*yin(  97)*zin( 116)+xin( 173)*yin( 145)*zin( 164))
          eri_value(  262)=eri_value(  262)+d23bra( 30)*d11ket(  1)*(xin(  28)*yin(   5)*zin(  17)+xin(  76)*yin(  53)*zin(  65)+xin( 124)*yin( 101)*zin( 113)+xin( 172)*yin( 149)*zin( 161))
          eri_value(  263)=eri_value(  263)+d23bra( 30)*d11ket(  2)*(xin(  27)*yin(   6)*zin(  17)+xin(  75)*yin(  54)*zin(  65)+xin( 123)*yin( 102)*zin( 113)+xin( 171)*yin( 150)*zin( 161))
          eri_value(  264)=eri_value(  264)+d23bra( 30)*d11ket(  3)*(xin(  27)*yin(   5)*zin(  18)+xin(  75)*yin(  53)*zin(  66)+xin( 123)*yin( 101)*zin( 114)+xin( 171)*yin( 149)*zin( 162))
          eri_value(  265)=eri_value(  265)+d23bra( 30)*d11ket(  4)*(xin(  26)*yin(   7)*zin(  17)+xin(  74)*yin(  55)*zin(  65)+xin( 122)*yin( 103)*zin( 113)+xin( 170)*yin( 151)*zin( 161))
          eri_value(  266)=eri_value(  266)+d23bra( 30)*d11ket(  5)*(xin(  25)*yin(   8)*zin(  17)+xin(  73)*yin(  56)*zin(  65)+xin( 121)*yin( 104)*zin( 113)+xin( 169)*yin( 152)*zin( 161))
          eri_value(  267)=eri_value(  267)+d23bra( 30)*d11ket(  6)*(xin(  25)*yin(   7)*zin(  18)+xin(  73)*yin(  55)*zin(  66)+xin( 121)*yin( 103)*zin( 114)+xin( 169)*yin( 151)*zin( 162))
          eri_value(  268)=eri_value(  268)+d23bra( 30)*d11ket(  7)*(xin(  26)*yin(   5)*zin(  19)+xin(  74)*yin(  53)*zin(  67)+xin( 122)*yin( 101)*zin( 115)+xin( 170)*yin( 149)*zin( 163))
          eri_value(  269)=eri_value(  269)+d23bra( 30)*d11ket(  8)*(xin(  25)*yin(   6)*zin(  19)+xin(  73)*yin(  54)*zin(  67)+xin( 121)*yin( 102)*zin( 115)+xin( 169)*yin( 150)*zin( 163))
          eri_value(  270)=eri_value(  270)+d23bra( 30)*d11ket(  9)*(xin(  25)*yin(   5)*zin(  20)+xin(  73)*yin(  53)*zin(  68)+xin( 121)*yin( 101)*zin( 116)+xin( 169)*yin( 149)*zin( 164))
          eri_value(  271)=eri_value(  271)+d23bra( 31)*d11ket(  1)*(xin(  24)*yin(  25)*zin(   1)+xin(  72)*yin(  73)*zin(  49)+xin( 120)*yin( 121)*zin(  97)+xin( 168)*yin( 169)*zin( 145))
          eri_value(  272)=eri_value(  272)+d23bra( 31)*d11ket(  2)*(xin(  23)*yin(  26)*zin(   1)+xin(  71)*yin(  74)*zin(  49)+xin( 119)*yin( 122)*zin(  97)+xin( 167)*yin( 170)*zin( 145))
          eri_value(  273)=eri_value(  273)+d23bra( 31)*d11ket(  3)*(xin(  23)*yin(  25)*zin(   2)+xin(  71)*yin(  73)*zin(  50)+xin( 119)*yin( 121)*zin(  98)+xin( 167)*yin( 169)*zin( 146))
          eri_value(  274)=eri_value(  274)+d23bra( 31)*d11ket(  4)*(xin(  22)*yin(  27)*zin(   1)+xin(  70)*yin(  75)*zin(  49)+xin( 118)*yin( 123)*zin(  97)+xin( 166)*yin( 171)*zin( 145))
          eri_value(  275)=eri_value(  275)+d23bra( 31)*d11ket(  5)*(xin(  21)*yin(  28)*zin(   1)+xin(  69)*yin(  76)*zin(  49)+xin( 117)*yin( 124)*zin(  97)+xin( 165)*yin( 172)*zin( 145))
          eri_value(  276)=eri_value(  276)+d23bra( 31)*d11ket(  6)*(xin(  21)*yin(  27)*zin(   2)+xin(  69)*yin(  75)*zin(  50)+xin( 117)*yin( 123)*zin(  98)+xin( 165)*yin( 171)*zin( 146))
          eri_value(  277)=eri_value(  277)+d23bra( 31)*d11ket(  7)*(xin(  22)*yin(  25)*zin(   3)+xin(  70)*yin(  73)*zin(  51)+xin( 118)*yin( 121)*zin(  99)+xin( 166)*yin( 169)*zin( 147))
          eri_value(  278)=eri_value(  278)+d23bra( 31)*d11ket(  8)*(xin(  21)*yin(  26)*zin(   3)+xin(  69)*yin(  74)*zin(  51)+xin( 117)*yin( 122)*zin(  99)+xin( 165)*yin( 170)*zin( 147))
          eri_value(  279)=eri_value(  279)+d23bra( 31)*d11ket(  9)*(xin(  21)*yin(  25)*zin(   4)+xin(  69)*yin(  73)*zin(  52)+xin( 117)*yin( 121)*zin( 100)+xin( 165)*yin( 169)*zin( 148))
          eri_value(  280)=eri_value(  280)+d23bra( 32)*d11ket(  1)*(xin(  16)*yin(  33)*zin(   1)+xin(  64)*yin(  81)*zin(  49)+xin( 112)*yin( 129)*zin(  97)+xin( 160)*yin( 177)*zin( 145))
          eri_value(  281)=eri_value(  281)+d23bra( 32)*d11ket(  2)*(xin(  15)*yin(  34)*zin(   1)+xin(  63)*yin(  82)*zin(  49)+xin( 111)*yin( 130)*zin(  97)+xin( 159)*yin( 178)*zin( 145))
          eri_value(  282)=eri_value(  282)+d23bra( 32)*d11ket(  3)*(xin(  15)*yin(  33)*zin(   2)+xin(  63)*yin(  81)*zin(  50)+xin( 111)*yin( 129)*zin(  98)+xin( 159)*yin( 177)*zin( 146))
          eri_value(  283)=eri_value(  283)+d23bra( 32)*d11ket(  4)*(xin(  14)*yin(  35)*zin(   1)+xin(  62)*yin(  83)*zin(  49)+xin( 110)*yin( 131)*zin(  97)+xin( 158)*yin( 179)*zin( 145))
          eri_value(  284)=eri_value(  284)+d23bra( 32)*d11ket(  5)*(xin(  13)*yin(  36)*zin(   1)+xin(  61)*yin(  84)*zin(  49)+xin( 109)*yin( 132)*zin(  97)+xin( 157)*yin( 180)*zin( 145))
          eri_value(  285)=eri_value(  285)+d23bra( 32)*d11ket(  6)*(xin(  13)*yin(  35)*zin(   2)+xin(  61)*yin(  83)*zin(  50)+xin( 109)*yin( 131)*zin(  98)+xin( 157)*yin( 179)*zin( 146))
          eri_value(  286)=eri_value(  286)+d23bra( 32)*d11ket(  7)*(xin(  14)*yin(  33)*zin(   3)+xin(  62)*yin(  81)*zin(  51)+xin( 110)*yin( 129)*zin(  99)+xin( 158)*yin( 177)*zin( 147))
          eri_value(  287)=eri_value(  287)+d23bra( 32)*d11ket(  8)*(xin(  13)*yin(  34)*zin(   3)+xin(  61)*yin(  82)*zin(  51)+xin( 109)*yin( 130)*zin(  99)+xin( 157)*yin( 178)*zin( 147))
          eri_value(  288)=eri_value(  288)+d23bra( 32)*d11ket(  9)*(xin(  13)*yin(  33)*zin(   4)+xin(  61)*yin(  81)*zin(  52)+xin( 109)*yin( 129)*zin( 100)+xin( 157)*yin( 177)*zin( 148))
          eri_value(  289)=eri_value(  289)+d23bra( 33)*d11ket(  1)*(xin(  16)*yin(  25)*zin(   9)+xin(  64)*yin(  73)*zin(  57)+xin( 112)*yin( 121)*zin( 105)+xin( 160)*yin( 169)*zin( 153))
          eri_value(  290)=eri_value(  290)+d23bra( 33)*d11ket(  2)*(xin(  15)*yin(  26)*zin(   9)+xin(  63)*yin(  74)*zin(  57)+xin( 111)*yin( 122)*zin( 105)+xin( 159)*yin( 170)*zin( 153))
          eri_value(  291)=eri_value(  291)+d23bra( 33)*d11ket(  3)*(xin(  15)*yin(  25)*zin(  10)+xin(  63)*yin(  73)*zin(  58)+xin( 111)*yin( 121)*zin( 106)+xin( 159)*yin( 169)*zin( 154))
          eri_value(  292)=eri_value(  292)+d23bra( 33)*d11ket(  4)*(xin(  14)*yin(  27)*zin(   9)+xin(  62)*yin(  75)*zin(  57)+xin( 110)*yin( 123)*zin( 105)+xin( 158)*yin( 171)*zin( 153))
          eri_value(  293)=eri_value(  293)+d23bra( 33)*d11ket(  5)*(xin(  13)*yin(  28)*zin(   9)+xin(  61)*yin(  76)*zin(  57)+xin( 109)*yin( 124)*zin( 105)+xin( 157)*yin( 172)*zin( 153))
          eri_value(  294)=eri_value(  294)+d23bra( 33)*d11ket(  6)*(xin(  13)*yin(  27)*zin(  10)+xin(  61)*yin(  75)*zin(  58)+xin( 109)*yin( 123)*zin( 106)+xin( 157)*yin( 171)*zin( 154))
          eri_value(  295)=eri_value(  295)+d23bra( 33)*d11ket(  7)*(xin(  14)*yin(  25)*zin(  11)+xin(  62)*yin(  73)*zin(  59)+xin( 110)*yin( 121)*zin( 107)+xin( 158)*yin( 169)*zin( 155))
          eri_value(  296)=eri_value(  296)+d23bra( 33)*d11ket(  8)*(xin(  13)*yin(  26)*zin(  11)+xin(  61)*yin(  74)*zin(  59)+xin( 109)*yin( 122)*zin( 107)+xin( 157)*yin( 170)*zin( 155))
          eri_value(  297)=eri_value(  297)+d23bra( 33)*d11ket(  9)*(xin(  13)*yin(  25)*zin(  12)+xin(  61)*yin(  73)*zin(  60)+xin( 109)*yin( 121)*zin( 108)+xin( 157)*yin( 169)*zin( 156))
          eri_value(  298)=eri_value(  298)+d23bra( 34)*d11ket(  1)*(xin(  20)*yin(  29)*zin(   1)+xin(  68)*yin(  77)*zin(  49)+xin( 116)*yin( 125)*zin(  97)+xin( 164)*yin( 173)*zin( 145))
          eri_value(  299)=eri_value(  299)+d23bra( 34)*d11ket(  2)*(xin(  19)*yin(  30)*zin(   1)+xin(  67)*yin(  78)*zin(  49)+xin( 115)*yin( 126)*zin(  97)+xin( 163)*yin( 174)*zin( 145))
          eri_value(  300)=eri_value(  300)+d23bra( 34)*d11ket(  3)*(xin(  19)*yin(  29)*zin(   2)+xin(  67)*yin(  77)*zin(  50)+xin( 115)*yin( 125)*zin(  98)+xin( 163)*yin( 173)*zin( 146))
          eri_value(  301)=eri_value(  301)+d23bra( 34)*d11ket(  4)*(xin(  18)*yin(  31)*zin(   1)+xin(  66)*yin(  79)*zin(  49)+xin( 114)*yin( 127)*zin(  97)+xin( 162)*yin( 175)*zin( 145))
          eri_value(  302)=eri_value(  302)+d23bra( 34)*d11ket(  5)*(xin(  17)*yin(  32)*zin(   1)+xin(  65)*yin(  80)*zin(  49)+xin( 113)*yin( 128)*zin(  97)+xin( 161)*yin( 176)*zin( 145))
          eri_value(  303)=eri_value(  303)+d23bra( 34)*d11ket(  6)*(xin(  17)*yin(  31)*zin(   2)+xin(  65)*yin(  79)*zin(  50)+xin( 113)*yin( 127)*zin(  98)+xin( 161)*yin( 175)*zin( 146))
          eri_value(  304)=eri_value(  304)+d23bra( 34)*d11ket(  7)*(xin(  18)*yin(  29)*zin(   3)+xin(  66)*yin(  77)*zin(  51)+xin( 114)*yin( 125)*zin(  99)+xin( 162)*yin( 173)*zin( 147))
          eri_value(  305)=eri_value(  305)+d23bra( 34)*d11ket(  8)*(xin(  17)*yin(  30)*zin(   3)+xin(  65)*yin(  78)*zin(  51)+xin( 113)*yin( 126)*zin(  99)+xin( 161)*yin( 174)*zin( 147))
          eri_value(  306)=eri_value(  306)+d23bra( 34)*d11ket(  9)*(xin(  17)*yin(  29)*zin(   4)+xin(  65)*yin(  77)*zin(  52)+xin( 113)*yin( 125)*zin( 100)+xin( 161)*yin( 173)*zin( 148))
          eri_value(  307)=eri_value(  307)+d23bra( 35)*d11ket(  1)*(xin(  20)*yin(  25)*zin(   5)+xin(  68)*yin(  73)*zin(  53)+xin( 116)*yin( 121)*zin( 101)+xin( 164)*yin( 169)*zin( 149))
          eri_value(  308)=eri_value(  308)+d23bra( 35)*d11ket(  2)*(xin(  19)*yin(  26)*zin(   5)+xin(  67)*yin(  74)*zin(  53)+xin( 115)*yin( 122)*zin( 101)+xin( 163)*yin( 170)*zin( 149))
          eri_value(  309)=eri_value(  309)+d23bra( 35)*d11ket(  3)*(xin(  19)*yin(  25)*zin(   6)+xin(  67)*yin(  73)*zin(  54)+xin( 115)*yin( 121)*zin( 102)+xin( 163)*yin( 169)*zin( 150))
          eri_value(  310)=eri_value(  310)+d23bra( 35)*d11ket(  4)*(xin(  18)*yin(  27)*zin(   5)+xin(  66)*yin(  75)*zin(  53)+xin( 114)*yin( 123)*zin( 101)+xin( 162)*yin( 171)*zin( 149))
          eri_value(  311)=eri_value(  311)+d23bra( 35)*d11ket(  5)*(xin(  17)*yin(  28)*zin(   5)+xin(  65)*yin(  76)*zin(  53)+xin( 113)*yin( 124)*zin( 101)+xin( 161)*yin( 172)*zin( 149))
          eri_value(  312)=eri_value(  312)+d23bra( 35)*d11ket(  6)*(xin(  17)*yin(  27)*zin(   6)+xin(  65)*yin(  75)*zin(  54)+xin( 113)*yin( 123)*zin( 102)+xin( 161)*yin( 171)*zin( 150))
          eri_value(  313)=eri_value(  313)+d23bra( 35)*d11ket(  7)*(xin(  18)*yin(  25)*zin(   7)+xin(  66)*yin(  73)*zin(  55)+xin( 114)*yin( 121)*zin( 103)+xin( 162)*yin( 169)*zin( 151))
          eri_value(  314)=eri_value(  314)+d23bra( 35)*d11ket(  8)*(xin(  17)*yin(  26)*zin(   7)+xin(  65)*yin(  74)*zin(  55)+xin( 113)*yin( 122)*zin( 103)+xin( 161)*yin( 170)*zin( 151))
          eri_value(  315)=eri_value(  315)+d23bra( 35)*d11ket(  9)*(xin(  17)*yin(  25)*zin(   8)+xin(  65)*yin(  73)*zin(  56)+xin( 113)*yin( 121)*zin( 104)+xin( 161)*yin( 169)*zin( 152))
          eri_value(  316)=eri_value(  316)+d23bra( 36)*d11ket(  1)*(xin(  16)*yin(  29)*zin(   5)+xin(  64)*yin(  77)*zin(  53)+xin( 112)*yin( 125)*zin( 101)+xin( 160)*yin( 173)*zin( 149))
          eri_value(  317)=eri_value(  317)+d23bra( 36)*d11ket(  2)*(xin(  15)*yin(  30)*zin(   5)+xin(  63)*yin(  78)*zin(  53)+xin( 111)*yin( 126)*zin( 101)+xin( 159)*yin( 174)*zin( 149))
          eri_value(  318)=eri_value(  318)+d23bra( 36)*d11ket(  3)*(xin(  15)*yin(  29)*zin(   6)+xin(  63)*yin(  77)*zin(  54)+xin( 111)*yin( 125)*zin( 102)+xin( 159)*yin( 173)*zin( 150))
          eri_value(  319)=eri_value(  319)+d23bra( 36)*d11ket(  4)*(xin(  14)*yin(  31)*zin(   5)+xin(  62)*yin(  79)*zin(  53)+xin( 110)*yin( 127)*zin( 101)+xin( 158)*yin( 175)*zin( 149))
          eri_value(  320)=eri_value(  320)+d23bra( 36)*d11ket(  5)*(xin(  13)*yin(  32)*zin(   5)+xin(  61)*yin(  80)*zin(  53)+xin( 109)*yin( 128)*zin( 101)+xin( 157)*yin( 176)*zin( 149))
          eri_value(  321)=eri_value(  321)+d23bra( 36)*d11ket(  6)*(xin(  13)*yin(  31)*zin(   6)+xin(  61)*yin(  79)*zin(  54)+xin( 109)*yin( 127)*zin( 102)+xin( 157)*yin( 175)*zin( 150))
          eri_value(  322)=eri_value(  322)+d23bra( 36)*d11ket(  7)*(xin(  14)*yin(  29)*zin(   7)+xin(  62)*yin(  77)*zin(  55)+xin( 110)*yin( 125)*zin( 103)+xin( 158)*yin( 173)*zin( 151))
          eri_value(  323)=eri_value(  323)+d23bra( 36)*d11ket(  8)*(xin(  13)*yin(  30)*zin(   7)+xin(  61)*yin(  78)*zin(  55)+xin( 109)*yin( 126)*zin( 103)+xin( 157)*yin( 174)*zin( 151))
          eri_value(  324)=eri_value(  324)+d23bra( 36)*d11ket(  9)*(xin(  13)*yin(  29)*zin(   8)+xin(  61)*yin(  77)*zin(  56)+xin( 109)*yin( 125)*zin( 104)+xin( 157)*yin( 173)*zin( 152))
          eri_value(  325)=eri_value(  325)+d23bra( 37)*d11ket(  1)*(xin(  12)*yin(  25)*zin(  13)+xin(  60)*yin(  73)*zin(  61)+xin( 108)*yin( 121)*zin( 109)+xin( 156)*yin( 169)*zin( 157))
          eri_value(  326)=eri_value(  326)+d23bra( 37)*d11ket(  2)*(xin(  11)*yin(  26)*zin(  13)+xin(  59)*yin(  74)*zin(  61)+xin( 107)*yin( 122)*zin( 109)+xin( 155)*yin( 170)*zin( 157))
          eri_value(  327)=eri_value(  327)+d23bra( 37)*d11ket(  3)*(xin(  11)*yin(  25)*zin(  14)+xin(  59)*yin(  73)*zin(  62)+xin( 107)*yin( 121)*zin( 110)+xin( 155)*yin( 169)*zin( 158))
          eri_value(  328)=eri_value(  328)+d23bra( 37)*d11ket(  4)*(xin(  10)*yin(  27)*zin(  13)+xin(  58)*yin(  75)*zin(  61)+xin( 106)*yin( 123)*zin( 109)+xin( 154)*yin( 171)*zin( 157))
          eri_value(  329)=eri_value(  329)+d23bra( 37)*d11ket(  5)*(xin(   9)*yin(  28)*zin(  13)+xin(  57)*yin(  76)*zin(  61)+xin( 105)*yin( 124)*zin( 109)+xin( 153)*yin( 172)*zin( 157))
          eri_value(  330)=eri_value(  330)+d23bra( 37)*d11ket(  6)*(xin(   9)*yin(  27)*zin(  14)+xin(  57)*yin(  75)*zin(  62)+xin( 105)*yin( 123)*zin( 110)+xin( 153)*yin( 171)*zin( 158))
          eri_value(  331)=eri_value(  331)+d23bra( 37)*d11ket(  7)*(xin(  10)*yin(  25)*zin(  15)+xin(  58)*yin(  73)*zin(  63)+xin( 106)*yin( 121)*zin( 111)+xin( 154)*yin( 169)*zin( 159))
          eri_value(  332)=eri_value(  332)+d23bra( 37)*d11ket(  8)*(xin(   9)*yin(  26)*zin(  15)+xin(  57)*yin(  74)*zin(  63)+xin( 105)*yin( 122)*zin( 111)+xin( 153)*yin( 170)*zin( 159))
          eri_value(  333)=eri_value(  333)+d23bra( 37)*d11ket(  9)*(xin(   9)*yin(  25)*zin(  16)+xin(  57)*yin(  73)*zin(  64)+xin( 105)*yin( 121)*zin( 112)+xin( 153)*yin( 169)*zin( 160))
          eri_value(  334)=eri_value(  334)+d23bra( 38)*d11ket(  1)*(xin(   4)*yin(  33)*zin(  13)+xin(  52)*yin(  81)*zin(  61)+xin( 100)*yin( 129)*zin( 109)+xin( 148)*yin( 177)*zin( 157))
          eri_value(  335)=eri_value(  335)+d23bra( 38)*d11ket(  2)*(xin(   3)*yin(  34)*zin(  13)+xin(  51)*yin(  82)*zin(  61)+xin(  99)*yin( 130)*zin( 109)+xin( 147)*yin( 178)*zin( 157))
          eri_value(  336)=eri_value(  336)+d23bra( 38)*d11ket(  3)*(xin(   3)*yin(  33)*zin(  14)+xin(  51)*yin(  81)*zin(  62)+xin(  99)*yin( 129)*zin( 110)+xin( 147)*yin( 177)*zin( 158))
          eri_value(  337)=eri_value(  337)+d23bra( 38)*d11ket(  4)*(xin(   2)*yin(  35)*zin(  13)+xin(  50)*yin(  83)*zin(  61)+xin(  98)*yin( 131)*zin( 109)+xin( 146)*yin( 179)*zin( 157))
          eri_value(  338)=eri_value(  338)+d23bra( 38)*d11ket(  5)*(xin(   1)*yin(  36)*zin(  13)+xin(  49)*yin(  84)*zin(  61)+xin(  97)*yin( 132)*zin( 109)+xin( 145)*yin( 180)*zin( 157))
          eri_value(  339)=eri_value(  339)+d23bra( 38)*d11ket(  6)*(xin(   1)*yin(  35)*zin(  14)+xin(  49)*yin(  83)*zin(  62)+xin(  97)*yin( 131)*zin( 110)+xin( 145)*yin( 179)*zin( 158))
          eri_value(  340)=eri_value(  340)+d23bra( 38)*d11ket(  7)*(xin(   2)*yin(  33)*zin(  15)+xin(  50)*yin(  81)*zin(  63)+xin(  98)*yin( 129)*zin( 111)+xin( 146)*yin( 177)*zin( 159))
          eri_value(  341)=eri_value(  341)+d23bra( 38)*d11ket(  8)*(xin(   1)*yin(  34)*zin(  15)+xin(  49)*yin(  82)*zin(  63)+xin(  97)*yin( 130)*zin( 111)+xin( 145)*yin( 178)*zin( 159))
          eri_value(  342)=eri_value(  342)+d23bra( 38)*d11ket(  9)*(xin(   1)*yin(  33)*zin(  16)+xin(  49)*yin(  81)*zin(  64)+xin(  97)*yin( 129)*zin( 112)+xin( 145)*yin( 177)*zin( 160))
          eri_value(  343)=eri_value(  343)+d23bra( 39)*d11ket(  1)*(xin(   4)*yin(  25)*zin(  21)+xin(  52)*yin(  73)*zin(  69)+xin( 100)*yin( 121)*zin( 117)+xin( 148)*yin( 169)*zin( 165))
          eri_value(  344)=eri_value(  344)+d23bra( 39)*d11ket(  2)*(xin(   3)*yin(  26)*zin(  21)+xin(  51)*yin(  74)*zin(  69)+xin(  99)*yin( 122)*zin( 117)+xin( 147)*yin( 170)*zin( 165))
          eri_value(  345)=eri_value(  345)+d23bra( 39)*d11ket(  3)*(xin(   3)*yin(  25)*zin(  22)+xin(  51)*yin(  73)*zin(  70)+xin(  99)*yin( 121)*zin( 118)+xin( 147)*yin( 169)*zin( 166))
          eri_value(  346)=eri_value(  346)+d23bra( 39)*d11ket(  4)*(xin(   2)*yin(  27)*zin(  21)+xin(  50)*yin(  75)*zin(  69)+xin(  98)*yin( 123)*zin( 117)+xin( 146)*yin( 171)*zin( 165))
          eri_value(  347)=eri_value(  347)+d23bra( 39)*d11ket(  5)*(xin(   1)*yin(  28)*zin(  21)+xin(  49)*yin(  76)*zin(  69)+xin(  97)*yin( 124)*zin( 117)+xin( 145)*yin( 172)*zin( 165))
          eri_value(  348)=eri_value(  348)+d23bra( 39)*d11ket(  6)*(xin(   1)*yin(  27)*zin(  22)+xin(  49)*yin(  75)*zin(  70)+xin(  97)*yin( 123)*zin( 118)+xin( 145)*yin( 171)*zin( 166))
          eri_value(  349)=eri_value(  349)+d23bra( 39)*d11ket(  7)*(xin(   2)*yin(  25)*zin(  23)+xin(  50)*yin(  73)*zin(  71)+xin(  98)*yin( 121)*zin( 119)+xin( 146)*yin( 169)*zin( 167))
          eri_value(  350)=eri_value(  350)+d23bra( 39)*d11ket(  8)*(xin(   1)*yin(  26)*zin(  23)+xin(  49)*yin(  74)*zin(  71)+xin(  97)*yin( 122)*zin( 119)+xin( 145)*yin( 170)*zin( 167))
          eri_value(  351)=eri_value(  351)+d23bra( 39)*d11ket(  9)*(xin(   1)*yin(  25)*zin(  24)+xin(  49)*yin(  73)*zin(  72)+xin(  97)*yin( 121)*zin( 120)+xin( 145)*yin( 169)*zin( 168))
          eri_value(  352)=eri_value(  352)+d23bra( 40)*d11ket(  1)*(xin(   8)*yin(  29)*zin(  13)+xin(  56)*yin(  77)*zin(  61)+xin( 104)*yin( 125)*zin( 109)+xin( 152)*yin( 173)*zin( 157))
          eri_value(  353)=eri_value(  353)+d23bra( 40)*d11ket(  2)*(xin(   7)*yin(  30)*zin(  13)+xin(  55)*yin(  78)*zin(  61)+xin( 103)*yin( 126)*zin( 109)+xin( 151)*yin( 174)*zin( 157))
          eri_value(  354)=eri_value(  354)+d23bra( 40)*d11ket(  3)*(xin(   7)*yin(  29)*zin(  14)+xin(  55)*yin(  77)*zin(  62)+xin( 103)*yin( 125)*zin( 110)+xin( 151)*yin( 173)*zin( 158))
          eri_value(  355)=eri_value(  355)+d23bra( 40)*d11ket(  4)*(xin(   6)*yin(  31)*zin(  13)+xin(  54)*yin(  79)*zin(  61)+xin( 102)*yin( 127)*zin( 109)+xin( 150)*yin( 175)*zin( 157))
          eri_value(  356)=eri_value(  356)+d23bra( 40)*d11ket(  5)*(xin(   5)*yin(  32)*zin(  13)+xin(  53)*yin(  80)*zin(  61)+xin( 101)*yin( 128)*zin( 109)+xin( 149)*yin( 176)*zin( 157))
          eri_value(  357)=eri_value(  357)+d23bra( 40)*d11ket(  6)*(xin(   5)*yin(  31)*zin(  14)+xin(  53)*yin(  79)*zin(  62)+xin( 101)*yin( 127)*zin( 110)+xin( 149)*yin( 175)*zin( 158))
          eri_value(  358)=eri_value(  358)+d23bra( 40)*d11ket(  7)*(xin(   6)*yin(  29)*zin(  15)+xin(  54)*yin(  77)*zin(  63)+xin( 102)*yin( 125)*zin( 111)+xin( 150)*yin( 173)*zin( 159))
          eri_value(  359)=eri_value(  359)+d23bra( 40)*d11ket(  8)*(xin(   5)*yin(  30)*zin(  15)+xin(  53)*yin(  78)*zin(  63)+xin( 101)*yin( 126)*zin( 111)+xin( 149)*yin( 174)*zin( 159))
          eri_value(  360)=eri_value(  360)+d23bra( 40)*d11ket(  9)*(xin(   5)*yin(  29)*zin(  16)+xin(  53)*yin(  77)*zin(  64)+xin( 101)*yin( 125)*zin( 112)+xin( 149)*yin( 173)*zin( 160))
          eri_value(  361)=eri_value(  361)+d23bra( 41)*d11ket(  1)*(xin(   8)*yin(  25)*zin(  17)+xin(  56)*yin(  73)*zin(  65)+xin( 104)*yin( 121)*zin( 113)+xin( 152)*yin( 169)*zin( 161))
          eri_value(  362)=eri_value(  362)+d23bra( 41)*d11ket(  2)*(xin(   7)*yin(  26)*zin(  17)+xin(  55)*yin(  74)*zin(  65)+xin( 103)*yin( 122)*zin( 113)+xin( 151)*yin( 170)*zin( 161))
          eri_value(  363)=eri_value(  363)+d23bra( 41)*d11ket(  3)*(xin(   7)*yin(  25)*zin(  18)+xin(  55)*yin(  73)*zin(  66)+xin( 103)*yin( 121)*zin( 114)+xin( 151)*yin( 169)*zin( 162))
          eri_value(  364)=eri_value(  364)+d23bra( 41)*d11ket(  4)*(xin(   6)*yin(  27)*zin(  17)+xin(  54)*yin(  75)*zin(  65)+xin( 102)*yin( 123)*zin( 113)+xin( 150)*yin( 171)*zin( 161))
          eri_value(  365)=eri_value(  365)+d23bra( 41)*d11ket(  5)*(xin(   5)*yin(  28)*zin(  17)+xin(  53)*yin(  76)*zin(  65)+xin( 101)*yin( 124)*zin( 113)+xin( 149)*yin( 172)*zin( 161))
          eri_value(  366)=eri_value(  366)+d23bra( 41)*d11ket(  6)*(xin(   5)*yin(  27)*zin(  18)+xin(  53)*yin(  75)*zin(  66)+xin( 101)*yin( 123)*zin( 114)+xin( 149)*yin( 171)*zin( 162))
          eri_value(  367)=eri_value(  367)+d23bra( 41)*d11ket(  7)*(xin(   6)*yin(  25)*zin(  19)+xin(  54)*yin(  73)*zin(  67)+xin( 102)*yin( 121)*zin( 115)+xin( 150)*yin( 169)*zin( 163))
          eri_value(  368)=eri_value(  368)+d23bra( 41)*d11ket(  8)*(xin(   5)*yin(  26)*zin(  19)+xin(  53)*yin(  74)*zin(  67)+xin( 101)*yin( 122)*zin( 115)+xin( 149)*yin( 170)*zin( 163))
          eri_value(  369)=eri_value(  369)+d23bra( 41)*d11ket(  9)*(xin(   5)*yin(  25)*zin(  20)+xin(  53)*yin(  73)*zin(  68)+xin( 101)*yin( 121)*zin( 116)+xin( 149)*yin( 169)*zin( 164))
          eri_value(  370)=eri_value(  370)+d23bra( 42)*d11ket(  1)*(xin(   4)*yin(  29)*zin(  17)+xin(  52)*yin(  77)*zin(  65)+xin( 100)*yin( 125)*zin( 113)+xin( 148)*yin( 173)*zin( 161))
          eri_value(  371)=eri_value(  371)+d23bra( 42)*d11ket(  2)*(xin(   3)*yin(  30)*zin(  17)+xin(  51)*yin(  78)*zin(  65)+xin(  99)*yin( 126)*zin( 113)+xin( 147)*yin( 174)*zin( 161))
          eri_value(  372)=eri_value(  372)+d23bra( 42)*d11ket(  3)*(xin(   3)*yin(  29)*zin(  18)+xin(  51)*yin(  77)*zin(  66)+xin(  99)*yin( 125)*zin( 114)+xin( 147)*yin( 173)*zin( 162))
          eri_value(  373)=eri_value(  373)+d23bra( 42)*d11ket(  4)*(xin(   2)*yin(  31)*zin(  17)+xin(  50)*yin(  79)*zin(  65)+xin(  98)*yin( 127)*zin( 113)+xin( 146)*yin( 175)*zin( 161))
          eri_value(  374)=eri_value(  374)+d23bra( 42)*d11ket(  5)*(xin(   1)*yin(  32)*zin(  17)+xin(  49)*yin(  80)*zin(  65)+xin(  97)*yin( 128)*zin( 113)+xin( 145)*yin( 176)*zin( 161))
          eri_value(  375)=eri_value(  375)+d23bra( 42)*d11ket(  6)*(xin(   1)*yin(  31)*zin(  18)+xin(  49)*yin(  79)*zin(  66)+xin(  97)*yin( 127)*zin( 114)+xin( 145)*yin( 175)*zin( 162))
          eri_value(  376)=eri_value(  376)+d23bra( 42)*d11ket(  7)*(xin(   2)*yin(  29)*zin(  19)+xin(  50)*yin(  77)*zin(  67)+xin(  98)*yin( 125)*zin( 115)+xin( 146)*yin( 173)*zin( 163))
          eri_value(  377)=eri_value(  377)+d23bra( 42)*d11ket(  8)*(xin(   1)*yin(  30)*zin(  19)+xin(  49)*yin(  78)*zin(  67)+xin(  97)*yin( 126)*zin( 115)+xin( 145)*yin( 174)*zin( 163))
          eri_value(  378)=eri_value(  378)+d23bra( 42)*d11ket(  9)*(xin(   1)*yin(  29)*zin(  20)+xin(  49)*yin(  77)*zin(  68)+xin(  97)*yin( 125)*zin( 116)+xin( 145)*yin( 173)*zin( 164))
          eri_value(  379)=eri_value(  379)+d23bra( 43)*d11ket(  1)*(xin(  24)*yin(   1)*zin(  25)+xin(  72)*yin(  49)*zin(  73)+xin( 120)*yin(  97)*zin( 121)+xin( 168)*yin( 145)*zin( 169))
          eri_value(  380)=eri_value(  380)+d23bra( 43)*d11ket(  2)*(xin(  23)*yin(   2)*zin(  25)+xin(  71)*yin(  50)*zin(  73)+xin( 119)*yin(  98)*zin( 121)+xin( 167)*yin( 146)*zin( 169))
          eri_value(  381)=eri_value(  381)+d23bra( 43)*d11ket(  3)*(xin(  23)*yin(   1)*zin(  26)+xin(  71)*yin(  49)*zin(  74)+xin( 119)*yin(  97)*zin( 122)+xin( 167)*yin( 145)*zin( 170))
          eri_value(  382)=eri_value(  382)+d23bra( 43)*d11ket(  4)*(xin(  22)*yin(   3)*zin(  25)+xin(  70)*yin(  51)*zin(  73)+xin( 118)*yin(  99)*zin( 121)+xin( 166)*yin( 147)*zin( 169))
          eri_value(  383)=eri_value(  383)+d23bra( 43)*d11ket(  5)*(xin(  21)*yin(   4)*zin(  25)+xin(  69)*yin(  52)*zin(  73)+xin( 117)*yin( 100)*zin( 121)+xin( 165)*yin( 148)*zin( 169))
          eri_value(  384)=eri_value(  384)+d23bra( 43)*d11ket(  6)*(xin(  21)*yin(   3)*zin(  26)+xin(  69)*yin(  51)*zin(  74)+xin( 117)*yin(  99)*zin( 122)+xin( 165)*yin( 147)*zin( 170))
          eri_value(  385)=eri_value(  385)+d23bra( 43)*d11ket(  7)*(xin(  22)*yin(   1)*zin(  27)+xin(  70)*yin(  49)*zin(  75)+xin( 118)*yin(  97)*zin( 123)+xin( 166)*yin( 145)*zin( 171))
          eri_value(  386)=eri_value(  386)+d23bra( 43)*d11ket(  8)*(xin(  21)*yin(   2)*zin(  27)+xin(  69)*yin(  50)*zin(  75)+xin( 117)*yin(  98)*zin( 123)+xin( 165)*yin( 146)*zin( 171))
          eri_value(  387)=eri_value(  387)+d23bra( 43)*d11ket(  9)*(xin(  21)*yin(   1)*zin(  28)+xin(  69)*yin(  49)*zin(  76)+xin( 117)*yin(  97)*zin( 124)+xin( 165)*yin( 145)*zin( 172))
          eri_value(  388)=eri_value(  388)+d23bra( 44)*d11ket(  1)*(xin(  16)*yin(   9)*zin(  25)+xin(  64)*yin(  57)*zin(  73)+xin( 112)*yin( 105)*zin( 121)+xin( 160)*yin( 153)*zin( 169))
          eri_value(  389)=eri_value(  389)+d23bra( 44)*d11ket(  2)*(xin(  15)*yin(  10)*zin(  25)+xin(  63)*yin(  58)*zin(  73)+xin( 111)*yin( 106)*zin( 121)+xin( 159)*yin( 154)*zin( 169))
          eri_value(  390)=eri_value(  390)+d23bra( 44)*d11ket(  3)*(xin(  15)*yin(   9)*zin(  26)+xin(  63)*yin(  57)*zin(  74)+xin( 111)*yin( 105)*zin( 122)+xin( 159)*yin( 153)*zin( 170))
          eri_value(  391)=eri_value(  391)+d23bra( 44)*d11ket(  4)*(xin(  14)*yin(  11)*zin(  25)+xin(  62)*yin(  59)*zin(  73)+xin( 110)*yin( 107)*zin( 121)+xin( 158)*yin( 155)*zin( 169))
          eri_value(  392)=eri_value(  392)+d23bra( 44)*d11ket(  5)*(xin(  13)*yin(  12)*zin(  25)+xin(  61)*yin(  60)*zin(  73)+xin( 109)*yin( 108)*zin( 121)+xin( 157)*yin( 156)*zin( 169))
          eri_value(  393)=eri_value(  393)+d23bra( 44)*d11ket(  6)*(xin(  13)*yin(  11)*zin(  26)+xin(  61)*yin(  59)*zin(  74)+xin( 109)*yin( 107)*zin( 122)+xin( 157)*yin( 155)*zin( 170))
          eri_value(  394)=eri_value(  394)+d23bra( 44)*d11ket(  7)*(xin(  14)*yin(   9)*zin(  27)+xin(  62)*yin(  57)*zin(  75)+xin( 110)*yin( 105)*zin( 123)+xin( 158)*yin( 153)*zin( 171))
          eri_value(  395)=eri_value(  395)+d23bra( 44)*d11ket(  8)*(xin(  13)*yin(  10)*zin(  27)+xin(  61)*yin(  58)*zin(  75)+xin( 109)*yin( 106)*zin( 123)+xin( 157)*yin( 154)*zin( 171))
          eri_value(  396)=eri_value(  396)+d23bra( 44)*d11ket(  9)*(xin(  13)*yin(   9)*zin(  28)+xin(  61)*yin(  57)*zin(  76)+xin( 109)*yin( 105)*zin( 124)+xin( 157)*yin( 153)*zin( 172))
          eri_value(  397)=eri_value(  397)+d23bra( 45)*d11ket(  1)*(xin(  16)*yin(   1)*zin(  33)+xin(  64)*yin(  49)*zin(  81)+xin( 112)*yin(  97)*zin( 129)+xin( 160)*yin( 145)*zin( 177))
          eri_value(  398)=eri_value(  398)+d23bra( 45)*d11ket(  2)*(xin(  15)*yin(   2)*zin(  33)+xin(  63)*yin(  50)*zin(  81)+xin( 111)*yin(  98)*zin( 129)+xin( 159)*yin( 146)*zin( 177))
          eri_value(  399)=eri_value(  399)+d23bra( 45)*d11ket(  3)*(xin(  15)*yin(   1)*zin(  34)+xin(  63)*yin(  49)*zin(  82)+xin( 111)*yin(  97)*zin( 130)+xin( 159)*yin( 145)*zin( 178))
          eri_value(  400)=eri_value(  400)+d23bra( 45)*d11ket(  4)*(xin(  14)*yin(   3)*zin(  33)+xin(  62)*yin(  51)*zin(  81)+xin( 110)*yin(  99)*zin( 129)+xin( 158)*yin( 147)*zin( 177))
          eri_value(  401)=eri_value(  401)+d23bra( 45)*d11ket(  5)*(xin(  13)*yin(   4)*zin(  33)+xin(  61)*yin(  52)*zin(  81)+xin( 109)*yin( 100)*zin( 129)+xin( 157)*yin( 148)*zin( 177))
          eri_value(  402)=eri_value(  402)+d23bra( 45)*d11ket(  6)*(xin(  13)*yin(   3)*zin(  34)+xin(  61)*yin(  51)*zin(  82)+xin( 109)*yin(  99)*zin( 130)+xin( 157)*yin( 147)*zin( 178))
          eri_value(  403)=eri_value(  403)+d23bra( 45)*d11ket(  7)*(xin(  14)*yin(   1)*zin(  35)+xin(  62)*yin(  49)*zin(  83)+xin( 110)*yin(  97)*zin( 131)+xin( 158)*yin( 145)*zin( 179))
          eri_value(  404)=eri_value(  404)+d23bra( 45)*d11ket(  8)*(xin(  13)*yin(   2)*zin(  35)+xin(  61)*yin(  50)*zin(  83)+xin( 109)*yin(  98)*zin( 131)+xin( 157)*yin( 146)*zin( 179))
          eri_value(  405)=eri_value(  405)+d23bra( 45)*d11ket(  9)*(xin(  13)*yin(   1)*zin(  36)+xin(  61)*yin(  49)*zin(  84)+xin( 109)*yin(  97)*zin( 132)+xin( 157)*yin( 145)*zin( 180))
          eri_value(  406)=eri_value(  406)+d23bra( 46)*d11ket(  1)*(xin(  20)*yin(   5)*zin(  25)+xin(  68)*yin(  53)*zin(  73)+xin( 116)*yin( 101)*zin( 121)+xin( 164)*yin( 149)*zin( 169))
          eri_value(  407)=eri_value(  407)+d23bra( 46)*d11ket(  2)*(xin(  19)*yin(   6)*zin(  25)+xin(  67)*yin(  54)*zin(  73)+xin( 115)*yin( 102)*zin( 121)+xin( 163)*yin( 150)*zin( 169))
          eri_value(  408)=eri_value(  408)+d23bra( 46)*d11ket(  3)*(xin(  19)*yin(   5)*zin(  26)+xin(  67)*yin(  53)*zin(  74)+xin( 115)*yin( 101)*zin( 122)+xin( 163)*yin( 149)*zin( 170))
          eri_value(  409)=eri_value(  409)+d23bra( 46)*d11ket(  4)*(xin(  18)*yin(   7)*zin(  25)+xin(  66)*yin(  55)*zin(  73)+xin( 114)*yin( 103)*zin( 121)+xin( 162)*yin( 151)*zin( 169))
          eri_value(  410)=eri_value(  410)+d23bra( 46)*d11ket(  5)*(xin(  17)*yin(   8)*zin(  25)+xin(  65)*yin(  56)*zin(  73)+xin( 113)*yin( 104)*zin( 121)+xin( 161)*yin( 152)*zin( 169))
          eri_value(  411)=eri_value(  411)+d23bra( 46)*d11ket(  6)*(xin(  17)*yin(   7)*zin(  26)+xin(  65)*yin(  55)*zin(  74)+xin( 113)*yin( 103)*zin( 122)+xin( 161)*yin( 151)*zin( 170))
          eri_value(  412)=eri_value(  412)+d23bra( 46)*d11ket(  7)*(xin(  18)*yin(   5)*zin(  27)+xin(  66)*yin(  53)*zin(  75)+xin( 114)*yin( 101)*zin( 123)+xin( 162)*yin( 149)*zin( 171))
          eri_value(  413)=eri_value(  413)+d23bra( 46)*d11ket(  8)*(xin(  17)*yin(   6)*zin(  27)+xin(  65)*yin(  54)*zin(  75)+xin( 113)*yin( 102)*zin( 123)+xin( 161)*yin( 150)*zin( 171))
          eri_value(  414)=eri_value(  414)+d23bra( 46)*d11ket(  9)*(xin(  17)*yin(   5)*zin(  28)+xin(  65)*yin(  53)*zin(  76)+xin( 113)*yin( 101)*zin( 124)+xin( 161)*yin( 149)*zin( 172))
          eri_value(  415)=eri_value(  415)+d23bra( 47)*d11ket(  1)*(xin(  20)*yin(   1)*zin(  29)+xin(  68)*yin(  49)*zin(  77)+xin( 116)*yin(  97)*zin( 125)+xin( 164)*yin( 145)*zin( 173))
          eri_value(  416)=eri_value(  416)+d23bra( 47)*d11ket(  2)*(xin(  19)*yin(   2)*zin(  29)+xin(  67)*yin(  50)*zin(  77)+xin( 115)*yin(  98)*zin( 125)+xin( 163)*yin( 146)*zin( 173))
          eri_value(  417)=eri_value(  417)+d23bra( 47)*d11ket(  3)*(xin(  19)*yin(   1)*zin(  30)+xin(  67)*yin(  49)*zin(  78)+xin( 115)*yin(  97)*zin( 126)+xin( 163)*yin( 145)*zin( 174))
          eri_value(  418)=eri_value(  418)+d23bra( 47)*d11ket(  4)*(xin(  18)*yin(   3)*zin(  29)+xin(  66)*yin(  51)*zin(  77)+xin( 114)*yin(  99)*zin( 125)+xin( 162)*yin( 147)*zin( 173))
          eri_value(  419)=eri_value(  419)+d23bra( 47)*d11ket(  5)*(xin(  17)*yin(   4)*zin(  29)+xin(  65)*yin(  52)*zin(  77)+xin( 113)*yin( 100)*zin( 125)+xin( 161)*yin( 148)*zin( 173))
          eri_value(  420)=eri_value(  420)+d23bra( 47)*d11ket(  6)*(xin(  17)*yin(   3)*zin(  30)+xin(  65)*yin(  51)*zin(  78)+xin( 113)*yin(  99)*zin( 126)+xin( 161)*yin( 147)*zin( 174))
          eri_value(  421)=eri_value(  421)+d23bra( 47)*d11ket(  7)*(xin(  18)*yin(   1)*zin(  31)+xin(  66)*yin(  49)*zin(  79)+xin( 114)*yin(  97)*zin( 127)+xin( 162)*yin( 145)*zin( 175))
          eri_value(  422)=eri_value(  422)+d23bra( 47)*d11ket(  8)*(xin(  17)*yin(   2)*zin(  31)+xin(  65)*yin(  50)*zin(  79)+xin( 113)*yin(  98)*zin( 127)+xin( 161)*yin( 146)*zin( 175))
          eri_value(  423)=eri_value(  423)+d23bra( 47)*d11ket(  9)*(xin(  17)*yin(   1)*zin(  32)+xin(  65)*yin(  49)*zin(  80)+xin( 113)*yin(  97)*zin( 128)+xin( 161)*yin( 145)*zin( 176))
          eri_value(  424)=eri_value(  424)+d23bra( 48)*d11ket(  1)*(xin(  16)*yin(   5)*zin(  29)+xin(  64)*yin(  53)*zin(  77)+xin( 112)*yin( 101)*zin( 125)+xin( 160)*yin( 149)*zin( 173))
          eri_value(  425)=eri_value(  425)+d23bra( 48)*d11ket(  2)*(xin(  15)*yin(   6)*zin(  29)+xin(  63)*yin(  54)*zin(  77)+xin( 111)*yin( 102)*zin( 125)+xin( 159)*yin( 150)*zin( 173))
          eri_value(  426)=eri_value(  426)+d23bra( 48)*d11ket(  3)*(xin(  15)*yin(   5)*zin(  30)+xin(  63)*yin(  53)*zin(  78)+xin( 111)*yin( 101)*zin( 126)+xin( 159)*yin( 149)*zin( 174))
          eri_value(  427)=eri_value(  427)+d23bra( 48)*d11ket(  4)*(xin(  14)*yin(   7)*zin(  29)+xin(  62)*yin(  55)*zin(  77)+xin( 110)*yin( 103)*zin( 125)+xin( 158)*yin( 151)*zin( 173))
          eri_value(  428)=eri_value(  428)+d23bra( 48)*d11ket(  5)*(xin(  13)*yin(   8)*zin(  29)+xin(  61)*yin(  56)*zin(  77)+xin( 109)*yin( 104)*zin( 125)+xin( 157)*yin( 152)*zin( 173))
          eri_value(  429)=eri_value(  429)+d23bra( 48)*d11ket(  6)*(xin(  13)*yin(   7)*zin(  30)+xin(  61)*yin(  55)*zin(  78)+xin( 109)*yin( 103)*zin( 126)+xin( 157)*yin( 151)*zin( 174))
          eri_value(  430)=eri_value(  430)+d23bra( 48)*d11ket(  7)*(xin(  14)*yin(   5)*zin(  31)+xin(  62)*yin(  53)*zin(  79)+xin( 110)*yin( 101)*zin( 127)+xin( 158)*yin( 149)*zin( 175))
          eri_value(  431)=eri_value(  431)+d23bra( 48)*d11ket(  8)*(xin(  13)*yin(   6)*zin(  31)+xin(  61)*yin(  54)*zin(  79)+xin( 109)*yin( 102)*zin( 127)+xin( 157)*yin( 150)*zin( 175))
          eri_value(  432)=eri_value(  432)+d23bra( 48)*d11ket(  9)*(xin(  13)*yin(   5)*zin(  32)+xin(  61)*yin(  53)*zin(  80)+xin( 109)*yin( 101)*zin( 128)+xin( 157)*yin( 149)*zin( 176))
          eri_value(  433)=eri_value(  433)+d23bra( 49)*d11ket(  1)*(xin(  12)*yin(  13)*zin(  25)+xin(  60)*yin(  61)*zin(  73)+xin( 108)*yin( 109)*zin( 121)+xin( 156)*yin( 157)*zin( 169))
          eri_value(  434)=eri_value(  434)+d23bra( 49)*d11ket(  2)*(xin(  11)*yin(  14)*zin(  25)+xin(  59)*yin(  62)*zin(  73)+xin( 107)*yin( 110)*zin( 121)+xin( 155)*yin( 158)*zin( 169))
          eri_value(  435)=eri_value(  435)+d23bra( 49)*d11ket(  3)*(xin(  11)*yin(  13)*zin(  26)+xin(  59)*yin(  61)*zin(  74)+xin( 107)*yin( 109)*zin( 122)+xin( 155)*yin( 157)*zin( 170))
          eri_value(  436)=eri_value(  436)+d23bra( 49)*d11ket(  4)*(xin(  10)*yin(  15)*zin(  25)+xin(  58)*yin(  63)*zin(  73)+xin( 106)*yin( 111)*zin( 121)+xin( 154)*yin( 159)*zin( 169))
          eri_value(  437)=eri_value(  437)+d23bra( 49)*d11ket(  5)*(xin(   9)*yin(  16)*zin(  25)+xin(  57)*yin(  64)*zin(  73)+xin( 105)*yin( 112)*zin( 121)+xin( 153)*yin( 160)*zin( 169))
          eri_value(  438)=eri_value(  438)+d23bra( 49)*d11ket(  6)*(xin(   9)*yin(  15)*zin(  26)+xin(  57)*yin(  63)*zin(  74)+xin( 105)*yin( 111)*zin( 122)+xin( 153)*yin( 159)*zin( 170))
          eri_value(  439)=eri_value(  439)+d23bra( 49)*d11ket(  7)*(xin(  10)*yin(  13)*zin(  27)+xin(  58)*yin(  61)*zin(  75)+xin( 106)*yin( 109)*zin( 123)+xin( 154)*yin( 157)*zin( 171))
          eri_value(  440)=eri_value(  440)+d23bra( 49)*d11ket(  8)*(xin(   9)*yin(  14)*zin(  27)+xin(  57)*yin(  62)*zin(  75)+xin( 105)*yin( 110)*zin( 123)+xin( 153)*yin( 158)*zin( 171))
          eri_value(  441)=eri_value(  441)+d23bra( 49)*d11ket(  9)*(xin(   9)*yin(  13)*zin(  28)+xin(  57)*yin(  61)*zin(  76)+xin( 105)*yin( 109)*zin( 124)+xin( 153)*yin( 157)*zin( 172))
          eri_value(  442)=eri_value(  442)+d23bra( 50)*d11ket(  1)*(xin(   4)*yin(  21)*zin(  25)+xin(  52)*yin(  69)*zin(  73)+xin( 100)*yin( 117)*zin( 121)+xin( 148)*yin( 165)*zin( 169))
          eri_value(  443)=eri_value(  443)+d23bra( 50)*d11ket(  2)*(xin(   3)*yin(  22)*zin(  25)+xin(  51)*yin(  70)*zin(  73)+xin(  99)*yin( 118)*zin( 121)+xin( 147)*yin( 166)*zin( 169))
          eri_value(  444)=eri_value(  444)+d23bra( 50)*d11ket(  3)*(xin(   3)*yin(  21)*zin(  26)+xin(  51)*yin(  69)*zin(  74)+xin(  99)*yin( 117)*zin( 122)+xin( 147)*yin( 165)*zin( 170))
          eri_value(  445)=eri_value(  445)+d23bra( 50)*d11ket(  4)*(xin(   2)*yin(  23)*zin(  25)+xin(  50)*yin(  71)*zin(  73)+xin(  98)*yin( 119)*zin( 121)+xin( 146)*yin( 167)*zin( 169))
          eri_value(  446)=eri_value(  446)+d23bra( 50)*d11ket(  5)*(xin(   1)*yin(  24)*zin(  25)+xin(  49)*yin(  72)*zin(  73)+xin(  97)*yin( 120)*zin( 121)+xin( 145)*yin( 168)*zin( 169))
          eri_value(  447)=eri_value(  447)+d23bra( 50)*d11ket(  6)*(xin(   1)*yin(  23)*zin(  26)+xin(  49)*yin(  71)*zin(  74)+xin(  97)*yin( 119)*zin( 122)+xin( 145)*yin( 167)*zin( 170))
          eri_value(  448)=eri_value(  448)+d23bra( 50)*d11ket(  7)*(xin(   2)*yin(  21)*zin(  27)+xin(  50)*yin(  69)*zin(  75)+xin(  98)*yin( 117)*zin( 123)+xin( 146)*yin( 165)*zin( 171))
          eri_value(  449)=eri_value(  449)+d23bra( 50)*d11ket(  8)*(xin(   1)*yin(  22)*zin(  27)+xin(  49)*yin(  70)*zin(  75)+xin(  97)*yin( 118)*zin( 123)+xin( 145)*yin( 166)*zin( 171))
          eri_value(  450)=eri_value(  450)+d23bra( 50)*d11ket(  9)*(xin(   1)*yin(  21)*zin(  28)+xin(  49)*yin(  69)*zin(  76)+xin(  97)*yin( 117)*zin( 124)+xin( 145)*yin( 165)*zin( 172))
          eri_value(  451)=eri_value(  451)+d23bra( 51)*d11ket(  1)*(xin(   4)*yin(  13)*zin(  33)+xin(  52)*yin(  61)*zin(  81)+xin( 100)*yin( 109)*zin( 129)+xin( 148)*yin( 157)*zin( 177))
          eri_value(  452)=eri_value(  452)+d23bra( 51)*d11ket(  2)*(xin(   3)*yin(  14)*zin(  33)+xin(  51)*yin(  62)*zin(  81)+xin(  99)*yin( 110)*zin( 129)+xin( 147)*yin( 158)*zin( 177))
          eri_value(  453)=eri_value(  453)+d23bra( 51)*d11ket(  3)*(xin(   3)*yin(  13)*zin(  34)+xin(  51)*yin(  61)*zin(  82)+xin(  99)*yin( 109)*zin( 130)+xin( 147)*yin( 157)*zin( 178))
          eri_value(  454)=eri_value(  454)+d23bra( 51)*d11ket(  4)*(xin(   2)*yin(  15)*zin(  33)+xin(  50)*yin(  63)*zin(  81)+xin(  98)*yin( 111)*zin( 129)+xin( 146)*yin( 159)*zin( 177))
          eri_value(  455)=eri_value(  455)+d23bra( 51)*d11ket(  5)*(xin(   1)*yin(  16)*zin(  33)+xin(  49)*yin(  64)*zin(  81)+xin(  97)*yin( 112)*zin( 129)+xin( 145)*yin( 160)*zin( 177))
          eri_value(  456)=eri_value(  456)+d23bra( 51)*d11ket(  6)*(xin(   1)*yin(  15)*zin(  34)+xin(  49)*yin(  63)*zin(  82)+xin(  97)*yin( 111)*zin( 130)+xin( 145)*yin( 159)*zin( 178))
          eri_value(  457)=eri_value(  457)+d23bra( 51)*d11ket(  7)*(xin(   2)*yin(  13)*zin(  35)+xin(  50)*yin(  61)*zin(  83)+xin(  98)*yin( 109)*zin( 131)+xin( 146)*yin( 157)*zin( 179))
          eri_value(  458)=eri_value(  458)+d23bra( 51)*d11ket(  8)*(xin(   1)*yin(  14)*zin(  35)+xin(  49)*yin(  62)*zin(  83)+xin(  97)*yin( 110)*zin( 131)+xin( 145)*yin( 158)*zin( 179))
          eri_value(  459)=eri_value(  459)+d23bra( 51)*d11ket(  9)*(xin(   1)*yin(  13)*zin(  36)+xin(  49)*yin(  61)*zin(  84)+xin(  97)*yin( 109)*zin( 132)+xin( 145)*yin( 157)*zin( 180))
          eri_value(  460)=eri_value(  460)+d23bra( 52)*d11ket(  1)*(xin(   8)*yin(  17)*zin(  25)+xin(  56)*yin(  65)*zin(  73)+xin( 104)*yin( 113)*zin( 121)+xin( 152)*yin( 161)*zin( 169))
          eri_value(  461)=eri_value(  461)+d23bra( 52)*d11ket(  2)*(xin(   7)*yin(  18)*zin(  25)+xin(  55)*yin(  66)*zin(  73)+xin( 103)*yin( 114)*zin( 121)+xin( 151)*yin( 162)*zin( 169))
          eri_value(  462)=eri_value(  462)+d23bra( 52)*d11ket(  3)*(xin(   7)*yin(  17)*zin(  26)+xin(  55)*yin(  65)*zin(  74)+xin( 103)*yin( 113)*zin( 122)+xin( 151)*yin( 161)*zin( 170))
          eri_value(  463)=eri_value(  463)+d23bra( 52)*d11ket(  4)*(xin(   6)*yin(  19)*zin(  25)+xin(  54)*yin(  67)*zin(  73)+xin( 102)*yin( 115)*zin( 121)+xin( 150)*yin( 163)*zin( 169))
          eri_value(  464)=eri_value(  464)+d23bra( 52)*d11ket(  5)*(xin(   5)*yin(  20)*zin(  25)+xin(  53)*yin(  68)*zin(  73)+xin( 101)*yin( 116)*zin( 121)+xin( 149)*yin( 164)*zin( 169))
          eri_value(  465)=eri_value(  465)+d23bra( 52)*d11ket(  6)*(xin(   5)*yin(  19)*zin(  26)+xin(  53)*yin(  67)*zin(  74)+xin( 101)*yin( 115)*zin( 122)+xin( 149)*yin( 163)*zin( 170))
          eri_value(  466)=eri_value(  466)+d23bra( 52)*d11ket(  7)*(xin(   6)*yin(  17)*zin(  27)+xin(  54)*yin(  65)*zin(  75)+xin( 102)*yin( 113)*zin( 123)+xin( 150)*yin( 161)*zin( 171))
          eri_value(  467)=eri_value(  467)+d23bra( 52)*d11ket(  8)*(xin(   5)*yin(  18)*zin(  27)+xin(  53)*yin(  66)*zin(  75)+xin( 101)*yin( 114)*zin( 123)+xin( 149)*yin( 162)*zin( 171))
          eri_value(  468)=eri_value(  468)+d23bra( 52)*d11ket(  9)*(xin(   5)*yin(  17)*zin(  28)+xin(  53)*yin(  65)*zin(  76)+xin( 101)*yin( 113)*zin( 124)+xin( 149)*yin( 161)*zin( 172))
          eri_value(  469)=eri_value(  469)+d23bra( 53)*d11ket(  1)*(xin(   8)*yin(  13)*zin(  29)+xin(  56)*yin(  61)*zin(  77)+xin( 104)*yin( 109)*zin( 125)+xin( 152)*yin( 157)*zin( 173))
          eri_value(  470)=eri_value(  470)+d23bra( 53)*d11ket(  2)*(xin(   7)*yin(  14)*zin(  29)+xin(  55)*yin(  62)*zin(  77)+xin( 103)*yin( 110)*zin( 125)+xin( 151)*yin( 158)*zin( 173))
          eri_value(  471)=eri_value(  471)+d23bra( 53)*d11ket(  3)*(xin(   7)*yin(  13)*zin(  30)+xin(  55)*yin(  61)*zin(  78)+xin( 103)*yin( 109)*zin( 126)+xin( 151)*yin( 157)*zin( 174))
          eri_value(  472)=eri_value(  472)+d23bra( 53)*d11ket(  4)*(xin(   6)*yin(  15)*zin(  29)+xin(  54)*yin(  63)*zin(  77)+xin( 102)*yin( 111)*zin( 125)+xin( 150)*yin( 159)*zin( 173))
          eri_value(  473)=eri_value(  473)+d23bra( 53)*d11ket(  5)*(xin(   5)*yin(  16)*zin(  29)+xin(  53)*yin(  64)*zin(  77)+xin( 101)*yin( 112)*zin( 125)+xin( 149)*yin( 160)*zin( 173))
          eri_value(  474)=eri_value(  474)+d23bra( 53)*d11ket(  6)*(xin(   5)*yin(  15)*zin(  30)+xin(  53)*yin(  63)*zin(  78)+xin( 101)*yin( 111)*zin( 126)+xin( 149)*yin( 159)*zin( 174))
          eri_value(  475)=eri_value(  475)+d23bra( 53)*d11ket(  7)*(xin(   6)*yin(  13)*zin(  31)+xin(  54)*yin(  61)*zin(  79)+xin( 102)*yin( 109)*zin( 127)+xin( 150)*yin( 157)*zin( 175))
          eri_value(  476)=eri_value(  476)+d23bra( 53)*d11ket(  8)*(xin(   5)*yin(  14)*zin(  31)+xin(  53)*yin(  62)*zin(  79)+xin( 101)*yin( 110)*zin( 127)+xin( 149)*yin( 158)*zin( 175))
          eri_value(  477)=eri_value(  477)+d23bra( 53)*d11ket(  9)*(xin(   5)*yin(  13)*zin(  32)+xin(  53)*yin(  61)*zin(  80)+xin( 101)*yin( 109)*zin( 128)+xin( 149)*yin( 157)*zin( 176))
          eri_value(  478)=eri_value(  478)+d23bra( 54)*d11ket(  1)*(xin(   4)*yin(  17)*zin(  29)+xin(  52)*yin(  65)*zin(  77)+xin( 100)*yin( 113)*zin( 125)+xin( 148)*yin( 161)*zin( 173))
          eri_value(  479)=eri_value(  479)+d23bra( 54)*d11ket(  2)*(xin(   3)*yin(  18)*zin(  29)+xin(  51)*yin(  66)*zin(  77)+xin(  99)*yin( 114)*zin( 125)+xin( 147)*yin( 162)*zin( 173))
          eri_value(  480)=eri_value(  480)+d23bra( 54)*d11ket(  3)*(xin(   3)*yin(  17)*zin(  30)+xin(  51)*yin(  65)*zin(  78)+xin(  99)*yin( 113)*zin( 126)+xin( 147)*yin( 161)*zin( 174))
          eri_value(  481)=eri_value(  481)+d23bra( 54)*d11ket(  4)*(xin(   2)*yin(  19)*zin(  29)+xin(  50)*yin(  67)*zin(  77)+xin(  98)*yin( 115)*zin( 125)+xin( 146)*yin( 163)*zin( 173))
          eri_value(  482)=eri_value(  482)+d23bra( 54)*d11ket(  5)*(xin(   1)*yin(  20)*zin(  29)+xin(  49)*yin(  68)*zin(  77)+xin(  97)*yin( 116)*zin( 125)+xin( 145)*yin( 164)*zin( 173))
          eri_value(  483)=eri_value(  483)+d23bra( 54)*d11ket(  6)*(xin(   1)*yin(  19)*zin(  30)+xin(  49)*yin(  67)*zin(  78)+xin(  97)*yin( 115)*zin( 126)+xin( 145)*yin( 163)*zin( 174))
          eri_value(  484)=eri_value(  484)+d23bra( 54)*d11ket(  7)*(xin(   2)*yin(  17)*zin(  31)+xin(  50)*yin(  65)*zin(  79)+xin(  98)*yin( 113)*zin( 127)+xin( 146)*yin( 161)*zin( 175))
          eri_value(  485)=eri_value(  485)+d23bra( 54)*d11ket(  8)*(xin(   1)*yin(  18)*zin(  31)+xin(  49)*yin(  66)*zin(  79)+xin(  97)*yin( 114)*zin( 127)+xin( 145)*yin( 162)*zin( 175))
          eri_value(  486)=eri_value(  486)+d23bra( 54)*d11ket(  9)*(xin(   1)*yin(  17)*zin(  32)+xin(  49)*yin(  65)*zin(  80)+xin(  97)*yin( 113)*zin( 128)+xin( 145)*yin( 161)*zin( 176))
          eri_value(  487)=eri_value(  487)+d23bra( 55)*d11ket(  1)*(xin(  24)*yin(  13)*zin(  13)+xin(  72)*yin(  61)*zin(  61)+xin( 120)*yin( 109)*zin( 109)+xin( 168)*yin( 157)*zin( 157))
          eri_value(  488)=eri_value(  488)+d23bra( 55)*d11ket(  2)*(xin(  23)*yin(  14)*zin(  13)+xin(  71)*yin(  62)*zin(  61)+xin( 119)*yin( 110)*zin( 109)+xin( 167)*yin( 158)*zin( 157))
          eri_value(  489)=eri_value(  489)+d23bra( 55)*d11ket(  3)*(xin(  23)*yin(  13)*zin(  14)+xin(  71)*yin(  61)*zin(  62)+xin( 119)*yin( 109)*zin( 110)+xin( 167)*yin( 157)*zin( 158))
          eri_value(  490)=eri_value(  490)+d23bra( 55)*d11ket(  4)*(xin(  22)*yin(  15)*zin(  13)+xin(  70)*yin(  63)*zin(  61)+xin( 118)*yin( 111)*zin( 109)+xin( 166)*yin( 159)*zin( 157))
          eri_value(  491)=eri_value(  491)+d23bra( 55)*d11ket(  5)*(xin(  21)*yin(  16)*zin(  13)+xin(  69)*yin(  64)*zin(  61)+xin( 117)*yin( 112)*zin( 109)+xin( 165)*yin( 160)*zin( 157))
          eri_value(  492)=eri_value(  492)+d23bra( 55)*d11ket(  6)*(xin(  21)*yin(  15)*zin(  14)+xin(  69)*yin(  63)*zin(  62)+xin( 117)*yin( 111)*zin( 110)+xin( 165)*yin( 159)*zin( 158))
          eri_value(  493)=eri_value(  493)+d23bra( 55)*d11ket(  7)*(xin(  22)*yin(  13)*zin(  15)+xin(  70)*yin(  61)*zin(  63)+xin( 118)*yin( 109)*zin( 111)+xin( 166)*yin( 157)*zin( 159))
          eri_value(  494)=eri_value(  494)+d23bra( 55)*d11ket(  8)*(xin(  21)*yin(  14)*zin(  15)+xin(  69)*yin(  62)*zin(  63)+xin( 117)*yin( 110)*zin( 111)+xin( 165)*yin( 158)*zin( 159))
          eri_value(  495)=eri_value(  495)+d23bra( 55)*d11ket(  9)*(xin(  21)*yin(  13)*zin(  16)+xin(  69)*yin(  61)*zin(  64)+xin( 117)*yin( 109)*zin( 112)+xin( 165)*yin( 157)*zin( 160))
          eri_value(  496)=eri_value(  496)+d23bra( 56)*d11ket(  1)*(xin(  16)*yin(  21)*zin(  13)+xin(  64)*yin(  69)*zin(  61)+xin( 112)*yin( 117)*zin( 109)+xin( 160)*yin( 165)*zin( 157))
          eri_value(  497)=eri_value(  497)+d23bra( 56)*d11ket(  2)*(xin(  15)*yin(  22)*zin(  13)+xin(  63)*yin(  70)*zin(  61)+xin( 111)*yin( 118)*zin( 109)+xin( 159)*yin( 166)*zin( 157))
          eri_value(  498)=eri_value(  498)+d23bra( 56)*d11ket(  3)*(xin(  15)*yin(  21)*zin(  14)+xin(  63)*yin(  69)*zin(  62)+xin( 111)*yin( 117)*zin( 110)+xin( 159)*yin( 165)*zin( 158))
          eri_value(  499)=eri_value(  499)+d23bra( 56)*d11ket(  4)*(xin(  14)*yin(  23)*zin(  13)+xin(  62)*yin(  71)*zin(  61)+xin( 110)*yin( 119)*zin( 109)+xin( 158)*yin( 167)*zin( 157))
          eri_value(  500)=eri_value(  500)+d23bra( 56)*d11ket(  5)*(xin(  13)*yin(  24)*zin(  13)+xin(  61)*yin(  72)*zin(  61)+xin( 109)*yin( 120)*zin( 109)+xin( 157)*yin( 168)*zin( 157))
          eri_value(  501)=eri_value(  501)+d23bra( 56)*d11ket(  6)*(xin(  13)*yin(  23)*zin(  14)+xin(  61)*yin(  71)*zin(  62)+xin( 109)*yin( 119)*zin( 110)+xin( 157)*yin( 167)*zin( 158))
          eri_value(  502)=eri_value(  502)+d23bra( 56)*d11ket(  7)*(xin(  14)*yin(  21)*zin(  15)+xin(  62)*yin(  69)*zin(  63)+xin( 110)*yin( 117)*zin( 111)+xin( 158)*yin( 165)*zin( 159))
          eri_value(  503)=eri_value(  503)+d23bra( 56)*d11ket(  8)*(xin(  13)*yin(  22)*zin(  15)+xin(  61)*yin(  70)*zin(  63)+xin( 109)*yin( 118)*zin( 111)+xin( 157)*yin( 166)*zin( 159))
          eri_value(  504)=eri_value(  504)+d23bra( 56)*d11ket(  9)*(xin(  13)*yin(  21)*zin(  16)+xin(  61)*yin(  69)*zin(  64)+xin( 109)*yin( 117)*zin( 112)+xin( 157)*yin( 165)*zin( 160))
          eri_value(  505)=eri_value(  505)+d23bra( 57)*d11ket(  1)*(xin(  16)*yin(  13)*zin(  21)+xin(  64)*yin(  61)*zin(  69)+xin( 112)*yin( 109)*zin( 117)+xin( 160)*yin( 157)*zin( 165))
          eri_value(  506)=eri_value(  506)+d23bra( 57)*d11ket(  2)*(xin(  15)*yin(  14)*zin(  21)+xin(  63)*yin(  62)*zin(  69)+xin( 111)*yin( 110)*zin( 117)+xin( 159)*yin( 158)*zin( 165))
          eri_value(  507)=eri_value(  507)+d23bra( 57)*d11ket(  3)*(xin(  15)*yin(  13)*zin(  22)+xin(  63)*yin(  61)*zin(  70)+xin( 111)*yin( 109)*zin( 118)+xin( 159)*yin( 157)*zin( 166))
          eri_value(  508)=eri_value(  508)+d23bra( 57)*d11ket(  4)*(xin(  14)*yin(  15)*zin(  21)+xin(  62)*yin(  63)*zin(  69)+xin( 110)*yin( 111)*zin( 117)+xin( 158)*yin( 159)*zin( 165))
          eri_value(  509)=eri_value(  509)+d23bra( 57)*d11ket(  5)*(xin(  13)*yin(  16)*zin(  21)+xin(  61)*yin(  64)*zin(  69)+xin( 109)*yin( 112)*zin( 117)+xin( 157)*yin( 160)*zin( 165))
          eri_value(  510)=eri_value(  510)+d23bra( 57)*d11ket(  6)*(xin(  13)*yin(  15)*zin(  22)+xin(  61)*yin(  63)*zin(  70)+xin( 109)*yin( 111)*zin( 118)+xin( 157)*yin( 159)*zin( 166))
          eri_value(  511)=eri_value(  511)+d23bra( 57)*d11ket(  7)*(xin(  14)*yin(  13)*zin(  23)+xin(  62)*yin(  61)*zin(  71)+xin( 110)*yin( 109)*zin( 119)+xin( 158)*yin( 157)*zin( 167))
          eri_value(  512)=eri_value(  512)+d23bra( 57)*d11ket(  8)*(xin(  13)*yin(  14)*zin(  23)+xin(  61)*yin(  62)*zin(  71)+xin( 109)*yin( 110)*zin( 119)+xin( 157)*yin( 158)*zin( 167))
          eri_value(  513)=eri_value(  513)+d23bra( 57)*d11ket(  9)*(xin(  13)*yin(  13)*zin(  24)+xin(  61)*yin(  61)*zin(  72)+xin( 109)*yin( 109)*zin( 120)+xin( 157)*yin( 157)*zin( 168))
          eri_value(  514)=eri_value(  514)+d23bra( 58)*d11ket(  1)*(xin(  20)*yin(  17)*zin(  13)+xin(  68)*yin(  65)*zin(  61)+xin( 116)*yin( 113)*zin( 109)+xin( 164)*yin( 161)*zin( 157))
          eri_value(  515)=eri_value(  515)+d23bra( 58)*d11ket(  2)*(xin(  19)*yin(  18)*zin(  13)+xin(  67)*yin(  66)*zin(  61)+xin( 115)*yin( 114)*zin( 109)+xin( 163)*yin( 162)*zin( 157))
          eri_value(  516)=eri_value(  516)+d23bra( 58)*d11ket(  3)*(xin(  19)*yin(  17)*zin(  14)+xin(  67)*yin(  65)*zin(  62)+xin( 115)*yin( 113)*zin( 110)+xin( 163)*yin( 161)*zin( 158))
          eri_value(  517)=eri_value(  517)+d23bra( 58)*d11ket(  4)*(xin(  18)*yin(  19)*zin(  13)+xin(  66)*yin(  67)*zin(  61)+xin( 114)*yin( 115)*zin( 109)+xin( 162)*yin( 163)*zin( 157))
          eri_value(  518)=eri_value(  518)+d23bra( 58)*d11ket(  5)*(xin(  17)*yin(  20)*zin(  13)+xin(  65)*yin(  68)*zin(  61)+xin( 113)*yin( 116)*zin( 109)+xin( 161)*yin( 164)*zin( 157))
          eri_value(  519)=eri_value(  519)+d23bra( 58)*d11ket(  6)*(xin(  17)*yin(  19)*zin(  14)+xin(  65)*yin(  67)*zin(  62)+xin( 113)*yin( 115)*zin( 110)+xin( 161)*yin( 163)*zin( 158))
          eri_value(  520)=eri_value(  520)+d23bra( 58)*d11ket(  7)*(xin(  18)*yin(  17)*zin(  15)+xin(  66)*yin(  65)*zin(  63)+xin( 114)*yin( 113)*zin( 111)+xin( 162)*yin( 161)*zin( 159))
          eri_value(  521)=eri_value(  521)+d23bra( 58)*d11ket(  8)*(xin(  17)*yin(  18)*zin(  15)+xin(  65)*yin(  66)*zin(  63)+xin( 113)*yin( 114)*zin( 111)+xin( 161)*yin( 162)*zin( 159))
          eri_value(  522)=eri_value(  522)+d23bra( 58)*d11ket(  9)*(xin(  17)*yin(  17)*zin(  16)+xin(  65)*yin(  65)*zin(  64)+xin( 113)*yin( 113)*zin( 112)+xin( 161)*yin( 161)*zin( 160))
          eri_value(  523)=eri_value(  523)+d23bra( 59)*d11ket(  1)*(xin(  20)*yin(  13)*zin(  17)+xin(  68)*yin(  61)*zin(  65)+xin( 116)*yin( 109)*zin( 113)+xin( 164)*yin( 157)*zin( 161))
          eri_value(  524)=eri_value(  524)+d23bra( 59)*d11ket(  2)*(xin(  19)*yin(  14)*zin(  17)+xin(  67)*yin(  62)*zin(  65)+xin( 115)*yin( 110)*zin( 113)+xin( 163)*yin( 158)*zin( 161))
          eri_value(  525)=eri_value(  525)+d23bra( 59)*d11ket(  3)*(xin(  19)*yin(  13)*zin(  18)+xin(  67)*yin(  61)*zin(  66)+xin( 115)*yin( 109)*zin( 114)+xin( 163)*yin( 157)*zin( 162))
          eri_value(  526)=eri_value(  526)+d23bra( 59)*d11ket(  4)*(xin(  18)*yin(  15)*zin(  17)+xin(  66)*yin(  63)*zin(  65)+xin( 114)*yin( 111)*zin( 113)+xin( 162)*yin( 159)*zin( 161))
          eri_value(  527)=eri_value(  527)+d23bra( 59)*d11ket(  5)*(xin(  17)*yin(  16)*zin(  17)+xin(  65)*yin(  64)*zin(  65)+xin( 113)*yin( 112)*zin( 113)+xin( 161)*yin( 160)*zin( 161))
          eri_value(  528)=eri_value(  528)+d23bra( 59)*d11ket(  6)*(xin(  17)*yin(  15)*zin(  18)+xin(  65)*yin(  63)*zin(  66)+xin( 113)*yin( 111)*zin( 114)+xin( 161)*yin( 159)*zin( 162))
          eri_value(  529)=eri_value(  529)+d23bra( 59)*d11ket(  7)*(xin(  18)*yin(  13)*zin(  19)+xin(  66)*yin(  61)*zin(  67)+xin( 114)*yin( 109)*zin( 115)+xin( 162)*yin( 157)*zin( 163))
          eri_value(  530)=eri_value(  530)+d23bra( 59)*d11ket(  8)*(xin(  17)*yin(  14)*zin(  19)+xin(  65)*yin(  62)*zin(  67)+xin( 113)*yin( 110)*zin( 115)+xin( 161)*yin( 158)*zin( 163))
          eri_value(  531)=eri_value(  531)+d23bra( 59)*d11ket(  9)*(xin(  17)*yin(  13)*zin(  20)+xin(  65)*yin(  61)*zin(  68)+xin( 113)*yin( 109)*zin( 116)+xin( 161)*yin( 157)*zin( 164))
          eri_value(  532)=eri_value(  532)+d23bra( 60)*d11ket(  1)*(xin(  16)*yin(  17)*zin(  17)+xin(  64)*yin(  65)*zin(  65)+xin( 112)*yin( 113)*zin( 113)+xin( 160)*yin( 161)*zin( 161))
          eri_value(  533)=eri_value(  533)+d23bra( 60)*d11ket(  2)*(xin(  15)*yin(  18)*zin(  17)+xin(  63)*yin(  66)*zin(  65)+xin( 111)*yin( 114)*zin( 113)+xin( 159)*yin( 162)*zin( 161))
          eri_value(  534)=eri_value(  534)+d23bra( 60)*d11ket(  3)*(xin(  15)*yin(  17)*zin(  18)+xin(  63)*yin(  65)*zin(  66)+xin( 111)*yin( 113)*zin( 114)+xin( 159)*yin( 161)*zin( 162))
          eri_value(  535)=eri_value(  535)+d23bra( 60)*d11ket(  4)*(xin(  14)*yin(  19)*zin(  17)+xin(  62)*yin(  67)*zin(  65)+xin( 110)*yin( 115)*zin( 113)+xin( 158)*yin( 163)*zin( 161))
          eri_value(  536)=eri_value(  536)+d23bra( 60)*d11ket(  5)*(xin(  13)*yin(  20)*zin(  17)+xin(  61)*yin(  68)*zin(  65)+xin( 109)*yin( 116)*zin( 113)+xin( 157)*yin( 164)*zin( 161))
          eri_value(  537)=eri_value(  537)+d23bra( 60)*d11ket(  6)*(xin(  13)*yin(  19)*zin(  18)+xin(  61)*yin(  67)*zin(  66)+xin( 109)*yin( 115)*zin( 114)+xin( 157)*yin( 163)*zin( 162))
          eri_value(  538)=eri_value(  538)+d23bra( 60)*d11ket(  7)*(xin(  14)*yin(  17)*zin(  19)+xin(  62)*yin(  65)*zin(  67)+xin( 110)*yin( 113)*zin( 115)+xin( 158)*yin( 161)*zin( 163))
          eri_value(  539)=eri_value(  539)+d23bra( 60)*d11ket(  8)*(xin(  13)*yin(  18)*zin(  19)+xin(  61)*yin(  66)*zin(  67)+xin( 109)*yin( 114)*zin( 115)+xin( 157)*yin( 162)*zin( 163))
          eri_value(  540)=eri_value(  540)+d23bra( 60)*d11ket(  9)*(xin(  13)*yin(  17)*zin(  20)+xin(  61)*yin(  65)*zin(  68)+xin( 109)*yin( 113)*zin( 116)+xin( 157)*yin( 161)*zin( 164))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxl = 3
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*54 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

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

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n11ket)
                              deallocate (xint11ket)

                              end subroutine int3211
                              end submodule
