! The total angular momentum of this class is:           8
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3221_impl
contains
  module subroutine int3221(df_pair, pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n12ket(:)
    real(dp), allocatable :: xint23bra(:), xint12ket(:)
    integer(kind=int64) :: ndfbra, npdket
    real(dp) :: scutdfbra, scutpdket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2
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
    real(dp) :: xin(360), yin(360), zin(360)
    real(dp) :: eri_value(1080)
    real(dp) :: d23bra(60), d12ket(18)
    integer(kind=int64) :: ix(10), jx(6), kx(6), lx(3)
    integer(kind=int64) :: iy(10), jy(6), ky(6), ly(3)
    integer(kind=int64) :: iz(10), jz(6), kz(6), lz(3)
    integer(kind=int64) :: in(6), in1(6), kn(4)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(18), kly(18), klz(18)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 19
    in1(3) = 37
    in1(4) = 55
    in1(5) = 61
    in1(6) = 67

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

    jx(1) = 12
    jx(2) = 0
    jx(3) = 0
    jx(4) = 6
    jx(5) = 6
    jx(6) = 0

    ix(1) = 55
    ix(2) = 1
    ix(3) = 1
    ix(4) = 37
    ix(5) = 37
    ix(6) = 19
    ix(7) = 1
    ix(8) = 19
    ix(9) = 1
    ix(10) = 19

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
    jy(2) = 12
    jy(3) = 0
    jy(4) = 6
    jy(5) = 0
    jy(6) = 6

    iy(1) = 1
    iy(2) = 55
    iy(3) = 1
    iy(4) = 19
    iy(5) = 1
    iy(6) = 37
    iy(7) = 37
    iy(8) = 1
    iy(9) = 19
    iy(10) = 19

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
    jz(3) = 12
    jz(4) = 0
    jz(5) = 6
    jz(6) = 6

    iz(1) = 1
    iz(2) = 1
    iz(3) = 55
    iz(4) = 1
    iz(5) = 19
    iz(6) = 1
    iz(7) = 19
    iz(8) = 37
    iz(9) = 37
    iz(10) = 19

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 67
    ijx(2) = 55
    ijx(3) = 55
    ijx(4) = 61
    ijx(5) = 61
    ijx(6) = 55
    ijx(7) = 13
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 7
    ijx(11) = 7
    ijx(12) = 1
    ijx(13) = 13
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 7
    ijx(17) = 7
    ijx(18) = 1
    ijx(19) = 49
    ijx(20) = 37
    ijx(21) = 37
    ijx(22) = 43
    ijx(23) = 43
    ijx(24) = 37
    ijx(25) = 49
    ijx(26) = 37
    ijx(27) = 37
    ijx(28) = 43
    ijx(29) = 43
    ijx(30) = 37
    ijx(31) = 31
    ijx(32) = 19
    ijx(33) = 19
    ijx(34) = 25
    ijx(35) = 25
    ijx(36) = 19
    ijx(37) = 13
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 7
    ijx(41) = 7
    ijx(42) = 1
    ijx(43) = 31
    ijx(44) = 19
    ijx(45) = 19
    ijx(46) = 25
    ijx(47) = 25
    ijx(48) = 19
    ijx(49) = 13
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 7
    ijx(53) = 7
    ijx(54) = 1
    ijx(55) = 31
    ijx(56) = 19
    ijx(57) = 19
    ijx(58) = 25
    ijx(59) = 25
    ijx(60) = 19

    ijy(1) = 1
    ijy(2) = 13
    ijy(3) = 1
    ijy(4) = 7
    ijy(5) = 1
    ijy(6) = 7
    ijy(7) = 55
    ijy(8) = 67
    ijy(9) = 55
    ijy(10) = 61
    ijy(11) = 55
    ijy(12) = 61
    ijy(13) = 1
    ijy(14) = 13
    ijy(15) = 1
    ijy(16) = 7
    ijy(17) = 1
    ijy(18) = 7
    ijy(19) = 19
    ijy(20) = 31
    ijy(21) = 19
    ijy(22) = 25
    ijy(23) = 19
    ijy(24) = 25
    ijy(25) = 1
    ijy(26) = 13
    ijy(27) = 1
    ijy(28) = 7
    ijy(29) = 1
    ijy(30) = 7
    ijy(31) = 37
    ijy(32) = 49
    ijy(33) = 37
    ijy(34) = 43
    ijy(35) = 37
    ijy(36) = 43
    ijy(37) = 37
    ijy(38) = 49
    ijy(39) = 37
    ijy(40) = 43
    ijy(41) = 37
    ijy(42) = 43
    ijy(43) = 1
    ijy(44) = 13
    ijy(45) = 1
    ijy(46) = 7
    ijy(47) = 1
    ijy(48) = 7
    ijy(49) = 19
    ijy(50) = 31
    ijy(51) = 19
    ijy(52) = 25
    ijy(53) = 19
    ijy(54) = 25
    ijy(55) = 19
    ijy(56) = 31
    ijy(57) = 19
    ijy(58) = 25
    ijy(59) = 19
    ijy(60) = 25

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 13
    ijz(4) = 1
    ijz(5) = 7
    ijz(6) = 7
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 13
    ijz(10) = 1
    ijz(11) = 7
    ijz(12) = 7
    ijz(13) = 55
    ijz(14) = 55
    ijz(15) = 67
    ijz(16) = 55
    ijz(17) = 61
    ijz(18) = 61
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 13
    ijz(22) = 1
    ijz(23) = 7
    ijz(24) = 7
    ijz(25) = 19
    ijz(26) = 19
    ijz(27) = 31
    ijz(28) = 19
    ijz(29) = 25
    ijz(30) = 25
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 13
    ijz(34) = 1
    ijz(35) = 7
    ijz(36) = 7
    ijz(37) = 19
    ijz(38) = 19
    ijz(39) = 31
    ijz(40) = 19
    ijz(41) = 25
    ijz(42) = 25
    ijz(43) = 37
    ijz(44) = 37
    ijz(45) = 49
    ijz(46) = 37
    ijz(47) = 43
    ijz(48) = 43
    ijz(49) = 37
    ijz(50) = 37
    ijz(51) = 49
    ijz(52) = 37
    ijz(53) = 43
    ijz(54) = 43
    ijz(55) = 19
    ijz(56) = 19
    ijz(57) = 31
    ijz(58) = 19
    ijz(59) = 25
    ijz(60) = 25

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

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n12ket(res%n_p_shl*res%n_d_shl))
    allocate (xint12ket(res%n_p_shl*res%n_d_shl))

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

    if ((ndfbra*npdket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*npdket
    ntile = int(ndfbra*npdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*npdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint12ket, n12ket, df_pair, pd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d12ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint12ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n12ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
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

                                      ! i2 = in(2) =   19
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(19) = xc00
                                      yin(19) = yc00
                                      zin(19) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   21
                                      ! i2 =   19

                                      xin(21) = xcp00*xin(19) + cp10
                                      yin(21) = ycp00*yin(19) + cp10
                                      zin(21) = zcp00*zin(19) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   19

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   37
                                      ! i3 =    1
                                      ! i4 =   19

                                      xin(37) = c10*xin(1) + xc00*xin(19)
                                      yin(37) = c10*yin(1) + yc00*yin(19)
                                      zin(37) = c10*zin(1) + zc00*zin(19)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   39
                                      ! i5 =   37
                                      ! i4 =   19

                                      xin(39) = xcp00*xin(37) + cp10*xin(19)
                                      yin(39) = ycp00*yin(37) + cp10*yin(19)
                                      zin(39) = zcp00*zin(37) + cp10*zin(19)

                                      ! ------------------

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   37

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   55
                                      ! i3 =   19
                                      ! i4 =   37

                                      xin(55) = c10*xin(19) + xc00*xin(37)
                                      yin(55) = c10*yin(19) + yc00*yin(37)
                                      zin(55) = c10*zin(19) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   57
                                      ! i5 =   55
                                      ! i4 =   37

                                      xin(57) = xcp00*xin(55) + cp10*xin(37)
                                      yin(57) = ycp00*yin(55) + cp10*yin(37)
                                      zin(57) = zcp00*zin(55) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   55

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   37
                                      ! i4 =   55

                                      xin(61) = c10*xin(37) + xc00*xin(55)
                                      yin(61) = c10*yin(37) + yc00*yin(55)
                                      zin(61) = c10*zin(37) + zc00*zin(55)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   63
                                      ! i5 =   61
                                      ! i4 =   55

                                      xin(63) = xcp00*xin(61) + cp10*xin(55)
                                      yin(63) = ycp00*yin(61) + cp10*yin(55)
                                      zin(63) = zcp00*zin(61) + cp10*zin(55)

                                      ! ------------------

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   61

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   67
                                      ! i3 =   55
                                      ! i4 =   61

                                      xin(67) = c10*xin(55) + xc00*xin(61)
                                      yin(67) = c10*yin(55) + yc00*yin(61)
                                      zin(67) = c10*zin(55) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   69
                                      ! i5 =   67
                                      ! i4 =   61

                                      xin(69) = xcp00*xin(67) + cp10*xin(61)
                                      yin(69) = ycp00*yin(67) + cp10*yin(61)
                                      zin(69) = zcp00*zin(67) + cp10*zin(61)

                                      ! ------------------

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   67

                                      ! n =    6

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

                                      ! i3 = i2+kn(n+1) =   23

                                      xin(23) = xc00*xin(5) + c01*xin(3)
                                      yin(23) = yc00*yin(5) + c01*yin(3)
                                      zin(23) = zc00*zin(5) + c01*zin(3)

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

                                      ! i3 = i2+kn(n+1) =   24

                                      xin(24) = xc00*xin(6) + c01*xin(5)
                                      yin(24) = yc00*yin(6) + c01*yin(5)
                                      zin(24) = zc00*zin(6) + c01*zin(5)

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
                                      ! i4 = i2 =   19

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   37

                                      xin(41) = c10*xin(5) + xc00*xin(23) + c01*xin(21)
                                      yin(41) = c10*yin(5) + yc00*yin(23) + c01*yin(21)
                                      zin(41) = c10*zin(5) + zc00*zin(23) + c01*zin(21)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   37

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   55

                                      xin(59) = c10*xin(23) + xc00*xin(41) + c01*xin(39)
                                      yin(59) = c10*yin(23) + yc00*yin(41) + c01*yin(39)
                                      zin(59) = c10*zin(23) + zc00*zin(41) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   55

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   61

                                      xin(65) = c10*xin(41) + xc00*xin(59) + c01*xin(57)
                                      yin(65) = c10*yin(41) + yc00*yin(59) + c01*yin(57)
                                      zin(65) = c10*zin(41) + zc00*zin(59) + c01*zin(57)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   61

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   67

                                      xin(71) = c10*xin(59) + xc00*xin(65) + c01*xin(63)
                                      yin(71) = c10*yin(59) + yc00*yin(65) + c01*yin(63)
                                      zin(71) = c10*zin(59) + zc00*zin(65) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   67

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   19

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   37

                                      xin(42) = c10*xin(6) + xc00*xin(24) + c01*xin(23)
                                      yin(42) = c10*yin(6) + yc00*yin(24) + c01*yin(23)
                                      zin(42) = c10*zin(6) + zc00*zin(24) + c01*zin(23)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   37

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   55

                                      xin(60) = c10*xin(24) + xc00*xin(42) + c01*xin(41)
                                      yin(60) = c10*yin(24) + yc00*yin(42) + c01*yin(41)
                                      zin(60) = c10*zin(24) + zc00*zin(42) + c01*zin(41)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   55

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   61

                                      xin(66) = c10*xin(42) + xc00*xin(60) + c01*xin(59)
                                      yin(66) = c10*yin(42) + yc00*yin(60) + c01*yin(59)
                                      zin(66) = c10*zin(42) + zc00*zin(60) + c01*zin(59)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   61

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   67

                                      xin(72) = c10*xin(60) + xc00*xin(66) + c01*xin(65)
                                      yin(72) = c10*yin(60) + yc00*yin(66) + c01*yin(65)
                                      zin(72) = c10*zin(60) + zc00*zin(66) + c01*zin(65)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   67

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   67

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   67

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   61

                                      xin(67) = xin(67) + dxij*xin(61)
                                      yin(67) = yin(67) + dyij*yin(61)
                                      zin(67) = zin(67) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   55

                                      xin(61) = xin(61) + dxij*xin(55)
                                      yin(61) = yin(61) + dyij*yin(55)
                                      zin(61) = zin(61) + dzij*zin(55)

                                      ! i3 = i4 =   55
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   67

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   61

                                      xin(67) = xin(67) + dxij*xin(61)
                                      yin(67) = yin(67) + dyij*yin(61)
                                      zin(67) = zin(67) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(19) + dxij*xin(1)
                                      yin(7) = yin(19) + dyij*yin(1)
                                      zin(7) = zin(19) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! ni =    2

                                      xin(25) = xin(37) + dxij*xin(19)
                                      yin(25) = yin(37) + dyij*yin(19)
                                      zin(25) = zin(37) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    3

                                      xin(43) = xin(55) + dxij*xin(37)
                                      yin(43) = yin(55) + dyij*yin(37)
                                      zin(43) = zin(55) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    2

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(25) + dxij*xin(7)
                                      yin(13) = yin(25) + dyij*yin(7)
                                      zin(13) = zin(25) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    2

                                      xin(31) = xin(43) + dxij*xin(25)
                                      yin(31) = yin(43) + dyij*yin(25)
                                      zin(31) = zin(43) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! ni =    3

                                      xin(49) = xin(61) + dxij*xin(43)
                                      yin(49) = yin(61) + dyij*yin(43)
                                      zin(49) = zin(61) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(69) = xin(69) + dxij*xin(63)
                                      yin(69) = yin(69) + dyij*yin(63)
                                      zin(69) = zin(69) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   57

                                      xin(63) = xin(63) + dxij*xin(57)
                                      yin(63) = yin(63) + dyij*yin(57)
                                      zin(63) = zin(63) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(69) = xin(69) + dxij*xin(63)
                                      yin(69) = yin(69) + dyij*yin(63)
                                      zin(69) = zin(69) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(21) + dxij*xin(3)
                                      yin(9) = yin(21) + dyij*yin(3)
                                      zin(9) = zin(21) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   27

                                      ! ni =    2

                                      xin(27) = xin(39) + dxij*xin(21)
                                      yin(27) = yin(39) + dyij*yin(21)
                                      zin(27) = zin(39) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    3

                                      xin(45) = xin(57) + dxij*xin(39)
                                      yin(45) = yin(57) + dyij*yin(39)
                                      zin(45) = zin(57) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    2

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    3

                                      xin(15) = xin(27) + dxij*xin(9)
                                      yin(15) = yin(27) + dyij*yin(9)
                                      zin(15) = zin(27) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    2

                                      xin(33) = xin(45) + dxij*xin(27)
                                      yin(33) = yin(45) + dyij*yin(27)
                                      zin(33) = zin(45) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   51

                                      ! ni =    3

                                      xin(51) = xin(63) + dxij*xin(45)
                                      yin(51) = yin(63) + dyij*yin(45)
                                      zin(51) = zin(63) + dzij*zin(45)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   65

                                      xin(71) = xin(71) + dxij*xin(65)
                                      yin(71) = yin(71) + dyij*yin(65)
                                      zin(71) = zin(71) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   59

                                      xin(65) = xin(65) + dxij*xin(59)
                                      yin(65) = yin(65) + dyij*yin(59)
                                      zin(65) = zin(65) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   65

                                      xin(71) = xin(71) + dxij*xin(65)
                                      yin(71) = yin(71) + dyij*yin(65)
                                      zin(71) = zin(71) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(23) + dxij*xin(5)
                                      yin(11) = yin(23) + dyij*yin(5)
                                      zin(11) = zin(23) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    2

                                      xin(29) = xin(41) + dxij*xin(23)
                                      yin(29) = yin(41) + dyij*yin(23)
                                      zin(29) = zin(41) + dzij*zin(23)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    3

                                      xin(47) = xin(59) + dxij*xin(41)
                                      yin(47) = yin(59) + dyij*yin(41)
                                      zin(47) = zin(59) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    3

                                      xin(17) = xin(29) + dxij*xin(11)
                                      yin(17) = yin(29) + dyij*yin(11)
                                      zin(17) = zin(29) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    2

                                      xin(35) = xin(47) + dxij*xin(29)
                                      yin(35) = yin(47) + dyij*yin(29)
                                      zin(35) = zin(47) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    3

                                      xin(53) = xin(65) + dxij*xin(47)
                                      yin(53) = yin(65) + dyij*yin(47)
                                      zin(53) = zin(65) + dzij*zin(47)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   23

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   66

                                      xin(72) = xin(72) + dxij*xin(66)
                                      yin(72) = yin(72) + dyij*yin(66)
                                      zin(72) = zin(72) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   60

                                      xin(66) = xin(66) + dxij*xin(60)
                                      yin(66) = yin(66) + dyij*yin(60)
                                      zin(66) = zin(66) + dzij*zin(60)

                                      ! i3 = i4 =   60
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   66

                                      xin(72) = xin(72) + dxij*xin(66)
                                      yin(72) = yin(72) + dyij*yin(66)
                                      zin(72) = zin(72) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   12

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(24) + dxij*xin(6)
                                      yin(12) = yin(24) + dyij*yin(6)
                                      zin(12) = zin(24) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    2

                                      xin(30) = xin(42) + dxij*xin(24)
                                      yin(30) = yin(42) + dyij*yin(24)
                                      zin(30) = zin(42) + dzij*zin(24)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    3

                                      xin(48) = xin(60) + dxij*xin(42)
                                      yin(48) = yin(60) + dyij*yin(42)
                                      zin(48) = zin(60) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   18

                                      ! nj =    2

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(30) + dxij*xin(12)
                                      yin(18) = yin(30) + dyij*yin(12)
                                      zin(18) = zin(30) + dzij*zin(12)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    2

                                      xin(36) = xin(48) + dxij*xin(30)
                                      yin(36) = yin(48) + dyij*yin(30)
                                      zin(36) = zin(48) + dzij*zin(30)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    3

                                      xin(54) = xin(66) + dxij*xin(48)
                                      yin(54) = yin(66) + dyij*yin(48)
                                      zin(54) = zin(66) + dzij*zin(48)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   24

                                      ! nj =    3

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   19

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   72

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

                                      ! i1 = in(1) =   73

                                      xin(73) = 1.0_dp
                                      yin(73) = 1.0_dp
                                      zin(73) = f00

                                      ! i2 = in(2) =   91
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(91) = xc00
                                      yin(91) = yc00
                                      zin(91) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   75

                                      xin(75) = xcp00
                                      yin(75) = ycp00
                                      zin(75) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   93
                                      ! i2 =   91

                                      xin(93) = xcp00*xin(91) + cp10
                                      yin(93) = ycp00*yin(91) + cp10
                                      zin(93) = zcp00*zin(91) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   91

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  109
                                      ! i3 =   73
                                      ! i4 =   91

                                      xin(109) = c10*xin(73) + xc00*xin(91)
                                      yin(109) = c10*yin(73) + yc00*yin(91)
                                      zin(109) = c10*zin(73) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  111
                                      ! i5 =  109
                                      ! i4 =   91

                                      xin(111) = xcp00*xin(109) + cp10*xin(91)
                                      yin(111) = ycp00*yin(109) + cp10*yin(91)
                                      zin(111) = zcp00*zin(109) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  109

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  127
                                      ! i3 =   91
                                      ! i4 =  109

                                      xin(127) = c10*xin(91) + xc00*xin(109)
                                      yin(127) = c10*yin(91) + yc00*yin(109)
                                      zin(127) = c10*zin(91) + zc00*zin(109)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  129
                                      ! i5 =  127
                                      ! i4 =  109

                                      xin(129) = xcp00*xin(127) + cp10*xin(109)
                                      yin(129) = ycp00*yin(127) + cp10*yin(109)
                                      zin(129) = zcp00*zin(127) + cp10*zin(109)

                                      ! ------------------

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  127

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  133
                                      ! i3 =  109
                                      ! i4 =  127

                                      xin(133) = c10*xin(109) + xc00*xin(127)
                                      yin(133) = c10*yin(109) + yc00*yin(127)
                                      zin(133) = c10*zin(109) + zc00*zin(127)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  135
                                      ! i5 =  133
                                      ! i4 =  127

                                      xin(135) = xcp00*xin(133) + cp10*xin(127)
                                      yin(135) = ycp00*yin(133) + cp10*yin(127)
                                      zin(135) = zcp00*zin(133) + cp10*zin(127)

                                      ! ------------------

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  133

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  139
                                      ! i3 =  127
                                      ! i4 =  133

                                      xin(139) = c10*xin(127) + xc00*xin(133)
                                      yin(139) = c10*yin(127) + yc00*yin(133)
                                      zin(139) = c10*zin(127) + zc00*zin(133)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  141
                                      ! i5 =  139
                                      ! i4 =  133

                                      xin(141) = xcp00*xin(139) + cp10*xin(133)
                                      yin(141) = ycp00*yin(139) + cp10*yin(133)
                                      zin(141) = zcp00*zin(139) + cp10*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  139

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   73
                                      ! i4 = i1+k2 =   75

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   77
                                      ! i3 =   73
                                      ! i4 =   75

                                      xin(77) = cp01*xin(73) + xcp00*xin(75)
                                      yin(77) = cp01*yin(73) + ycp00*yin(75)
                                      zin(77) = cp01*zin(73) + zcp00*zin(75)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   95

                                      xin(95) = xc00*xin(77) + c01*xin(75)
                                      yin(95) = yc00*yin(77) + c01*yin(75)
                                      zin(95) = zc00*zin(77) + c01*zin(75)

                                      ! ------------------

                                      ! i3 = i4 =   75
                                      ! i4 = i5 =   77

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   78
                                      ! i3 =   75
                                      ! i4 =   77

                                      xin(78) = cp01*xin(75) + xcp00*xin(77)
                                      yin(78) = cp01*yin(75) + ycp00*yin(77)
                                      zin(78) = cp01*zin(75) + zcp00*zin(77)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   96

                                      xin(96) = xc00*xin(78) + c01*xin(77)
                                      yin(96) = yc00*yin(78) + c01*yin(77)
                                      zin(96) = zc00*zin(78) + c01*zin(77)

                                      ! ------------------

                                      ! i3 = i4 =   77
                                      ! i4 = i5 =   78

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   91

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  109

                                      xin(113) = c10*xin(77) + xc00*xin(95) + c01*xin(93)
                                      yin(113) = c10*yin(77) + yc00*yin(95) + c01*yin(93)
                                      zin(113) = c10*zin(77) + zc00*zin(95) + c01*zin(93)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  109

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  127

                                      xin(131) = c10*xin(95) + xc00*xin(113) + c01*xin(111)
                                      yin(131) = c10*yin(95) + yc00*yin(113) + c01*yin(111)
                                      zin(131) = c10*zin(95) + zc00*zin(113) + c01*zin(111)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  127

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  133

                                      xin(137) = c10*xin(113) + xc00*xin(131) + c01*xin(129)
                                      yin(137) = c10*yin(113) + yc00*yin(131) + c01*yin(129)
                                      zin(137) = c10*zin(113) + zc00*zin(131) + c01*zin(129)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  133

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  139

                                      xin(143) = c10*xin(131) + xc00*xin(137) + c01*xin(135)
                                      yin(143) = c10*yin(131) + yc00*yin(137) + c01*yin(135)
                                      zin(143) = c10*zin(131) + zc00*zin(137) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  139

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   91

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  109

                                      xin(114) = c10*xin(78) + xc00*xin(96) + c01*xin(95)
                                      yin(114) = c10*yin(78) + yc00*yin(96) + c01*yin(95)
                                      zin(114) = c10*zin(78) + zc00*zin(96) + c01*zin(95)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  109

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  127

                                      xin(132) = c10*xin(96) + xc00*xin(114) + c01*xin(113)
                                      yin(132) = c10*yin(96) + yc00*yin(114) + c01*yin(113)
                                      zin(132) = c10*zin(96) + zc00*zin(114) + c01*zin(113)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  127

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  133

                                      xin(138) = c10*xin(114) + xc00*xin(132) + c01*xin(131)
                                      yin(138) = c10*yin(114) + yc00*yin(132) + c01*yin(131)
                                      zin(138) = c10*zin(114) + zc00*zin(132) + c01*zin(131)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  133

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  139

                                      xin(144) = c10*xin(132) + xc00*xin(138) + c01*xin(137)
                                      yin(144) = c10*yin(132) + yc00*yin(138) + c01*yin(137)
                                      zin(144) = c10*zin(132) + zc00*zin(138) + c01*zin(137)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  139

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  139

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(139) = xin(139) + dxij*xin(133)
                                      yin(139) = yin(139) + dyij*yin(133)
                                      zin(139) = zin(139) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  127

                                      xin(133) = xin(133) + dxij*xin(127)
                                      yin(133) = yin(133) + dyij*yin(127)
                                      zin(133) = zin(133) + dzij*zin(127)

                                      ! i3 = i4 =  127
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(139) = xin(139) + dxij*xin(133)
                                      yin(139) = yin(139) + dyij*yin(133)
                                      zin(139) = zin(139) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   79

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   79

                                      ! do ni = 1,    3

                                      xin(79) = xin(91) + dxij*xin(73)
                                      yin(79) = yin(91) + dyij*yin(73)
                                      zin(79) = zin(91) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! ni =    2

                                      xin(97) = xin(109) + dxij*xin(91)
                                      yin(97) = yin(109) + dyij*yin(91)
                                      zin(97) = zin(109) + dzij*zin(91)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    3

                                      xin(115) = xin(127) + dxij*xin(109)
                                      yin(115) = yin(127) + dyij*yin(109)
                                      zin(115) = zin(127) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   85

                                      ! nj =    2

                                      ! i4 = i3 =   85

                                      ! do ni = 1,    3

                                      xin(85) = xin(97) + dxij*xin(79)
                                      yin(85) = yin(97) + dyij*yin(79)
                                      zin(85) = zin(97) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    2

                                      xin(103) = xin(115) + dxij*xin(97)
                                      yin(103) = yin(115) + dyij*yin(97)
                                      zin(103) = zin(115) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! ni =    3

                                      xin(121) = xin(133) + dxij*xin(115)
                                      yin(121) = yin(133) + dyij*yin(115)
                                      zin(121) = zin(133) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   91

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(141) = xin(141) + dxij*xin(135)
                                      yin(141) = yin(141) + dyij*yin(135)
                                      zin(141) = zin(141) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  129

                                      xin(135) = xin(135) + dxij*xin(129)
                                      yin(135) = yin(135) + dyij*yin(129)
                                      zin(135) = zin(135) + dzij*zin(129)

                                      ! i3 = i4 =  129
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(141) = xin(141) + dxij*xin(135)
                                      yin(141) = yin(141) + dyij*yin(135)
                                      zin(141) = zin(141) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   81

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   81

                                      ! do ni = 1,    3

                                      xin(81) = xin(93) + dxij*xin(75)
                                      yin(81) = yin(93) + dyij*yin(75)
                                      zin(81) = zin(93) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   99

                                      ! ni =    2

                                      xin(99) = xin(111) + dxij*xin(93)
                                      yin(99) = yin(111) + dyij*yin(93)
                                      zin(99) = zin(111) + dzij*zin(93)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    3

                                      xin(117) = xin(129) + dxij*xin(111)
                                      yin(117) = yin(129) + dyij*yin(111)
                                      zin(117) = zin(129) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   87

                                      ! nj =    2

                                      ! i4 = i3 =   87

                                      ! do ni = 1,    3

                                      xin(87) = xin(99) + dxij*xin(81)
                                      yin(87) = yin(99) + dyij*yin(81)
                                      zin(87) = zin(99) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    2

                                      xin(105) = xin(117) + dxij*xin(99)
                                      yin(105) = yin(117) + dyij*yin(99)
                                      zin(105) = zin(117) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  123

                                      ! ni =    3

                                      xin(123) = xin(135) + dxij*xin(117)
                                      yin(123) = yin(135) + dyij*yin(117)
                                      zin(123) = zin(135) + dzij*zin(117)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   93

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(143) = xin(143) + dxij*xin(137)
                                      yin(143) = yin(143) + dyij*yin(137)
                                      zin(143) = zin(143) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  131

                                      xin(137) = xin(137) + dxij*xin(131)
                                      yin(137) = yin(137) + dyij*yin(131)
                                      zin(137) = zin(137) + dzij*zin(131)

                                      ! i3 = i4 =  131
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(143) = xin(143) + dxij*xin(137)
                                      yin(143) = yin(143) + dyij*yin(137)
                                      zin(143) = zin(143) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   83

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   83

                                      ! do ni = 1,    3

                                      xin(83) = xin(95) + dxij*xin(77)
                                      yin(83) = yin(95) + dyij*yin(77)
                                      zin(83) = zin(95) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  101

                                      ! ni =    2

                                      xin(101) = xin(113) + dxij*xin(95)
                                      yin(101) = yin(113) + dyij*yin(95)
                                      zin(101) = zin(113) + dzij*zin(95)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    3

                                      xin(119) = xin(131) + dxij*xin(113)
                                      yin(119) = yin(131) + dyij*yin(113)
                                      zin(119) = zin(131) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   89

                                      ! nj =    2

                                      ! i4 = i3 =   89

                                      ! do ni = 1,    3

                                      xin(89) = xin(101) + dxij*xin(83)
                                      yin(89) = yin(101) + dyij*yin(83)
                                      zin(89) = zin(101) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    2

                                      xin(107) = xin(119) + dxij*xin(101)
                                      yin(107) = yin(119) + dyij*yin(101)
                                      zin(107) = zin(119) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    3

                                      xin(125) = xin(137) + dxij*xin(119)
                                      yin(125) = yin(137) + dyij*yin(119)
                                      zin(125) = zin(137) + dzij*zin(119)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   95

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(144) = xin(144) + dxij*xin(138)
                                      yin(144) = yin(144) + dyij*yin(138)
                                      zin(144) = zin(144) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  132

                                      xin(138) = xin(138) + dxij*xin(132)
                                      yin(138) = yin(138) + dyij*yin(132)
                                      zin(138) = zin(138) + dzij*zin(132)

                                      ! i3 = i4 =  132
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(144) = xin(144) + dxij*xin(138)
                                      yin(144) = yin(144) + dyij*yin(138)
                                      zin(144) = zin(144) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   84

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   84

                                      ! do ni = 1,    3

                                      xin(84) = xin(96) + dxij*xin(78)
                                      yin(84) = yin(96) + dyij*yin(78)
                                      zin(84) = zin(96) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  102

                                      ! ni =    2

                                      xin(102) = xin(114) + dxij*xin(96)
                                      yin(102) = yin(114) + dyij*yin(96)
                                      zin(102) = zin(114) + dzij*zin(96)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    3

                                      xin(120) = xin(132) + dxij*xin(114)
                                      yin(120) = yin(132) + dyij*yin(114)
                                      zin(120) = zin(132) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   90

                                      ! nj =    2

                                      ! i4 = i3 =   90

                                      ! do ni = 1,    3

                                      xin(90) = xin(102) + dxij*xin(84)
                                      yin(90) = yin(102) + dyij*yin(84)
                                      zin(90) = zin(102) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    2

                                      xin(108) = xin(120) + dxij*xin(102)
                                      yin(108) = yin(120) + dyij*yin(102)
                                      zin(108) = zin(120) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    3

                                      xin(126) = xin(138) + dxij*xin(120)
                                      yin(126) = yin(138) + dyij*yin(120)
                                      zin(126) = zin(138) + dzij*zin(120)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   96

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =   73

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

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

                                      ! i1 = in(1) =  145

                                      xin(145) = 1.0_dp
                                      yin(145) = 1.0_dp
                                      zin(145) = f00

                                      ! i2 = in(2) =  163
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(163) = xc00
                                      yin(163) = yc00
                                      zin(163) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  147

                                      xin(147) = xcp00
                                      yin(147) = ycp00
                                      zin(147) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  165
                                      ! i2 =  163

                                      xin(165) = xcp00*xin(163) + cp10
                                      yin(165) = ycp00*yin(163) + cp10
                                      zin(165) = zcp00*zin(163) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  163

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  181
                                      ! i3 =  145
                                      ! i4 =  163

                                      xin(181) = c10*xin(145) + xc00*xin(163)
                                      yin(181) = c10*yin(145) + yc00*yin(163)
                                      zin(181) = c10*zin(145) + zc00*zin(163)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  183
                                      ! i5 =  181
                                      ! i4 =  163

                                      xin(183) = xcp00*xin(181) + cp10*xin(163)
                                      yin(183) = ycp00*yin(181) + cp10*yin(163)
                                      zin(183) = zcp00*zin(181) + cp10*zin(163)

                                      ! ------------------

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  181

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  199
                                      ! i3 =  163
                                      ! i4 =  181

                                      xin(199) = c10*xin(163) + xc00*xin(181)
                                      yin(199) = c10*yin(163) + yc00*yin(181)
                                      zin(199) = c10*zin(163) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  201
                                      ! i5 =  199
                                      ! i4 =  181

                                      xin(201) = xcp00*xin(199) + cp10*xin(181)
                                      yin(201) = ycp00*yin(199) + cp10*yin(181)
                                      zin(201) = zcp00*zin(199) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  199

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  205
                                      ! i3 =  181
                                      ! i4 =  199

                                      xin(205) = c10*xin(181) + xc00*xin(199)
                                      yin(205) = c10*yin(181) + yc00*yin(199)
                                      zin(205) = c10*zin(181) + zc00*zin(199)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  207
                                      ! i5 =  205
                                      ! i4 =  199

                                      xin(207) = xcp00*xin(205) + cp10*xin(199)
                                      yin(207) = ycp00*yin(205) + cp10*yin(199)
                                      zin(207) = zcp00*zin(205) + cp10*zin(199)

                                      ! ------------------

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  205

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  211
                                      ! i3 =  199
                                      ! i4 =  205

                                      xin(211) = c10*xin(199) + xc00*xin(205)
                                      yin(211) = c10*yin(199) + yc00*yin(205)
                                      zin(211) = c10*zin(199) + zc00*zin(205)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  213
                                      ! i5 =  211
                                      ! i4 =  205

                                      xin(213) = xcp00*xin(211) + cp10*xin(205)
                                      yin(213) = ycp00*yin(211) + cp10*yin(205)
                                      zin(213) = zcp00*zin(211) + cp10*zin(205)

                                      ! ------------------

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  211

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  147

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  149
                                      ! i3 =  145
                                      ! i4 =  147

                                      xin(149) = cp01*xin(145) + xcp00*xin(147)
                                      yin(149) = cp01*yin(145) + ycp00*yin(147)
                                      zin(149) = cp01*zin(145) + zcp00*zin(147)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  167

                                      xin(167) = xc00*xin(149) + c01*xin(147)
                                      yin(167) = yc00*yin(149) + c01*yin(147)
                                      zin(167) = zc00*zin(149) + c01*zin(147)

                                      ! ------------------

                                      ! i3 = i4 =  147
                                      ! i4 = i5 =  149

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  150
                                      ! i3 =  147
                                      ! i4 =  149

                                      xin(150) = cp01*xin(147) + xcp00*xin(149)
                                      yin(150) = cp01*yin(147) + ycp00*yin(149)
                                      zin(150) = cp01*zin(147) + zcp00*zin(149)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  168

                                      xin(168) = xc00*xin(150) + c01*xin(149)
                                      yin(168) = yc00*yin(150) + c01*yin(149)
                                      zin(168) = zc00*zin(150) + c01*zin(149)

                                      ! ------------------

                                      ! i3 = i4 =  149
                                      ! i4 = i5 =  150

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  163

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  181

                                      xin(185) = c10*xin(149) + xc00*xin(167) + c01*xin(165)
                                      yin(185) = c10*yin(149) + yc00*yin(167) + c01*yin(165)
                                      zin(185) = c10*zin(149) + zc00*zin(167) + c01*zin(165)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  181

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  199

                                      xin(203) = c10*xin(167) + xc00*xin(185) + c01*xin(183)
                                      yin(203) = c10*yin(167) + yc00*yin(185) + c01*yin(183)
                                      zin(203) = c10*zin(167) + zc00*zin(185) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  199

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  205

                                      xin(209) = c10*xin(185) + xc00*xin(203) + c01*xin(201)
                                      yin(209) = c10*yin(185) + yc00*yin(203) + c01*yin(201)
                                      zin(209) = c10*zin(185) + zc00*zin(203) + c01*zin(201)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  205

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  211

                                      xin(215) = c10*xin(203) + xc00*xin(209) + c01*xin(207)
                                      yin(215) = c10*yin(203) + yc00*yin(209) + c01*yin(207)
                                      zin(215) = c10*zin(203) + zc00*zin(209) + c01*zin(207)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  211

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  163

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  181

                                      xin(186) = c10*xin(150) + xc00*xin(168) + c01*xin(167)
                                      yin(186) = c10*yin(150) + yc00*yin(168) + c01*yin(167)
                                      zin(186) = c10*zin(150) + zc00*zin(168) + c01*zin(167)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  181

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  199

                                      xin(204) = c10*xin(168) + xc00*xin(186) + c01*xin(185)
                                      yin(204) = c10*yin(168) + yc00*yin(186) + c01*yin(185)
                                      zin(204) = c10*zin(168) + zc00*zin(186) + c01*zin(185)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  199

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  205

                                      xin(210) = c10*xin(186) + xc00*xin(204) + c01*xin(203)
                                      yin(210) = c10*yin(186) + yc00*yin(204) + c01*yin(203)
                                      zin(210) = c10*zin(186) + zc00*zin(204) + c01*zin(203)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  205

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  211

                                      xin(216) = c10*xin(204) + xc00*xin(210) + c01*xin(209)
                                      yin(216) = c10*yin(204) + yc00*yin(210) + c01*yin(209)
                                      zin(216) = c10*zin(204) + zc00*zin(210) + c01*zin(209)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  211

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  211

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(211) = xin(211) + dxij*xin(205)
                                      yin(211) = yin(211) + dyij*yin(205)
                                      zin(211) = zin(211) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  199

                                      xin(205) = xin(205) + dxij*xin(199)
                                      yin(205) = yin(205) + dyij*yin(199)
                                      zin(205) = zin(205) + dzij*zin(199)

                                      ! i3 = i4 =  199
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(211) = xin(211) + dxij*xin(205)
                                      yin(211) = yin(211) + dyij*yin(205)
                                      zin(211) = zin(211) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  151

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(163) + dxij*xin(145)
                                      yin(151) = yin(163) + dyij*yin(145)
                                      zin(151) = zin(163) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! ni =    2

                                      xin(169) = xin(181) + dxij*xin(163)
                                      yin(169) = yin(181) + dyij*yin(163)
                                      zin(169) = zin(181) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    3

                                      xin(187) = xin(199) + dxij*xin(181)
                                      yin(187) = yin(199) + dyij*yin(181)
                                      zin(187) = zin(199) + dzij*zin(181)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    2

                                      ! i4 = i3 =  157

                                      ! do ni = 1,    3

                                      xin(157) = xin(169) + dxij*xin(151)
                                      yin(157) = yin(169) + dyij*yin(151)
                                      zin(157) = zin(169) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    2

                                      xin(175) = xin(187) + dxij*xin(169)
                                      yin(175) = yin(187) + dyij*yin(169)
                                      zin(175) = zin(187) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! ni =    3

                                      xin(193) = xin(205) + dxij*xin(187)
                                      yin(193) = yin(205) + dyij*yin(187)
                                      zin(193) = zin(205) + dzij*zin(187)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  163

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(213) = xin(213) + dxij*xin(207)
                                      yin(213) = yin(213) + dyij*yin(207)
                                      zin(213) = zin(213) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  201

                                      xin(207) = xin(207) + dxij*xin(201)
                                      yin(207) = yin(207) + dyij*yin(201)
                                      zin(207) = zin(207) + dzij*zin(201)

                                      ! i3 = i4 =  201
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(213) = xin(213) + dxij*xin(207)
                                      yin(213) = yin(213) + dyij*yin(207)
                                      zin(213) = zin(213) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  153

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(165) + dxij*xin(147)
                                      yin(153) = yin(165) + dyij*yin(147)
                                      zin(153) = zin(165) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  171

                                      ! ni =    2

                                      xin(171) = xin(183) + dxij*xin(165)
                                      yin(171) = yin(183) + dyij*yin(165)
                                      zin(171) = zin(183) + dzij*zin(165)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    3

                                      xin(189) = xin(201) + dxij*xin(183)
                                      yin(189) = yin(201) + dyij*yin(183)
                                      zin(189) = zin(201) + dzij*zin(183)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    2

                                      ! i4 = i3 =  159

                                      ! do ni = 1,    3

                                      xin(159) = xin(171) + dxij*xin(153)
                                      yin(159) = yin(171) + dyij*yin(153)
                                      zin(159) = zin(171) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    2

                                      xin(177) = xin(189) + dxij*xin(171)
                                      yin(177) = yin(189) + dyij*yin(171)
                                      zin(177) = zin(189) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  195

                                      ! ni =    3

                                      xin(195) = xin(207) + dxij*xin(189)
                                      yin(195) = yin(207) + dyij*yin(189)
                                      zin(195) = zin(207) + dzij*zin(189)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  165

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  209

                                      xin(215) = xin(215) + dxij*xin(209)
                                      yin(215) = yin(215) + dyij*yin(209)
                                      zin(215) = zin(215) + dzij*zin(209)

                                      ! i3 = i4 =  209
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  203

                                      xin(209) = xin(209) + dxij*xin(203)
                                      yin(209) = yin(209) + dyij*yin(203)
                                      zin(209) = zin(209) + dzij*zin(203)

                                      ! i3 = i4 =  203
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  209

                                      xin(215) = xin(215) + dxij*xin(209)
                                      yin(215) = yin(215) + dyij*yin(209)
                                      zin(215) = zin(215) + dzij*zin(209)

                                      ! i3 = i4 =  209
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  155

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(167) + dxij*xin(149)
                                      yin(155) = yin(167) + dyij*yin(149)
                                      zin(155) = zin(167) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    2

                                      xin(173) = xin(185) + dxij*xin(167)
                                      yin(173) = yin(185) + dyij*yin(167)
                                      zin(173) = zin(185) + dzij*zin(167)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    3

                                      xin(191) = xin(203) + dxij*xin(185)
                                      yin(191) = yin(203) + dyij*yin(185)
                                      zin(191) = zin(203) + dzij*zin(185)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  161

                                      ! nj =    2

                                      ! i4 = i3 =  161

                                      ! do ni = 1,    3

                                      xin(161) = xin(173) + dxij*xin(155)
                                      yin(161) = yin(173) + dyij*yin(155)
                                      zin(161) = zin(173) + dzij*zin(155)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    2

                                      xin(179) = xin(191) + dxij*xin(173)
                                      yin(179) = yin(191) + dyij*yin(173)
                                      zin(179) = zin(191) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  197

                                      ! ni =    3

                                      xin(197) = xin(209) + dxij*xin(191)
                                      yin(197) = yin(209) + dyij*yin(191)
                                      zin(197) = zin(209) + dzij*zin(191)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  167

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  210

                                      xin(216) = xin(216) + dxij*xin(210)
                                      yin(216) = yin(216) + dyij*yin(210)
                                      zin(216) = zin(216) + dzij*zin(210)

                                      ! i3 = i4 =  210
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  204

                                      xin(210) = xin(210) + dxij*xin(204)
                                      yin(210) = yin(210) + dyij*yin(204)
                                      zin(210) = zin(210) + dzij*zin(204)

                                      ! i3 = i4 =  204
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  210

                                      xin(216) = xin(216) + dxij*xin(210)
                                      yin(216) = yin(216) + dyij*yin(210)
                                      zin(216) = zin(216) + dzij*zin(210)

                                      ! i3 = i4 =  210
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  156

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  156

                                      ! do ni = 1,    3

                                      xin(156) = xin(168) + dxij*xin(150)
                                      yin(156) = yin(168) + dyij*yin(150)
                                      zin(156) = zin(168) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  174

                                      ! ni =    2

                                      xin(174) = xin(186) + dxij*xin(168)
                                      yin(174) = yin(186) + dyij*yin(168)
                                      zin(174) = zin(186) + dzij*zin(168)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    3

                                      xin(192) = xin(204) + dxij*xin(186)
                                      yin(192) = yin(204) + dyij*yin(186)
                                      zin(192) = zin(204) + dzij*zin(186)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  210

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  162

                                      ! nj =    2

                                      ! i4 = i3 =  162

                                      ! do ni = 1,    3

                                      xin(162) = xin(174) + dxij*xin(156)
                                      yin(162) = yin(174) + dyij*yin(156)
                                      zin(162) = zin(174) + dzij*zin(156)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    2

                                      xin(180) = xin(192) + dxij*xin(174)
                                      yin(180) = yin(192) + dyij*yin(174)
                                      zin(180) = zin(192) + dzij*zin(174)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  198

                                      ! ni =    3

                                      xin(198) = xin(210) + dxij*xin(192)
                                      yin(198) = yin(210) + dyij*yin(192)
                                      zin(198) = zin(210) + dzij*zin(192)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  168

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  145

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  199

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  216

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

                                      ! i1 = in(1) =  217

                                      xin(217) = 1.0_dp
                                      yin(217) = 1.0_dp
                                      zin(217) = f00

                                      ! i2 = in(2) =  235
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(235) = xc00
                                      yin(235) = yc00
                                      zin(235) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  219

                                      xin(219) = xcp00
                                      yin(219) = ycp00
                                      zin(219) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  237
                                      ! i2 =  235

                                      xin(237) = xcp00*xin(235) + cp10
                                      yin(237) = ycp00*yin(235) + cp10
                                      zin(237) = zcp00*zin(235) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  235

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  253
                                      ! i3 =  217
                                      ! i4 =  235

                                      xin(253) = c10*xin(217) + xc00*xin(235)
                                      yin(253) = c10*yin(217) + yc00*yin(235)
                                      zin(253) = c10*zin(217) + zc00*zin(235)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  255
                                      ! i5 =  253
                                      ! i4 =  235

                                      xin(255) = xcp00*xin(253) + cp10*xin(235)
                                      yin(255) = ycp00*yin(253) + cp10*yin(235)
                                      zin(255) = zcp00*zin(253) + cp10*zin(235)

                                      ! ------------------

                                      ! i3 = i4 =  235
                                      ! i4 = i5 =  253

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  271
                                      ! i3 =  235
                                      ! i4 =  253

                                      xin(271) = c10*xin(235) + xc00*xin(253)
                                      yin(271) = c10*yin(235) + yc00*yin(253)
                                      zin(271) = c10*zin(235) + zc00*zin(253)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  273
                                      ! i5 =  271
                                      ! i4 =  253

                                      xin(273) = xcp00*xin(271) + cp10*xin(253)
                                      yin(273) = ycp00*yin(271) + cp10*yin(253)
                                      zin(273) = zcp00*zin(271) + cp10*zin(253)

                                      ! ------------------

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  271

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  277
                                      ! i3 =  253
                                      ! i4 =  271

                                      xin(277) = c10*xin(253) + xc00*xin(271)
                                      yin(277) = c10*yin(253) + yc00*yin(271)
                                      zin(277) = c10*zin(253) + zc00*zin(271)

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

                                      ! n =    5

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

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  217
                                      ! i4 = i1+k2 =  219

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  221
                                      ! i3 =  217
                                      ! i4 =  219

                                      xin(221) = cp01*xin(217) + xcp00*xin(219)
                                      yin(221) = cp01*yin(217) + ycp00*yin(219)
                                      zin(221) = cp01*zin(217) + zcp00*zin(219)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  239

                                      xin(239) = xc00*xin(221) + c01*xin(219)
                                      yin(239) = yc00*yin(221) + c01*yin(219)
                                      zin(239) = zc00*zin(221) + c01*zin(219)

                                      ! ------------------

                                      ! i3 = i4 =  219
                                      ! i4 = i5 =  221

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  222
                                      ! i3 =  219
                                      ! i4 =  221

                                      xin(222) = cp01*xin(219) + xcp00*xin(221)
                                      yin(222) = cp01*yin(219) + ycp00*yin(221)
                                      zin(222) = cp01*zin(219) + zcp00*zin(221)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  240

                                      xin(240) = xc00*xin(222) + c01*xin(221)
                                      yin(240) = yc00*yin(222) + c01*yin(221)
                                      zin(240) = zc00*zin(222) + c01*zin(221)

                                      ! ------------------

                                      ! i3 = i4 =  221
                                      ! i4 = i5 =  222

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  235

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  253

                                      xin(257) = c10*xin(221) + xc00*xin(239) + c01*xin(237)
                                      yin(257) = c10*yin(221) + yc00*yin(239) + c01*yin(237)
                                      zin(257) = c10*zin(221) + zc00*zin(239) + c01*zin(237)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  235
                                      ! i4 = i5 =  253

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  271

                                      xin(275) = c10*xin(239) + xc00*xin(257) + c01*xin(255)
                                      yin(275) = c10*yin(239) + yc00*yin(257) + c01*yin(255)
                                      zin(275) = c10*zin(239) + zc00*zin(257) + c01*zin(255)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  271

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  277

                                      xin(281) = c10*xin(257) + xc00*xin(275) + c01*xin(273)
                                      yin(281) = c10*yin(257) + yc00*yin(275) + c01*yin(273)
                                      zin(281) = c10*zin(257) + zc00*zin(275) + c01*zin(273)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  277

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  283

                                      xin(287) = c10*xin(275) + xc00*xin(281) + c01*xin(279)
                                      yin(287) = c10*yin(275) + yc00*yin(281) + c01*yin(279)
                                      zin(287) = c10*zin(275) + zc00*zin(281) + c01*zin(279)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  277
                                      ! i4 = i5 =  283

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  235

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  253

                                      xin(258) = c10*xin(222) + xc00*xin(240) + c01*xin(239)
                                      yin(258) = c10*yin(222) + yc00*yin(240) + c01*yin(239)
                                      zin(258) = c10*zin(222) + zc00*zin(240) + c01*zin(239)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  235
                                      ! i4 = i5 =  253

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  271

                                      xin(276) = c10*xin(240) + xc00*xin(258) + c01*xin(257)
                                      yin(276) = c10*yin(240) + yc00*yin(258) + c01*yin(257)
                                      zin(276) = c10*zin(240) + zc00*zin(258) + c01*zin(257)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  253
                                      ! i4 = i5 =  271

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  277

                                      xin(282) = c10*xin(258) + xc00*xin(276) + c01*xin(275)
                                      yin(282) = c10*yin(258) + yc00*yin(276) + c01*yin(275)
                                      zin(282) = c10*zin(258) + zc00*zin(276) + c01*zin(275)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  277

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  283

                                      xin(288) = c10*xin(276) + xc00*xin(282) + c01*xin(281)
                                      yin(288) = c10*yin(276) + yc00*yin(282) + c01*yin(281)
                                      zin(288) = c10*zin(276) + zc00*zin(282) + c01*zin(281)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  277
                                      ! i4 = i5 =  283

                                      ! nn =    6

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

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(283) = xin(283) + dxij*xin(277)
                                      yin(283) = yin(283) + dyij*yin(277)
                                      zin(283) = zin(283) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  271

                                      xin(277) = xin(277) + dxij*xin(271)
                                      yin(277) = yin(277) + dyij*yin(271)
                                      zin(277) = zin(277) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(283) = xin(283) + dxij*xin(277)
                                      yin(283) = yin(283) + dyij*yin(277)
                                      zin(283) = zin(283) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  223

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  223

                                      ! do ni = 1,    3

                                      xin(223) = xin(235) + dxij*xin(217)
                                      yin(223) = yin(235) + dyij*yin(217)
                                      zin(223) = zin(235) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! ni =    2

                                      xin(241) = xin(253) + dxij*xin(235)
                                      yin(241) = yin(253) + dyij*yin(235)
                                      zin(241) = zin(253) + dzij*zin(235)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  259

                                      ! ni =    3

                                      xin(259) = xin(271) + dxij*xin(253)
                                      yin(259) = yin(271) + dyij*yin(253)
                                      zin(259) = zin(271) + dzij*zin(253)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  229

                                      ! nj =    2

                                      ! i4 = i3 =  229

                                      ! do ni = 1,    3

                                      xin(229) = xin(241) + dxij*xin(223)
                                      yin(229) = yin(241) + dyij*yin(223)
                                      zin(229) = zin(241) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                      ! ni =    2

                                      xin(247) = xin(259) + dxij*xin(241)
                                      yin(247) = yin(259) + dyij*yin(241)
                                      zin(247) = zin(259) + dzij*zin(241)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  265

                                      ! ni =    3

                                      xin(265) = xin(277) + dxij*xin(259)
                                      yin(265) = yin(277) + dyij*yin(259)
                                      zin(265) = zin(277) + dzij*zin(259)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  235

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(285) = xin(285) + dxij*xin(279)
                                      yin(285) = yin(285) + dyij*yin(279)
                                      zin(285) = zin(285) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  273

                                      xin(279) = xin(279) + dxij*xin(273)
                                      yin(279) = yin(279) + dyij*yin(273)
                                      zin(279) = zin(279) + dzij*zin(273)

                                      ! i3 = i4 =  273
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(285) = xin(285) + dxij*xin(279)
                                      yin(285) = yin(285) + dyij*yin(279)
                                      zin(285) = zin(285) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  225

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  225

                                      ! do ni = 1,    3

                                      xin(225) = xin(237) + dxij*xin(219)
                                      yin(225) = yin(237) + dyij*yin(219)
                                      zin(225) = zin(237) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  243

                                      ! ni =    2

                                      xin(243) = xin(255) + dxij*xin(237)
                                      yin(243) = yin(255) + dyij*yin(237)
                                      zin(243) = zin(255) + dzij*zin(237)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  261

                                      ! ni =    3

                                      xin(261) = xin(273) + dxij*xin(255)
                                      yin(261) = yin(273) + dyij*yin(255)
                                      zin(261) = zin(273) + dzij*zin(255)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  279

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  231

                                      ! nj =    2

                                      ! i4 = i3 =  231

                                      ! do ni = 1,    3

                                      xin(231) = xin(243) + dxij*xin(225)
                                      yin(231) = yin(243) + dyij*yin(225)
                                      zin(231) = zin(243) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  249

                                      ! ni =    2

                                      xin(249) = xin(261) + dxij*xin(243)
                                      yin(249) = yin(261) + dyij*yin(243)
                                      zin(249) = zin(261) + dzij*zin(243)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  267

                                      ! ni =    3

                                      xin(267) = xin(279) + dxij*xin(261)
                                      yin(267) = yin(279) + dyij*yin(261)
                                      zin(267) = zin(279) + dzij*zin(261)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  285

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  237

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  281

                                      xin(287) = xin(287) + dxij*xin(281)
                                      yin(287) = yin(287) + dyij*yin(281)
                                      zin(287) = zin(287) + dzij*zin(281)

                                      ! i3 = i4 =  281
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  275

                                      xin(281) = xin(281) + dxij*xin(275)
                                      yin(281) = yin(281) + dyij*yin(275)
                                      zin(281) = zin(281) + dzij*zin(275)

                                      ! i3 = i4 =  275
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  281

                                      xin(287) = xin(287) + dxij*xin(281)
                                      yin(287) = yin(287) + dyij*yin(281)
                                      zin(287) = zin(287) + dzij*zin(281)

                                      ! i3 = i4 =  281
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  227

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  227

                                      ! do ni = 1,    3

                                      xin(227) = xin(239) + dxij*xin(221)
                                      yin(227) = yin(239) + dyij*yin(221)
                                      zin(227) = zin(239) + dzij*zin(221)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  245

                                      ! ni =    2

                                      xin(245) = xin(257) + dxij*xin(239)
                                      yin(245) = yin(257) + dyij*yin(239)
                                      zin(245) = zin(257) + dzij*zin(239)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  263

                                      ! ni =    3

                                      xin(263) = xin(275) + dxij*xin(257)
                                      yin(263) = yin(275) + dyij*yin(257)
                                      zin(263) = zin(275) + dzij*zin(257)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  281

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  233

                                      ! nj =    2

                                      ! i4 = i3 =  233

                                      ! do ni = 1,    3

                                      xin(233) = xin(245) + dxij*xin(227)
                                      yin(233) = yin(245) + dyij*yin(227)
                                      zin(233) = zin(245) + dzij*zin(227)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    2

                                      xin(251) = xin(263) + dxij*xin(245)
                                      yin(251) = yin(263) + dyij*yin(245)
                                      zin(251) = zin(263) + dzij*zin(245)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  269

                                      ! ni =    3

                                      xin(269) = xin(281) + dxij*xin(263)
                                      yin(269) = yin(281) + dyij*yin(263)
                                      zin(269) = zin(281) + dzij*zin(263)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  239

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  282

                                      xin(288) = xin(288) + dxij*xin(282)
                                      yin(288) = yin(288) + dyij*yin(282)
                                      zin(288) = zin(288) + dzij*zin(282)

                                      ! i3 = i4 =  282
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  276

                                      xin(282) = xin(282) + dxij*xin(276)
                                      yin(282) = yin(282) + dyij*yin(276)
                                      zin(282) = zin(282) + dzij*zin(276)

                                      ! i3 = i4 =  276
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  282

                                      xin(288) = xin(288) + dxij*xin(282)
                                      yin(288) = yin(288) + dyij*yin(282)
                                      zin(288) = zin(288) + dzij*zin(282)

                                      ! i3 = i4 =  282
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  228

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  228

                                      ! do ni = 1,    3

                                      xin(228) = xin(240) + dxij*xin(222)
                                      yin(228) = yin(240) + dyij*yin(222)
                                      zin(228) = zin(240) + dzij*zin(222)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  246

                                      ! ni =    2

                                      xin(246) = xin(258) + dxij*xin(240)
                                      yin(246) = yin(258) + dyij*yin(240)
                                      zin(246) = zin(258) + dzij*zin(240)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  264

                                      ! ni =    3

                                      xin(264) = xin(276) + dxij*xin(258)
                                      yin(264) = yin(276) + dyij*yin(258)
                                      zin(264) = zin(276) + dzij*zin(258)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  282

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  234

                                      ! nj =    2

                                      ! i4 = i3 =  234

                                      ! do ni = 1,    3

                                      xin(234) = xin(246) + dxij*xin(228)
                                      yin(234) = yin(246) + dyij*yin(228)
                                      zin(234) = zin(246) + dzij*zin(228)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  252

                                      ! ni =    2

                                      xin(252) = xin(264) + dxij*xin(246)
                                      yin(252) = yin(264) + dyij*yin(246)
                                      zin(252) = zin(264) + dzij*zin(246)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  270

                                      ! ni =    3

                                      xin(270) = xin(282) + dxij*xin(264)
                                      yin(270) = yin(282) + dyij*yin(264)
                                      zin(270) = zin(282) + dzij*zin(264)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  240

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  217

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  271

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  288

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

                                      ! i1 = in(1) =  289

                                      xin(289) = 1.0_dp
                                      yin(289) = 1.0_dp
                                      zin(289) = f00

                                      ! i2 = in(2) =  307
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(307) = xc00
                                      yin(307) = yc00
                                      zin(307) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  291

                                      xin(291) = xcp00
                                      yin(291) = ycp00
                                      zin(291) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  309
                                      ! i2 =  307

                                      xin(309) = xcp00*xin(307) + cp10
                                      yin(309) = ycp00*yin(307) + cp10
                                      zin(309) = zcp00*zin(307) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  307

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  325
                                      ! i3 =  289
                                      ! i4 =  307

                                      xin(325) = c10*xin(289) + xc00*xin(307)
                                      yin(325) = c10*yin(289) + yc00*yin(307)
                                      zin(325) = c10*zin(289) + zc00*zin(307)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  327
                                      ! i5 =  325
                                      ! i4 =  307

                                      xin(327) = xcp00*xin(325) + cp10*xin(307)
                                      yin(327) = ycp00*yin(325) + cp10*yin(307)
                                      zin(327) = zcp00*zin(325) + cp10*zin(307)

                                      ! ------------------

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  325

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  343
                                      ! i3 =  307
                                      ! i4 =  325

                                      xin(343) = c10*xin(307) + xc00*xin(325)
                                      yin(343) = c10*yin(307) + yc00*yin(325)
                                      zin(343) = c10*zin(307) + zc00*zin(325)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  345
                                      ! i5 =  343
                                      ! i4 =  325

                                      xin(345) = xcp00*xin(343) + cp10*xin(325)
                                      yin(345) = ycp00*yin(343) + cp10*yin(325)
                                      zin(345) = zcp00*zin(343) + cp10*zin(325)

                                      ! ------------------

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  343

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  349
                                      ! i3 =  325
                                      ! i4 =  343

                                      xin(349) = c10*xin(325) + xc00*xin(343)
                                      yin(349) = c10*yin(325) + yc00*yin(343)
                                      zin(349) = c10*zin(325) + zc00*zin(343)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  351
                                      ! i5 =  349
                                      ! i4 =  343

                                      xin(351) = xcp00*xin(349) + cp10*xin(343)
                                      yin(351) = ycp00*yin(349) + cp10*yin(343)
                                      zin(351) = zcp00*zin(349) + cp10*zin(343)

                                      ! ------------------

                                      ! i3 = i4 =  343
                                      ! i4 = i5 =  349

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  355
                                      ! i3 =  343
                                      ! i4 =  349

                                      xin(355) = c10*xin(343) + xc00*xin(349)
                                      yin(355) = c10*yin(343) + yc00*yin(349)
                                      zin(355) = c10*zin(343) + zc00*zin(349)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  357
                                      ! i5 =  355
                                      ! i4 =  349

                                      xin(357) = xcp00*xin(355) + cp10*xin(349)
                                      yin(357) = ycp00*yin(355) + cp10*yin(349)
                                      zin(357) = zcp00*zin(355) + cp10*zin(349)

                                      ! ------------------

                                      ! i3 = i4 =  349
                                      ! i4 = i5 =  355

                                      ! n =    6

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

                                      ! i3 = i2+kn(n+1) =  311

                                      xin(311) = xc00*xin(293) + c01*xin(291)
                                      yin(311) = yc00*yin(293) + c01*yin(291)
                                      zin(311) = zc00*zin(293) + c01*zin(291)

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

                                      ! i3 = i2+kn(n+1) =  312

                                      xin(312) = xc00*xin(294) + c01*xin(293)
                                      yin(312) = yc00*yin(294) + c01*yin(293)
                                      zin(312) = zc00*zin(294) + c01*zin(293)

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
                                      ! i4 = i2 =  307

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  325

                                      xin(329) = c10*xin(293) + xc00*xin(311) + c01*xin(309)
                                      yin(329) = c10*yin(293) + yc00*yin(311) + c01*yin(309)
                                      zin(329) = c10*zin(293) + zc00*zin(311) + c01*zin(309)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  325

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  343

                                      xin(347) = c10*xin(311) + xc00*xin(329) + c01*xin(327)
                                      yin(347) = c10*yin(311) + yc00*yin(329) + c01*yin(327)
                                      zin(347) = c10*zin(311) + zc00*zin(329) + c01*zin(327)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  343

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  349

                                      xin(353) = c10*xin(329) + xc00*xin(347) + c01*xin(345)
                                      yin(353) = c10*yin(329) + yc00*yin(347) + c01*yin(345)
                                      zin(353) = c10*zin(329) + zc00*zin(347) + c01*zin(345)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  343
                                      ! i4 = i5 =  349

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  355

                                      xin(359) = c10*xin(347) + xc00*xin(353) + c01*xin(351)
                                      yin(359) = c10*yin(347) + yc00*yin(353) + c01*yin(351)
                                      zin(359) = c10*zin(347) + zc00*zin(353) + c01*zin(351)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  349
                                      ! i4 = i5 =  355

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  307

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  325

                                      xin(330) = c10*xin(294) + xc00*xin(312) + c01*xin(311)
                                      yin(330) = c10*yin(294) + yc00*yin(312) + c01*yin(311)
                                      zin(330) = c10*zin(294) + zc00*zin(312) + c01*zin(311)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  325

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  343

                                      xin(348) = c10*xin(312) + xc00*xin(330) + c01*xin(329)
                                      yin(348) = c10*yin(312) + yc00*yin(330) + c01*yin(329)
                                      zin(348) = c10*zin(312) + zc00*zin(330) + c01*zin(329)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  325
                                      ! i4 = i5 =  343

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  349

                                      xin(354) = c10*xin(330) + xc00*xin(348) + c01*xin(347)
                                      yin(354) = c10*yin(330) + yc00*yin(348) + c01*yin(347)
                                      zin(354) = c10*zin(330) + zc00*zin(348) + c01*zin(347)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  343
                                      ! i4 = i5 =  349

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  355

                                      xin(360) = c10*xin(348) + xc00*xin(354) + c01*xin(353)
                                      yin(360) = c10*yin(348) + yc00*yin(354) + c01*yin(353)
                                      zin(360) = c10*zin(348) + zc00*zin(354) + c01*zin(353)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  349
                                      ! i4 = i5 =  355

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  355

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  355

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  349

                                      xin(355) = xin(355) + dxij*xin(349)
                                      yin(355) = yin(355) + dyij*yin(349)
                                      zin(355) = zin(355) + dzij*zin(349)

                                      ! i3 = i4 =  349
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  343

                                      xin(349) = xin(349) + dxij*xin(343)
                                      yin(349) = yin(349) + dyij*yin(343)
                                      zin(349) = zin(349) + dzij*zin(343)

                                      ! i3 = i4 =  343
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  355

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  349

                                      xin(355) = xin(355) + dxij*xin(349)
                                      yin(355) = yin(355) + dyij*yin(349)
                                      zin(355) = zin(355) + dzij*zin(349)

                                      ! i3 = i4 =  349
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  295

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  295

                                      ! do ni = 1,    3

                                      xin(295) = xin(307) + dxij*xin(289)
                                      yin(295) = yin(307) + dyij*yin(289)
                                      zin(295) = zin(307) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! ni =    2

                                      xin(313) = xin(325) + dxij*xin(307)
                                      yin(313) = yin(325) + dyij*yin(307)
                                      zin(313) = zin(325) + dzij*zin(307)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  331

                                      ! ni =    3

                                      xin(331) = xin(343) + dxij*xin(325)
                                      yin(331) = yin(343) + dyij*yin(325)
                                      zin(331) = zin(343) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  301

                                      ! nj =    2

                                      ! i4 = i3 =  301

                                      ! do ni = 1,    3

                                      xin(301) = xin(313) + dxij*xin(295)
                                      yin(301) = yin(313) + dyij*yin(295)
                                      zin(301) = zin(313) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  319

                                      ! ni =    2

                                      xin(319) = xin(331) + dxij*xin(313)
                                      yin(319) = yin(331) + dyij*yin(313)
                                      zin(319) = zin(331) + dzij*zin(313)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  337

                                      ! ni =    3

                                      xin(337) = xin(349) + dxij*xin(331)
                                      yin(337) = yin(349) + dyij*yin(331)
                                      zin(337) = zin(349) + dzij*zin(331)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  355

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  307

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  357

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  351

                                      xin(357) = xin(357) + dxij*xin(351)
                                      yin(357) = yin(357) + dyij*yin(351)
                                      zin(357) = zin(357) + dzij*zin(351)

                                      ! i3 = i4 =  351
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  345

                                      xin(351) = xin(351) + dxij*xin(345)
                                      yin(351) = yin(351) + dyij*yin(345)
                                      zin(351) = zin(351) + dzij*zin(345)

                                      ! i3 = i4 =  345
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  357

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  351

                                      xin(357) = xin(357) + dxij*xin(351)
                                      yin(357) = yin(357) + dyij*yin(351)
                                      zin(357) = zin(357) + dzij*zin(351)

                                      ! i3 = i4 =  351
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  297

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  297

                                      ! do ni = 1,    3

                                      xin(297) = xin(309) + dxij*xin(291)
                                      yin(297) = yin(309) + dyij*yin(291)
                                      zin(297) = zin(309) + dzij*zin(291)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  315

                                      ! ni =    2

                                      xin(315) = xin(327) + dxij*xin(309)
                                      yin(315) = yin(327) + dyij*yin(309)
                                      zin(315) = zin(327) + dzij*zin(309)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  333

                                      ! ni =    3

                                      xin(333) = xin(345) + dxij*xin(327)
                                      yin(333) = yin(345) + dyij*yin(327)
                                      zin(333) = zin(345) + dzij*zin(327)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  351

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  303

                                      ! nj =    2

                                      ! i4 = i3 =  303

                                      ! do ni = 1,    3

                                      xin(303) = xin(315) + dxij*xin(297)
                                      yin(303) = yin(315) + dyij*yin(297)
                                      zin(303) = zin(315) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  321

                                      ! ni =    2

                                      xin(321) = xin(333) + dxij*xin(315)
                                      yin(321) = yin(333) + dyij*yin(315)
                                      zin(321) = zin(333) + dzij*zin(315)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  339

                                      ! ni =    3

                                      xin(339) = xin(351) + dxij*xin(333)
                                      yin(339) = yin(351) + dyij*yin(333)
                                      zin(339) = zin(351) + dzij*zin(333)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  357

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  309

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  359

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  353

                                      xin(359) = xin(359) + dxij*xin(353)
                                      yin(359) = yin(359) + dyij*yin(353)
                                      zin(359) = zin(359) + dzij*zin(353)

                                      ! i3 = i4 =  353
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  347

                                      xin(353) = xin(353) + dxij*xin(347)
                                      yin(353) = yin(353) + dyij*yin(347)
                                      zin(353) = zin(353) + dzij*zin(347)

                                      ! i3 = i4 =  347
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  359

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  353

                                      xin(359) = xin(359) + dxij*xin(353)
                                      yin(359) = yin(359) + dyij*yin(353)
                                      zin(359) = zin(359) + dzij*zin(353)

                                      ! i3 = i4 =  353
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  299

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  299

                                      ! do ni = 1,    3

                                      xin(299) = xin(311) + dxij*xin(293)
                                      yin(299) = yin(311) + dyij*yin(293)
                                      zin(299) = zin(311) + dzij*zin(293)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  317

                                      ! ni =    2

                                      xin(317) = xin(329) + dxij*xin(311)
                                      yin(317) = yin(329) + dyij*yin(311)
                                      zin(317) = zin(329) + dzij*zin(311)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  335

                                      ! ni =    3

                                      xin(335) = xin(347) + dxij*xin(329)
                                      yin(335) = yin(347) + dyij*yin(329)
                                      zin(335) = zin(347) + dzij*zin(329)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  353

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  305

                                      ! nj =    2

                                      ! i4 = i3 =  305

                                      ! do ni = 1,    3

                                      xin(305) = xin(317) + dxij*xin(299)
                                      yin(305) = yin(317) + dyij*yin(299)
                                      zin(305) = zin(317) + dzij*zin(299)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  323

                                      ! ni =    2

                                      xin(323) = xin(335) + dxij*xin(317)
                                      yin(323) = yin(335) + dyij*yin(317)
                                      zin(323) = zin(335) + dzij*zin(317)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  341

                                      ! ni =    3

                                      xin(341) = xin(353) + dxij*xin(335)
                                      yin(341) = yin(353) + dyij*yin(335)
                                      zin(341) = zin(353) + dzij*zin(335)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  359

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  311

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  360

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  354

                                      xin(360) = xin(360) + dxij*xin(354)
                                      yin(360) = yin(360) + dyij*yin(354)
                                      zin(360) = zin(360) + dzij*zin(354)

                                      ! i3 = i4 =  354
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  348

                                      xin(354) = xin(354) + dxij*xin(348)
                                      yin(354) = yin(354) + dyij*yin(348)
                                      zin(354) = zin(354) + dzij*zin(348)

                                      ! i3 = i4 =  348
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  360

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  354

                                      xin(360) = xin(360) + dxij*xin(354)
                                      yin(360) = yin(360) + dyij*yin(354)
                                      zin(360) = zin(360) + dzij*zin(354)

                                      ! i3 = i4 =  354
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  300

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  300

                                      ! do ni = 1,    3

                                      xin(300) = xin(312) + dxij*xin(294)
                                      yin(300) = yin(312) + dyij*yin(294)
                                      zin(300) = zin(312) + dzij*zin(294)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  318

                                      ! ni =    2

                                      xin(318) = xin(330) + dxij*xin(312)
                                      yin(318) = yin(330) + dyij*yin(312)
                                      zin(318) = zin(330) + dzij*zin(312)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  336

                                      ! ni =    3

                                      xin(336) = xin(348) + dxij*xin(330)
                                      yin(336) = yin(348) + dyij*yin(330)
                                      zin(336) = zin(348) + dzij*zin(330)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  354

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  306

                                      ! nj =    2

                                      ! i4 = i3 =  306

                                      ! do ni = 1,    3

                                      xin(306) = xin(318) + dxij*xin(300)
                                      yin(306) = yin(318) + dyij*yin(300)
                                      zin(306) = zin(318) + dzij*zin(300)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  324

                                      ! ni =    2

                                      xin(324) = xin(336) + dxij*xin(318)
                                      yin(324) = yin(336) + dyij*yin(318)
                                      zin(324) = zin(336) + dzij*zin(318)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  342

                                      ! ni =    3

                                      xin(342) = xin(354) + dxij*xin(336)
                                      yin(342) = yin(354) + dyij*yin(336)
                                      zin(342) = zin(354) + dzij*zin(336)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  360

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  312

                                      ! nj =    3

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  307

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  325

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  343

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

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

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  360

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 1080! loop over all integrals

                                        l = n - 18*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d23bra(j)*d12ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 72)*yin(my + 72)*zin(mz + 72) & ! root  2
                                                        + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  3
                                                        + xin(mx + 216)*yin(my + 216)*zin(mz + 216) & ! root  4
                                                        + xin(mx + 288)*yin(my + 288)*zin(mz + 288)) ! root  5

                                        j = int(n/18) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*108 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

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

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n12ket)
                              deallocate (xint12ket)

                              end subroutine int3221
                              end submodule
