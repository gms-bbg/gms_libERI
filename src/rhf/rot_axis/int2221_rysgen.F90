! The total angular momentum of this class is:           7
! The algorithm chosen is: Rys quadrature
submodule(rot_axis_kernels) int2221_impl
contains
  module subroutine int2221(dd_pair, pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: dd_pair, pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n22bra(:), n12ket(:)
    real(dp), allocatable :: xint22bra(:), xint12ket(:)
    integer(kind=int64) :: nddbra, npdket
    real(dp) :: scutddbra, scutpdket, test
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(216), yin(216), zin(216)
    real(dp) :: eri_value(648)
    real(dp) :: d22bra(36), d12ket(18)
    integer(kind=int64) :: ix(6), jx(6), kx(6), lx(3)
    integer(kind=int64) :: iy(6), jy(6), ky(6), ly(3)
    integer(kind=int64) :: iz(6), jz(6), kz(6), lz(3)
    integer(kind=int64) :: in(5), in1(5), kn(4)
    integer(kind=int64) :: ijx(36), ijy(36), ijz(36)
    integer(kind=int64) :: klx(18), kly(18), klz(18)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 19
    in1(3) = 37
    in1(4) = 43
    in1(5) = 49

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

    ix(1) = 37
    ix(2) = 1
    ix(3) = 1
    ix(4) = 19
    ix(5) = 19
    ix(6) = 1

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
    iy(2) = 37
    iy(3) = 1
    iy(4) = 19
    iy(5) = 1
    iy(6) = 19

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
    iz(3) = 37
    iz(4) = 1
    iz(5) = 19
    iz(6) = 19

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 49
    ijx(2) = 37
    ijx(3) = 37
    ijx(4) = 43
    ijx(5) = 43
    ijx(6) = 37
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
    ijx(19) = 31
    ijx(20) = 19
    ijx(21) = 19
    ijx(22) = 25
    ijx(23) = 25
    ijx(24) = 19
    ijx(25) = 31
    ijx(26) = 19
    ijx(27) = 19
    ijx(28) = 25
    ijx(29) = 25
    ijx(30) = 19
    ijx(31) = 13
    ijx(32) = 1
    ijx(33) = 1
    ijx(34) = 7
    ijx(35) = 7
    ijx(36) = 1

    ijy(1) = 1
    ijy(2) = 13
    ijy(3) = 1
    ijy(4) = 7
    ijy(5) = 1
    ijy(6) = 7
    ijy(7) = 37
    ijy(8) = 49
    ijy(9) = 37
    ijy(10) = 43
    ijy(11) = 37
    ijy(12) = 43
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
    ijy(31) = 19
    ijy(32) = 31
    ijy(33) = 19
    ijy(34) = 25
    ijy(35) = 19
    ijy(36) = 25

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
    ijz(13) = 37
    ijz(14) = 37
    ijz(15) = 49
    ijz(16) = 37
    ijz(17) = 43
    ijz(18) = 43
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
    ijz(31) = 19
    ijz(32) = 19
    ijz(33) = 31
    ijz(34) = 19
    ijz(35) = 25
    ijz(36) = 25

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

    allocate (n22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (n12ket(res%n_p_shl*res%n_d_shl))
    allocate (xint12ket(res%n_p_shl*res%n_d_shl))

    ! Start screening

    scutddbra = cutoff_schwarz/maxval(dd_pair%xints)
    nddbra = 0
    do ij = 1, res%n_d_shl*(res%n_d_shl + 1)/2
      if (dd_pair%xints(ij) .ge. scutddbra) then
        nddbra = nddbra + 1
        xint22bra(nddbra) = dd_pair%xints(ij)
        n22bra(nddbra) = ij
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

    if ((nddbra*npdket) .le. nchunksize_int64) nchunksize_int64 = nddbra*npdket
    ntile = int(nddbra*npdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nddbra*npdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nddbra, xint22bra, n22bra, xint12ket, n12ket, dd_pair, pd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d12ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d22bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,iandj)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nddbra) + 1
              kl_tmp = (iquart - 1)/nddbra + 1

              test = xint22bra(ij_tmp)*xint12ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n22bra(ij_tmp)
                kl = n12ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_d_shl(ish_tmp)
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

                    t_expon_ab = dd_pair%t_expon_ab(dd_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = dd_pair%expon_a(dd_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = dd_pair%expon_b(dd_pair%pair_loc(ij) + bra_loop)
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

                    d22bra(1) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(2) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(3) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(4) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(5) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(6) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(7) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(8) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(9) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(10) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(11) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(12) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(13) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(14) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(15) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)
                    d22bra(16) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(17) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(18) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(19) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(20) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(21) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(22) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(23) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(24) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(25) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(26) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(27) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(28) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(29) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(30) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(31) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(32) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(33) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d22bra(34) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(35) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3
                    d22bra(36) = dd_pair%d_coeff_alt(dd_pair%pair_loc(ij) + bra_loop)*sqrt3*sqrt3

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

                                      ! do n = 2,   4

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

                                      ! i5 = in(n+1) =   43
                                      ! i3 =   19
                                      ! i4 =   37

                                      xin(43) = c10*xin(19) + xc00*xin(37)
                                      yin(43) = c10*yin(19) + yc00*yin(37)
                                      zin(43) = c10*zin(19) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   45
                                      ! i5 =   43
                                      ! i4 =   37

                                      xin(45) = xcp00*xin(43) + cp10*xin(37)
                                      yin(45) = ycp00*yin(43) + cp10*yin(37)
                                      zin(45) = zcp00*zin(43) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =   37
                                      ! i4 =   43

                                      xin(49) = c10*xin(37) + xc00*xin(43)
                                      yin(49) = c10*yin(37) + yc00*yin(43)
                                      zin(49) = c10*zin(37) + zc00*zin(43)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   51
                                      ! i5 =   49
                                      ! i4 =   43

                                      xin(51) = xcp00*xin(49) + cp10*xin(43)
                                      yin(51) = ycp00*yin(49) + cp10*yin(43)
                                      zin(51) = zcp00*zin(49) + cp10*zin(43)

                                      ! ------------------

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   49

                                      ! n =    5

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

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   37

                                      xin(41) = c10*xin(5) + xc00*xin(23) + c01*xin(21)
                                      yin(41) = c10*yin(5) + yc00*yin(23) + c01*yin(21)
                                      zin(41) = c10*zin(5) + zc00*zin(23) + c01*zin(21)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   37

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   43

                                      xin(47) = c10*xin(23) + xc00*xin(41) + c01*xin(39)
                                      yin(47) = c10*yin(23) + yc00*yin(41) + c01*yin(39)
                                      zin(47) = c10*zin(23) + zc00*zin(41) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   49

                                      xin(53) = c10*xin(41) + xc00*xin(47) + c01*xin(45)
                                      yin(53) = c10*yin(41) + yc00*yin(47) + c01*yin(45)
                                      zin(53) = c10*zin(41) + zc00*zin(47) + c01*zin(45)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   49

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   19

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   37

                                      xin(42) = c10*xin(6) + xc00*xin(24) + c01*xin(23)
                                      yin(42) = c10*yin(6) + yc00*yin(24) + c01*yin(23)
                                      zin(42) = c10*zin(6) + zc00*zin(24) + c01*zin(23)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   37

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   43

                                      xin(48) = c10*xin(24) + xc00*xin(42) + c01*xin(41)
                                      yin(48) = c10*yin(24) + yc00*yin(42) + c01*yin(41)
                                      zin(48) = c10*zin(24) + zc00*zin(42) + c01*zin(41)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   49

                                      xin(54) = c10*xin(42) + xc00*xin(48) + c01*xin(47)
                                      yin(54) = c10*yin(42) + yc00*yin(48) + c01*yin(47)
                                      zin(54) = c10*zin(42) + zc00*zin(48) + c01*zin(47)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   49

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   49

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   49

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(49) = xin(49) + dxij*xin(43)
                                      yin(49) = yin(49) + dyij*yin(43)
                                      zin(49) = zin(49) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   37

                                      xin(43) = xin(43) + dxij*xin(37)
                                      yin(43) = yin(43) + dyij*yin(37)
                                      zin(43) = zin(43) + dzij*zin(37)

                                      ! i3 = i4 =   37
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   49

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(49) = xin(49) + dxij*xin(43)
                                      yin(49) = yin(49) + dyij*yin(43)
                                      zin(49) = zin(49) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    2

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   51

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(51) = xin(51) + dxij*xin(45)
                                      yin(51) = yin(51) + dyij*yin(45)
                                      zin(51) = zin(51) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   39

                                      xin(45) = xin(45) + dxij*xin(39)
                                      yin(45) = yin(45) + dyij*yin(39)
                                      zin(45) = zin(45) + dzij*zin(39)

                                      ! i3 = i4 =   39
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   51

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   45

                                      xin(51) = xin(51) + dxij*xin(45)
                                      yin(51) = yin(51) + dyij*yin(45)
                                      zin(51) = zin(51) + dzij*zin(45)

                                      ! i3 = i4 =   45
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    2

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   53

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   47

                                      xin(53) = xin(53) + dxij*xin(47)
                                      yin(53) = yin(53) + dyij*yin(47)
                                      zin(53) = zin(53) + dzij*zin(47)

                                      ! i3 = i4 =   47
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   41

                                      xin(47) = xin(47) + dxij*xin(41)
                                      yin(47) = yin(47) + dyij*yin(41)
                                      zin(47) = zin(47) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   53

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   47

                                      xin(53) = xin(53) + dxij*xin(47)
                                      yin(53) = yin(53) + dyij*yin(47)
                                      zin(53) = zin(53) + dzij*zin(47)

                                      ! i3 = i4 =   47
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   23

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   54

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   48

                                      xin(54) = xin(54) + dxij*xin(48)
                                      yin(54) = yin(54) + dyij*yin(48)
                                      zin(54) = zin(54) + dzij*zin(48)

                                      ! i3 = i4 =   48
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   42

                                      xin(48) = xin(48) + dxij*xin(42)
                                      yin(48) = yin(48) + dyij*yin(42)
                                      zin(48) = zin(48) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   54

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   48

                                      xin(54) = xin(54) + dxij*xin(48)
                                      yin(54) = yin(54) + dyij*yin(48)
                                      zin(54) = zin(54) + dzij*zin(48)

                                      ! i3 = i4 =   48
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   12

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   18

                                      ! nj =    2

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   54

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

                                      ! i1 = in(1) =   55

                                      xin(55) = 1.0_dp
                                      yin(55) = 1.0_dp
                                      zin(55) = f00

                                      ! i2 = in(2) =   73
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(73) = xc00
                                      yin(73) = yc00
                                      zin(73) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   57

                                      xin(57) = xcp00
                                      yin(57) = ycp00
                                      zin(57) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   75
                                      ! i2 =   73

                                      xin(75) = xcp00*xin(73) + cp10
                                      yin(75) = ycp00*yin(73) + cp10
                                      zin(75) = zcp00*zin(73) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   55
                                      ! i4 = i2 =   73

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   55
                                      ! i4 =   73

                                      xin(91) = c10*xin(55) + xc00*xin(73)
                                      yin(91) = c10*yin(55) + yc00*yin(73)
                                      zin(91) = c10*zin(55) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   93
                                      ! i5 =   91
                                      ! i4 =   73

                                      xin(93) = xcp00*xin(91) + cp10*xin(73)
                                      yin(93) = ycp00*yin(91) + cp10*yin(73)
                                      zin(93) = zcp00*zin(91) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   91

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   97
                                      ! i3 =   73
                                      ! i4 =   91

                                      xin(97) = c10*xin(73) + xc00*xin(91)
                                      yin(97) = c10*yin(73) + yc00*yin(91)
                                      zin(97) = c10*zin(73) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   99
                                      ! i5 =   97
                                      ! i4 =   91

                                      xin(99) = xcp00*xin(97) + cp10*xin(91)
                                      yin(99) = ycp00*yin(97) + cp10*yin(91)
                                      zin(99) = zcp00*zin(97) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   97

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  103
                                      ! i3 =   91
                                      ! i4 =   97

                                      xin(103) = c10*xin(91) + xc00*xin(97)
                                      yin(103) = c10*yin(91) + yc00*yin(97)
                                      zin(103) = c10*zin(91) + zc00*zin(97)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  105
                                      ! i5 =  103
                                      ! i4 =   97

                                      xin(105) = xcp00*xin(103) + cp10*xin(97)
                                      yin(105) = ycp00*yin(103) + cp10*yin(97)
                                      zin(105) = zcp00*zin(103) + cp10*zin(97)

                                      ! ------------------

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  103

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   55
                                      ! i4 = i1+k2 =   57

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   59
                                      ! i3 =   55
                                      ! i4 =   57

                                      xin(59) = cp01*xin(55) + xcp00*xin(57)
                                      yin(59) = cp01*yin(55) + ycp00*yin(57)
                                      zin(59) = cp01*zin(55) + zcp00*zin(57)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   77

                                      xin(77) = xc00*xin(59) + c01*xin(57)
                                      yin(77) = yc00*yin(59) + c01*yin(57)
                                      zin(77) = zc00*zin(59) + c01*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   59

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   60
                                      ! i3 =   57
                                      ! i4 =   59

                                      xin(60) = cp01*xin(57) + xcp00*xin(59)
                                      yin(60) = cp01*yin(57) + ycp00*yin(59)
                                      zin(60) = cp01*zin(57) + zcp00*zin(59)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   78

                                      xin(78) = xc00*xin(60) + c01*xin(59)
                                      yin(78) = yc00*yin(60) + c01*yin(59)
                                      zin(78) = zc00*zin(60) + c01*zin(59)

                                      ! ------------------

                                      ! i3 = i4 =   59
                                      ! i4 = i5 =   60

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   55
                                      ! i4 = i2 =   73

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   91

                                      xin(95) = c10*xin(59) + xc00*xin(77) + c01*xin(75)
                                      yin(95) = c10*yin(59) + yc00*yin(77) + c01*yin(75)
                                      zin(95) = c10*zin(59) + zc00*zin(77) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   91

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   97

                                      xin(101) = c10*xin(77) + xc00*xin(95) + c01*xin(93)
                                      yin(101) = c10*yin(77) + yc00*yin(95) + c01*yin(93)
                                      zin(101) = c10*zin(77) + zc00*zin(95) + c01*zin(93)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   97

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  103

                                      xin(107) = c10*xin(95) + xc00*xin(101) + c01*xin(99)
                                      yin(107) = c10*yin(95) + yc00*yin(101) + c01*yin(99)
                                      zin(107) = c10*zin(95) + zc00*zin(101) + c01*zin(99)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  103

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =   55
                                      ! i4 = i2 =   73

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   91

                                      xin(96) = c10*xin(60) + xc00*xin(78) + c01*xin(77)
                                      yin(96) = c10*yin(60) + yc00*yin(78) + c01*yin(77)
                                      zin(96) = c10*zin(60) + zc00*zin(78) + c01*zin(77)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   91

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   97

                                      xin(102) = c10*xin(78) + xc00*xin(96) + c01*xin(95)
                                      yin(102) = c10*yin(78) + yc00*yin(96) + c01*yin(95)
                                      zin(102) = c10*zin(78) + zc00*zin(96) + c01*zin(95)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =   97

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  103

                                      xin(108) = c10*xin(96) + xc00*xin(102) + c01*xin(101)
                                      yin(108) = c10*yin(96) + yc00*yin(102) + c01*yin(101)
                                      zin(108) = c10*zin(96) + zc00*zin(102) + c01*zin(101)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  103

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  103

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  103

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   97

                                      xin(103) = xin(103) + dxij*xin(97)
                                      yin(103) = yin(103) + dyij*yin(97)
                                      zin(103) = zin(103) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   91

                                      xin(97) = xin(97) + dxij*xin(91)
                                      yin(97) = yin(97) + dyij*yin(91)
                                      zin(97) = zin(97) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  103

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   97

                                      xin(103) = xin(103) + dxij*xin(97)
                                      yin(103) = yin(103) + dyij*yin(97)
                                      zin(103) = zin(103) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   61

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   61

                                      ! do ni = 1,    2

                                      xin(61) = xin(73) + dxij*xin(55)
                                      yin(61) = yin(73) + dyij*yin(55)
                                      zin(61) = zin(73) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    2

                                      xin(79) = xin(91) + dxij*xin(73)
                                      yin(79) = yin(91) + dyij*yin(73)
                                      zin(79) = zin(91) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   67

                                      ! nj =    2

                                      ! i4 = i3 =   67

                                      ! do ni = 1,    2

                                      xin(67) = xin(79) + dxij*xin(61)
                                      yin(67) = yin(79) + dyij*yin(61)
                                      zin(67) = zin(79) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    2

                                      xin(85) = xin(97) + dxij*xin(79)
                                      yin(85) = yin(97) + dyij*yin(79)
                                      zin(85) = zin(97) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   73

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  105

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   99

                                      xin(105) = xin(105) + dxij*xin(99)
                                      yin(105) = yin(105) + dyij*yin(99)
                                      zin(105) = zin(105) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   93

                                      xin(99) = xin(99) + dxij*xin(93)
                                      yin(99) = yin(99) + dyij*yin(93)
                                      zin(99) = zin(99) + dzij*zin(93)

                                      ! i3 = i4 =   93
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  105

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   99

                                      xin(105) = xin(105) + dxij*xin(99)
                                      yin(105) = yin(105) + dyij*yin(99)
                                      zin(105) = zin(105) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   63

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   63

                                      ! do ni = 1,    2

                                      xin(63) = xin(75) + dxij*xin(57)
                                      yin(63) = yin(75) + dyij*yin(57)
                                      zin(63) = zin(75) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    2

                                      xin(81) = xin(93) + dxij*xin(75)
                                      yin(81) = yin(93) + dyij*yin(75)
                                      zin(81) = zin(93) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   99

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   69

                                      ! nj =    2

                                      ! i4 = i3 =   69

                                      ! do ni = 1,    2

                                      xin(69) = xin(81) + dxij*xin(63)
                                      yin(69) = yin(81) + dyij*yin(63)
                                      zin(69) = zin(81) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    2

                                      xin(87) = xin(99) + dxij*xin(81)
                                      yin(87) = yin(99) + dyij*yin(81)
                                      zin(87) = zin(99) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   75

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  101

                                      xin(107) = xin(107) + dxij*xin(101)
                                      yin(107) = yin(107) + dyij*yin(101)
                                      zin(107) = zin(107) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   95

                                      xin(101) = xin(101) + dxij*xin(95)
                                      yin(101) = yin(101) + dyij*yin(95)
                                      zin(101) = zin(101) + dzij*zin(95)

                                      ! i3 = i4 =   95
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  101

                                      xin(107) = xin(107) + dxij*xin(101)
                                      yin(107) = yin(107) + dyij*yin(101)
                                      zin(107) = zin(107) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   65

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   65

                                      ! do ni = 1,    2

                                      xin(65) = xin(77) + dxij*xin(59)
                                      yin(65) = yin(77) + dyij*yin(59)
                                      zin(65) = zin(77) + dzij*zin(59)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    2

                                      xin(83) = xin(95) + dxij*xin(77)
                                      yin(83) = yin(95) + dyij*yin(77)
                                      zin(83) = zin(95) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  101

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   71

                                      ! nj =    2

                                      ! i4 = i3 =   71

                                      ! do ni = 1,    2

                                      xin(71) = xin(83) + dxij*xin(65)
                                      yin(71) = yin(83) + dyij*yin(65)
                                      zin(71) = zin(83) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    2

                                      xin(89) = xin(101) + dxij*xin(83)
                                      yin(89) = yin(101) + dyij*yin(83)
                                      zin(89) = zin(101) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   77

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  102

                                      xin(108) = xin(108) + dxij*xin(102)
                                      yin(108) = yin(108) + dyij*yin(102)
                                      zin(108) = zin(108) + dzij*zin(102)

                                      ! i3 = i4 =  102
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   96

                                      xin(102) = xin(102) + dxij*xin(96)
                                      yin(102) = yin(102) + dyij*yin(96)
                                      zin(102) = zin(102) + dzij*zin(96)

                                      ! i3 = i4 =   96
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  102

                                      xin(108) = xin(108) + dxij*xin(102)
                                      yin(108) = yin(108) + dyij*yin(102)
                                      zin(108) = zin(108) + dzij*zin(102)

                                      ! i3 = i4 =  102
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   66

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   66

                                      ! do ni = 1,    2

                                      xin(66) = xin(78) + dxij*xin(60)
                                      yin(66) = yin(78) + dyij*yin(60)
                                      zin(66) = zin(78) + dzij*zin(60)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    2

                                      xin(84) = xin(96) + dxij*xin(78)
                                      yin(84) = yin(96) + dyij*yin(78)
                                      zin(84) = zin(96) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  102

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   72

                                      ! nj =    2

                                      ! i4 = i3 =   72

                                      ! do ni = 1,    2

                                      xin(72) = xin(84) + dxij*xin(66)
                                      yin(72) = yin(84) + dyij*yin(66)
                                      zin(72) = zin(84) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    2

                                      xin(90) = xin(102) + dxij*xin(84)
                                      yin(90) = yin(102) + dyij*yin(84)
                                      zin(90) = zin(102) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   78

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =   55

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! ni = ni + 1 =    1

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

                                      ! end do

                                      ! ni = ni + 1 =    2

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

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  108

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

                                      ! i1 = in(1) =  109

                                      xin(109) = 1.0_dp
                                      yin(109) = 1.0_dp
                                      zin(109) = f00

                                      ! i2 = in(2) =  127
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(127) = xc00
                                      yin(127) = yc00
                                      zin(127) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  111

                                      xin(111) = xcp00
                                      yin(111) = ycp00
                                      zin(111) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  129
                                      ! i2 =  127

                                      xin(129) = xcp00*xin(127) + cp10
                                      yin(129) = ycp00*yin(127) + cp10
                                      zin(129) = zcp00*zin(127) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  127

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  145
                                      ! i3 =  109
                                      ! i4 =  127

                                      xin(145) = c10*xin(109) + xc00*xin(127)
                                      yin(145) = c10*yin(109) + yc00*yin(127)
                                      zin(145) = c10*zin(109) + zc00*zin(127)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  147
                                      ! i5 =  145
                                      ! i4 =  127

                                      xin(147) = xcp00*xin(145) + cp10*xin(127)
                                      yin(147) = ycp00*yin(145) + cp10*yin(127)
                                      zin(147) = zcp00*zin(145) + cp10*zin(127)

                                      ! ------------------

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  145

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  151
                                      ! i3 =  127
                                      ! i4 =  145

                                      xin(151) = c10*xin(127) + xc00*xin(145)
                                      yin(151) = c10*yin(127) + yc00*yin(145)
                                      zin(151) = c10*zin(127) + zc00*zin(145)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  153
                                      ! i5 =  151
                                      ! i4 =  145

                                      xin(153) = xcp00*xin(151) + cp10*xin(145)
                                      yin(153) = ycp00*yin(151) + cp10*yin(145)
                                      zin(153) = zcp00*zin(151) + cp10*zin(145)

                                      ! ------------------

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  151

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  157
                                      ! i3 =  145
                                      ! i4 =  151

                                      xin(157) = c10*xin(145) + xc00*xin(151)
                                      yin(157) = c10*yin(145) + yc00*yin(151)
                                      zin(157) = c10*zin(145) + zc00*zin(151)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  159
                                      ! i5 =  157
                                      ! i4 =  151

                                      xin(159) = xcp00*xin(157) + cp10*xin(151)
                                      yin(159) = ycp00*yin(157) + cp10*yin(151)
                                      zin(159) = zcp00*zin(157) + cp10*zin(151)

                                      ! ------------------

                                      ! i3 = i4 =  151
                                      ! i4 = i5 =  157

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  109
                                      ! i4 = i1+k2 =  111

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  113
                                      ! i3 =  109
                                      ! i4 =  111

                                      xin(113) = cp01*xin(109) + xcp00*xin(111)
                                      yin(113) = cp01*yin(109) + ycp00*yin(111)
                                      zin(113) = cp01*zin(109) + zcp00*zin(111)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  131

                                      xin(131) = xc00*xin(113) + c01*xin(111)
                                      yin(131) = yc00*yin(113) + c01*yin(111)
                                      zin(131) = zc00*zin(113) + c01*zin(111)

                                      ! ------------------

                                      ! i3 = i4 =  111
                                      ! i4 = i5 =  113

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  114
                                      ! i3 =  111
                                      ! i4 =  113

                                      xin(114) = cp01*xin(111) + xcp00*xin(113)
                                      yin(114) = cp01*yin(111) + ycp00*yin(113)
                                      zin(114) = cp01*zin(111) + zcp00*zin(113)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  132

                                      xin(132) = xc00*xin(114) + c01*xin(113)
                                      yin(132) = yc00*yin(114) + c01*yin(113)
                                      zin(132) = zc00*zin(114) + c01*zin(113)

                                      ! ------------------

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  114

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  127

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  145

                                      xin(149) = c10*xin(113) + xc00*xin(131) + c01*xin(129)
                                      yin(149) = c10*yin(113) + yc00*yin(131) + c01*yin(129)
                                      zin(149) = c10*zin(113) + zc00*zin(131) + c01*zin(129)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  145

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  151

                                      xin(155) = c10*xin(131) + xc00*xin(149) + c01*xin(147)
                                      yin(155) = c10*yin(131) + yc00*yin(149) + c01*yin(147)
                                      zin(155) = c10*zin(131) + zc00*zin(149) + c01*zin(147)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  151

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  157

                                      xin(161) = c10*xin(149) + xc00*xin(155) + c01*xin(153)
                                      yin(161) = c10*yin(149) + yc00*yin(155) + c01*yin(153)
                                      zin(161) = c10*zin(149) + zc00*zin(155) + c01*zin(153)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  151
                                      ! i4 = i5 =  157

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  127

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  145

                                      xin(150) = c10*xin(114) + xc00*xin(132) + c01*xin(131)
                                      yin(150) = c10*yin(114) + yc00*yin(132) + c01*yin(131)
                                      zin(150) = c10*zin(114) + zc00*zin(132) + c01*zin(131)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  127
                                      ! i4 = i5 =  145

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  151

                                      xin(156) = c10*xin(132) + xc00*xin(150) + c01*xin(149)
                                      yin(156) = c10*yin(132) + yc00*yin(150) + c01*yin(149)
                                      zin(156) = c10*zin(132) + zc00*zin(150) + c01*zin(149)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  145
                                      ! i4 = i5 =  151

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  157

                                      xin(162) = c10*xin(150) + xc00*xin(156) + c01*xin(155)
                                      yin(162) = c10*yin(150) + yc00*yin(156) + c01*yin(155)
                                      zin(162) = c10*zin(150) + zc00*zin(156) + c01*zin(155)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  151
                                      ! i4 = i5 =  157

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  157

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  157

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  151

                                      xin(157) = xin(157) + dxij*xin(151)
                                      yin(157) = yin(157) + dyij*yin(151)
                                      zin(157) = zin(157) + dzij*zin(151)

                                      ! i3 = i4 =  151
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  145

                                      xin(151) = xin(151) + dxij*xin(145)
                                      yin(151) = yin(151) + dyij*yin(145)
                                      zin(151) = zin(151) + dzij*zin(145)

                                      ! i3 = i4 =  145
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  157

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  151

                                      xin(157) = xin(157) + dxij*xin(151)
                                      yin(157) = yin(157) + dyij*yin(151)
                                      zin(157) = zin(157) + dzij*zin(151)

                                      ! i3 = i4 =  151
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  115

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  115

                                      ! do ni = 1,    2

                                      xin(115) = xin(127) + dxij*xin(109)
                                      yin(115) = yin(127) + dyij*yin(109)
                                      zin(115) = zin(127) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    2

                                      xin(133) = xin(145) + dxij*xin(127)
                                      yin(133) = yin(145) + dyij*yin(127)
                                      zin(133) = zin(145) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  151

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  121

                                      ! nj =    2

                                      ! i4 = i3 =  121

                                      ! do ni = 1,    2

                                      xin(121) = xin(133) + dxij*xin(115)
                                      yin(121) = yin(133) + dyij*yin(115)
                                      zin(121) = zin(133) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    2

                                      xin(139) = xin(151) + dxij*xin(133)
                                      yin(139) = yin(151) + dyij*yin(133)
                                      zin(139) = zin(151) + dzij*zin(133)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  127

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  159

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  153

                                      xin(159) = xin(159) + dxij*xin(153)
                                      yin(159) = yin(159) + dyij*yin(153)
                                      zin(159) = zin(159) + dzij*zin(153)

                                      ! i3 = i4 =  153
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  147

                                      xin(153) = xin(153) + dxij*xin(147)
                                      yin(153) = yin(153) + dyij*yin(147)
                                      zin(153) = zin(153) + dzij*zin(147)

                                      ! i3 = i4 =  147
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  159

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  153

                                      xin(159) = xin(159) + dxij*xin(153)
                                      yin(159) = yin(159) + dyij*yin(153)
                                      zin(159) = zin(159) + dzij*zin(153)

                                      ! i3 = i4 =  153
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  117

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  117

                                      ! do ni = 1,    2

                                      xin(117) = xin(129) + dxij*xin(111)
                                      yin(117) = yin(129) + dyij*yin(111)
                                      zin(117) = zin(129) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    2

                                      xin(135) = xin(147) + dxij*xin(129)
                                      yin(135) = yin(147) + dyij*yin(129)
                                      zin(135) = zin(147) + dzij*zin(129)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  153

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  123

                                      ! nj =    2

                                      ! i4 = i3 =  123

                                      ! do ni = 1,    2

                                      xin(123) = xin(135) + dxij*xin(117)
                                      yin(123) = yin(135) + dyij*yin(117)
                                      zin(123) = zin(135) + dzij*zin(117)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    2

                                      xin(141) = xin(153) + dxij*xin(135)
                                      yin(141) = yin(153) + dyij*yin(135)
                                      zin(141) = zin(153) + dzij*zin(135)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  159

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  129

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  161

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  155

                                      xin(161) = xin(161) + dxij*xin(155)
                                      yin(161) = yin(161) + dyij*yin(155)
                                      zin(161) = zin(161) + dzij*zin(155)

                                      ! i3 = i4 =  155
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  149

                                      xin(155) = xin(155) + dxij*xin(149)
                                      yin(155) = yin(155) + dyij*yin(149)
                                      zin(155) = zin(155) + dzij*zin(149)

                                      ! i3 = i4 =  149
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  161

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  155

                                      xin(161) = xin(161) + dxij*xin(155)
                                      yin(161) = yin(161) + dyij*yin(155)
                                      zin(161) = zin(161) + dzij*zin(155)

                                      ! i3 = i4 =  155
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  119

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  119

                                      ! do ni = 1,    2

                                      xin(119) = xin(131) + dxij*xin(113)
                                      yin(119) = yin(131) + dyij*yin(113)
                                      zin(119) = zin(131) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    2

                                      xin(137) = xin(149) + dxij*xin(131)
                                      yin(137) = yin(149) + dyij*yin(131)
                                      zin(137) = zin(149) + dzij*zin(131)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  155

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  125

                                      ! nj =    2

                                      ! i4 = i3 =  125

                                      ! do ni = 1,    2

                                      xin(125) = xin(137) + dxij*xin(119)
                                      yin(125) = yin(137) + dyij*yin(119)
                                      zin(125) = zin(137) + dzij*zin(119)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    2

                                      xin(143) = xin(155) + dxij*xin(137)
                                      yin(143) = yin(155) + dyij*yin(137)
                                      zin(143) = zin(155) + dzij*zin(137)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  131

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  162

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  156

                                      xin(162) = xin(162) + dxij*xin(156)
                                      yin(162) = yin(162) + dyij*yin(156)
                                      zin(162) = zin(162) + dzij*zin(156)

                                      ! i3 = i4 =  156
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  150

                                      xin(156) = xin(156) + dxij*xin(150)
                                      yin(156) = yin(156) + dyij*yin(150)
                                      zin(156) = zin(156) + dzij*zin(150)

                                      ! i3 = i4 =  150
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  162

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  156

                                      xin(162) = xin(162) + dxij*xin(156)
                                      yin(162) = yin(162) + dyij*yin(156)
                                      zin(162) = zin(162) + dzij*zin(156)

                                      ! i3 = i4 =  156
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  120

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  120

                                      ! do ni = 1,    2

                                      xin(120) = xin(132) + dxij*xin(114)
                                      yin(120) = yin(132) + dyij*yin(114)
                                      zin(120) = zin(132) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    2

                                      xin(138) = xin(150) + dxij*xin(132)
                                      yin(138) = yin(150) + dyij*yin(132)
                                      zin(138) = zin(150) + dzij*zin(132)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  156

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  126

                                      ! nj =    2

                                      ! i4 = i3 =  126

                                      ! do ni = 1,    2

                                      xin(126) = xin(138) + dxij*xin(120)
                                      yin(126) = yin(138) + dyij*yin(120)
                                      zin(126) = zin(138) + dzij*zin(120)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    2

                                      xin(144) = xin(156) + dxij*xin(138)
                                      yin(144) = yin(156) + dyij*yin(138)
                                      zin(144) = zin(156) + dzij*zin(138)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  162

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  132

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  109

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! ni = ni + 1 =    1

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

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  162

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

                                      ! i1 = in(1) =  163

                                      xin(163) = 1.0_dp
                                      yin(163) = 1.0_dp
                                      zin(163) = f00

                                      ! i2 = in(2) =  181
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(181) = xc00
                                      yin(181) = yc00
                                      zin(181) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  165

                                      xin(165) = xcp00
                                      yin(165) = ycp00
                                      zin(165) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  183
                                      ! i2 =  181

                                      xin(183) = xcp00*xin(181) + cp10
                                      yin(183) = ycp00*yin(181) + cp10
                                      zin(183) = zcp00*zin(181) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  163
                                      ! i4 = i2 =  181

                                      ! do n = 2,   4

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

                                      ! n =    3

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

                                      ! n =    4

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

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  163
                                      ! i4 = i1+k2 =  165

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  167
                                      ! i3 =  163
                                      ! i4 =  165

                                      xin(167) = cp01*xin(163) + xcp00*xin(165)
                                      yin(167) = cp01*yin(163) + ycp00*yin(165)
                                      zin(167) = cp01*zin(163) + zcp00*zin(165)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  185

                                      xin(185) = xc00*xin(167) + c01*xin(165)
                                      yin(185) = yc00*yin(167) + c01*yin(165)
                                      zin(185) = zc00*zin(167) + c01*zin(165)

                                      ! ------------------

                                      ! i3 = i4 =  165
                                      ! i4 = i5 =  167

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  168
                                      ! i3 =  165
                                      ! i4 =  167

                                      xin(168) = cp01*xin(165) + xcp00*xin(167)
                                      yin(168) = cp01*yin(165) + ycp00*yin(167)
                                      zin(168) = cp01*zin(165) + zcp00*zin(167)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  186

                                      xin(186) = xc00*xin(168) + c01*xin(167)
                                      yin(186) = yc00*yin(168) + c01*yin(167)
                                      zin(186) = zc00*zin(168) + c01*zin(167)

                                      ! ------------------

                                      ! i3 = i4 =  167
                                      ! i4 = i5 =  168

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  163
                                      ! i4 = i2 =  181

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  199

                                      xin(203) = c10*xin(167) + xc00*xin(185) + c01*xin(183)
                                      yin(203) = c10*yin(167) + yc00*yin(185) + c01*yin(183)
                                      zin(203) = c10*zin(167) + zc00*zin(185) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  199

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  205

                                      xin(209) = c10*xin(185) + xc00*xin(203) + c01*xin(201)
                                      yin(209) = c10*yin(185) + yc00*yin(203) + c01*yin(201)
                                      zin(209) = c10*zin(185) + zc00*zin(203) + c01*zin(201)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  205

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  211

                                      xin(215) = c10*xin(203) + xc00*xin(209) + c01*xin(207)
                                      yin(215) = c10*yin(203) + yc00*yin(209) + c01*yin(207)
                                      zin(215) = c10*zin(203) + zc00*zin(209) + c01*zin(207)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  211

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  163
                                      ! i4 = i2 =  181

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  199

                                      xin(204) = c10*xin(168) + xc00*xin(186) + c01*xin(185)
                                      yin(204) = c10*yin(168) + yc00*yin(186) + c01*yin(185)
                                      zin(204) = c10*zin(168) + zc00*zin(186) + c01*zin(185)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  199

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  205

                                      xin(210) = c10*xin(186) + xc00*xin(204) + c01*xin(203)
                                      yin(210) = c10*yin(186) + yc00*yin(204) + c01*yin(203)
                                      zin(210) = c10*zin(186) + zc00*zin(204) + c01*zin(203)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  205

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  211

                                      xin(216) = c10*xin(204) + xc00*xin(210) + c01*xin(209)
                                      yin(216) = c10*yin(204) + yc00*yin(210) + c01*yin(209)
                                      zin(216) = c10*zin(204) + zc00*zin(210) + c01*zin(209)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  211

                                      ! nn =    5

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

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(211) = xin(211) + dxij*xin(205)
                                      yin(211) = yin(211) + dyij*yin(205)
                                      zin(211) = zin(211) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  199

                                      xin(205) = xin(205) + dxij*xin(199)
                                      yin(205) = yin(205) + dyij*yin(199)
                                      zin(205) = zin(205) + dzij*zin(199)

                                      ! i3 = i4 =  199
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(211) = xin(211) + dxij*xin(205)
                                      yin(211) = yin(211) + dyij*yin(205)
                                      zin(211) = zin(211) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  169

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  169

                                      ! do ni = 1,    2

                                      xin(169) = xin(181) + dxij*xin(163)
                                      yin(169) = yin(181) + dyij*yin(163)
                                      zin(169) = zin(181) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    2

                                      xin(187) = xin(199) + dxij*xin(181)
                                      yin(187) = yin(199) + dyij*yin(181)
                                      zin(187) = zin(199) + dzij*zin(181)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  175

                                      ! nj =    2

                                      ! i4 = i3 =  175

                                      ! do ni = 1,    2

                                      xin(175) = xin(187) + dxij*xin(169)
                                      yin(175) = yin(187) + dyij*yin(169)
                                      zin(175) = zin(187) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! ni =    2

                                      xin(193) = xin(205) + dxij*xin(187)
                                      yin(193) = yin(205) + dyij*yin(187)
                                      zin(193) = zin(205) + dzij*zin(187)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  181

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(213) = xin(213) + dxij*xin(207)
                                      yin(213) = yin(213) + dyij*yin(207)
                                      zin(213) = zin(213) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  201

                                      xin(207) = xin(207) + dxij*xin(201)
                                      yin(207) = yin(207) + dyij*yin(201)
                                      zin(207) = zin(207) + dzij*zin(201)

                                      ! i3 = i4 =  201
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(213) = xin(213) + dxij*xin(207)
                                      yin(213) = yin(213) + dyij*yin(207)
                                      zin(213) = zin(213) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  171

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  171

                                      ! do ni = 1,    2

                                      xin(171) = xin(183) + dxij*xin(165)
                                      yin(171) = yin(183) + dyij*yin(165)
                                      zin(171) = zin(183) + dzij*zin(165)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    2

                                      xin(189) = xin(201) + dxij*xin(183)
                                      yin(189) = yin(201) + dyij*yin(183)
                                      zin(189) = zin(201) + dzij*zin(183)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  177

                                      ! nj =    2

                                      ! i4 = i3 =  177

                                      ! do ni = 1,    2

                                      xin(177) = xin(189) + dxij*xin(171)
                                      yin(177) = yin(189) + dyij*yin(171)
                                      zin(177) = zin(189) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  195

                                      ! ni =    2

                                      xin(195) = xin(207) + dxij*xin(189)
                                      yin(195) = yin(207) + dyij*yin(189)
                                      zin(195) = zin(207) + dzij*zin(189)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  183

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  209

                                      xin(215) = xin(215) + dxij*xin(209)
                                      yin(215) = yin(215) + dyij*yin(209)
                                      zin(215) = zin(215) + dzij*zin(209)

                                      ! i3 = i4 =  209
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  203

                                      xin(209) = xin(209) + dxij*xin(203)
                                      yin(209) = yin(209) + dyij*yin(203)
                                      zin(209) = zin(209) + dzij*zin(203)

                                      ! i3 = i4 =  203
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  209

                                      xin(215) = xin(215) + dxij*xin(209)
                                      yin(215) = yin(215) + dyij*yin(209)
                                      zin(215) = zin(215) + dzij*zin(209)

                                      ! i3 = i4 =  209
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  173

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  173

                                      ! do ni = 1,    2

                                      xin(173) = xin(185) + dxij*xin(167)
                                      yin(173) = yin(185) + dyij*yin(167)
                                      zin(173) = zin(185) + dzij*zin(167)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    2

                                      xin(191) = xin(203) + dxij*xin(185)
                                      yin(191) = yin(203) + dyij*yin(185)
                                      zin(191) = zin(203) + dzij*zin(185)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  179

                                      ! nj =    2

                                      ! i4 = i3 =  179

                                      ! do ni = 1,    2

                                      xin(179) = xin(191) + dxij*xin(173)
                                      yin(179) = yin(191) + dyij*yin(173)
                                      zin(179) = zin(191) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  197

                                      ! ni =    2

                                      xin(197) = xin(209) + dxij*xin(191)
                                      yin(197) = yin(209) + dyij*yin(191)
                                      zin(197) = zin(209) + dzij*zin(191)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  185

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  210

                                      xin(216) = xin(216) + dxij*xin(210)
                                      yin(216) = yin(216) + dyij*yin(210)
                                      zin(216) = zin(216) + dzij*zin(210)

                                      ! i3 = i4 =  210
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  204

                                      xin(210) = xin(210) + dxij*xin(204)
                                      yin(210) = yin(210) + dyij*yin(204)
                                      zin(210) = zin(210) + dzij*zin(204)

                                      ! i3 = i4 =  204
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  210

                                      xin(216) = xin(216) + dxij*xin(210)
                                      yin(216) = yin(216) + dyij*yin(210)
                                      zin(216) = zin(216) + dzij*zin(210)

                                      ! i3 = i4 =  210
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  174

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  174

                                      ! do ni = 1,    2

                                      xin(174) = xin(186) + dxij*xin(168)
                                      yin(174) = yin(186) + dyij*yin(168)
                                      zin(174) = zin(186) + dzij*zin(168)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    2

                                      xin(192) = xin(204) + dxij*xin(186)
                                      yin(192) = yin(204) + dyij*yin(186)
                                      zin(192) = zin(204) + dzij*zin(186)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  210

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  180

                                      ! nj =    2

                                      ! i4 = i3 =  180

                                      ! do ni = 1,    2

                                      xin(180) = xin(192) + dxij*xin(174)
                                      yin(180) = yin(192) + dyij*yin(174)
                                      zin(180) = zin(192) + dzij*zin(174)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  198

                                      ! ni =    2

                                      xin(198) = xin(210) + dxij*xin(192)
                                      yin(198) = yin(210) + dyij*yin(192)
                                      zin(198) = zin(210) + dzij*zin(192)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  186

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =  163

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! ni = ni + 1 =    1

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

                                      ! ni = ni + 1 =    2

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

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  216

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 648! loop over all integrals

                                        l = n - 18*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d22bra(j)*d12ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 54)*yin(my + 54)*zin(mz + 54) & ! root  2
                                                        + xin(mx + 108)*yin(my + 108)*zin(mz + 108) & ! root  3
                                                        + xin(mx + 162)*yin(my + 162)*zin(mz + 162)) ! root  4

                                        j = int(n/18) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxj2 = 6
                                  iandj = ish .eq. jsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 6 ! # of cartesians in i

                                    if (iandj) maxj2 = i

                                    ii1 = i + loci
                                    ip = (i - 1)*108 ! Stride between functions in i

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

                              deallocate (n22bra)
                              deallocate (xint22bra)
                              deallocate (n12ket)
                              deallocate (xint12ket)

                              end subroutine int2221
                              end submodule
