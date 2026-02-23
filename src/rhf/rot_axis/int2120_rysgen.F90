! The total angular momentum of this class is:           5
! The algorithm chosen is: Rys quadrature
submodule(rot_axis_kernels) int2120_impl
contains
  module subroutine int2120(pd_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pd_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n12bra(:), n02ket(:)
    real(dp), allocatable :: xint12bra(:), xint02ket(:)
    integer(kind=int64) :: npdbra, nsdket
    real(dp) :: scutpdbra, scutsdket, test
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
    real(dp) :: roots(3), wghts(3)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(3), wts(3), alpha(3), beta(3), wrk(3)
    real(dp) :: xin(54), yin(54), zin(54)
    real(dp) :: eri_value(108)
    real(dp) :: d12bra(18), d02ket(6)
    integer(kind=int64) :: ix(6), jx(3), kx(6), lx(1)
    integer(kind=int64) :: iy(6), jy(3), ky(6), ly(1)
    integer(kind=int64) :: iz(6), jz(3), kz(6), lz(1)
    integer(kind=int64) :: in(4), in1(4), kn(3)
    integer(kind=int64) :: ijx(18), ijy(18), ijz(18)
    integer(kind=int64) :: klx(6), kly(6), klz(6)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 7
    in1(3) = 13
    in1(4) = 16

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

    jx(1) = 3
    jx(2) = 0
    jx(3) = 0

    ix(1) = 13
    ix(2) = 1
    ix(3) = 1
    ix(4) = 7
    ix(5) = 7
    ix(6) = 1

    ! y-arrays

    ly(1) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 1

    jy(1) = 0
    jy(2) = 3
    jy(3) = 0

    iy(1) = 1
    iy(2) = 13
    iy(3) = 1
    iy(4) = 7
    iy(5) = 1
    iy(6) = 7

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
    jz(3) = 3

    iz(1) = 1
    iz(2) = 1
    iz(3) = 13
    iz(4) = 1
    iz(5) = 7
    iz(6) = 7

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 16
    ijx(2) = 13
    ijx(3) = 13
    ijx(4) = 4
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 4
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 10
    ijx(11) = 7
    ijx(12) = 7
    ijx(13) = 10
    ijx(14) = 7
    ijx(15) = 7
    ijx(16) = 4
    ijx(17) = 1
    ijx(18) = 1

    ijy(1) = 1
    ijy(2) = 4
    ijy(3) = 1
    ijy(4) = 13
    ijy(5) = 16
    ijy(6) = 13
    ijy(7) = 1
    ijy(8) = 4
    ijy(9) = 1
    ijy(10) = 7
    ijy(11) = 10
    ijy(12) = 7
    ijy(13) = 1
    ijy(14) = 4
    ijy(15) = 1
    ijy(16) = 7
    ijy(17) = 10
    ijy(18) = 7

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 4
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 4
    ijz(7) = 13
    ijz(8) = 13
    ijz(9) = 16
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 4
    ijz(13) = 7
    ijz(14) = 7
    ijz(15) = 10
    ijz(16) = 7
    ijz(17) = 7
    ijz(18) = 10

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

    allocate (n12bra(res%n_p_shl*res%n_d_shl))
    allocate (xint12bra(res%n_p_shl*res%n_d_shl))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))

    ! Start screening

    scutpdbra = cutoff_schwarz/maxval(pd_pair%xints)
    npdbra = 0
    do ij = 1, res%n_p_shl*res%n_d_shl
      if (pd_pair%xints(ij) .ge. scutpdbra) then
        npdbra = npdbra + 1
        xint12bra(npdbra) = pd_pair%xints(ij)
        n12bra(npdbra) = ij
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

    if ((npdbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = npdbra*nsdket
    ntile = int(npdbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = npdbra*nsdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npdbra, xint12bra, n12bra, xint02ket, n02ket, pd_pair, sd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d02ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d12bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, npdbra) + 1
              kl_tmp = (iquart - 1)/npdbra + 1

              test = xint12bra(ij_tmp)*xint02ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n12bra(ij_tmp)
                kl = n02ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_d_shl) + 1
                jsh_tmp = (ij - 1)/res%n_d_shl + 1
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_d_shl(ish_tmp)
                jsh = res%i_p_shl(jsh_tmp)
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

                    t_expon_ab = pd_pair%t_expon_ab(pd_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = pd_pair%expon_b(pd_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = pd_pair%expon_a(pd_pair%pair_loc(ij) + bra_loop)
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

                    d12bra(1) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(2) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(3) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(4) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(5) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(6) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(7) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(8) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(9) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)
                    d12bra(10) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(11) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(12) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(13) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(14) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(15) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(16) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(17) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d12bra(18) = pd_pair%d_coeff_alt(pd_pair%pair_loc(ij) + bra_loop)*sqrt3

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

                                      ! do n = 2,   3

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

                                      ! i5 = in(n+1) =   16
                                      ! i3 =    7
                                      ! i4 =   13

                                      xin(16) = c10*xin(7) + xc00*xin(13)
                                      yin(16) = c10*yin(7) + yc00*yin(13)
                                      zin(16) = c10*zin(7) + zc00*zin(13)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   17
                                      ! i5 =   16
                                      ! i4 =   13

                                      xin(17) = xcp00*xin(16) + cp10*xin(13)
                                      yin(17) = ycp00*yin(16) + cp10*yin(13)
                                      zin(17) = zcp00*zin(16) + cp10*zin(13)

                                      ! ------------------

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   16

                                      ! n =    4

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

                                      ! i3 = i2+kn(n+1) =    9

                                      xin(9) = xc00*xin(3) + c01*xin(2)
                                      yin(9) = yc00*yin(3) + c01*yin(2)
                                      zin(9) = zc00*zin(3) + c01*zin(2)

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
                                      ! i4 = i2 =    7

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   13

                                      xin(15) = c10*xin(3) + xc00*xin(9) + c01*xin(8)
                                      yin(15) = c10*yin(3) + yc00*yin(9) + c01*yin(8)
                                      zin(15) = c10*zin(3) + zc00*zin(9) + c01*zin(8)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =   13

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   16

                                      xin(18) = c10*xin(9) + xc00*xin(15) + c01*xin(14)
                                      yin(18) = c10*yin(9) + yc00*yin(15) + c01*yin(14)
                                      zin(18) = c10*zin(9) + zc00*zin(15) + c01*zin(14)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   16

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   16

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   16

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   13

                                      xin(16) = xin(16) + dxij*xin(13)
                                      yin(16) = yin(16) + dyij*yin(13)
                                      zin(16) = zin(16) + dzij*zin(13)

                                      ! i3 = i4 =   13
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    4

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    4

                                      ! do ni = 1,    2

                                      xin(4) = xin(7) + dxij*xin(1)
                                      yin(4) = yin(7) + dyij*yin(1)
                                      zin(4) = zin(7) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   10

                                      ! ni =    2

                                      xin(10) = xin(13) + dxij*xin(7)
                                      yin(10) = yin(13) + dyij*yin(7)
                                      zin(10) = zin(13) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    7

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   17

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   14

                                      xin(17) = xin(17) + dxij*xin(14)
                                      yin(17) = yin(17) + dyij*yin(14)
                                      zin(17) = zin(17) + dzij*zin(14)

                                      ! i3 = i4 =   14
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    2

                                      xin(5) = xin(8) + dxij*xin(2)
                                      yin(5) = yin(8) + dyij*yin(2)
                                      zin(5) = zin(8) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   11

                                      ! ni =    2

                                      xin(11) = xin(14) + dxij*xin(8)
                                      yin(11) = yin(14) + dyij*yin(8)
                                      zin(11) = zin(14) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    8

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   18

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   15

                                      xin(18) = xin(18) + dxij*xin(15)
                                      yin(18) = yin(18) + dyij*yin(15)
                                      zin(18) = zin(18) + dzij*zin(15)

                                      ! i3 = i4 =   15
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    2

                                      xin(6) = xin(9) + dxij*xin(3)
                                      yin(6) = yin(9) + dyij*yin(3)
                                      zin(6) = zin(9) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   12

                                      ! ni =    2

                                      xin(12) = xin(15) + dxij*xin(9)
                                      yin(12) = yin(15) + dyij*yin(9)
                                      zin(12) = zin(15) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   18

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

                                      ! i1 = in(1) =   19

                                      xin(19) = 1.0_dp
                                      yin(19) = 1.0_dp
                                      zin(19) = f00

                                      ! i2 = in(2) =   25
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(25) = xc00
                                      yin(25) = yc00
                                      zin(25) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   20

                                      xin(20) = xcp00
                                      yin(20) = ycp00
                                      zin(20) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   26
                                      ! i2 =   25

                                      xin(26) = xcp00*xin(25) + cp10
                                      yin(26) = ycp00*yin(25) + cp10
                                      zin(26) = zcp00*zin(25) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   19
                                      ! i4 = i2 =   25

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   31
                                      ! i3 =   19
                                      ! i4 =   25

                                      xin(31) = c10*xin(19) + xc00*xin(25)
                                      yin(31) = c10*yin(19) + yc00*yin(25)
                                      zin(31) = c10*zin(19) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   32
                                      ! i5 =   31
                                      ! i4 =   25

                                      xin(32) = xcp00*xin(31) + cp10*xin(25)
                                      yin(32) = ycp00*yin(31) + cp10*yin(25)
                                      zin(32) = zcp00*zin(31) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   31

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   34
                                      ! i3 =   25
                                      ! i4 =   31

                                      xin(34) = c10*xin(25) + xc00*xin(31)
                                      yin(34) = c10*yin(25) + yc00*yin(31)
                                      zin(34) = c10*zin(25) + zc00*zin(31)

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

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   19
                                      ! i4 = i1+k2 =   20

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   21
                                      ! i3 =   19
                                      ! i4 =   20

                                      xin(21) = cp01*xin(19) + xcp00*xin(20)
                                      yin(21) = cp01*yin(19) + ycp00*yin(20)
                                      zin(21) = cp01*zin(19) + zcp00*zin(20)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   27

                                      xin(27) = xc00*xin(21) + c01*xin(20)
                                      yin(27) = yc00*yin(21) + c01*yin(20)
                                      zin(27) = zc00*zin(21) + c01*zin(20)

                                      ! ------------------

                                      ! i3 = i4 =   20
                                      ! i4 = i5 =   21

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   19
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   31

                                      xin(33) = c10*xin(21) + xc00*xin(27) + c01*xin(26)
                                      yin(33) = c10*yin(21) + yc00*yin(27) + c01*yin(26)
                                      zin(33) = c10*zin(21) + zc00*zin(27) + c01*zin(26)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   31

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   34

                                      xin(36) = c10*xin(27) + xc00*xin(33) + c01*xin(32)
                                      yin(36) = c10*yin(27) + yc00*yin(33) + c01*yin(32)
                                      zin(36) = c10*zin(27) + zc00*zin(33) + c01*zin(32)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   34

                                      ! nn =    4

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

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   34

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(34) = xin(34) + dxij*xin(31)
                                      yin(34) = yin(34) + dyij*yin(31)
                                      zin(34) = zin(34) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   22

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   22

                                      ! do ni = 1,    2

                                      xin(22) = xin(25) + dxij*xin(19)
                                      yin(22) = yin(25) + dyij*yin(19)
                                      zin(22) = zin(25) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   28

                                      ! ni =    2

                                      xin(28) = xin(31) + dxij*xin(25)
                                      yin(28) = yin(31) + dyij*yin(25)
                                      zin(28) = zin(31) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   35

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   32

                                      xin(35) = xin(35) + dxij*xin(32)
                                      yin(35) = yin(35) + dyij*yin(32)
                                      zin(35) = zin(35) + dzij*zin(32)

                                      ! i3 = i4 =   32
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   23

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   23

                                      ! do ni = 1,    2

                                      xin(23) = xin(26) + dxij*xin(20)
                                      yin(23) = yin(26) + dyij*yin(20)
                                      zin(23) = zin(26) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    2

                                      xin(29) = xin(32) + dxij*xin(26)
                                      yin(29) = yin(32) + dyij*yin(26)
                                      zin(29) = zin(32) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   26

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   33

                                      xin(36) = xin(36) + dxij*xin(33)
                                      yin(36) = yin(36) + dyij*yin(33)
                                      zin(36) = zin(36) + dzij*zin(33)

                                      ! i3 = i4 =   33
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   24

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   24

                                      ! do ni = 1,    2

                                      xin(24) = xin(27) + dxij*xin(21)
                                      yin(24) = yin(27) + dyij*yin(21)
                                      zin(24) = zin(27) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    2

                                      xin(30) = xin(33) + dxij*xin(27)
                                      yin(30) = yin(33) + dyij*yin(27)
                                      zin(30) = zin(33) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   36

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

                                      ! i1 = in(1) =   37

                                      xin(37) = 1.0_dp
                                      yin(37) = 1.0_dp
                                      zin(37) = f00

                                      ! i2 = in(2) =   43
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(43) = xc00
                                      yin(43) = yc00
                                      zin(43) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   38

                                      xin(38) = xcp00
                                      yin(38) = ycp00
                                      zin(38) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   44
                                      ! i2 =   43

                                      xin(44) = xcp00*xin(43) + cp10
                                      yin(44) = ycp00*yin(43) + cp10
                                      zin(44) = zcp00*zin(43) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   43

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =   37
                                      ! i4 =   43

                                      xin(49) = c10*xin(37) + xc00*xin(43)
                                      yin(49) = c10*yin(37) + yc00*yin(43)
                                      zin(49) = c10*zin(37) + zc00*zin(43)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   50
                                      ! i5 =   49
                                      ! i4 =   43

                                      xin(50) = xcp00*xin(49) + cp10*xin(43)
                                      yin(50) = ycp00*yin(49) + cp10*yin(43)
                                      zin(50) = zcp00*zin(49) + cp10*zin(43)

                                      ! ------------------

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   49

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   52
                                      ! i3 =   43
                                      ! i4 =   49

                                      xin(52) = c10*xin(43) + xc00*xin(49)
                                      yin(52) = c10*yin(43) + yc00*yin(49)
                                      zin(52) = c10*zin(43) + zc00*zin(49)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   53
                                      ! i5 =   52
                                      ! i4 =   49

                                      xin(53) = xcp00*xin(52) + cp10*xin(49)
                                      yin(53) = ycp00*yin(52) + cp10*yin(49)
                                      zin(53) = zcp00*zin(52) + cp10*zin(49)

                                      ! ------------------

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   52

                                      ! n =    4

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

                                      ! i3 = i2+kn(n+1) =   45

                                      xin(45) = xc00*xin(39) + c01*xin(38)
                                      yin(45) = yc00*yin(39) + c01*yin(38)
                                      zin(45) = zc00*zin(39) + c01*zin(38)

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
                                      ! i4 = i2 =   43

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   49

                                      xin(51) = c10*xin(39) + xc00*xin(45) + c01*xin(44)
                                      yin(51) = c10*yin(39) + yc00*yin(45) + c01*yin(44)
                                      zin(51) = c10*zin(39) + zc00*zin(45) + c01*zin(44)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   43
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   52

                                      xin(54) = c10*xin(45) + xc00*xin(51) + c01*xin(50)
                                      yin(54) = c10*yin(45) + yc00*yin(51) + c01*yin(50)
                                      zin(54) = c10*zin(45) + zc00*zin(51) + c01*zin(50)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   52

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   52

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   52

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   49

                                      xin(52) = xin(52) + dxij*xin(49)
                                      yin(52) = yin(52) + dyij*yin(49)
                                      zin(52) = zin(52) + dzij*zin(49)

                                      ! i3 = i4 =   49
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   40

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   40

                                      ! do ni = 1,    2

                                      xin(40) = xin(43) + dxij*xin(37)
                                      yin(40) = yin(43) + dyij*yin(37)
                                      zin(40) = zin(43) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    2

                                      xin(46) = xin(49) + dxij*xin(43)
                                      yin(46) = yin(49) + dyij*yin(43)
                                      zin(46) = zin(49) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   52

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   43

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   53

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   50

                                      xin(53) = xin(53) + dxij*xin(50)
                                      yin(53) = yin(53) + dyij*yin(50)
                                      zin(53) = zin(53) + dzij*zin(50)

                                      ! i3 = i4 =   50
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   41

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   41

                                      ! do ni = 1,    2

                                      xin(41) = xin(44) + dxij*xin(38)
                                      yin(41) = yin(44) + dyij*yin(38)
                                      zin(41) = zin(44) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    2

                                      xin(47) = xin(50) + dxij*xin(44)
                                      yin(47) = yin(50) + dyij*yin(44)
                                      zin(47) = zin(50) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   44

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   54

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   51

                                      xin(54) = xin(54) + dxij*xin(51)
                                      yin(54) = yin(54) + dyij*yin(51)
                                      zin(54) = zin(54) + dzij*zin(51)

                                      ! i3 = i4 =   51
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   42

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   42

                                      ! do ni = 1,    2

                                      xin(42) = xin(45) + dxij*xin(39)
                                      yin(42) = yin(45) + dyij*yin(39)
                                      zin(42) = zin(45) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    2

                                      xin(48) = xin(51) + dxij*xin(45)
                                      yin(48) = yin(51) + dyij*yin(45)
                                      zin(48) = zin(51) + dzij*zin(45)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   45

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   54

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

       eri_value(1) = eri_value(1) + d12bra(1)*d02ket(1)*(xin(18)*yin(1)*zin(1) + xin(36)*yin(19)*zin(19) + xin(54)*yin(37)*zin(37))
       eri_value(2) = eri_value(2) + d12bra(1)*d02ket(2)*(xin(16)*yin(3)*zin(1) + xin(34)*yin(21)*zin(19) + xin(52)*yin(39)*zin(37))
       eri_value(3) = eri_value(3) + d12bra(1)*d02ket(3)*(xin(16)*yin(1)*zin(3) + xin(34)*yin(19)*zin(21) + xin(52)*yin(37)*zin(39))
       eri_value(4) = eri_value(4) + d12bra(1)*d02ket(4)*(xin(17)*yin(2)*zin(1) + xin(35)*yin(20)*zin(19) + xin(53)*yin(38)*zin(37))
       eri_value(5) = eri_value(5) + d12bra(1)*d02ket(5)*(xin(17)*yin(1)*zin(2) + xin(35)*yin(19)*zin(20) + xin(53)*yin(37)*zin(38))
       eri_value(6) = eri_value(6) + d12bra(1)*d02ket(6)*(xin(16)*yin(2)*zin(2) + xin(34)*yin(20)*zin(20) + xin(52)*yin(38)*zin(38))
       eri_value(7) = eri_value(7) + d12bra(2)*d02ket(1)*(xin(15)*yin(4)*zin(1) + xin(33)*yin(22)*zin(19) + xin(51)*yin(40)*zin(37))
       eri_value(8) = eri_value(8) + d12bra(2)*d02ket(2)*(xin(13)*yin(6)*zin(1) + xin(31)*yin(24)*zin(19) + xin(49)*yin(42)*zin(37))
       eri_value(9) = eri_value(9) + d12bra(2)*d02ket(3)*(xin(13)*yin(4)*zin(3) + xin(31)*yin(22)*zin(21) + xin(49)*yin(40)*zin(39))
     eri_value(10) = eri_value(10) + d12bra(2)*d02ket(4)*(xin(14)*yin(5)*zin(1) + xin(32)*yin(23)*zin(19) + xin(50)*yin(41)*zin(37))
     eri_value(11) = eri_value(11) + d12bra(2)*d02ket(5)*(xin(14)*yin(4)*zin(2) + xin(32)*yin(22)*zin(20) + xin(50)*yin(40)*zin(38))
     eri_value(12) = eri_value(12) + d12bra(2)*d02ket(6)*(xin(13)*yin(5)*zin(2) + xin(31)*yin(23)*zin(20) + xin(49)*yin(41)*zin(38))
     eri_value(13) = eri_value(13) + d12bra(3)*d02ket(1)*(xin(15)*yin(1)*zin(4) + xin(33)*yin(19)*zin(22) + xin(51)*yin(37)*zin(40))
     eri_value(14) = eri_value(14) + d12bra(3)*d02ket(2)*(xin(13)*yin(3)*zin(4) + xin(31)*yin(21)*zin(22) + xin(49)*yin(39)*zin(40))
     eri_value(15) = eri_value(15) + d12bra(3)*d02ket(3)*(xin(13)*yin(1)*zin(6) + xin(31)*yin(19)*zin(24) + xin(49)*yin(37)*zin(42))
     eri_value(16) = eri_value(16) + d12bra(3)*d02ket(4)*(xin(14)*yin(2)*zin(4) + xin(32)*yin(20)*zin(22) + xin(50)*yin(38)*zin(40))
     eri_value(17) = eri_value(17) + d12bra(3)*d02ket(5)*(xin(14)*yin(1)*zin(5) + xin(32)*yin(19)*zin(23) + xin(50)*yin(37)*zin(41))
     eri_value(18) = eri_value(18) + d12bra(3)*d02ket(6)*(xin(13)*yin(2)*zin(5) + xin(31)*yin(20)*zin(23) + xin(49)*yin(38)*zin(41))
     eri_value(19) = eri_value(19) + d12bra(4)*d02ket(1)*(xin(6)*yin(13)*zin(1) + xin(24)*yin(31)*zin(19) + xin(42)*yin(49)*zin(37))
     eri_value(20) = eri_value(20) + d12bra(4)*d02ket(2)*(xin(4)*yin(15)*zin(1) + xin(22)*yin(33)*zin(19) + xin(40)*yin(51)*zin(37))
     eri_value(21) = eri_value(21) + d12bra(4)*d02ket(3)*(xin(4)*yin(13)*zin(3) + xin(22)*yin(31)*zin(21) + xin(40)*yin(49)*zin(39))
     eri_value(22) = eri_value(22) + d12bra(4)*d02ket(4)*(xin(5)*yin(14)*zin(1) + xin(23)*yin(32)*zin(19) + xin(41)*yin(50)*zin(37))
     eri_value(23) = eri_value(23) + d12bra(4)*d02ket(5)*(xin(5)*yin(13)*zin(2) + xin(23)*yin(31)*zin(20) + xin(41)*yin(49)*zin(38))
     eri_value(24) = eri_value(24) + d12bra(4)*d02ket(6)*(xin(4)*yin(14)*zin(2) + xin(22)*yin(32)*zin(20) + xin(40)*yin(50)*zin(38))
     eri_value(25) = eri_value(25) + d12bra(5)*d02ket(1)*(xin(3)*yin(16)*zin(1) + xin(21)*yin(34)*zin(19) + xin(39)*yin(52)*zin(37))
     eri_value(26) = eri_value(26) + d12bra(5)*d02ket(2)*(xin(1)*yin(18)*zin(1) + xin(19)*yin(36)*zin(19) + xin(37)*yin(54)*zin(37))
     eri_value(27) = eri_value(27) + d12bra(5)*d02ket(3)*(xin(1)*yin(16)*zin(3) + xin(19)*yin(34)*zin(21) + xin(37)*yin(52)*zin(39))
     eri_value(28) = eri_value(28) + d12bra(5)*d02ket(4)*(xin(2)*yin(17)*zin(1) + xin(20)*yin(35)*zin(19) + xin(38)*yin(53)*zin(37))
     eri_value(29) = eri_value(29) + d12bra(5)*d02ket(5)*(xin(2)*yin(16)*zin(2) + xin(20)*yin(34)*zin(20) + xin(38)*yin(52)*zin(38))
     eri_value(30) = eri_value(30) + d12bra(5)*d02ket(6)*(xin(1)*yin(17)*zin(2) + xin(19)*yin(35)*zin(20) + xin(37)*yin(53)*zin(38))
     eri_value(31) = eri_value(31) + d12bra(6)*d02ket(1)*(xin(3)*yin(13)*zin(4) + xin(21)*yin(31)*zin(22) + xin(39)*yin(49)*zin(40))
     eri_value(32) = eri_value(32) + d12bra(6)*d02ket(2)*(xin(1)*yin(15)*zin(4) + xin(19)*yin(33)*zin(22) + xin(37)*yin(51)*zin(40))
     eri_value(33) = eri_value(33) + d12bra(6)*d02ket(3)*(xin(1)*yin(13)*zin(6) + xin(19)*yin(31)*zin(24) + xin(37)*yin(49)*zin(42))
     eri_value(34) = eri_value(34) + d12bra(6)*d02ket(4)*(xin(2)*yin(14)*zin(4) + xin(20)*yin(32)*zin(22) + xin(38)*yin(50)*zin(40))
     eri_value(35) = eri_value(35) + d12bra(6)*d02ket(5)*(xin(2)*yin(13)*zin(5) + xin(20)*yin(31)*zin(23) + xin(38)*yin(49)*zin(41))
     eri_value(36) = eri_value(36) + d12bra(6)*d02ket(6)*(xin(1)*yin(14)*zin(5) + xin(19)*yin(32)*zin(23) + xin(37)*yin(50)*zin(41))
     eri_value(37) = eri_value(37) + d12bra(7)*d02ket(1)*(xin(6)*yin(1)*zin(13) + xin(24)*yin(19)*zin(31) + xin(42)*yin(37)*zin(49))
     eri_value(38) = eri_value(38) + d12bra(7)*d02ket(2)*(xin(4)*yin(3)*zin(13) + xin(22)*yin(21)*zin(31) + xin(40)*yin(39)*zin(49))
     eri_value(39) = eri_value(39) + d12bra(7)*d02ket(3)*(xin(4)*yin(1)*zin(15) + xin(22)*yin(19)*zin(33) + xin(40)*yin(37)*zin(51))
     eri_value(40) = eri_value(40) + d12bra(7)*d02ket(4)*(xin(5)*yin(2)*zin(13) + xin(23)*yin(20)*zin(31) + xin(41)*yin(38)*zin(49))
     eri_value(41) = eri_value(41) + d12bra(7)*d02ket(5)*(xin(5)*yin(1)*zin(14) + xin(23)*yin(19)*zin(32) + xin(41)*yin(37)*zin(50))
     eri_value(42) = eri_value(42) + d12bra(7)*d02ket(6)*(xin(4)*yin(2)*zin(14) + xin(22)*yin(20)*zin(32) + xin(40)*yin(38)*zin(50))
     eri_value(43) = eri_value(43) + d12bra(8)*d02ket(1)*(xin(3)*yin(4)*zin(13) + xin(21)*yin(22)*zin(31) + xin(39)*yin(40)*zin(49))
     eri_value(44) = eri_value(44) + d12bra(8)*d02ket(2)*(xin(1)*yin(6)*zin(13) + xin(19)*yin(24)*zin(31) + xin(37)*yin(42)*zin(49))
     eri_value(45) = eri_value(45) + d12bra(8)*d02ket(3)*(xin(1)*yin(4)*zin(15) + xin(19)*yin(22)*zin(33) + xin(37)*yin(40)*zin(51))
     eri_value(46) = eri_value(46) + d12bra(8)*d02ket(4)*(xin(2)*yin(5)*zin(13) + xin(20)*yin(23)*zin(31) + xin(38)*yin(41)*zin(49))
     eri_value(47) = eri_value(47) + d12bra(8)*d02ket(5)*(xin(2)*yin(4)*zin(14) + xin(20)*yin(22)*zin(32) + xin(38)*yin(40)*zin(50))
     eri_value(48) = eri_value(48) + d12bra(8)*d02ket(6)*(xin(1)*yin(5)*zin(14) + xin(19)*yin(23)*zin(32) + xin(37)*yin(41)*zin(50))
     eri_value(49) = eri_value(49) + d12bra(9)*d02ket(1)*(xin(3)*yin(1)*zin(16) + xin(21)*yin(19)*zin(34) + xin(39)*yin(37)*zin(52))
     eri_value(50) = eri_value(50) + d12bra(9)*d02ket(2)*(xin(1)*yin(3)*zin(16) + xin(19)*yin(21)*zin(34) + xin(37)*yin(39)*zin(52))
     eri_value(51) = eri_value(51) + d12bra(9)*d02ket(3)*(xin(1)*yin(1)*zin(18) + xin(19)*yin(19)*zin(36) + xin(37)*yin(37)*zin(54))
     eri_value(52) = eri_value(52) + d12bra(9)*d02ket(4)*(xin(2)*yin(2)*zin(16) + xin(20)*yin(20)*zin(34) + xin(38)*yin(38)*zin(52))
     eri_value(53) = eri_value(53) + d12bra(9)*d02ket(5)*(xin(2)*yin(1)*zin(17) + xin(20)*yin(19)*zin(35) + xin(38)*yin(37)*zin(53))
     eri_value(54) = eri_value(54) + d12bra(9)*d02ket(6)*(xin(1)*yin(2)*zin(17) + xin(19)*yin(20)*zin(35) + xin(37)*yin(38)*zin(53))
    eri_value(55) = eri_value(55) + d12bra(10)*d02ket(1)*(xin(12)*yin(7)*zin(1) + xin(30)*yin(25)*zin(19) + xin(48)*yin(43)*zin(37))
    eri_value(56) = eri_value(56) + d12bra(10)*d02ket(2)*(xin(10)*yin(9)*zin(1) + xin(28)*yin(27)*zin(19) + xin(46)*yin(45)*zin(37))
    eri_value(57) = eri_value(57) + d12bra(10)*d02ket(3)*(xin(10)*yin(7)*zin(3) + xin(28)*yin(25)*zin(21) + xin(46)*yin(43)*zin(39))
    eri_value(58) = eri_value(58) + d12bra(10)*d02ket(4)*(xin(11)*yin(8)*zin(1) + xin(29)*yin(26)*zin(19) + xin(47)*yin(44)*zin(37))
    eri_value(59) = eri_value(59) + d12bra(10)*d02ket(5)*(xin(11)*yin(7)*zin(2) + xin(29)*yin(25)*zin(20) + xin(47)*yin(43)*zin(38))
    eri_value(60) = eri_value(60) + d12bra(10)*d02ket(6)*(xin(10)*yin(8)*zin(2) + xin(28)*yin(26)*zin(20) + xin(46)*yin(44)*zin(38))
    eri_value(61) = eri_value(61) + d12bra(11)*d02ket(1)*(xin(9)*yin(10)*zin(1) + xin(27)*yin(28)*zin(19) + xin(45)*yin(46)*zin(37))
    eri_value(62) = eri_value(62) + d12bra(11)*d02ket(2)*(xin(7)*yin(12)*zin(1) + xin(25)*yin(30)*zin(19) + xin(43)*yin(48)*zin(37))
    eri_value(63) = eri_value(63) + d12bra(11)*d02ket(3)*(xin(7)*yin(10)*zin(3) + xin(25)*yin(28)*zin(21) + xin(43)*yin(46)*zin(39))
    eri_value(64) = eri_value(64) + d12bra(11)*d02ket(4)*(xin(8)*yin(11)*zin(1) + xin(26)*yin(29)*zin(19) + xin(44)*yin(47)*zin(37))
    eri_value(65) = eri_value(65) + d12bra(11)*d02ket(5)*(xin(8)*yin(10)*zin(2) + xin(26)*yin(28)*zin(20) + xin(44)*yin(46)*zin(38))
    eri_value(66) = eri_value(66) + d12bra(11)*d02ket(6)*(xin(7)*yin(11)*zin(2) + xin(25)*yin(29)*zin(20) + xin(43)*yin(47)*zin(38))
     eri_value(67) = eri_value(67) + d12bra(12)*d02ket(1)*(xin(9)*yin(7)*zin(4) + xin(27)*yin(25)*zin(22) + xin(45)*yin(43)*zin(40))
     eri_value(68) = eri_value(68) + d12bra(12)*d02ket(2)*(xin(7)*yin(9)*zin(4) + xin(25)*yin(27)*zin(22) + xin(43)*yin(45)*zin(40))
     eri_value(69) = eri_value(69) + d12bra(12)*d02ket(3)*(xin(7)*yin(7)*zin(6) + xin(25)*yin(25)*zin(24) + xin(43)*yin(43)*zin(42))
     eri_value(70) = eri_value(70) + d12bra(12)*d02ket(4)*(xin(8)*yin(8)*zin(4) + xin(26)*yin(26)*zin(22) + xin(44)*yin(44)*zin(40))
     eri_value(71) = eri_value(71) + d12bra(12)*d02ket(5)*(xin(8)*yin(7)*zin(5) + xin(26)*yin(25)*zin(23) + xin(44)*yin(43)*zin(41))
     eri_value(72) = eri_value(72) + d12bra(12)*d02ket(6)*(xin(7)*yin(8)*zin(5) + xin(25)*yin(26)*zin(23) + xin(43)*yin(44)*zin(41))
    eri_value(73) = eri_value(73) + d12bra(13)*d02ket(1)*(xin(12)*yin(1)*zin(7) + xin(30)*yin(19)*zin(25) + xin(48)*yin(37)*zin(43))
    eri_value(74) = eri_value(74) + d12bra(13)*d02ket(2)*(xin(10)*yin(3)*zin(7) + xin(28)*yin(21)*zin(25) + xin(46)*yin(39)*zin(43))
    eri_value(75) = eri_value(75) + d12bra(13)*d02ket(3)*(xin(10)*yin(1)*zin(9) + xin(28)*yin(19)*zin(27) + xin(46)*yin(37)*zin(45))
    eri_value(76) = eri_value(76) + d12bra(13)*d02ket(4)*(xin(11)*yin(2)*zin(7) + xin(29)*yin(20)*zin(25) + xin(47)*yin(38)*zin(43))
    eri_value(77) = eri_value(77) + d12bra(13)*d02ket(5)*(xin(11)*yin(1)*zin(8) + xin(29)*yin(19)*zin(26) + xin(47)*yin(37)*zin(44))
    eri_value(78) = eri_value(78) + d12bra(13)*d02ket(6)*(xin(10)*yin(2)*zin(8) + xin(28)*yin(20)*zin(26) + xin(46)*yin(38)*zin(44))
     eri_value(79) = eri_value(79) + d12bra(14)*d02ket(1)*(xin(9)*yin(4)*zin(7) + xin(27)*yin(22)*zin(25) + xin(45)*yin(40)*zin(43))
     eri_value(80) = eri_value(80) + d12bra(14)*d02ket(2)*(xin(7)*yin(6)*zin(7) + xin(25)*yin(24)*zin(25) + xin(43)*yin(42)*zin(43))
     eri_value(81) = eri_value(81) + d12bra(14)*d02ket(3)*(xin(7)*yin(4)*zin(9) + xin(25)*yin(22)*zin(27) + xin(43)*yin(40)*zin(45))
     eri_value(82) = eri_value(82) + d12bra(14)*d02ket(4)*(xin(8)*yin(5)*zin(7) + xin(26)*yin(23)*zin(25) + xin(44)*yin(41)*zin(43))
     eri_value(83) = eri_value(83) + d12bra(14)*d02ket(5)*(xin(8)*yin(4)*zin(8) + xin(26)*yin(22)*zin(26) + xin(44)*yin(40)*zin(44))
     eri_value(84) = eri_value(84) + d12bra(14)*d02ket(6)*(xin(7)*yin(5)*zin(8) + xin(25)*yin(23)*zin(26) + xin(43)*yin(41)*zin(44))
    eri_value(85) = eri_value(85) + d12bra(15)*d02ket(1)*(xin(9)*yin(1)*zin(10) + xin(27)*yin(19)*zin(28) + xin(45)*yin(37)*zin(46))
    eri_value(86) = eri_value(86) + d12bra(15)*d02ket(2)*(xin(7)*yin(3)*zin(10) + xin(25)*yin(21)*zin(28) + xin(43)*yin(39)*zin(46))
    eri_value(87) = eri_value(87) + d12bra(15)*d02ket(3)*(xin(7)*yin(1)*zin(12) + xin(25)*yin(19)*zin(30) + xin(43)*yin(37)*zin(48))
    eri_value(88) = eri_value(88) + d12bra(15)*d02ket(4)*(xin(8)*yin(2)*zin(10) + xin(26)*yin(20)*zin(28) + xin(44)*yin(38)*zin(46))
    eri_value(89) = eri_value(89) + d12bra(15)*d02ket(5)*(xin(8)*yin(1)*zin(11) + xin(26)*yin(19)*zin(29) + xin(44)*yin(37)*zin(47))
    eri_value(90) = eri_value(90) + d12bra(15)*d02ket(6)*(xin(7)*yin(2)*zin(11) + xin(25)*yin(20)*zin(29) + xin(43)*yin(38)*zin(47))
     eri_value(91) = eri_value(91) + d12bra(16)*d02ket(1)*(xin(6)*yin(7)*zin(7) + xin(24)*yin(25)*zin(25) + xin(42)*yin(43)*zin(43))
     eri_value(92) = eri_value(92) + d12bra(16)*d02ket(2)*(xin(4)*yin(9)*zin(7) + xin(22)*yin(27)*zin(25) + xin(40)*yin(45)*zin(43))
     eri_value(93) = eri_value(93) + d12bra(16)*d02ket(3)*(xin(4)*yin(7)*zin(9) + xin(22)*yin(25)*zin(27) + xin(40)*yin(43)*zin(45))
     eri_value(94) = eri_value(94) + d12bra(16)*d02ket(4)*(xin(5)*yin(8)*zin(7) + xin(23)*yin(26)*zin(25) + xin(41)*yin(44)*zin(43))
     eri_value(95) = eri_value(95) + d12bra(16)*d02ket(5)*(xin(5)*yin(7)*zin(8) + xin(23)*yin(25)*zin(26) + xin(41)*yin(43)*zin(44))
     eri_value(96) = eri_value(96) + d12bra(16)*d02ket(6)*(xin(4)*yin(8)*zin(8) + xin(22)*yin(26)*zin(26) + xin(40)*yin(44)*zin(44))
    eri_value(97) = eri_value(97) + d12bra(17)*d02ket(1)*(xin(3)*yin(10)*zin(7) + xin(21)*yin(28)*zin(25) + xin(39)*yin(46)*zin(43))
    eri_value(98) = eri_value(98) + d12bra(17)*d02ket(2)*(xin(1)*yin(12)*zin(7) + xin(19)*yin(30)*zin(25) + xin(37)*yin(48)*zin(43))
    eri_value(99) = eri_value(99) + d12bra(17)*d02ket(3)*(xin(1)*yin(10)*zin(9) + xin(19)*yin(28)*zin(27) + xin(37)*yin(46)*zin(45))
  eri_value(100) = eri_value(100) + d12bra(17)*d02ket(4)*(xin(2)*yin(11)*zin(7) + xin(20)*yin(29)*zin(25) + xin(38)*yin(47)*zin(43))
  eri_value(101) = eri_value(101) + d12bra(17)*d02ket(5)*(xin(2)*yin(10)*zin(8) + xin(20)*yin(28)*zin(26) + xin(38)*yin(46)*zin(44))
  eri_value(102) = eri_value(102) + d12bra(17)*d02ket(6)*(xin(1)*yin(11)*zin(8) + xin(19)*yin(29)*zin(26) + xin(37)*yin(47)*zin(44))
  eri_value(103) = eri_value(103) + d12bra(18)*d02ket(1)*(xin(3)*yin(7)*zin(10) + xin(21)*yin(25)*zin(28) + xin(39)*yin(43)*zin(46))
  eri_value(104) = eri_value(104) + d12bra(18)*d02ket(2)*(xin(1)*yin(9)*zin(10) + xin(19)*yin(27)*zin(28) + xin(37)*yin(45)*zin(46))
  eri_value(105) = eri_value(105) + d12bra(18)*d02ket(3)*(xin(1)*yin(7)*zin(12) + xin(19)*yin(25)*zin(30) + xin(37)*yin(43)*zin(48))
  eri_value(106) = eri_value(106) + d12bra(18)*d02ket(4)*(xin(2)*yin(8)*zin(10) + xin(20)*yin(26)*zin(28) + xin(38)*yin(44)*zin(46))
  eri_value(107) = eri_value(107) + d12bra(18)*d02ket(5)*(xin(2)*yin(7)*zin(11) + xin(20)*yin(25)*zin(29) + xin(38)*yin(43)*zin(47))
  eri_value(108) = eri_value(108) + d12bra(18)*d02ket(6)*(xin(1)*yin(8)*zin(11) + xin(19)*yin(26)*zin(29) + xin(37)*yin(44)*zin(47))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 6 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*18 ! Stride between functions in i

                                    do j = 1, 3 ! # of cartesians in j

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

                              deallocate (n12bra)
                              deallocate (xint12bra)
                              deallocate (n02ket)
                              deallocate (xint02ket)

                              end subroutine int2120
                              end submodule
