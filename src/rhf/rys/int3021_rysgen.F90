! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int2130_impl
contains
  module subroutine int2130(pd_pair, sf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pd_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n12bra(:), n03ket(:)
    real(dp), allocatable :: xint12bra(:), xint03ket(:)
    integer(kind=int64) :: npdbra, nsfket
    real(dp) :: scutpdbra, scutsfket, test
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
    real(dp) :: xin(96), yin(96), zin(96)
    real(dp) :: eri_value(180)
    real(dp) :: d12bra(18), d03ket(10)
    integer(kind=int64) :: ix(6), jx(3), kx(10), lx(1)
    integer(kind=int64) :: iy(6), jy(3), ky(10), ly(1)
    integer(kind=int64) :: iz(6), jz(3), kz(10), lz(1)
    integer(kind=int64) :: in(4), in1(4), kn(4)
    integer(kind=int64) :: ijx(18), ijy(18), ijz(18)
    integer(kind=int64) :: klx(10), kly(10), klz(10)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 9
    in1(3) = 17
    in1(4) = 21

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

    jx(1) = 4
    jx(2) = 0
    jx(3) = 0

    ix(1) = 17
    ix(2) = 1
    ix(3) = 1
    ix(4) = 9
    ix(5) = 9
    ix(6) = 1

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
    jy(2) = 4
    jy(3) = 0

    iy(1) = 1
    iy(2) = 17
    iy(3) = 1
    iy(4) = 9
    iy(5) = 1
    iy(6) = 9

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
    jz(3) = 4

    iz(1) = 1
    iz(2) = 1
    iz(3) = 17
    iz(4) = 1
    iz(5) = 9
    iz(6) = 9

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 21
    ijx(2) = 17
    ijx(3) = 17
    ijx(4) = 5
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 5
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 13
    ijx(11) = 9
    ijx(12) = 9
    ijx(13) = 13
    ijx(14) = 9
    ijx(15) = 9
    ijx(16) = 5
    ijx(17) = 1
    ijx(18) = 1

    ijy(1) = 1
    ijy(2) = 5
    ijy(3) = 1
    ijy(4) = 17
    ijy(5) = 21
    ijy(6) = 17
    ijy(7) = 1
    ijy(8) = 5
    ijy(9) = 1
    ijy(10) = 9
    ijy(11) = 13
    ijy(12) = 9
    ijy(13) = 1
    ijy(14) = 5
    ijy(15) = 1
    ijy(16) = 9
    ijy(17) = 13
    ijy(18) = 9

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 5
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 5
    ijz(7) = 17
    ijz(8) = 17
    ijz(9) = 21
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 5
    ijz(13) = 9
    ijz(14) = 9
    ijz(15) = 13
    ijz(16) = 9
    ijz(17) = 9
    ijz(18) = 13

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

    allocate (n12bra(res%n_p_shl*res%n_d_shl))
    allocate (xint12bra(res%n_p_shl*res%n_d_shl))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

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

    if ((npdbra*nsfket) .le. nchunksize_int64) nchunksize_int64 = npdbra*nsfket
    ntile = int(npdbra*nsfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = npdbra*nsfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npdbra, xint12bra, n12bra, xint03ket, n03ket, pd_pair, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d03ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
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

              test = xint12bra(ij_tmp)*xint03ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n12bra(ij_tmp)
                kl = n03ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_d_shl) + 1
                jsh_tmp = (ij - 1)/res%n_d_shl + 1
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_d_shl(ish_tmp)
                jsh = res%i_p_shl(jsh_tmp)
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

                                      ! i2 = in(2) =    9
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(9) = xc00
                                      yin(9) = yc00
                                      zin(9) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   10
                                      ! i2 =    9

                                      xin(10) = xcp00*xin(9) + cp10
                                      yin(10) = ycp00*yin(9) + cp10
                                      zin(10) = zcp00*zin(9) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    9

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   17
                                      ! i3 =    1
                                      ! i4 =    9

                                      xin(17) = c10*xin(1) + xc00*xin(9)
                                      yin(17) = c10*yin(1) + yc00*yin(9)
                                      zin(17) = c10*zin(1) + zc00*zin(9)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   18
                                      ! i5 =   17
                                      ! i4 =    9

                                      xin(18) = xcp00*xin(17) + cp10*xin(9)
                                      yin(18) = ycp00*yin(17) + cp10*yin(9)
                                      zin(18) = zcp00*zin(17) + cp10*zin(9)

                                      ! ------------------

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   17

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   21
                                      ! i3 =    9
                                      ! i4 =   17

                                      xin(21) = c10*xin(9) + xc00*xin(17)
                                      yin(21) = c10*yin(9) + yc00*yin(17)
                                      zin(21) = c10*zin(9) + zc00*zin(17)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   22
                                      ! i5 =   21
                                      ! i4 =   17

                                      xin(22) = xcp00*xin(21) + cp10*xin(17)
                                      yin(22) = ycp00*yin(21) + cp10*yin(17)
                                      zin(22) = zcp00*zin(21) + cp10*zin(17)

                                      ! ------------------

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   21

                                      ! n =    4

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

                                      ! i3 = i2+kn(n+1) =   11

                                      xin(11) = xc00*xin(3) + c01*xin(2)
                                      yin(11) = yc00*yin(3) + c01*yin(2)
                                      zin(11) = zc00*zin(3) + c01*zin(2)

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

                                      ! i3 = i2+kn(n+1) =   12

                                      xin(12) = xc00*xin(4) + c01*xin(3)
                                      yin(12) = yc00*yin(4) + c01*yin(3)
                                      zin(12) = zc00*zin(4) + c01*zin(3)

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
                                      ! i4 = i2 =    9

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   17

                                      xin(19) = c10*xin(3) + xc00*xin(11) + c01*xin(10)
                                      yin(19) = c10*yin(3) + yc00*yin(11) + c01*yin(10)
                                      zin(19) = c10*zin(3) + zc00*zin(11) + c01*zin(10)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   17

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   21

                                      xin(23) = c10*xin(11) + xc00*xin(19) + c01*xin(18)
                                      yin(23) = c10*yin(11) + yc00*yin(19) + c01*yin(18)
                                      zin(23) = c10*zin(11) + zc00*zin(19) + c01*zin(18)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   21

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    9

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   17

                                      xin(20) = c10*xin(4) + xc00*xin(12) + c01*xin(11)
                                      yin(20) = c10*yin(4) + yc00*yin(12) + c01*yin(11)
                                      zin(20) = c10*zin(4) + zc00*zin(12) + c01*zin(11)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   17

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   21

                                      xin(24) = c10*xin(12) + xc00*xin(20) + c01*xin(19)
                                      yin(24) = c10*yin(12) + yc00*yin(20) + c01*yin(19)
                                      zin(24) = c10*zin(12) + zc00*zin(20) + c01*zin(19)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   21

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   21

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   21

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   17

                                      xin(21) = xin(21) + dxij*xin(17)
                                      yin(21) = yin(21) + dyij*yin(17)
                                      zin(21) = zin(21) + dzij*zin(17)

                                      ! i3 = i4 =   17
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    2

                                      xin(5) = xin(9) + dxij*xin(1)
                                      yin(5) = yin(9) + dyij*yin(1)
                                      zin(5) = zin(9) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   13

                                      ! ni =    2

                                      xin(13) = xin(17) + dxij*xin(9)
                                      yin(13) = yin(17) + dyij*yin(9)
                                      zin(13) = zin(17) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   22

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   18

                                      xin(22) = xin(22) + dxij*xin(18)
                                      yin(22) = yin(22) + dyij*yin(18)
                                      zin(22) = zin(22) + dzij*zin(18)

                                      ! i3 = i4 =   18
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    2

                                      xin(6) = xin(10) + dxij*xin(2)
                                      yin(6) = yin(10) + dyij*yin(2)
                                      zin(6) = zin(10) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   14

                                      ! ni =    2

                                      xin(14) = xin(18) + dxij*xin(10)
                                      yin(14) = yin(18) + dyij*yin(10)
                                      zin(14) = zin(18) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   23

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   19

                                      xin(23) = xin(23) + dxij*xin(19)
                                      yin(23) = yin(23) + dyij*yin(19)
                                      zin(23) = zin(23) + dzij*zin(19)

                                      ! i3 = i4 =   19
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    2

                                      xin(7) = xin(11) + dxij*xin(3)
                                      yin(7) = yin(11) + dyij*yin(3)
                                      zin(7) = zin(11) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   15

                                      ! ni =    2

                                      xin(15) = xin(19) + dxij*xin(11)
                                      yin(15) = yin(19) + dyij*yin(11)
                                      zin(15) = zin(19) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   20

                                      xin(24) = xin(24) + dxij*xin(20)
                                      yin(24) = yin(24) + dyij*yin(20)
                                      zin(24) = zin(24) + dzij*zin(20)

                                      ! i3 = i4 =   20
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    2

                                      xin(8) = xin(12) + dxij*xin(4)
                                      yin(8) = yin(12) + dyij*yin(4)
                                      zin(8) = zin(12) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    2

                                      xin(16) = xin(20) + dxij*xin(12)
                                      yin(16) = yin(20) + dyij*yin(12)
                                      zin(16) = zin(20) + dzij*zin(12)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =   33
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(33) = xc00
                                      yin(33) = yc00
                                      zin(33) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   26

                                      xin(26) = xcp00
                                      yin(26) = ycp00
                                      zin(26) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   34
                                      ! i2 =   33

                                      xin(34) = xcp00*xin(33) + cp10
                                      yin(34) = ycp00*yin(33) + cp10
                                      zin(34) = zcp00*zin(33) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   33

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   41
                                      ! i3 =   25
                                      ! i4 =   33

                                      xin(41) = c10*xin(25) + xc00*xin(33)
                                      yin(41) = c10*yin(25) + yc00*yin(33)
                                      zin(41) = c10*zin(25) + zc00*zin(33)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   42
                                      ! i5 =   41
                                      ! i4 =   33

                                      xin(42) = xcp00*xin(41) + cp10*xin(33)
                                      yin(42) = ycp00*yin(41) + cp10*yin(33)
                                      zin(42) = zcp00*zin(41) + cp10*zin(33)

                                      ! ------------------

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   41

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   45
                                      ! i3 =   33
                                      ! i4 =   41

                                      xin(45) = c10*xin(33) + xc00*xin(41)
                                      yin(45) = c10*yin(33) + yc00*yin(41)
                                      zin(45) = c10*zin(33) + zc00*zin(41)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   46
                                      ! i5 =   45
                                      ! i4 =   41

                                      xin(46) = xcp00*xin(45) + cp10*xin(41)
                                      yin(46) = ycp00*yin(45) + cp10*yin(41)
                                      zin(46) = zcp00*zin(45) + cp10*zin(41)

                                      ! ------------------

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   25
                                      ! i4 = i1+k2 =   26

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   27
                                      ! i3 =   25
                                      ! i4 =   26

                                      xin(27) = cp01*xin(25) + xcp00*xin(26)
                                      yin(27) = cp01*yin(25) + ycp00*yin(26)
                                      zin(27) = cp01*zin(25) + zcp00*zin(26)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   35

                                      xin(35) = xc00*xin(27) + c01*xin(26)
                                      yin(35) = yc00*yin(27) + c01*yin(26)
                                      zin(35) = zc00*zin(27) + c01*zin(26)

                                      ! ------------------

                                      ! i3 = i4 =   26
                                      ! i4 = i5 =   27

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   28
                                      ! i3 =   26
                                      ! i4 =   27

                                      xin(28) = cp01*xin(26) + xcp00*xin(27)
                                      yin(28) = cp01*yin(26) + ycp00*yin(27)
                                      zin(28) = cp01*zin(26) + zcp00*zin(27)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   36

                                      xin(36) = xc00*xin(28) + c01*xin(27)
                                      yin(36) = yc00*yin(28) + c01*yin(27)
                                      zin(36) = zc00*zin(28) + c01*zin(27)

                                      ! ------------------

                                      ! i3 = i4 =   27
                                      ! i4 = i5 =   28

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   33

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   41

                                      xin(43) = c10*xin(27) + xc00*xin(35) + c01*xin(34)
                                      yin(43) = c10*yin(27) + yc00*yin(35) + c01*yin(34)
                                      zin(43) = c10*zin(27) + zc00*zin(35) + c01*zin(34)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   41

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   45

                                      xin(47) = c10*xin(35) + xc00*xin(43) + c01*xin(42)
                                      yin(47) = c10*yin(35) + yc00*yin(43) + c01*yin(42)
                                      zin(47) = c10*zin(35) + zc00*zin(43) + c01*zin(42)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   33

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   41

                                      xin(44) = c10*xin(28) + xc00*xin(36) + c01*xin(35)
                                      yin(44) = c10*yin(28) + yc00*yin(36) + c01*yin(35)
                                      zin(44) = c10*zin(28) + zc00*zin(36) + c01*zin(35)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   33
                                      ! i4 = i5 =   41

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   45

                                      xin(48) = c10*xin(36) + xc00*xin(44) + c01*xin(43)
                                      yin(48) = c10*yin(36) + yc00*yin(44) + c01*yin(43)
                                      zin(48) = c10*zin(36) + zc00*zin(44) + c01*zin(43)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   45

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(45) = xin(45) + dxij*xin(41)
                                      yin(45) = yin(45) + dyij*yin(41)
                                      zin(45) = zin(45) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   29

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   29

                                      ! do ni = 1,    2

                                      xin(29) = xin(33) + dxij*xin(25)
                                      yin(29) = yin(33) + dyij*yin(25)
                                      zin(29) = zin(33) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! ni =    2

                                      xin(37) = xin(41) + dxij*xin(33)
                                      yin(37) = yin(41) + dyij*yin(33)
                                      zin(37) = zin(41) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   33

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   42

                                      xin(46) = xin(46) + dxij*xin(42)
                                      yin(46) = yin(46) + dyij*yin(42)
                                      zin(46) = zin(46) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   30

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   30

                                      ! do ni = 1,    2

                                      xin(30) = xin(34) + dxij*xin(26)
                                      yin(30) = yin(34) + dyij*yin(26)
                                      zin(30) = zin(34) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   38

                                      ! ni =    2

                                      xin(38) = xin(42) + dxij*xin(34)
                                      yin(38) = yin(42) + dyij*yin(34)
                                      zin(38) = zin(42) + dzij*zin(34)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   34

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(47) = xin(47) + dxij*xin(43)
                                      yin(47) = yin(47) + dyij*yin(43)
                                      zin(47) = zin(47) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   31

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   31

                                      ! do ni = 1,    2

                                      xin(31) = xin(35) + dxij*xin(27)
                                      yin(31) = yin(35) + dyij*yin(27)
                                      zin(31) = zin(35) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    2

                                      xin(39) = xin(43) + dxij*xin(35)
                                      yin(39) = yin(43) + dyij*yin(35)
                                      zin(39) = zin(43) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   35

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(48) = xin(48) + dxij*xin(44)
                                      yin(48) = yin(48) + dyij*yin(44)
                                      zin(48) = zin(48) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   32

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   32

                                      ! do ni = 1,    2

                                      xin(32) = xin(36) + dxij*xin(28)
                                      yin(32) = yin(36) + dyij*yin(28)
                                      zin(32) = zin(36) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    2

                                      xin(40) = xin(44) + dxij*xin(36)
                                      yin(40) = yin(44) + dyij*yin(36)
                                      zin(40) = zin(44) + dzij*zin(36)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   36

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =   57
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(57) = xc00
                                      yin(57) = yc00
                                      zin(57) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   50

                                      xin(50) = xcp00
                                      yin(50) = ycp00
                                      zin(50) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   58
                                      ! i2 =   57

                                      xin(58) = xcp00*xin(57) + cp10
                                      yin(58) = ycp00*yin(57) + cp10
                                      zin(58) = zcp00*zin(57) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   57

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   65
                                      ! i3 =   49
                                      ! i4 =   57

                                      xin(65) = c10*xin(49) + xc00*xin(57)
                                      yin(65) = c10*yin(49) + yc00*yin(57)
                                      zin(65) = c10*zin(49) + zc00*zin(57)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   66
                                      ! i5 =   65
                                      ! i4 =   57

                                      xin(66) = xcp00*xin(65) + cp10*xin(57)
                                      yin(66) = ycp00*yin(65) + cp10*yin(57)
                                      zin(66) = zcp00*zin(65) + cp10*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   69
                                      ! i3 =   57
                                      ! i4 =   65

                                      xin(69) = c10*xin(57) + xc00*xin(65)
                                      yin(69) = c10*yin(57) + yc00*yin(65)
                                      zin(69) = c10*zin(57) + zc00*zin(65)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   70
                                      ! i5 =   69
                                      ! i4 =   65

                                      xin(70) = xcp00*xin(69) + cp10*xin(65)
                                      yin(70) = ycp00*yin(69) + cp10*yin(65)
                                      zin(70) = zcp00*zin(69) + cp10*zin(65)

                                      ! ------------------

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   69

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   49
                                      ! i4 = i1+k2 =   50

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   51
                                      ! i3 =   49
                                      ! i4 =   50

                                      xin(51) = cp01*xin(49) + xcp00*xin(50)
                                      yin(51) = cp01*yin(49) + ycp00*yin(50)
                                      zin(51) = cp01*zin(49) + zcp00*zin(50)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   59

                                      xin(59) = xc00*xin(51) + c01*xin(50)
                                      yin(59) = yc00*yin(51) + c01*yin(50)
                                      zin(59) = zc00*zin(51) + c01*zin(50)

                                      ! ------------------

                                      ! i3 = i4 =   50
                                      ! i4 = i5 =   51

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   52
                                      ! i3 =   50
                                      ! i4 =   51

                                      xin(52) = cp01*xin(50) + xcp00*xin(51)
                                      yin(52) = cp01*yin(50) + ycp00*yin(51)
                                      zin(52) = cp01*zin(50) + zcp00*zin(51)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   60

                                      xin(60) = xc00*xin(52) + c01*xin(51)
                                      yin(60) = yc00*yin(52) + c01*yin(51)
                                      zin(60) = zc00*zin(52) + c01*zin(51)

                                      ! ------------------

                                      ! i3 = i4 =   51
                                      ! i4 = i5 =   52

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   57

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   65

                                      xin(67) = c10*xin(51) + xc00*xin(59) + c01*xin(58)
                                      yin(67) = c10*yin(51) + yc00*yin(59) + c01*yin(58)
                                      zin(67) = c10*zin(51) + zc00*zin(59) + c01*zin(58)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   69

                                      xin(71) = c10*xin(59) + xc00*xin(67) + c01*xin(66)
                                      yin(71) = c10*yin(59) + yc00*yin(67) + c01*yin(66)
                                      zin(71) = c10*zin(59) + zc00*zin(67) + c01*zin(66)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   69

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   57

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   65

                                      xin(68) = c10*xin(52) + xc00*xin(60) + c01*xin(59)
                                      yin(68) = c10*yin(52) + yc00*yin(60) + c01*yin(59)
                                      zin(68) = c10*zin(52) + zc00*zin(60) + c01*zin(59)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   69

                                      xin(72) = c10*xin(60) + xc00*xin(68) + c01*xin(67)
                                      yin(72) = c10*yin(60) + yc00*yin(68) + c01*yin(67)
                                      zin(72) = c10*zin(60) + zc00*zin(68) + c01*zin(67)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   69

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   69

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   65

                                      xin(69) = xin(69) + dxij*xin(65)
                                      yin(69) = yin(69) + dyij*yin(65)
                                      zin(69) = zin(69) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   53

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   53

                                      ! do ni = 1,    2

                                      xin(53) = xin(57) + dxij*xin(49)
                                      yin(53) = yin(57) + dyij*yin(49)
                                      zin(53) = zin(57) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    2

                                      xin(61) = xin(65) + dxij*xin(57)
                                      yin(61) = yin(65) + dyij*yin(57)
                                      zin(61) = zin(65) + dzij*zin(57)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   57

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   70

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   66

                                      xin(70) = xin(70) + dxij*xin(66)
                                      yin(70) = yin(70) + dyij*yin(66)
                                      zin(70) = zin(70) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   54

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   54

                                      ! do ni = 1,    2

                                      xin(54) = xin(58) + dxij*xin(50)
                                      yin(54) = yin(58) + dyij*yin(50)
                                      zin(54) = zin(58) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   62

                                      ! ni =    2

                                      xin(62) = xin(66) + dxij*xin(58)
                                      yin(62) = yin(66) + dyij*yin(58)
                                      zin(62) = zin(66) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   58

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   67

                                      xin(71) = xin(71) + dxij*xin(67)
                                      yin(71) = yin(71) + dyij*yin(67)
                                      zin(71) = zin(71) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   55

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   55

                                      ! do ni = 1,    2

                                      xin(55) = xin(59) + dxij*xin(51)
                                      yin(55) = yin(59) + dyij*yin(51)
                                      zin(55) = zin(59) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    2

                                      xin(63) = xin(67) + dxij*xin(59)
                                      yin(63) = yin(67) + dyij*yin(59)
                                      zin(63) = zin(67) + dzij*zin(59)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   59

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   68

                                      xin(72) = xin(72) + dxij*xin(68)
                                      yin(72) = yin(72) + dyij*yin(68)
                                      zin(72) = zin(72) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   56

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   56

                                      ! do ni = 1,    2

                                      xin(56) = xin(60) + dxij*xin(52)
                                      yin(56) = yin(60) + dyij*yin(52)
                                      zin(56) = zin(60) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    2

                                      xin(64) = xin(68) + dxij*xin(60)
                                      yin(64) = yin(68) + dyij*yin(60)
                                      zin(64) = zin(68) + dzij*zin(60)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   60

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =   81
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(81) = xc00
                                      yin(81) = yc00
                                      zin(81) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   74

                                      xin(74) = xcp00
                                      yin(74) = ycp00
                                      zin(74) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   82
                                      ! i2 =   81

                                      xin(82) = xcp00*xin(81) + cp10
                                      yin(82) = ycp00*yin(81) + cp10
                                      zin(82) = zcp00*zin(81) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   81

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   89
                                      ! i3 =   73
                                      ! i4 =   81

                                      xin(89) = c10*xin(73) + xc00*xin(81)
                                      yin(89) = c10*yin(73) + yc00*yin(81)
                                      zin(89) = c10*zin(73) + zc00*zin(81)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   90
                                      ! i5 =   89
                                      ! i4 =   81

                                      xin(90) = xcp00*xin(89) + cp10*xin(81)
                                      yin(90) = ycp00*yin(89) + cp10*yin(81)
                                      zin(90) = zcp00*zin(89) + cp10*zin(81)

                                      ! ------------------

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   89

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   93
                                      ! i3 =   81
                                      ! i4 =   89

                                      xin(93) = c10*xin(81) + xc00*xin(89)
                                      yin(93) = c10*yin(81) + yc00*yin(89)
                                      zin(93) = c10*zin(81) + zc00*zin(89)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   94
                                      ! i5 =   93
                                      ! i4 =   89

                                      xin(94) = xcp00*xin(93) + cp10*xin(89)
                                      yin(94) = ycp00*yin(93) + cp10*yin(89)
                                      zin(94) = zcp00*zin(93) + cp10*zin(89)

                                      ! ------------------

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   73
                                      ! i4 = i1+k2 =   74

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   75
                                      ! i3 =   73
                                      ! i4 =   74

                                      xin(75) = cp01*xin(73) + xcp00*xin(74)
                                      yin(75) = cp01*yin(73) + ycp00*yin(74)
                                      zin(75) = cp01*zin(73) + zcp00*zin(74)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   83

                                      xin(83) = xc00*xin(75) + c01*xin(74)
                                      yin(83) = yc00*yin(75) + c01*yin(74)
                                      zin(83) = zc00*zin(75) + c01*zin(74)

                                      ! ------------------

                                      ! i3 = i4 =   74
                                      ! i4 = i5 =   75

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   76
                                      ! i3 =   74
                                      ! i4 =   75

                                      xin(76) = cp01*xin(74) + xcp00*xin(75)
                                      yin(76) = cp01*yin(74) + ycp00*yin(75)
                                      zin(76) = cp01*zin(74) + zcp00*zin(75)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   84

                                      xin(84) = xc00*xin(76) + c01*xin(75)
                                      yin(84) = yc00*yin(76) + c01*yin(75)
                                      zin(84) = zc00*zin(76) + c01*zin(75)

                                      ! ------------------

                                      ! i3 = i4 =   75
                                      ! i4 = i5 =   76

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   81

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   89

                                      xin(91) = c10*xin(75) + xc00*xin(83) + c01*xin(82)
                                      yin(91) = c10*yin(75) + yc00*yin(83) + c01*yin(82)
                                      zin(91) = c10*zin(75) + zc00*zin(83) + c01*zin(82)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   89

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   93

                                      xin(95) = c10*xin(83) + xc00*xin(91) + c01*xin(90)
                                      yin(95) = c10*yin(83) + yc00*yin(91) + c01*yin(90)
                                      zin(95) = c10*zin(83) + zc00*zin(91) + c01*zin(90)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   81

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   89

                                      xin(92) = c10*xin(76) + xc00*xin(84) + c01*xin(83)
                                      yin(92) = c10*yin(76) + yc00*yin(84) + c01*yin(83)
                                      zin(92) = c10*zin(76) + zc00*zin(84) + c01*zin(83)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   89

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   93

                                      xin(96) = c10*xin(84) + xc00*xin(92) + c01*xin(91)
                                      yin(96) = c10*yin(84) + yc00*yin(92) + c01*yin(91)
                                      zin(96) = c10*zin(84) + zc00*zin(92) + c01*zin(91)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   93

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   77

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   77

                                      ! do ni = 1,    2

                                      xin(77) = xin(81) + dxij*xin(73)
                                      yin(77) = yin(81) + dyij*yin(73)
                                      zin(77) = zin(81) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    2

                                      xin(85) = xin(89) + dxij*xin(81)
                                      yin(85) = yin(89) + dyij*yin(81)
                                      zin(85) = zin(89) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   81

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(94) = xin(94) + dxij*xin(90)
                                      yin(94) = yin(94) + dyij*yin(90)
                                      zin(94) = zin(94) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   78

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   78

                                      ! do ni = 1,    2

                                      xin(78) = xin(82) + dxij*xin(74)
                                      yin(78) = yin(82) + dyij*yin(74)
                                      zin(78) = zin(82) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   86

                                      ! ni =    2

                                      xin(86) = xin(90) + dxij*xin(82)
                                      yin(86) = yin(90) + dyij*yin(82)
                                      zin(86) = zin(90) + dzij*zin(82)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   82

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   79

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   79

                                      ! do ni = 1,    2

                                      xin(79) = xin(83) + dxij*xin(75)
                                      yin(79) = yin(83) + dyij*yin(75)
                                      zin(79) = zin(83) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    2

                                      xin(87) = xin(91) + dxij*xin(83)
                                      yin(87) = yin(91) + dyij*yin(83)
                                      zin(87) = zin(91) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   83

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    3

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   80

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   80

                                      ! do ni = 1,    2

                                      xin(80) = xin(84) + dxij*xin(76)
                                      yin(80) = yin(84) + dyij*yin(76)
                                      zin(80) = zin(84) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    2

                                      xin(88) = xin(92) + dxij*xin(84)
                                      yin(88) = yin(92) + dyij*yin(84)
                                      zin(88) = zin(92) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   84

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d12bra(  1)*d03ket(  1)*(xin(  24)*yin(   1)*zin(   1)+xin(  48)*yin(  25)*zin(  25)+xin(  72)*yin(  49)*zin(  49)+xin(  96)*yin(  73)*zin(  73))
          eri_value(    2)=eri_value(    2)+d12bra(  1)*d03ket(  2)*(xin(  21)*yin(   4)*zin(   1)+xin(  45)*yin(  28)*zin(  25)+xin(  69)*yin(  52)*zin(  49)+xin(  93)*yin(  76)*zin(  73))
          eri_value(    3)=eri_value(    3)+d12bra(  1)*d03ket(  3)*(xin(  21)*yin(   1)*zin(   4)+xin(  45)*yin(  25)*zin(  28)+xin(  69)*yin(  49)*zin(  52)+xin(  93)*yin(  73)*zin(  76))
          eri_value(    4)=eri_value(    4)+d12bra(  1)*d03ket(  4)*(xin(  23)*yin(   2)*zin(   1)+xin(  47)*yin(  26)*zin(  25)+xin(  71)*yin(  50)*zin(  49)+xin(  95)*yin(  74)*zin(  73))
          eri_value(    5)=eri_value(    5)+d12bra(  1)*d03ket(  5)*(xin(  23)*yin(   1)*zin(   2)+xin(  47)*yin(  25)*zin(  26)+xin(  71)*yin(  49)*zin(  50)+xin(  95)*yin(  73)*zin(  74))
          eri_value(    6)=eri_value(    6)+d12bra(  1)*d03ket(  6)*(xin(  22)*yin(   3)*zin(   1)+xin(  46)*yin(  27)*zin(  25)+xin(  70)*yin(  51)*zin(  49)+xin(  94)*yin(  75)*zin(  73))
          eri_value(    7)=eri_value(    7)+d12bra(  1)*d03ket(  7)*(xin(  21)*yin(   3)*zin(   2)+xin(  45)*yin(  27)*zin(  26)+xin(  69)*yin(  51)*zin(  50)+xin(  93)*yin(  75)*zin(  74))
          eri_value(    8)=eri_value(    8)+d12bra(  1)*d03ket(  8)*(xin(  22)*yin(   1)*zin(   3)+xin(  46)*yin(  25)*zin(  27)+xin(  70)*yin(  49)*zin(  51)+xin(  94)*yin(  73)*zin(  75))
          eri_value(    9)=eri_value(    9)+d12bra(  1)*d03ket(  9)*(xin(  21)*yin(   2)*zin(   3)+xin(  45)*yin(  26)*zin(  27)+xin(  69)*yin(  50)*zin(  51)+xin(  93)*yin(  74)*zin(  75))
          eri_value(   10)=eri_value(   10)+d12bra(  1)*d03ket( 10)*(xin(  22)*yin(   2)*zin(   2)+xin(  46)*yin(  26)*zin(  26)+xin(  70)*yin(  50)*zin(  50)+xin(  94)*yin(  74)*zin(  74))
          eri_value(   11)=eri_value(   11)+d12bra(  2)*d03ket(  1)*(xin(  20)*yin(   5)*zin(   1)+xin(  44)*yin(  29)*zin(  25)+xin(  68)*yin(  53)*zin(  49)+xin(  92)*yin(  77)*zin(  73))
          eri_value(   12)=eri_value(   12)+d12bra(  2)*d03ket(  2)*(xin(  17)*yin(   8)*zin(   1)+xin(  41)*yin(  32)*zin(  25)+xin(  65)*yin(  56)*zin(  49)+xin(  89)*yin(  80)*zin(  73))
          eri_value(   13)=eri_value(   13)+d12bra(  2)*d03ket(  3)*(xin(  17)*yin(   5)*zin(   4)+xin(  41)*yin(  29)*zin(  28)+xin(  65)*yin(  53)*zin(  52)+xin(  89)*yin(  77)*zin(  76))
          eri_value(   14)=eri_value(   14)+d12bra(  2)*d03ket(  4)*(xin(  19)*yin(   6)*zin(   1)+xin(  43)*yin(  30)*zin(  25)+xin(  67)*yin(  54)*zin(  49)+xin(  91)*yin(  78)*zin(  73))
          eri_value(   15)=eri_value(   15)+d12bra(  2)*d03ket(  5)*(xin(  19)*yin(   5)*zin(   2)+xin(  43)*yin(  29)*zin(  26)+xin(  67)*yin(  53)*zin(  50)+xin(  91)*yin(  77)*zin(  74))
          eri_value(   16)=eri_value(   16)+d12bra(  2)*d03ket(  6)*(xin(  18)*yin(   7)*zin(   1)+xin(  42)*yin(  31)*zin(  25)+xin(  66)*yin(  55)*zin(  49)+xin(  90)*yin(  79)*zin(  73))
          eri_value(   17)=eri_value(   17)+d12bra(  2)*d03ket(  7)*(xin(  17)*yin(   7)*zin(   2)+xin(  41)*yin(  31)*zin(  26)+xin(  65)*yin(  55)*zin(  50)+xin(  89)*yin(  79)*zin(  74))
          eri_value(   18)=eri_value(   18)+d12bra(  2)*d03ket(  8)*(xin(  18)*yin(   5)*zin(   3)+xin(  42)*yin(  29)*zin(  27)+xin(  66)*yin(  53)*zin(  51)+xin(  90)*yin(  77)*zin(  75))
          eri_value(   19)=eri_value(   19)+d12bra(  2)*d03ket(  9)*(xin(  17)*yin(   6)*zin(   3)+xin(  41)*yin(  30)*zin(  27)+xin(  65)*yin(  54)*zin(  51)+xin(  89)*yin(  78)*zin(  75))
          eri_value(   20)=eri_value(   20)+d12bra(  2)*d03ket( 10)*(xin(  18)*yin(   6)*zin(   2)+xin(  42)*yin(  30)*zin(  26)+xin(  66)*yin(  54)*zin(  50)+xin(  90)*yin(  78)*zin(  74))
          eri_value(   21)=eri_value(   21)+d12bra(  3)*d03ket(  1)*(xin(  20)*yin(   1)*zin(   5)+xin(  44)*yin(  25)*zin(  29)+xin(  68)*yin(  49)*zin(  53)+xin(  92)*yin(  73)*zin(  77))
          eri_value(   22)=eri_value(   22)+d12bra(  3)*d03ket(  2)*(xin(  17)*yin(   4)*zin(   5)+xin(  41)*yin(  28)*zin(  29)+xin(  65)*yin(  52)*zin(  53)+xin(  89)*yin(  76)*zin(  77))
          eri_value(   23)=eri_value(   23)+d12bra(  3)*d03ket(  3)*(xin(  17)*yin(   1)*zin(   8)+xin(  41)*yin(  25)*zin(  32)+xin(  65)*yin(  49)*zin(  56)+xin(  89)*yin(  73)*zin(  80))
          eri_value(   24)=eri_value(   24)+d12bra(  3)*d03ket(  4)*(xin(  19)*yin(   2)*zin(   5)+xin(  43)*yin(  26)*zin(  29)+xin(  67)*yin(  50)*zin(  53)+xin(  91)*yin(  74)*zin(  77))
          eri_value(   25)=eri_value(   25)+d12bra(  3)*d03ket(  5)*(xin(  19)*yin(   1)*zin(   6)+xin(  43)*yin(  25)*zin(  30)+xin(  67)*yin(  49)*zin(  54)+xin(  91)*yin(  73)*zin(  78))
          eri_value(   26)=eri_value(   26)+d12bra(  3)*d03ket(  6)*(xin(  18)*yin(   3)*zin(   5)+xin(  42)*yin(  27)*zin(  29)+xin(  66)*yin(  51)*zin(  53)+xin(  90)*yin(  75)*zin(  77))
          eri_value(   27)=eri_value(   27)+d12bra(  3)*d03ket(  7)*(xin(  17)*yin(   3)*zin(   6)+xin(  41)*yin(  27)*zin(  30)+xin(  65)*yin(  51)*zin(  54)+xin(  89)*yin(  75)*zin(  78))
          eri_value(   28)=eri_value(   28)+d12bra(  3)*d03ket(  8)*(xin(  18)*yin(   1)*zin(   7)+xin(  42)*yin(  25)*zin(  31)+xin(  66)*yin(  49)*zin(  55)+xin(  90)*yin(  73)*zin(  79))
          eri_value(   29)=eri_value(   29)+d12bra(  3)*d03ket(  9)*(xin(  17)*yin(   2)*zin(   7)+xin(  41)*yin(  26)*zin(  31)+xin(  65)*yin(  50)*zin(  55)+xin(  89)*yin(  74)*zin(  79))
          eri_value(   30)=eri_value(   30)+d12bra(  3)*d03ket( 10)*(xin(  18)*yin(   2)*zin(   6)+xin(  42)*yin(  26)*zin(  30)+xin(  66)*yin(  50)*zin(  54)+xin(  90)*yin(  74)*zin(  78))
          eri_value(   31)=eri_value(   31)+d12bra(  4)*d03ket(  1)*(xin(   8)*yin(  17)*zin(   1)+xin(  32)*yin(  41)*zin(  25)+xin(  56)*yin(  65)*zin(  49)+xin(  80)*yin(  89)*zin(  73))
          eri_value(   32)=eri_value(   32)+d12bra(  4)*d03ket(  2)*(xin(   5)*yin(  20)*zin(   1)+xin(  29)*yin(  44)*zin(  25)+xin(  53)*yin(  68)*zin(  49)+xin(  77)*yin(  92)*zin(  73))
          eri_value(   33)=eri_value(   33)+d12bra(  4)*d03ket(  3)*(xin(   5)*yin(  17)*zin(   4)+xin(  29)*yin(  41)*zin(  28)+xin(  53)*yin(  65)*zin(  52)+xin(  77)*yin(  89)*zin(  76))
          eri_value(   34)=eri_value(   34)+d12bra(  4)*d03ket(  4)*(xin(   7)*yin(  18)*zin(   1)+xin(  31)*yin(  42)*zin(  25)+xin(  55)*yin(  66)*zin(  49)+xin(  79)*yin(  90)*zin(  73))
          eri_value(   35)=eri_value(   35)+d12bra(  4)*d03ket(  5)*(xin(   7)*yin(  17)*zin(   2)+xin(  31)*yin(  41)*zin(  26)+xin(  55)*yin(  65)*zin(  50)+xin(  79)*yin(  89)*zin(  74))
          eri_value(   36)=eri_value(   36)+d12bra(  4)*d03ket(  6)*(xin(   6)*yin(  19)*zin(   1)+xin(  30)*yin(  43)*zin(  25)+xin(  54)*yin(  67)*zin(  49)+xin(  78)*yin(  91)*zin(  73))
          eri_value(   37)=eri_value(   37)+d12bra(  4)*d03ket(  7)*(xin(   5)*yin(  19)*zin(   2)+xin(  29)*yin(  43)*zin(  26)+xin(  53)*yin(  67)*zin(  50)+xin(  77)*yin(  91)*zin(  74))
          eri_value(   38)=eri_value(   38)+d12bra(  4)*d03ket(  8)*(xin(   6)*yin(  17)*zin(   3)+xin(  30)*yin(  41)*zin(  27)+xin(  54)*yin(  65)*zin(  51)+xin(  78)*yin(  89)*zin(  75))
          eri_value(   39)=eri_value(   39)+d12bra(  4)*d03ket(  9)*(xin(   5)*yin(  18)*zin(   3)+xin(  29)*yin(  42)*zin(  27)+xin(  53)*yin(  66)*zin(  51)+xin(  77)*yin(  90)*zin(  75))
          eri_value(   40)=eri_value(   40)+d12bra(  4)*d03ket( 10)*(xin(   6)*yin(  18)*zin(   2)+xin(  30)*yin(  42)*zin(  26)+xin(  54)*yin(  66)*zin(  50)+xin(  78)*yin(  90)*zin(  74))
          eri_value(   41)=eri_value(   41)+d12bra(  5)*d03ket(  1)*(xin(   4)*yin(  21)*zin(   1)+xin(  28)*yin(  45)*zin(  25)+xin(  52)*yin(  69)*zin(  49)+xin(  76)*yin(  93)*zin(  73))
          eri_value(   42)=eri_value(   42)+d12bra(  5)*d03ket(  2)*(xin(   1)*yin(  24)*zin(   1)+xin(  25)*yin(  48)*zin(  25)+xin(  49)*yin(  72)*zin(  49)+xin(  73)*yin(  96)*zin(  73))
          eri_value(   43)=eri_value(   43)+d12bra(  5)*d03ket(  3)*(xin(   1)*yin(  21)*zin(   4)+xin(  25)*yin(  45)*zin(  28)+xin(  49)*yin(  69)*zin(  52)+xin(  73)*yin(  93)*zin(  76))
          eri_value(   44)=eri_value(   44)+d12bra(  5)*d03ket(  4)*(xin(   3)*yin(  22)*zin(   1)+xin(  27)*yin(  46)*zin(  25)+xin(  51)*yin(  70)*zin(  49)+xin(  75)*yin(  94)*zin(  73))
          eri_value(   45)=eri_value(   45)+d12bra(  5)*d03ket(  5)*(xin(   3)*yin(  21)*zin(   2)+xin(  27)*yin(  45)*zin(  26)+xin(  51)*yin(  69)*zin(  50)+xin(  75)*yin(  93)*zin(  74))
          eri_value(   46)=eri_value(   46)+d12bra(  5)*d03ket(  6)*(xin(   2)*yin(  23)*zin(   1)+xin(  26)*yin(  47)*zin(  25)+xin(  50)*yin(  71)*zin(  49)+xin(  74)*yin(  95)*zin(  73))
          eri_value(   47)=eri_value(   47)+d12bra(  5)*d03ket(  7)*(xin(   1)*yin(  23)*zin(   2)+xin(  25)*yin(  47)*zin(  26)+xin(  49)*yin(  71)*zin(  50)+xin(  73)*yin(  95)*zin(  74))
          eri_value(   48)=eri_value(   48)+d12bra(  5)*d03ket(  8)*(xin(   2)*yin(  21)*zin(   3)+xin(  26)*yin(  45)*zin(  27)+xin(  50)*yin(  69)*zin(  51)+xin(  74)*yin(  93)*zin(  75))
          eri_value(   49)=eri_value(   49)+d12bra(  5)*d03ket(  9)*(xin(   1)*yin(  22)*zin(   3)+xin(  25)*yin(  46)*zin(  27)+xin(  49)*yin(  70)*zin(  51)+xin(  73)*yin(  94)*zin(  75))
          eri_value(   50)=eri_value(   50)+d12bra(  5)*d03ket( 10)*(xin(   2)*yin(  22)*zin(   2)+xin(  26)*yin(  46)*zin(  26)+xin(  50)*yin(  70)*zin(  50)+xin(  74)*yin(  94)*zin(  74))
          eri_value(   51)=eri_value(   51)+d12bra(  6)*d03ket(  1)*(xin(   4)*yin(  17)*zin(   5)+xin(  28)*yin(  41)*zin(  29)+xin(  52)*yin(  65)*zin(  53)+xin(  76)*yin(  89)*zin(  77))
          eri_value(   52)=eri_value(   52)+d12bra(  6)*d03ket(  2)*(xin(   1)*yin(  20)*zin(   5)+xin(  25)*yin(  44)*zin(  29)+xin(  49)*yin(  68)*zin(  53)+xin(  73)*yin(  92)*zin(  77))
          eri_value(   53)=eri_value(   53)+d12bra(  6)*d03ket(  3)*(xin(   1)*yin(  17)*zin(   8)+xin(  25)*yin(  41)*zin(  32)+xin(  49)*yin(  65)*zin(  56)+xin(  73)*yin(  89)*zin(  80))
          eri_value(   54)=eri_value(   54)+d12bra(  6)*d03ket(  4)*(xin(   3)*yin(  18)*zin(   5)+xin(  27)*yin(  42)*zin(  29)+xin(  51)*yin(  66)*zin(  53)+xin(  75)*yin(  90)*zin(  77))
          eri_value(   55)=eri_value(   55)+d12bra(  6)*d03ket(  5)*(xin(   3)*yin(  17)*zin(   6)+xin(  27)*yin(  41)*zin(  30)+xin(  51)*yin(  65)*zin(  54)+xin(  75)*yin(  89)*zin(  78))
          eri_value(   56)=eri_value(   56)+d12bra(  6)*d03ket(  6)*(xin(   2)*yin(  19)*zin(   5)+xin(  26)*yin(  43)*zin(  29)+xin(  50)*yin(  67)*zin(  53)+xin(  74)*yin(  91)*zin(  77))
          eri_value(   57)=eri_value(   57)+d12bra(  6)*d03ket(  7)*(xin(   1)*yin(  19)*zin(   6)+xin(  25)*yin(  43)*zin(  30)+xin(  49)*yin(  67)*zin(  54)+xin(  73)*yin(  91)*zin(  78))
          eri_value(   58)=eri_value(   58)+d12bra(  6)*d03ket(  8)*(xin(   2)*yin(  17)*zin(   7)+xin(  26)*yin(  41)*zin(  31)+xin(  50)*yin(  65)*zin(  55)+xin(  74)*yin(  89)*zin(  79))
          eri_value(   59)=eri_value(   59)+d12bra(  6)*d03ket(  9)*(xin(   1)*yin(  18)*zin(   7)+xin(  25)*yin(  42)*zin(  31)+xin(  49)*yin(  66)*zin(  55)+xin(  73)*yin(  90)*zin(  79))
          eri_value(   60)=eri_value(   60)+d12bra(  6)*d03ket( 10)*(xin(   2)*yin(  18)*zin(   6)+xin(  26)*yin(  42)*zin(  30)+xin(  50)*yin(  66)*zin(  54)+xin(  74)*yin(  90)*zin(  78))
          eri_value(   61)=eri_value(   61)+d12bra(  7)*d03ket(  1)*(xin(   8)*yin(   1)*zin(  17)+xin(  32)*yin(  25)*zin(  41)+xin(  56)*yin(  49)*zin(  65)+xin(  80)*yin(  73)*zin(  89))
          eri_value(   62)=eri_value(   62)+d12bra(  7)*d03ket(  2)*(xin(   5)*yin(   4)*zin(  17)+xin(  29)*yin(  28)*zin(  41)+xin(  53)*yin(  52)*zin(  65)+xin(  77)*yin(  76)*zin(  89))
          eri_value(   63)=eri_value(   63)+d12bra(  7)*d03ket(  3)*(xin(   5)*yin(   1)*zin(  20)+xin(  29)*yin(  25)*zin(  44)+xin(  53)*yin(  49)*zin(  68)+xin(  77)*yin(  73)*zin(  92))
          eri_value(   64)=eri_value(   64)+d12bra(  7)*d03ket(  4)*(xin(   7)*yin(   2)*zin(  17)+xin(  31)*yin(  26)*zin(  41)+xin(  55)*yin(  50)*zin(  65)+xin(  79)*yin(  74)*zin(  89))
          eri_value(   65)=eri_value(   65)+d12bra(  7)*d03ket(  5)*(xin(   7)*yin(   1)*zin(  18)+xin(  31)*yin(  25)*zin(  42)+xin(  55)*yin(  49)*zin(  66)+xin(  79)*yin(  73)*zin(  90))
          eri_value(   66)=eri_value(   66)+d12bra(  7)*d03ket(  6)*(xin(   6)*yin(   3)*zin(  17)+xin(  30)*yin(  27)*zin(  41)+xin(  54)*yin(  51)*zin(  65)+xin(  78)*yin(  75)*zin(  89))
          eri_value(   67)=eri_value(   67)+d12bra(  7)*d03ket(  7)*(xin(   5)*yin(   3)*zin(  18)+xin(  29)*yin(  27)*zin(  42)+xin(  53)*yin(  51)*zin(  66)+xin(  77)*yin(  75)*zin(  90))
          eri_value(   68)=eri_value(   68)+d12bra(  7)*d03ket(  8)*(xin(   6)*yin(   1)*zin(  19)+xin(  30)*yin(  25)*zin(  43)+xin(  54)*yin(  49)*zin(  67)+xin(  78)*yin(  73)*zin(  91))
          eri_value(   69)=eri_value(   69)+d12bra(  7)*d03ket(  9)*(xin(   5)*yin(   2)*zin(  19)+xin(  29)*yin(  26)*zin(  43)+xin(  53)*yin(  50)*zin(  67)+xin(  77)*yin(  74)*zin(  91))
          eri_value(   70)=eri_value(   70)+d12bra(  7)*d03ket( 10)*(xin(   6)*yin(   2)*zin(  18)+xin(  30)*yin(  26)*zin(  42)+xin(  54)*yin(  50)*zin(  66)+xin(  78)*yin(  74)*zin(  90))
          eri_value(   71)=eri_value(   71)+d12bra(  8)*d03ket(  1)*(xin(   4)*yin(   5)*zin(  17)+xin(  28)*yin(  29)*zin(  41)+xin(  52)*yin(  53)*zin(  65)+xin(  76)*yin(  77)*zin(  89))
          eri_value(   72)=eri_value(   72)+d12bra(  8)*d03ket(  2)*(xin(   1)*yin(   8)*zin(  17)+xin(  25)*yin(  32)*zin(  41)+xin(  49)*yin(  56)*zin(  65)+xin(  73)*yin(  80)*zin(  89))
          eri_value(   73)=eri_value(   73)+d12bra(  8)*d03ket(  3)*(xin(   1)*yin(   5)*zin(  20)+xin(  25)*yin(  29)*zin(  44)+xin(  49)*yin(  53)*zin(  68)+xin(  73)*yin(  77)*zin(  92))
          eri_value(   74)=eri_value(   74)+d12bra(  8)*d03ket(  4)*(xin(   3)*yin(   6)*zin(  17)+xin(  27)*yin(  30)*zin(  41)+xin(  51)*yin(  54)*zin(  65)+xin(  75)*yin(  78)*zin(  89))
          eri_value(   75)=eri_value(   75)+d12bra(  8)*d03ket(  5)*(xin(   3)*yin(   5)*zin(  18)+xin(  27)*yin(  29)*zin(  42)+xin(  51)*yin(  53)*zin(  66)+xin(  75)*yin(  77)*zin(  90))
          eri_value(   76)=eri_value(   76)+d12bra(  8)*d03ket(  6)*(xin(   2)*yin(   7)*zin(  17)+xin(  26)*yin(  31)*zin(  41)+xin(  50)*yin(  55)*zin(  65)+xin(  74)*yin(  79)*zin(  89))
          eri_value(   77)=eri_value(   77)+d12bra(  8)*d03ket(  7)*(xin(   1)*yin(   7)*zin(  18)+xin(  25)*yin(  31)*zin(  42)+xin(  49)*yin(  55)*zin(  66)+xin(  73)*yin(  79)*zin(  90))
          eri_value(   78)=eri_value(   78)+d12bra(  8)*d03ket(  8)*(xin(   2)*yin(   5)*zin(  19)+xin(  26)*yin(  29)*zin(  43)+xin(  50)*yin(  53)*zin(  67)+xin(  74)*yin(  77)*zin(  91))
          eri_value(   79)=eri_value(   79)+d12bra(  8)*d03ket(  9)*(xin(   1)*yin(   6)*zin(  19)+xin(  25)*yin(  30)*zin(  43)+xin(  49)*yin(  54)*zin(  67)+xin(  73)*yin(  78)*zin(  91))
          eri_value(   80)=eri_value(   80)+d12bra(  8)*d03ket( 10)*(xin(   2)*yin(   6)*zin(  18)+xin(  26)*yin(  30)*zin(  42)+xin(  50)*yin(  54)*zin(  66)+xin(  74)*yin(  78)*zin(  90))
          eri_value(   81)=eri_value(   81)+d12bra(  9)*d03ket(  1)*(xin(   4)*yin(   1)*zin(  21)+xin(  28)*yin(  25)*zin(  45)+xin(  52)*yin(  49)*zin(  69)+xin(  76)*yin(  73)*zin(  93))
          eri_value(   82)=eri_value(   82)+d12bra(  9)*d03ket(  2)*(xin(   1)*yin(   4)*zin(  21)+xin(  25)*yin(  28)*zin(  45)+xin(  49)*yin(  52)*zin(  69)+xin(  73)*yin(  76)*zin(  93))
          eri_value(   83)=eri_value(   83)+d12bra(  9)*d03ket(  3)*(xin(   1)*yin(   1)*zin(  24)+xin(  25)*yin(  25)*zin(  48)+xin(  49)*yin(  49)*zin(  72)+xin(  73)*yin(  73)*zin(  96))
          eri_value(   84)=eri_value(   84)+d12bra(  9)*d03ket(  4)*(xin(   3)*yin(   2)*zin(  21)+xin(  27)*yin(  26)*zin(  45)+xin(  51)*yin(  50)*zin(  69)+xin(  75)*yin(  74)*zin(  93))
          eri_value(   85)=eri_value(   85)+d12bra(  9)*d03ket(  5)*(xin(   3)*yin(   1)*zin(  22)+xin(  27)*yin(  25)*zin(  46)+xin(  51)*yin(  49)*zin(  70)+xin(  75)*yin(  73)*zin(  94))
          eri_value(   86)=eri_value(   86)+d12bra(  9)*d03ket(  6)*(xin(   2)*yin(   3)*zin(  21)+xin(  26)*yin(  27)*zin(  45)+xin(  50)*yin(  51)*zin(  69)+xin(  74)*yin(  75)*zin(  93))
          eri_value(   87)=eri_value(   87)+d12bra(  9)*d03ket(  7)*(xin(   1)*yin(   3)*zin(  22)+xin(  25)*yin(  27)*zin(  46)+xin(  49)*yin(  51)*zin(  70)+xin(  73)*yin(  75)*zin(  94))
          eri_value(   88)=eri_value(   88)+d12bra(  9)*d03ket(  8)*(xin(   2)*yin(   1)*zin(  23)+xin(  26)*yin(  25)*zin(  47)+xin(  50)*yin(  49)*zin(  71)+xin(  74)*yin(  73)*zin(  95))
          eri_value(   89)=eri_value(   89)+d12bra(  9)*d03ket(  9)*(xin(   1)*yin(   2)*zin(  23)+xin(  25)*yin(  26)*zin(  47)+xin(  49)*yin(  50)*zin(  71)+xin(  73)*yin(  74)*zin(  95))
          eri_value(   90)=eri_value(   90)+d12bra(  9)*d03ket( 10)*(xin(   2)*yin(   2)*zin(  22)+xin(  26)*yin(  26)*zin(  46)+xin(  50)*yin(  50)*zin(  70)+xin(  74)*yin(  74)*zin(  94))
          eri_value(   91)=eri_value(   91)+d12bra( 10)*d03ket(  1)*(xin(  16)*yin(   9)*zin(   1)+xin(  40)*yin(  33)*zin(  25)+xin(  64)*yin(  57)*zin(  49)+xin(  88)*yin(  81)*zin(  73))
          eri_value(   92)=eri_value(   92)+d12bra( 10)*d03ket(  2)*(xin(  13)*yin(  12)*zin(   1)+xin(  37)*yin(  36)*zin(  25)+xin(  61)*yin(  60)*zin(  49)+xin(  85)*yin(  84)*zin(  73))
          eri_value(   93)=eri_value(   93)+d12bra( 10)*d03ket(  3)*(xin(  13)*yin(   9)*zin(   4)+xin(  37)*yin(  33)*zin(  28)+xin(  61)*yin(  57)*zin(  52)+xin(  85)*yin(  81)*zin(  76))
          eri_value(   94)=eri_value(   94)+d12bra( 10)*d03ket(  4)*(xin(  15)*yin(  10)*zin(   1)+xin(  39)*yin(  34)*zin(  25)+xin(  63)*yin(  58)*zin(  49)+xin(  87)*yin(  82)*zin(  73))
          eri_value(   95)=eri_value(   95)+d12bra( 10)*d03ket(  5)*(xin(  15)*yin(   9)*zin(   2)+xin(  39)*yin(  33)*zin(  26)+xin(  63)*yin(  57)*zin(  50)+xin(  87)*yin(  81)*zin(  74))
          eri_value(   96)=eri_value(   96)+d12bra( 10)*d03ket(  6)*(xin(  14)*yin(  11)*zin(   1)+xin(  38)*yin(  35)*zin(  25)+xin(  62)*yin(  59)*zin(  49)+xin(  86)*yin(  83)*zin(  73))
          eri_value(   97)=eri_value(   97)+d12bra( 10)*d03ket(  7)*(xin(  13)*yin(  11)*zin(   2)+xin(  37)*yin(  35)*zin(  26)+xin(  61)*yin(  59)*zin(  50)+xin(  85)*yin(  83)*zin(  74))
          eri_value(   98)=eri_value(   98)+d12bra( 10)*d03ket(  8)*(xin(  14)*yin(   9)*zin(   3)+xin(  38)*yin(  33)*zin(  27)+xin(  62)*yin(  57)*zin(  51)+xin(  86)*yin(  81)*zin(  75))
          eri_value(   99)=eri_value(   99)+d12bra( 10)*d03ket(  9)*(xin(  13)*yin(  10)*zin(   3)+xin(  37)*yin(  34)*zin(  27)+xin(  61)*yin(  58)*zin(  51)+xin(  85)*yin(  82)*zin(  75))
          eri_value(  100)=eri_value(  100)+d12bra( 10)*d03ket( 10)*(xin(  14)*yin(  10)*zin(   2)+xin(  38)*yin(  34)*zin(  26)+xin(  62)*yin(  58)*zin(  50)+xin(  86)*yin(  82)*zin(  74))
          eri_value(  101)=eri_value(  101)+d12bra( 11)*d03ket(  1)*(xin(  12)*yin(  13)*zin(   1)+xin(  36)*yin(  37)*zin(  25)+xin(  60)*yin(  61)*zin(  49)+xin(  84)*yin(  85)*zin(  73))
          eri_value(  102)=eri_value(  102)+d12bra( 11)*d03ket(  2)*(xin(   9)*yin(  16)*zin(   1)+xin(  33)*yin(  40)*zin(  25)+xin(  57)*yin(  64)*zin(  49)+xin(  81)*yin(  88)*zin(  73))
          eri_value(  103)=eri_value(  103)+d12bra( 11)*d03ket(  3)*(xin(   9)*yin(  13)*zin(   4)+xin(  33)*yin(  37)*zin(  28)+xin(  57)*yin(  61)*zin(  52)+xin(  81)*yin(  85)*zin(  76))
          eri_value(  104)=eri_value(  104)+d12bra( 11)*d03ket(  4)*(xin(  11)*yin(  14)*zin(   1)+xin(  35)*yin(  38)*zin(  25)+xin(  59)*yin(  62)*zin(  49)+xin(  83)*yin(  86)*zin(  73))
          eri_value(  105)=eri_value(  105)+d12bra( 11)*d03ket(  5)*(xin(  11)*yin(  13)*zin(   2)+xin(  35)*yin(  37)*zin(  26)+xin(  59)*yin(  61)*zin(  50)+xin(  83)*yin(  85)*zin(  74))
          eri_value(  106)=eri_value(  106)+d12bra( 11)*d03ket(  6)*(xin(  10)*yin(  15)*zin(   1)+xin(  34)*yin(  39)*zin(  25)+xin(  58)*yin(  63)*zin(  49)+xin(  82)*yin(  87)*zin(  73))
          eri_value(  107)=eri_value(  107)+d12bra( 11)*d03ket(  7)*(xin(   9)*yin(  15)*zin(   2)+xin(  33)*yin(  39)*zin(  26)+xin(  57)*yin(  63)*zin(  50)+xin(  81)*yin(  87)*zin(  74))
          eri_value(  108)=eri_value(  108)+d12bra( 11)*d03ket(  8)*(xin(  10)*yin(  13)*zin(   3)+xin(  34)*yin(  37)*zin(  27)+xin(  58)*yin(  61)*zin(  51)+xin(  82)*yin(  85)*zin(  75))
          eri_value(  109)=eri_value(  109)+d12bra( 11)*d03ket(  9)*(xin(   9)*yin(  14)*zin(   3)+xin(  33)*yin(  38)*zin(  27)+xin(  57)*yin(  62)*zin(  51)+xin(  81)*yin(  86)*zin(  75))
          eri_value(  110)=eri_value(  110)+d12bra( 11)*d03ket( 10)*(xin(  10)*yin(  14)*zin(   2)+xin(  34)*yin(  38)*zin(  26)+xin(  58)*yin(  62)*zin(  50)+xin(  82)*yin(  86)*zin(  74))
          eri_value(  111)=eri_value(  111)+d12bra( 12)*d03ket(  1)*(xin(  12)*yin(   9)*zin(   5)+xin(  36)*yin(  33)*zin(  29)+xin(  60)*yin(  57)*zin(  53)+xin(  84)*yin(  81)*zin(  77))
          eri_value(  112)=eri_value(  112)+d12bra( 12)*d03ket(  2)*(xin(   9)*yin(  12)*zin(   5)+xin(  33)*yin(  36)*zin(  29)+xin(  57)*yin(  60)*zin(  53)+xin(  81)*yin(  84)*zin(  77))
          eri_value(  113)=eri_value(  113)+d12bra( 12)*d03ket(  3)*(xin(   9)*yin(   9)*zin(   8)+xin(  33)*yin(  33)*zin(  32)+xin(  57)*yin(  57)*zin(  56)+xin(  81)*yin(  81)*zin(  80))
          eri_value(  114)=eri_value(  114)+d12bra( 12)*d03ket(  4)*(xin(  11)*yin(  10)*zin(   5)+xin(  35)*yin(  34)*zin(  29)+xin(  59)*yin(  58)*zin(  53)+xin(  83)*yin(  82)*zin(  77))
          eri_value(  115)=eri_value(  115)+d12bra( 12)*d03ket(  5)*(xin(  11)*yin(   9)*zin(   6)+xin(  35)*yin(  33)*zin(  30)+xin(  59)*yin(  57)*zin(  54)+xin(  83)*yin(  81)*zin(  78))
          eri_value(  116)=eri_value(  116)+d12bra( 12)*d03ket(  6)*(xin(  10)*yin(  11)*zin(   5)+xin(  34)*yin(  35)*zin(  29)+xin(  58)*yin(  59)*zin(  53)+xin(  82)*yin(  83)*zin(  77))
          eri_value(  117)=eri_value(  117)+d12bra( 12)*d03ket(  7)*(xin(   9)*yin(  11)*zin(   6)+xin(  33)*yin(  35)*zin(  30)+xin(  57)*yin(  59)*zin(  54)+xin(  81)*yin(  83)*zin(  78))
          eri_value(  118)=eri_value(  118)+d12bra( 12)*d03ket(  8)*(xin(  10)*yin(   9)*zin(   7)+xin(  34)*yin(  33)*zin(  31)+xin(  58)*yin(  57)*zin(  55)+xin(  82)*yin(  81)*zin(  79))
          eri_value(  119)=eri_value(  119)+d12bra( 12)*d03ket(  9)*(xin(   9)*yin(  10)*zin(   7)+xin(  33)*yin(  34)*zin(  31)+xin(  57)*yin(  58)*zin(  55)+xin(  81)*yin(  82)*zin(  79))
          eri_value(  120)=eri_value(  120)+d12bra( 12)*d03ket( 10)*(xin(  10)*yin(  10)*zin(   6)+xin(  34)*yin(  34)*zin(  30)+xin(  58)*yin(  58)*zin(  54)+xin(  82)*yin(  82)*zin(  78))
          eri_value(  121)=eri_value(  121)+d12bra( 13)*d03ket(  1)*(xin(  16)*yin(   1)*zin(   9)+xin(  40)*yin(  25)*zin(  33)+xin(  64)*yin(  49)*zin(  57)+xin(  88)*yin(  73)*zin(  81))
          eri_value(  122)=eri_value(  122)+d12bra( 13)*d03ket(  2)*(xin(  13)*yin(   4)*zin(   9)+xin(  37)*yin(  28)*zin(  33)+xin(  61)*yin(  52)*zin(  57)+xin(  85)*yin(  76)*zin(  81))
          eri_value(  123)=eri_value(  123)+d12bra( 13)*d03ket(  3)*(xin(  13)*yin(   1)*zin(  12)+xin(  37)*yin(  25)*zin(  36)+xin(  61)*yin(  49)*zin(  60)+xin(  85)*yin(  73)*zin(  84))
          eri_value(  124)=eri_value(  124)+d12bra( 13)*d03ket(  4)*(xin(  15)*yin(   2)*zin(   9)+xin(  39)*yin(  26)*zin(  33)+xin(  63)*yin(  50)*zin(  57)+xin(  87)*yin(  74)*zin(  81))
          eri_value(  125)=eri_value(  125)+d12bra( 13)*d03ket(  5)*(xin(  15)*yin(   1)*zin(  10)+xin(  39)*yin(  25)*zin(  34)+xin(  63)*yin(  49)*zin(  58)+xin(  87)*yin(  73)*zin(  82))
          eri_value(  126)=eri_value(  126)+d12bra( 13)*d03ket(  6)*(xin(  14)*yin(   3)*zin(   9)+xin(  38)*yin(  27)*zin(  33)+xin(  62)*yin(  51)*zin(  57)+xin(  86)*yin(  75)*zin(  81))
          eri_value(  127)=eri_value(  127)+d12bra( 13)*d03ket(  7)*(xin(  13)*yin(   3)*zin(  10)+xin(  37)*yin(  27)*zin(  34)+xin(  61)*yin(  51)*zin(  58)+xin(  85)*yin(  75)*zin(  82))
          eri_value(  128)=eri_value(  128)+d12bra( 13)*d03ket(  8)*(xin(  14)*yin(   1)*zin(  11)+xin(  38)*yin(  25)*zin(  35)+xin(  62)*yin(  49)*zin(  59)+xin(  86)*yin(  73)*zin(  83))
          eri_value(  129)=eri_value(  129)+d12bra( 13)*d03ket(  9)*(xin(  13)*yin(   2)*zin(  11)+xin(  37)*yin(  26)*zin(  35)+xin(  61)*yin(  50)*zin(  59)+xin(  85)*yin(  74)*zin(  83))
          eri_value(  130)=eri_value(  130)+d12bra( 13)*d03ket( 10)*(xin(  14)*yin(   2)*zin(  10)+xin(  38)*yin(  26)*zin(  34)+xin(  62)*yin(  50)*zin(  58)+xin(  86)*yin(  74)*zin(  82))
          eri_value(  131)=eri_value(  131)+d12bra( 14)*d03ket(  1)*(xin(  12)*yin(   5)*zin(   9)+xin(  36)*yin(  29)*zin(  33)+xin(  60)*yin(  53)*zin(  57)+xin(  84)*yin(  77)*zin(  81))
          eri_value(  132)=eri_value(  132)+d12bra( 14)*d03ket(  2)*(xin(   9)*yin(   8)*zin(   9)+xin(  33)*yin(  32)*zin(  33)+xin(  57)*yin(  56)*zin(  57)+xin(  81)*yin(  80)*zin(  81))
          eri_value(  133)=eri_value(  133)+d12bra( 14)*d03ket(  3)*(xin(   9)*yin(   5)*zin(  12)+xin(  33)*yin(  29)*zin(  36)+xin(  57)*yin(  53)*zin(  60)+xin(  81)*yin(  77)*zin(  84))
          eri_value(  134)=eri_value(  134)+d12bra( 14)*d03ket(  4)*(xin(  11)*yin(   6)*zin(   9)+xin(  35)*yin(  30)*zin(  33)+xin(  59)*yin(  54)*zin(  57)+xin(  83)*yin(  78)*zin(  81))
          eri_value(  135)=eri_value(  135)+d12bra( 14)*d03ket(  5)*(xin(  11)*yin(   5)*zin(  10)+xin(  35)*yin(  29)*zin(  34)+xin(  59)*yin(  53)*zin(  58)+xin(  83)*yin(  77)*zin(  82))
          eri_value(  136)=eri_value(  136)+d12bra( 14)*d03ket(  6)*(xin(  10)*yin(   7)*zin(   9)+xin(  34)*yin(  31)*zin(  33)+xin(  58)*yin(  55)*zin(  57)+xin(  82)*yin(  79)*zin(  81))
          eri_value(  137)=eri_value(  137)+d12bra( 14)*d03ket(  7)*(xin(   9)*yin(   7)*zin(  10)+xin(  33)*yin(  31)*zin(  34)+xin(  57)*yin(  55)*zin(  58)+xin(  81)*yin(  79)*zin(  82))
          eri_value(  138)=eri_value(  138)+d12bra( 14)*d03ket(  8)*(xin(  10)*yin(   5)*zin(  11)+xin(  34)*yin(  29)*zin(  35)+xin(  58)*yin(  53)*zin(  59)+xin(  82)*yin(  77)*zin(  83))
          eri_value(  139)=eri_value(  139)+d12bra( 14)*d03ket(  9)*(xin(   9)*yin(   6)*zin(  11)+xin(  33)*yin(  30)*zin(  35)+xin(  57)*yin(  54)*zin(  59)+xin(  81)*yin(  78)*zin(  83))
          eri_value(  140)=eri_value(  140)+d12bra( 14)*d03ket( 10)*(xin(  10)*yin(   6)*zin(  10)+xin(  34)*yin(  30)*zin(  34)+xin(  58)*yin(  54)*zin(  58)+xin(  82)*yin(  78)*zin(  82))
          eri_value(  141)=eri_value(  141)+d12bra( 15)*d03ket(  1)*(xin(  12)*yin(   1)*zin(  13)+xin(  36)*yin(  25)*zin(  37)+xin(  60)*yin(  49)*zin(  61)+xin(  84)*yin(  73)*zin(  85))
          eri_value(  142)=eri_value(  142)+d12bra( 15)*d03ket(  2)*(xin(   9)*yin(   4)*zin(  13)+xin(  33)*yin(  28)*zin(  37)+xin(  57)*yin(  52)*zin(  61)+xin(  81)*yin(  76)*zin(  85))
          eri_value(  143)=eri_value(  143)+d12bra( 15)*d03ket(  3)*(xin(   9)*yin(   1)*zin(  16)+xin(  33)*yin(  25)*zin(  40)+xin(  57)*yin(  49)*zin(  64)+xin(  81)*yin(  73)*zin(  88))
          eri_value(  144)=eri_value(  144)+d12bra( 15)*d03ket(  4)*(xin(  11)*yin(   2)*zin(  13)+xin(  35)*yin(  26)*zin(  37)+xin(  59)*yin(  50)*zin(  61)+xin(  83)*yin(  74)*zin(  85))
          eri_value(  145)=eri_value(  145)+d12bra( 15)*d03ket(  5)*(xin(  11)*yin(   1)*zin(  14)+xin(  35)*yin(  25)*zin(  38)+xin(  59)*yin(  49)*zin(  62)+xin(  83)*yin(  73)*zin(  86))
          eri_value(  146)=eri_value(  146)+d12bra( 15)*d03ket(  6)*(xin(  10)*yin(   3)*zin(  13)+xin(  34)*yin(  27)*zin(  37)+xin(  58)*yin(  51)*zin(  61)+xin(  82)*yin(  75)*zin(  85))
          eri_value(  147)=eri_value(  147)+d12bra( 15)*d03ket(  7)*(xin(   9)*yin(   3)*zin(  14)+xin(  33)*yin(  27)*zin(  38)+xin(  57)*yin(  51)*zin(  62)+xin(  81)*yin(  75)*zin(  86))
          eri_value(  148)=eri_value(  148)+d12bra( 15)*d03ket(  8)*(xin(  10)*yin(   1)*zin(  15)+xin(  34)*yin(  25)*zin(  39)+xin(  58)*yin(  49)*zin(  63)+xin(  82)*yin(  73)*zin(  87))
          eri_value(  149)=eri_value(  149)+d12bra( 15)*d03ket(  9)*(xin(   9)*yin(   2)*zin(  15)+xin(  33)*yin(  26)*zin(  39)+xin(  57)*yin(  50)*zin(  63)+xin(  81)*yin(  74)*zin(  87))
          eri_value(  150)=eri_value(  150)+d12bra( 15)*d03ket( 10)*(xin(  10)*yin(   2)*zin(  14)+xin(  34)*yin(  26)*zin(  38)+xin(  58)*yin(  50)*zin(  62)+xin(  82)*yin(  74)*zin(  86))
          eri_value(  151)=eri_value(  151)+d12bra( 16)*d03ket(  1)*(xin(   8)*yin(   9)*zin(   9)+xin(  32)*yin(  33)*zin(  33)+xin(  56)*yin(  57)*zin(  57)+xin(  80)*yin(  81)*zin(  81))
          eri_value(  152)=eri_value(  152)+d12bra( 16)*d03ket(  2)*(xin(   5)*yin(  12)*zin(   9)+xin(  29)*yin(  36)*zin(  33)+xin(  53)*yin(  60)*zin(  57)+xin(  77)*yin(  84)*zin(  81))
          eri_value(  153)=eri_value(  153)+d12bra( 16)*d03ket(  3)*(xin(   5)*yin(   9)*zin(  12)+xin(  29)*yin(  33)*zin(  36)+xin(  53)*yin(  57)*zin(  60)+xin(  77)*yin(  81)*zin(  84))
          eri_value(  154)=eri_value(  154)+d12bra( 16)*d03ket(  4)*(xin(   7)*yin(  10)*zin(   9)+xin(  31)*yin(  34)*zin(  33)+xin(  55)*yin(  58)*zin(  57)+xin(  79)*yin(  82)*zin(  81))
          eri_value(  155)=eri_value(  155)+d12bra( 16)*d03ket(  5)*(xin(   7)*yin(   9)*zin(  10)+xin(  31)*yin(  33)*zin(  34)+xin(  55)*yin(  57)*zin(  58)+xin(  79)*yin(  81)*zin(  82))
          eri_value(  156)=eri_value(  156)+d12bra( 16)*d03ket(  6)*(xin(   6)*yin(  11)*zin(   9)+xin(  30)*yin(  35)*zin(  33)+xin(  54)*yin(  59)*zin(  57)+xin(  78)*yin(  83)*zin(  81))
          eri_value(  157)=eri_value(  157)+d12bra( 16)*d03ket(  7)*(xin(   5)*yin(  11)*zin(  10)+xin(  29)*yin(  35)*zin(  34)+xin(  53)*yin(  59)*zin(  58)+xin(  77)*yin(  83)*zin(  82))
          eri_value(  158)=eri_value(  158)+d12bra( 16)*d03ket(  8)*(xin(   6)*yin(   9)*zin(  11)+xin(  30)*yin(  33)*zin(  35)+xin(  54)*yin(  57)*zin(  59)+xin(  78)*yin(  81)*zin(  83))
          eri_value(  159)=eri_value(  159)+d12bra( 16)*d03ket(  9)*(xin(   5)*yin(  10)*zin(  11)+xin(  29)*yin(  34)*zin(  35)+xin(  53)*yin(  58)*zin(  59)+xin(  77)*yin(  82)*zin(  83))
          eri_value(  160)=eri_value(  160)+d12bra( 16)*d03ket( 10)*(xin(   6)*yin(  10)*zin(  10)+xin(  30)*yin(  34)*zin(  34)+xin(  54)*yin(  58)*zin(  58)+xin(  78)*yin(  82)*zin(  82))
          eri_value(  161)=eri_value(  161)+d12bra( 17)*d03ket(  1)*(xin(   4)*yin(  13)*zin(   9)+xin(  28)*yin(  37)*zin(  33)+xin(  52)*yin(  61)*zin(  57)+xin(  76)*yin(  85)*zin(  81))
          eri_value(  162)=eri_value(  162)+d12bra( 17)*d03ket(  2)*(xin(   1)*yin(  16)*zin(   9)+xin(  25)*yin(  40)*zin(  33)+xin(  49)*yin(  64)*zin(  57)+xin(  73)*yin(  88)*zin(  81))
          eri_value(  163)=eri_value(  163)+d12bra( 17)*d03ket(  3)*(xin(   1)*yin(  13)*zin(  12)+xin(  25)*yin(  37)*zin(  36)+xin(  49)*yin(  61)*zin(  60)+xin(  73)*yin(  85)*zin(  84))
          eri_value(  164)=eri_value(  164)+d12bra( 17)*d03ket(  4)*(xin(   3)*yin(  14)*zin(   9)+xin(  27)*yin(  38)*zin(  33)+xin(  51)*yin(  62)*zin(  57)+xin(  75)*yin(  86)*zin(  81))
          eri_value(  165)=eri_value(  165)+d12bra( 17)*d03ket(  5)*(xin(   3)*yin(  13)*zin(  10)+xin(  27)*yin(  37)*zin(  34)+xin(  51)*yin(  61)*zin(  58)+xin(  75)*yin(  85)*zin(  82))
          eri_value(  166)=eri_value(  166)+d12bra( 17)*d03ket(  6)*(xin(   2)*yin(  15)*zin(   9)+xin(  26)*yin(  39)*zin(  33)+xin(  50)*yin(  63)*zin(  57)+xin(  74)*yin(  87)*zin(  81))
          eri_value(  167)=eri_value(  167)+d12bra( 17)*d03ket(  7)*(xin(   1)*yin(  15)*zin(  10)+xin(  25)*yin(  39)*zin(  34)+xin(  49)*yin(  63)*zin(  58)+xin(  73)*yin(  87)*zin(  82))
          eri_value(  168)=eri_value(  168)+d12bra( 17)*d03ket(  8)*(xin(   2)*yin(  13)*zin(  11)+xin(  26)*yin(  37)*zin(  35)+xin(  50)*yin(  61)*zin(  59)+xin(  74)*yin(  85)*zin(  83))
          eri_value(  169)=eri_value(  169)+d12bra( 17)*d03ket(  9)*(xin(   1)*yin(  14)*zin(  11)+xin(  25)*yin(  38)*zin(  35)+xin(  49)*yin(  62)*zin(  59)+xin(  73)*yin(  86)*zin(  83))
          eri_value(  170)=eri_value(  170)+d12bra( 17)*d03ket( 10)*(xin(   2)*yin(  14)*zin(  10)+xin(  26)*yin(  38)*zin(  34)+xin(  50)*yin(  62)*zin(  58)+xin(  74)*yin(  86)*zin(  82))
          eri_value(  171)=eri_value(  171)+d12bra( 18)*d03ket(  1)*(xin(   4)*yin(   9)*zin(  13)+xin(  28)*yin(  33)*zin(  37)+xin(  52)*yin(  57)*zin(  61)+xin(  76)*yin(  81)*zin(  85))
          eri_value(  172)=eri_value(  172)+d12bra( 18)*d03ket(  2)*(xin(   1)*yin(  12)*zin(  13)+xin(  25)*yin(  36)*zin(  37)+xin(  49)*yin(  60)*zin(  61)+xin(  73)*yin(  84)*zin(  85))
          eri_value(  173)=eri_value(  173)+d12bra( 18)*d03ket(  3)*(xin(   1)*yin(   9)*zin(  16)+xin(  25)*yin(  33)*zin(  40)+xin(  49)*yin(  57)*zin(  64)+xin(  73)*yin(  81)*zin(  88))
          eri_value(  174)=eri_value(  174)+d12bra( 18)*d03ket(  4)*(xin(   3)*yin(  10)*zin(  13)+xin(  27)*yin(  34)*zin(  37)+xin(  51)*yin(  58)*zin(  61)+xin(  75)*yin(  82)*zin(  85))
          eri_value(  175)=eri_value(  175)+d12bra( 18)*d03ket(  5)*(xin(   3)*yin(   9)*zin(  14)+xin(  27)*yin(  33)*zin(  38)+xin(  51)*yin(  57)*zin(  62)+xin(  75)*yin(  81)*zin(  86))
          eri_value(  176)=eri_value(  176)+d12bra( 18)*d03ket(  6)*(xin(   2)*yin(  11)*zin(  13)+xin(  26)*yin(  35)*zin(  37)+xin(  50)*yin(  59)*zin(  61)+xin(  74)*yin(  83)*zin(  85))
          eri_value(  177)=eri_value(  177)+d12bra( 18)*d03ket(  7)*(xin(   1)*yin(  11)*zin(  14)+xin(  25)*yin(  35)*zin(  38)+xin(  49)*yin(  59)*zin(  62)+xin(  73)*yin(  83)*zin(  86))
          eri_value(  178)=eri_value(  178)+d12bra( 18)*d03ket(  8)*(xin(   2)*yin(   9)*zin(  15)+xin(  26)*yin(  33)*zin(  39)+xin(  50)*yin(  57)*zin(  63)+xin(  74)*yin(  81)*zin(  87))
          eri_value(  179)=eri_value(  179)+d12bra( 18)*d03ket(  9)*(xin(   1)*yin(  10)*zin(  15)+xin(  25)*yin(  34)*zin(  39)+xin(  49)*yin(  58)*zin(  63)+xin(  73)*yin(  82)*zin(  87))
          eri_value(  180)=eri_value(  180)+d12bra( 18)*d03ket( 10)*(xin(   2)*yin(  10)*zin(  14)+xin(  26)*yin(  34)*zin(  38)+xin(  50)*yin(  58)*zin(  62)+xin(  74)*yin(  82)*zin(  86))

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
                                    ip = (i - 1)*30 ! Stride between functions in i

                                    do j = 1, 3 ! # of cartesians in j

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

                              deallocate (n12bra)
                              deallocate (xint12bra)
                              deallocate (n03ket)
                              deallocate (xint03ket)

                              end subroutine int2130
                              end submodule
