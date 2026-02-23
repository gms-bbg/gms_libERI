! The total angular momentum of this class is:           7
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int2230_impl
contains
  module subroutine int2230(dd_pair, sf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: dd_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n22bra(:), n03ket(:)
    real(dp), allocatable :: xint22bra(:), xint03ket(:)
    integer(kind=int64) :: nddbra, nsfket
    real(dp) :: scutddbra, scutsfket, test
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(144), yin(144), zin(144)
    real(dp) :: eri_value(360)
    real(dp) :: d22bra(36), d03ket(10)
    integer(kind=int64) :: ix(6), jx(6), kx(10), lx(1)
    integer(kind=int64) :: iy(6), jy(6), ky(10), ly(1)
    integer(kind=int64) :: iz(6), jz(6), kz(10), lz(1)
    integer(kind=int64) :: in(5), in1(5), kn(4)
    integer(kind=int64) :: ijx(36), ijy(36), ijz(36)
    integer(kind=int64) :: klx(10), kly(10), klz(10)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 29
    in1(5) = 33

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

    jx(1) = 8
    jx(2) = 0
    jx(3) = 0
    jx(4) = 4
    jx(5) = 4
    jx(6) = 0

    ix(1) = 25
    ix(2) = 1
    ix(3) = 1
    ix(4) = 13
    ix(5) = 13
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
    jy(2) = 8
    jy(3) = 0
    jy(4) = 4
    jy(5) = 0
    jy(6) = 4

    iy(1) = 1
    iy(2) = 25
    iy(3) = 1
    iy(4) = 13
    iy(5) = 1
    iy(6) = 13

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
    jz(3) = 8
    jz(4) = 0
    jz(5) = 4
    jz(6) = 4

    iz(1) = 1
    iz(2) = 1
    iz(3) = 25
    iz(4) = 1
    iz(5) = 13
    iz(6) = 13

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 33
    ijx(2) = 25
    ijx(3) = 25
    ijx(4) = 29
    ijx(5) = 29
    ijx(6) = 25
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
    ijx(19) = 21
    ijx(20) = 13
    ijx(21) = 13
    ijx(22) = 17
    ijx(23) = 17
    ijx(24) = 13
    ijx(25) = 21
    ijx(26) = 13
    ijx(27) = 13
    ijx(28) = 17
    ijx(29) = 17
    ijx(30) = 13
    ijx(31) = 9
    ijx(32) = 1
    ijx(33) = 1
    ijx(34) = 5
    ijx(35) = 5
    ijx(36) = 1

    ijy(1) = 1
    ijy(2) = 9
    ijy(3) = 1
    ijy(4) = 5
    ijy(5) = 1
    ijy(6) = 5
    ijy(7) = 25
    ijy(8) = 33
    ijy(9) = 25
    ijy(10) = 29
    ijy(11) = 25
    ijy(12) = 29
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
    ijy(31) = 13
    ijy(32) = 21
    ijy(33) = 13
    ijy(34) = 17
    ijy(35) = 13
    ijy(36) = 17

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
    ijz(13) = 25
    ijz(14) = 25
    ijz(15) = 33
    ijz(16) = 25
    ijz(17) = 29
    ijz(18) = 29
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
    ijz(31) = 13
    ijz(32) = 13
    ijz(33) = 21
    ijz(34) = 13
    ijz(35) = 17
    ijz(36) = 17

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

    allocate (n22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

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

    if ((nddbra*nsfket) .le. nchunksize_int64) nchunksize_int64 = nddbra*nsfket
    ntile = int(nddbra*nsfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nddbra*nsfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nddbra, xint22bra, n22bra, xint03ket, n03ket, dd_pair, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d03ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d22bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,iandj)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nddbra) + 1
              kl_tmp = (iquart - 1)/nddbra + 1

              test = xint22bra(ij_tmp)*xint03ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n22bra(ij_tmp)
                kl = n03ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_d_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
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

                                      ! do n = 2,   4

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

                                      ! i5 = in(n+1) =   29
                                      ! i3 =   13
                                      ! i4 =   25

                                      xin(29) = c10*xin(13) + xc00*xin(25)
                                      yin(29) = c10*yin(13) + yc00*yin(25)
                                      zin(29) = c10*zin(13) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   30
                                      ! i5 =   29
                                      ! i4 =   25

                                      xin(30) = xcp00*xin(29) + cp10*xin(25)
                                      yin(30) = ycp00*yin(29) + cp10*yin(25)
                                      zin(30) = zcp00*zin(29) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   33
                                      ! i3 =   25
                                      ! i4 =   29

                                      xin(33) = c10*xin(25) + xc00*xin(29)
                                      yin(33) = c10*yin(25) + yc00*yin(29)
                                      zin(33) = c10*zin(25) + zc00*zin(29)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   34
                                      ! i5 =   33
                                      ! i4 =   29

                                      xin(34) = xcp00*xin(33) + cp10*xin(29)
                                      yin(34) = ycp00*yin(33) + cp10*yin(29)
                                      zin(34) = zcp00*zin(33) + cp10*zin(29)

                                      ! ------------------

                                      ! i3 = i4 =   29
                                      ! i4 = i5 =   33

                                      ! n =    5

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

                                      ! i3 = i2+kn(n+1) =   15

                                      xin(15) = xc00*xin(3) + c01*xin(2)
                                      yin(15) = yc00*yin(3) + c01*yin(2)
                                      zin(15) = zc00*zin(3) + c01*zin(2)

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

                                      ! i3 = i2+kn(n+1) =   16

                                      xin(16) = xc00*xin(4) + c01*xin(3)
                                      yin(16) = yc00*yin(4) + c01*yin(3)
                                      zin(16) = zc00*zin(4) + c01*zin(3)

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
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   25

                                      xin(27) = c10*xin(3) + xc00*xin(15) + c01*xin(14)
                                      yin(27) = c10*yin(3) + yc00*yin(15) + c01*yin(14)
                                      zin(27) = c10*zin(3) + zc00*zin(15) + c01*zin(14)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   29

                                      xin(31) = c10*xin(15) + xc00*xin(27) + c01*xin(26)
                                      yin(31) = c10*yin(15) + yc00*yin(27) + c01*yin(26)
                                      zin(31) = c10*zin(15) + zc00*zin(27) + c01*zin(26)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   33

                                      xin(35) = c10*xin(27) + xc00*xin(31) + c01*xin(30)
                                      yin(35) = c10*yin(27) + yc00*yin(31) + c01*yin(30)
                                      zin(35) = c10*zin(27) + zc00*zin(31) + c01*zin(30)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   29
                                      ! i4 = i5 =   33

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   25

                                      xin(28) = c10*xin(4) + xc00*xin(16) + c01*xin(15)
                                      yin(28) = c10*yin(4) + yc00*yin(16) + c01*yin(15)
                                      zin(28) = c10*zin(4) + zc00*zin(16) + c01*zin(15)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   29

                                      xin(32) = c10*xin(16) + xc00*xin(28) + c01*xin(27)
                                      yin(32) = c10*yin(16) + yc00*yin(28) + c01*yin(27)
                                      zin(32) = c10*zin(16) + zc00*zin(28) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   33

                                      xin(36) = c10*xin(28) + xc00*xin(32) + c01*xin(31)
                                      yin(36) = c10*yin(28) + yc00*yin(32) + c01*yin(31)
                                      zin(36) = c10*zin(28) + zc00*zin(32) + c01*zin(31)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   29
                                      ! i4 = i5 =   33

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   33

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   33

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   29

                                      xin(33) = xin(33) + dxij*xin(29)
                                      yin(33) = yin(33) + dyij*yin(29)
                                      zin(33) = zin(33) + dzij*zin(29)

                                      ! i3 = i4 =   29
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   25

                                      xin(29) = xin(29) + dxij*xin(25)
                                      yin(29) = yin(29) + dyij*yin(25)
                                      zin(29) = zin(29) + dzij*zin(25)

                                      ! i3 = i4 =   25
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   33

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   29

                                      xin(33) = xin(33) + dxij*xin(29)
                                      yin(33) = yin(33) + dyij*yin(29)
                                      zin(33) = zin(33) + dzij*zin(29)

                                      ! i3 = i4 =   29
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   34

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   30

                                      xin(34) = xin(34) + dxij*xin(30)
                                      yin(34) = yin(34) + dyij*yin(30)
                                      zin(34) = zin(34) + dzij*zin(30)

                                      ! i3 = i4 =   30
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   26

                                      xin(30) = xin(30) + dxij*xin(26)
                                      yin(30) = yin(30) + dyij*yin(26)
                                      zin(30) = zin(30) + dzij*zin(26)

                                      ! i3 = i4 =   26
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   34

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   30

                                      xin(34) = xin(34) + dxij*xin(30)
                                      yin(34) = yin(34) + dyij*yin(30)
                                      zin(34) = zin(34) + dzij*zin(30)

                                      ! i3 = i4 =   30
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    2

                                      xin(6) = xin(14) + dxij*xin(2)
                                      yin(6) = yin(14) + dyij*yin(2)
                                      zin(6) = zin(14) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    2

                                      xin(18) = xin(26) + dxij*xin(14)
                                      yin(18) = yin(26) + dyij*yin(14)
                                      zin(18) = zin(26) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    2

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    2

                                      xin(10) = xin(18) + dxij*xin(6)
                                      yin(10) = yin(18) + dyij*yin(6)
                                      zin(10) = zin(18) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    2

                                      xin(22) = xin(30) + dxij*xin(18)
                                      yin(22) = yin(30) + dyij*yin(18)
                                      zin(22) = zin(30) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   14

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   35

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(35) = xin(35) + dxij*xin(31)
                                      yin(35) = yin(35) + dyij*yin(31)
                                      zin(35) = zin(35) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   27

                                      xin(31) = xin(31) + dxij*xin(27)
                                      yin(31) = yin(31) + dyij*yin(27)
                                      zin(31) = zin(31) + dzij*zin(27)

                                      ! i3 = i4 =   27
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   35

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   31

                                      xin(35) = xin(35) + dxij*xin(31)
                                      yin(35) = yin(35) + dyij*yin(31)
                                      zin(35) = zin(35) + dzij*zin(31)

                                      ! i3 = i4 =   31
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   32

                                      xin(36) = xin(36) + dxij*xin(32)
                                      yin(36) = yin(36) + dyij*yin(32)
                                      zin(36) = zin(36) + dzij*zin(32)

                                      ! i3 = i4 =   32
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   28

                                      xin(32) = xin(32) + dxij*xin(28)
                                      yin(32) = yin(32) + dyij*yin(28)
                                      zin(32) = zin(32) + dzij*zin(28)

                                      ! i3 = i4 =   28
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   36

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   32

                                      xin(36) = xin(36) + dxij*xin(32)
                                      yin(36) = yin(36) + dyij*yin(32)
                                      zin(36) = zin(36) + dzij*zin(32)

                                      ! i3 = i4 =   32
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    2

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

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   16

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =   49
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(49) = xc00
                                      yin(49) = yc00
                                      zin(49) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   38

                                      xin(38) = xcp00
                                      yin(38) = ycp00
                                      zin(38) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   50
                                      ! i2 =   49

                                      xin(50) = xcp00*xin(49) + cp10
                                      yin(50) = ycp00*yin(49) + cp10
                                      zin(50) = zcp00*zin(49) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   49

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   37
                                      ! i4 =   49

                                      xin(61) = c10*xin(37) + xc00*xin(49)
                                      yin(61) = c10*yin(37) + yc00*yin(49)
                                      zin(61) = c10*zin(37) + zc00*zin(49)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   62
                                      ! i5 =   61
                                      ! i4 =   49

                                      xin(62) = xcp00*xin(61) + cp10*xin(49)
                                      yin(62) = ycp00*yin(61) + cp10*yin(49)
                                      zin(62) = zcp00*zin(61) + cp10*zin(49)

                                      ! ------------------

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   61

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   65
                                      ! i3 =   49
                                      ! i4 =   61

                                      xin(65) = c10*xin(49) + xc00*xin(61)
                                      yin(65) = c10*yin(49) + yc00*yin(61)
                                      zin(65) = c10*zin(49) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   66
                                      ! i5 =   65
                                      ! i4 =   61

                                      xin(66) = xcp00*xin(65) + cp10*xin(61)
                                      yin(66) = ycp00*yin(65) + cp10*yin(61)
                                      zin(66) = zcp00*zin(65) + cp10*zin(61)

                                      ! ------------------

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   65

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   69
                                      ! i3 =   61
                                      ! i4 =   65

                                      xin(69) = c10*xin(61) + xc00*xin(65)
                                      yin(69) = c10*yin(61) + yc00*yin(65)
                                      zin(69) = c10*zin(61) + zc00*zin(65)

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

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   37
                                      ! i4 = i1+k2 =   38

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   39
                                      ! i3 =   37
                                      ! i4 =   38

                                      xin(39) = cp01*xin(37) + xcp00*xin(38)
                                      yin(39) = cp01*yin(37) + ycp00*yin(38)
                                      zin(39) = cp01*zin(37) + zcp00*zin(38)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   51

                                      xin(51) = xc00*xin(39) + c01*xin(38)
                                      yin(51) = yc00*yin(39) + c01*yin(38)
                                      zin(51) = zc00*zin(39) + c01*zin(38)

                                      ! ------------------

                                      ! i3 = i4 =   38
                                      ! i4 = i5 =   39

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   40
                                      ! i3 =   38
                                      ! i4 =   39

                                      xin(40) = cp01*xin(38) + xcp00*xin(39)
                                      yin(40) = cp01*yin(38) + ycp00*yin(39)
                                      zin(40) = cp01*zin(38) + zcp00*zin(39)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   52

                                      xin(52) = xc00*xin(40) + c01*xin(39)
                                      yin(52) = yc00*yin(40) + c01*yin(39)
                                      zin(52) = zc00*zin(40) + c01*zin(39)

                                      ! ------------------

                                      ! i3 = i4 =   39
                                      ! i4 = i5 =   40

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   49

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   61

                                      xin(63) = c10*xin(39) + xc00*xin(51) + c01*xin(50)
                                      yin(63) = c10*yin(39) + yc00*yin(51) + c01*yin(50)
                                      zin(63) = c10*zin(39) + zc00*zin(51) + c01*zin(50)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   61

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   65

                                      xin(67) = c10*xin(51) + xc00*xin(63) + c01*xin(62)
                                      yin(67) = c10*yin(51) + yc00*yin(63) + c01*yin(62)
                                      zin(67) = c10*zin(51) + zc00*zin(63) + c01*zin(62)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   65

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   69

                                      xin(71) = c10*xin(63) + xc00*xin(67) + c01*xin(66)
                                      yin(71) = c10*yin(63) + yc00*yin(67) + c01*yin(66)
                                      zin(71) = c10*zin(63) + zc00*zin(67) + c01*zin(66)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   69

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   37
                                      ! i4 = i2 =   49

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   61

                                      xin(64) = c10*xin(40) + xc00*xin(52) + c01*xin(51)
                                      yin(64) = c10*yin(40) + yc00*yin(52) + c01*yin(51)
                                      zin(64) = c10*zin(40) + zc00*zin(52) + c01*zin(51)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   61

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   65

                                      xin(68) = c10*xin(52) + xc00*xin(64) + c01*xin(63)
                                      yin(68) = c10*yin(52) + yc00*yin(64) + c01*yin(63)
                                      zin(68) = c10*zin(52) + zc00*zin(64) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   65

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   69

                                      xin(72) = c10*xin(64) + xc00*xin(68) + c01*xin(67)
                                      yin(72) = c10*yin(64) + yc00*yin(68) + c01*yin(67)
                                      zin(72) = c10*zin(64) + zc00*zin(68) + c01*zin(67)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   65
                                      ! i4 = i5 =   69

                                      ! nn =    5

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

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   65

                                      xin(69) = xin(69) + dxij*xin(65)
                                      yin(69) = yin(69) + dyij*yin(65)
                                      zin(69) = zin(69) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   61

                                      xin(65) = xin(65) + dxij*xin(61)
                                      yin(65) = yin(65) + dyij*yin(61)
                                      zin(65) = zin(65) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   65

                                      xin(69) = xin(69) + dxij*xin(65)
                                      yin(69) = yin(69) + dyij*yin(65)
                                      zin(69) = zin(69) + dzij*zin(65)

                                      ! i3 = i4 =   65
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   41

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   41

                                      ! do ni = 1,    2

                                      xin(41) = xin(49) + dxij*xin(37)
                                      yin(41) = yin(49) + dyij*yin(37)
                                      zin(41) = zin(49) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    2

                                      xin(53) = xin(61) + dxij*xin(49)
                                      yin(53) = yin(61) + dyij*yin(49)
                                      zin(53) = zin(61) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   45

                                      ! nj =    2

                                      ! i4 = i3 =   45

                                      ! do ni = 1,    2

                                      xin(45) = xin(53) + dxij*xin(41)
                                      yin(45) = yin(53) + dyij*yin(41)
                                      zin(45) = zin(53) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    2

                                      xin(57) = xin(65) + dxij*xin(53)
                                      yin(57) = yin(65) + dyij*yin(53)
                                      zin(57) = zin(65) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   49

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   70

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   66

                                      xin(70) = xin(70) + dxij*xin(66)
                                      yin(70) = yin(70) + dyij*yin(66)
                                      zin(70) = zin(70) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   62

                                      xin(66) = xin(66) + dxij*xin(62)
                                      yin(66) = yin(66) + dyij*yin(62)
                                      zin(66) = zin(66) + dzij*zin(62)

                                      ! i3 = i4 =   62
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   70

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   66

                                      xin(70) = xin(70) + dxij*xin(66)
                                      yin(70) = yin(70) + dyij*yin(66)
                                      zin(70) = zin(70) + dzij*zin(66)

                                      ! i3 = i4 =   66
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   42

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   42

                                      ! do ni = 1,    2

                                      xin(42) = xin(50) + dxij*xin(38)
                                      yin(42) = yin(50) + dyij*yin(38)
                                      zin(42) = zin(50) + dzij*zin(38)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    2

                                      xin(54) = xin(62) + dxij*xin(50)
                                      yin(54) = yin(62) + dyij*yin(50)
                                      zin(54) = zin(62) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   46

                                      ! nj =    2

                                      ! i4 = i3 =   46

                                      ! do ni = 1,    2

                                      xin(46) = xin(54) + dxij*xin(42)
                                      yin(46) = yin(54) + dyij*yin(42)
                                      zin(46) = zin(54) + dzij*zin(42)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                      ! ni =    2

                                      xin(58) = xin(66) + dxij*xin(54)
                                      yin(58) = yin(66) + dyij*yin(54)
                                      zin(58) = zin(66) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   50

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   67

                                      xin(71) = xin(71) + dxij*xin(67)
                                      yin(71) = yin(71) + dyij*yin(67)
                                      zin(71) = zin(71) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   63

                                      xin(67) = xin(67) + dxij*xin(63)
                                      yin(67) = yin(67) + dyij*yin(63)
                                      zin(67) = zin(67) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   67

                                      xin(71) = xin(71) + dxij*xin(67)
                                      yin(71) = yin(71) + dyij*yin(67)
                                      zin(71) = zin(71) + dzij*zin(67)

                                      ! i3 = i4 =   67
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   43

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   43

                                      ! do ni = 1,    2

                                      xin(43) = xin(51) + dxij*xin(39)
                                      yin(43) = yin(51) + dyij*yin(39)
                                      zin(43) = zin(51) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    2

                                      xin(55) = xin(63) + dxij*xin(51)
                                      yin(55) = yin(63) + dyij*yin(51)
                                      zin(55) = zin(63) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   47

                                      ! nj =    2

                                      ! i4 = i3 =   47

                                      ! do ni = 1,    2

                                      xin(47) = xin(55) + dxij*xin(43)
                                      yin(47) = yin(55) + dyij*yin(43)
                                      zin(47) = zin(55) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    2

                                      xin(59) = xin(67) + dxij*xin(55)
                                      yin(59) = yin(67) + dyij*yin(55)
                                      zin(59) = zin(67) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   51

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   68

                                      xin(72) = xin(72) + dxij*xin(68)
                                      yin(72) = yin(72) + dyij*yin(68)
                                      zin(72) = zin(72) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   64

                                      xin(68) = xin(68) + dxij*xin(64)
                                      yin(68) = yin(68) + dyij*yin(64)
                                      zin(68) = zin(68) + dzij*zin(64)

                                      ! i3 = i4 =   64
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   68

                                      xin(72) = xin(72) + dxij*xin(68)
                                      yin(72) = yin(72) + dyij*yin(68)
                                      zin(72) = zin(72) + dzij*zin(68)

                                      ! i3 = i4 =   68
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   44

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   44

                                      ! do ni = 1,    2

                                      xin(44) = xin(52) + dxij*xin(40)
                                      yin(44) = yin(52) + dyij*yin(40)
                                      zin(44) = zin(52) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   56

                                      ! ni =    2

                                      xin(56) = xin(64) + dxij*xin(52)
                                      yin(56) = yin(64) + dyij*yin(52)
                                      zin(56) = zin(64) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   68

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   48

                                      ! nj =    2

                                      ! i4 = i3 =   48

                                      ! do ni = 1,    2

                                      xin(48) = xin(56) + dxij*xin(44)
                                      yin(48) = yin(56) + dyij*yin(44)
                                      zin(48) = zin(56) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                      ! ni =    2

                                      xin(60) = xin(68) + dxij*xin(56)
                                      yin(60) = yin(68) + dyij*yin(56)
                                      zin(60) = zin(68) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   52

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =   85
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(85) = xc00
                                      yin(85) = yc00
                                      zin(85) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   74

                                      xin(74) = xcp00
                                      yin(74) = ycp00
                                      zin(74) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   86
                                      ! i2 =   85

                                      xin(86) = xcp00*xin(85) + cp10
                                      yin(86) = ycp00*yin(85) + cp10
                                      zin(86) = zcp00*zin(85) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   85

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   97
                                      ! i3 =   73
                                      ! i4 =   85

                                      xin(97) = c10*xin(73) + xc00*xin(85)
                                      yin(97) = c10*yin(73) + yc00*yin(85)
                                      zin(97) = c10*zin(73) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   98
                                      ! i5 =   97
                                      ! i4 =   85

                                      xin(98) = xcp00*xin(97) + cp10*xin(85)
                                      yin(98) = ycp00*yin(97) + cp10*yin(85)
                                      zin(98) = zcp00*zin(97) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   97

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  101
                                      ! i3 =   85
                                      ! i4 =   97

                                      xin(101) = c10*xin(85) + xc00*xin(97)
                                      yin(101) = c10*yin(85) + yc00*yin(97)
                                      zin(101) = c10*zin(85) + zc00*zin(97)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  102
                                      ! i5 =  101
                                      ! i4 =   97

                                      xin(102) = xcp00*xin(101) + cp10*xin(97)
                                      yin(102) = ycp00*yin(101) + cp10*yin(97)
                                      zin(102) = zcp00*zin(101) + cp10*zin(97)

                                      ! ------------------

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  101

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  105
                                      ! i3 =   97
                                      ! i4 =  101

                                      xin(105) = c10*xin(97) + xc00*xin(101)
                                      yin(105) = c10*yin(97) + yc00*yin(101)
                                      zin(105) = c10*zin(97) + zc00*zin(101)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  106
                                      ! i5 =  105
                                      ! i4 =  101

                                      xin(106) = xcp00*xin(105) + cp10*xin(101)
                                      yin(106) = ycp00*yin(105) + cp10*yin(101)
                                      zin(106) = zcp00*zin(105) + cp10*zin(101)

                                      ! ------------------

                                      ! i3 = i4 =  101
                                      ! i4 = i5 =  105

                                      ! n =    5

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

                                      ! i3 = i2+kn(n+1) =   87

                                      xin(87) = xc00*xin(75) + c01*xin(74)
                                      yin(87) = yc00*yin(75) + c01*yin(74)
                                      zin(87) = zc00*zin(75) + c01*zin(74)

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

                                      ! i3 = i2+kn(n+1) =   88

                                      xin(88) = xc00*xin(76) + c01*xin(75)
                                      yin(88) = yc00*yin(76) + c01*yin(75)
                                      zin(88) = zc00*zin(76) + c01*zin(75)

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
                                      ! i4 = i2 =   85

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   97

                                      xin(99) = c10*xin(75) + xc00*xin(87) + c01*xin(86)
                                      yin(99) = c10*yin(75) + yc00*yin(87) + c01*yin(86)
                                      zin(99) = c10*zin(75) + zc00*zin(87) + c01*zin(86)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   97

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  101

                                      xin(103) = c10*xin(87) + xc00*xin(99) + c01*xin(98)
                                      yin(103) = c10*yin(87) + yc00*yin(99) + c01*yin(98)
                                      zin(103) = c10*zin(87) + zc00*zin(99) + c01*zin(98)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  101

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  105

                                      xin(107) = c10*xin(99) + xc00*xin(103) + c01*xin(102)
                                      yin(107) = c10*yin(99) + yc00*yin(103) + c01*yin(102)
                                      zin(107) = c10*zin(99) + zc00*zin(103) + c01*zin(102)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  101
                                      ! i4 = i5 =  105

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   85

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   97

                                      xin(100) = c10*xin(76) + xc00*xin(88) + c01*xin(87)
                                      yin(100) = c10*yin(76) + yc00*yin(88) + c01*yin(87)
                                      zin(100) = c10*zin(76) + zc00*zin(88) + c01*zin(87)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   97

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  101

                                      xin(104) = c10*xin(88) + xc00*xin(100) + c01*xin(99)
                                      yin(104) = c10*yin(88) + yc00*yin(100) + c01*yin(99)
                                      zin(104) = c10*zin(88) + zc00*zin(100) + c01*zin(99)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  101

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  105

                                      xin(108) = c10*xin(100) + xc00*xin(104) + c01*xin(103)
                                      yin(108) = c10*yin(100) + yc00*yin(104) + c01*yin(103)
                                      zin(108) = c10*zin(100) + zc00*zin(104) + c01*zin(103)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  101
                                      ! i4 = i5 =  105

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  105

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  105

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  101

                                      xin(105) = xin(105) + dxij*xin(101)
                                      yin(105) = yin(105) + dyij*yin(101)
                                      zin(105) = zin(105) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   97

                                      xin(101) = xin(101) + dxij*xin(97)
                                      yin(101) = yin(101) + dyij*yin(97)
                                      zin(101) = zin(101) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  105

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  101

                                      xin(105) = xin(105) + dxij*xin(101)
                                      yin(105) = yin(105) + dyij*yin(101)
                                      zin(105) = zin(105) + dzij*zin(101)

                                      ! i3 = i4 =  101
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   77

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   77

                                      ! do ni = 1,    2

                                      xin(77) = xin(85) + dxij*xin(73)
                                      yin(77) = yin(85) + dyij*yin(73)
                                      zin(77) = zin(85) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    2

                                      xin(89) = xin(97) + dxij*xin(85)
                                      yin(89) = yin(97) + dyij*yin(85)
                                      zin(89) = zin(97) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  101

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   81

                                      ! nj =    2

                                      ! i4 = i3 =   81

                                      ! do ni = 1,    2

                                      xin(81) = xin(89) + dxij*xin(77)
                                      yin(81) = yin(89) + dyij*yin(77)
                                      zin(81) = zin(89) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    2

                                      xin(93) = xin(101) + dxij*xin(89)
                                      yin(93) = yin(101) + dyij*yin(89)
                                      zin(93) = zin(101) + dzij*zin(89)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   85

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  102

                                      xin(106) = xin(106) + dxij*xin(102)
                                      yin(106) = yin(106) + dyij*yin(102)
                                      zin(106) = zin(106) + dzij*zin(102)

                                      ! i3 = i4 =  102
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   98

                                      xin(102) = xin(102) + dxij*xin(98)
                                      yin(102) = yin(102) + dyij*yin(98)
                                      zin(102) = zin(102) + dzij*zin(98)

                                      ! i3 = i4 =   98
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  102

                                      xin(106) = xin(106) + dxij*xin(102)
                                      yin(106) = yin(106) + dyij*yin(102)
                                      zin(106) = zin(106) + dzij*zin(102)

                                      ! i3 = i4 =  102
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   78

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   78

                                      ! do ni = 1,    2

                                      xin(78) = xin(86) + dxij*xin(74)
                                      yin(78) = yin(86) + dyij*yin(74)
                                      zin(78) = zin(86) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    2

                                      xin(90) = xin(98) + dxij*xin(86)
                                      yin(90) = yin(98) + dyij*yin(86)
                                      zin(90) = zin(98) + dzij*zin(86)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  102

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   82

                                      ! nj =    2

                                      ! i4 = i3 =   82

                                      ! do ni = 1,    2

                                      xin(82) = xin(90) + dxij*xin(78)
                                      yin(82) = yin(90) + dyij*yin(78)
                                      zin(82) = zin(90) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    2

                                      xin(94) = xin(102) + dxij*xin(90)
                                      yin(94) = yin(102) + dyij*yin(90)
                                      zin(94) = zin(102) + dzij*zin(90)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   86

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  103

                                      xin(107) = xin(107) + dxij*xin(103)
                                      yin(107) = yin(107) + dyij*yin(103)
                                      zin(107) = zin(107) + dzij*zin(103)

                                      ! i3 = i4 =  103
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   99

                                      xin(103) = xin(103) + dxij*xin(99)
                                      yin(103) = yin(103) + dyij*yin(99)
                                      zin(103) = zin(103) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  103

                                      xin(107) = xin(107) + dxij*xin(103)
                                      yin(107) = yin(107) + dyij*yin(103)
                                      zin(107) = zin(107) + dzij*zin(103)

                                      ! i3 = i4 =  103
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   79

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   79

                                      ! do ni = 1,    2

                                      xin(79) = xin(87) + dxij*xin(75)
                                      yin(79) = yin(87) + dyij*yin(75)
                                      zin(79) = zin(87) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    2

                                      xin(91) = xin(99) + dxij*xin(87)
                                      yin(91) = yin(99) + dyij*yin(87)
                                      zin(91) = zin(99) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   83

                                      ! nj =    2

                                      ! i4 = i3 =   83

                                      ! do ni = 1,    2

                                      xin(83) = xin(91) + dxij*xin(79)
                                      yin(83) = yin(91) + dyij*yin(79)
                                      zin(83) = zin(91) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    2

                                      xin(95) = xin(103) + dxij*xin(91)
                                      yin(95) = yin(103) + dyij*yin(91)
                                      zin(95) = zin(103) + dzij*zin(91)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   87

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  104

                                      xin(108) = xin(108) + dxij*xin(104)
                                      yin(108) = yin(108) + dyij*yin(104)
                                      zin(108) = zin(108) + dzij*zin(104)

                                      ! i3 = i4 =  104
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  100

                                      xin(104) = xin(104) + dxij*xin(100)
                                      yin(104) = yin(104) + dyij*yin(100)
                                      zin(104) = zin(104) + dzij*zin(100)

                                      ! i3 = i4 =  100
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  104

                                      xin(108) = xin(108) + dxij*xin(104)
                                      yin(108) = yin(108) + dyij*yin(104)
                                      zin(108) = zin(108) + dzij*zin(104)

                                      ! i3 = i4 =  104
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   80

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   80

                                      ! do ni = 1,    2

                                      xin(80) = xin(88) + dxij*xin(76)
                                      yin(80) = yin(88) + dyij*yin(76)
                                      zin(80) = zin(88) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   92

                                      ! ni =    2

                                      xin(92) = xin(100) + dxij*xin(88)
                                      yin(92) = yin(100) + dyij*yin(88)
                                      zin(92) = zin(100) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  104

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   84

                                      ! nj =    2

                                      ! i4 = i3 =   84

                                      ! do ni = 1,    2

                                      xin(84) = xin(92) + dxij*xin(80)
                                      yin(84) = yin(92) + dyij*yin(80)
                                      zin(84) = zin(92) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    2

                                      xin(96) = xin(104) + dxij*xin(92)
                                      yin(96) = yin(104) + dyij*yin(92)
                                      zin(96) = zin(104) + dzij*zin(92)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   88

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

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

                                      ! i2 = in(2) =  121
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(121) = xc00
                                      yin(121) = yc00
                                      zin(121) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  110

                                      xin(110) = xcp00
                                      yin(110) = ycp00
                                      zin(110) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  122
                                      ! i2 =  121

                                      xin(122) = xcp00*xin(121) + cp10
                                      yin(122) = ycp00*yin(121) + cp10
                                      zin(122) = zcp00*zin(121) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  121

                                      ! do n = 2,   4

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

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  137
                                      ! i3 =  121
                                      ! i4 =  133

                                      xin(137) = c10*xin(121) + xc00*xin(133)
                                      yin(137) = c10*yin(121) + yc00*yin(133)
                                      zin(137) = c10*zin(121) + zc00*zin(133)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  138
                                      ! i5 =  137
                                      ! i4 =  133

                                      xin(138) = xcp00*xin(137) + cp10*xin(133)
                                      yin(138) = ycp00*yin(137) + cp10*yin(133)
                                      zin(138) = zcp00*zin(137) + cp10*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  141
                                      ! i3 =  133
                                      ! i4 =  137

                                      xin(141) = c10*xin(133) + xc00*xin(137)
                                      yin(141) = c10*yin(133) + yc00*yin(137)
                                      zin(141) = c10*zin(133) + zc00*zin(137)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  142
                                      ! i5 =  141
                                      ! i4 =  137

                                      xin(142) = xcp00*xin(141) + cp10*xin(137)
                                      yin(142) = ycp00*yin(141) + cp10*yin(137)
                                      zin(142) = zcp00*zin(141) + cp10*zin(137)

                                      ! ------------------

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  109
                                      ! i4 = i1+k2 =  110

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  111
                                      ! i3 =  109
                                      ! i4 =  110

                                      xin(111) = cp01*xin(109) + xcp00*xin(110)
                                      yin(111) = cp01*yin(109) + ycp00*yin(110)
                                      zin(111) = cp01*zin(109) + zcp00*zin(110)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  123

                                      xin(123) = xc00*xin(111) + c01*xin(110)
                                      yin(123) = yc00*yin(111) + c01*yin(110)
                                      zin(123) = zc00*zin(111) + c01*zin(110)

                                      ! ------------------

                                      ! i3 = i4 =  110
                                      ! i4 = i5 =  111

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  112
                                      ! i3 =  110
                                      ! i4 =  111

                                      xin(112) = cp01*xin(110) + xcp00*xin(111)
                                      yin(112) = cp01*yin(110) + ycp00*yin(111)
                                      zin(112) = cp01*zin(110) + zcp00*zin(111)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  124

                                      xin(124) = xc00*xin(112) + c01*xin(111)
                                      yin(124) = yc00*yin(112) + c01*yin(111)
                                      zin(124) = zc00*zin(112) + c01*zin(111)

                                      ! ------------------

                                      ! i3 = i4 =  111
                                      ! i4 = i5 =  112

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  121

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  133

                                      xin(135) = c10*xin(111) + xc00*xin(123) + c01*xin(122)
                                      yin(135) = c10*yin(111) + yc00*yin(123) + c01*yin(122)
                                      zin(135) = c10*zin(111) + zc00*zin(123) + c01*zin(122)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  137

                                      xin(139) = c10*xin(123) + xc00*xin(135) + c01*xin(134)
                                      yin(139) = c10*yin(123) + yc00*yin(135) + c01*yin(134)
                                      zin(139) = c10*zin(123) + zc00*zin(135) + c01*zin(134)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  141

                                      xin(143) = c10*xin(135) + xc00*xin(139) + c01*xin(138)
                                      yin(143) = c10*yin(135) + yc00*yin(139) + c01*yin(138)
                                      zin(143) = c10*zin(135) + zc00*zin(139) + c01*zin(138)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  121

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  133

                                      xin(136) = c10*xin(112) + xc00*xin(124) + c01*xin(123)
                                      yin(136) = c10*yin(112) + yc00*yin(124) + c01*yin(123)
                                      zin(136) = c10*zin(112) + zc00*zin(124) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  137

                                      xin(140) = c10*xin(124) + xc00*xin(136) + c01*xin(135)
                                      yin(140) = c10*yin(124) + yc00*yin(136) + c01*yin(135)
                                      zin(140) = c10*zin(124) + zc00*zin(136) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  141

                                      xin(144) = c10*xin(136) + xc00*xin(140) + c01*xin(139)
                                      yin(144) = c10*yin(136) + yc00*yin(140) + c01*yin(139)
                                      zin(144) = c10*zin(136) + zc00*zin(140) + c01*zin(139)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  141

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  133

                                      xin(137) = xin(137) + dxij*xin(133)
                                      yin(137) = yin(137) + dyij*yin(133)
                                      zin(137) = zin(137) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  113

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  113

                                      ! do ni = 1,    2

                                      xin(113) = xin(121) + dxij*xin(109)
                                      yin(113) = yin(121) + dyij*yin(109)
                                      zin(113) = zin(121) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    2

                                      xin(125) = xin(133) + dxij*xin(121)
                                      yin(125) = yin(133) + dyij*yin(121)
                                      zin(125) = zin(133) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  117

                                      ! nj =    2

                                      ! i4 = i3 =  117

                                      ! do ni = 1,    2

                                      xin(117) = xin(125) + dxij*xin(113)
                                      yin(117) = yin(125) + dyij*yin(113)
                                      zin(117) = zin(125) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    2

                                      xin(129) = xin(137) + dxij*xin(125)
                                      yin(129) = yin(137) + dyij*yin(125)
                                      zin(129) = zin(137) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  121

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(142) = xin(142) + dxij*xin(138)
                                      yin(142) = yin(142) + dyij*yin(138)
                                      zin(142) = zin(142) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  134

                                      xin(138) = xin(138) + dxij*xin(134)
                                      yin(138) = yin(138) + dyij*yin(134)
                                      zin(138) = zin(138) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(142) = xin(142) + dxij*xin(138)
                                      yin(142) = yin(142) + dyij*yin(138)
                                      zin(142) = zin(142) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  114

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  114

                                      ! do ni = 1,    2

                                      xin(114) = xin(122) + dxij*xin(110)
                                      yin(114) = yin(122) + dyij*yin(110)
                                      zin(114) = zin(122) + dzij*zin(110)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    2

                                      xin(126) = xin(134) + dxij*xin(122)
                                      yin(126) = yin(134) + dyij*yin(122)
                                      zin(126) = zin(134) + dzij*zin(122)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  118

                                      ! nj =    2

                                      ! i4 = i3 =  118

                                      ! do ni = 1,    2

                                      xin(118) = xin(126) + dxij*xin(114)
                                      yin(118) = yin(126) + dyij*yin(114)
                                      zin(118) = zin(126) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                      ! ni =    2

                                      xin(130) = xin(138) + dxij*xin(126)
                                      yin(130) = yin(138) + dyij*yin(126)
                                      zin(130) = zin(138) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  122

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  135

                                      xin(139) = xin(139) + dxij*xin(135)
                                      yin(139) = yin(139) + dyij*yin(135)
                                      zin(139) = zin(139) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  115

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  115

                                      ! do ni = 1,    2

                                      xin(115) = xin(123) + dxij*xin(111)
                                      yin(115) = yin(123) + dyij*yin(111)
                                      zin(115) = zin(123) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    2

                                      xin(127) = xin(135) + dxij*xin(123)
                                      yin(127) = yin(135) + dyij*yin(123)
                                      zin(127) = zin(135) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  119

                                      ! nj =    2

                                      ! i4 = i3 =  119

                                      ! do ni = 1,    2

                                      xin(119) = xin(127) + dxij*xin(115)
                                      yin(119) = yin(127) + dyij*yin(115)
                                      zin(119) = zin(127) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    2

                                      xin(131) = xin(139) + dxij*xin(127)
                                      yin(131) = yin(139) + dyij*yin(127)
                                      zin(131) = zin(139) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  123

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  136

                                      xin(140) = xin(140) + dxij*xin(136)
                                      yin(140) = yin(140) + dyij*yin(136)
                                      zin(140) = zin(140) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  116

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  116

                                      ! do ni = 1,    2

                                      xin(116) = xin(124) + dxij*xin(112)
                                      yin(116) = yin(124) + dyij*yin(112)
                                      zin(116) = zin(124) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    2

                                      xin(128) = xin(136) + dxij*xin(124)
                                      yin(128) = yin(136) + dyij*yin(124)
                                      zin(128) = zin(136) + dzij*zin(124)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  140

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  120

                                      ! nj =    2

                                      ! i4 = i3 =  120

                                      ! do ni = 1,    2

                                      xin(120) = xin(128) + dxij*xin(116)
                                      yin(120) = yin(128) + dyij*yin(116)
                                      zin(120) = zin(128) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    2

                                      xin(132) = xin(140) + dxij*xin(128)
                                      yin(132) = yin(140) + dyij*yin(128)
                                      zin(132) = zin(140) + dzij*zin(128)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  124

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 360! loop over all integrals

                                        l = n - 10*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d22bra(j)*d03ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 36)*yin(my + 36)*zin(mz + 36) & ! root  2
                                                        + xin(mx + 72)*yin(my + 72)*zin(mz + 72) & ! root  3
                                                        + xin(mx + 108)*yin(my + 108)*zin(mz + 108)) ! root  4

                                        j = int(n/10) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*60 ! Stride between functions in i

                                    do j = 1, maxj2

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

                              deallocate (n22bra)
                              deallocate (xint22bra)
                              deallocate (n03ket)
                              deallocate (xint03ket)

                              end subroutine int2230
                              end submodule
