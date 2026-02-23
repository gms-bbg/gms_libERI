! The total angular momentum of this class is:           5
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3020_impl
contains
  module subroutine int3020(sf_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sf_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n03bra(:), n02ket(:)
    real(dp), allocatable :: xint03bra(:), xint02ket(:)
    integer(kind=int64) :: nsfbra, nsdket
    real(dp) :: scutsfbra, scutsdket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2
    integer(kind=int64) :: n, i1, i3, i4, i5, k3, k4, nn
    real(dp) :: cp10, c10, cp01, c01
    integer(kind=int64) :: nx, ny, nz, mx, my, mz
    integer(kind=int64) :: bra_loop, ket_loop, ijtop, kltop
    real(dp) :: t_expon_ab, t_expon_cd, t_inverse_expon_ab, t_inverse_expon_cd
    real(dp) :: t_expon_abcd_inverse, rho, expe, dum, rab, rcd
    real(dp) :: brrk, akxk, akyk, akzk, t_expon_c, t_expon_d, t_expon_a, t_expon_b
    real(dp) :: xa, ya, za, axak, ayak, azak, axai, ayai, azai, bbrrk, xb, yb, zb, bxbk
    real(dp) :: bybk, bzbk, bxbi, bybi, bzbi, xx, c1x, c2x, c3x, c4x, c1y, c2y, c3y, c4y
    real(dp) :: c1z, c2z, c3z, c4z, f00, u2, duminv, dm2inv, bp01, b00, b10, xcp00, xc00
    real(dp) :: ycp00, zcp00, zc00, yc00, dij
    real(dp) :: buff(9)
    real(dp) :: roots(3), wghts(3)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(3), wts(3), alpha(3), beta(3), wrk(3)
    real(dp) :: xin(36), yin(36), zin(36)
    real(dp) :: eri_value(60)
    real(dp) :: d03bra(10), d02ket(6)
    integer(kind=int64) :: ix(10), jx(1), kx(6), lx(1)
    integer(kind=int64) :: iy(10), jy(1), ky(6), ly(1)
    integer(kind=int64) :: iz(10), jz(1), kz(6), lz(1)
    integer(kind=int64) :: in(4), in1(4), kn(3)
    integer(kind=int64) :: ijx(10), ijy(10), ijz(10)
    integer(kind=int64) :: klx(6), kly(6), klz(6)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 4
    in1(3) = 7
    in1(4) = 10

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

    jx(1) = 0

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
    ky(2) = 2
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 1

    jy(1) = 0

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
    kz(2) = 0
    kz(3) = 2
    kz(4) = 0
    kz(5) = 1
    kz(6) = 1

    jz(1) = 0

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

    ijx(1) = 10
    ijx(2) = 1
    ijx(3) = 1
    ijx(4) = 7
    ijx(5) = 7
    ijx(6) = 4
    ijx(7) = 1
    ijx(8) = 4
    ijx(9) = 1
    ijx(10) = 4

    ijy(1) = 1
    ijy(2) = 10
    ijy(3) = 1
    ijy(4) = 4
    ijy(5) = 1
    ijy(6) = 7
    ijy(7) = 7
    ijy(8) = 1
    ijy(9) = 4
    ijy(10) = 4

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 10
    ijz(4) = 1
    ijz(5) = 4
    ijz(6) = 1
    ijz(7) = 4
    ijz(8) = 7
    ijz(9) = 7
    ijz(10) = 4

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

    allocate (n03bra(res%n_s_shl*res%n_f_shl))
    allocate (xint03bra(res%n_s_shl*res%n_f_shl))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))

    ! Start screening

    scutsfbra = cutoff_schwarz/maxval(sf_pair%xints)
    nsfbra = 0
    do ij = 1, res%n_s_shl*res%n_f_shl
      if (sf_pair%xints(ij) .ge. scutsfbra) then
        nsfbra = nsfbra + 1
        xint03bra(nsfbra) = sf_pair%xints(ij)
        n03bra(nsfbra) = ij
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

    if ((nsfbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = nsfbra*nsdket
    ntile = int(nsfbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nsfbra*nsdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nsfbra, xint03bra, n03bra, xint02ket, n02ket, sd_pair, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d02ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d03bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nsfbra) + 1
              kl_tmp = (iquart - 1)/nsfbra + 1

              test = xint03bra(ij_tmp)*xint02ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n03bra(ij_tmp)
                kl = n02ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_s_shl(jsh_tmp)
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

                    t_expon_ab = sf_pair%t_expon_ab(sf_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = sf_pair%expon_b(sf_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = sf_pair%expon_a(sf_pair%pair_loc(ij) + bra_loop)
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

                    d03bra(1) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)
                    d03bra(2) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)
                    d03bra(3) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)
                    d03bra(4) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(5) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(6) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(7) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(8) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(9) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d03bra(10) = sf_pair%d_coeff_alt(sf_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(4) = xc00
                                      yin(4) = yc00
                                      zin(4) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =    5
                                      ! i2 =    4

                                      xin(5) = xcp00*xin(4) + cp10
                                      yin(5) = ycp00*yin(4) + cp10
                                      zin(5) = zcp00*zin(4) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    4

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    7
                                      ! i3 =    1
                                      ! i4 =    4

                                      xin(7) = c10*xin(1) + xc00*xin(4)
                                      yin(7) = c10*yin(1) + yc00*yin(4)
                                      zin(7) = c10*zin(1) + zc00*zin(4)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =    8
                                      ! i5 =    7
                                      ! i4 =    4

                                      xin(8) = xcp00*xin(7) + cp10*xin(4)
                                      yin(8) = ycp00*yin(7) + cp10*yin(4)
                                      zin(8) = zcp00*zin(7) + cp10*zin(4)

                                      ! ------------------

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

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   11
                                      ! i5 =   10
                                      ! i4 =    7

                                      xin(11) = xcp00*xin(10) + cp10*xin(7)
                                      yin(11) = ycp00*yin(10) + cp10*yin(7)
                                      zin(11) = zcp00*zin(10) + cp10*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =   10

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

                                      ! i3 = i2+kn(n+1) =    6

                                      xin(6) = xc00*xin(3) + c01*xin(2)
                                      yin(6) = yc00*yin(3) + c01*yin(2)
                                      zin(6) = zc00*zin(3) + c01*zin(2)

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
                                      ! i4 = i2 =    4

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =    7

                                      xin(9) = c10*xin(3) + xc00*xin(6) + c01*xin(5)
                                      yin(9) = c10*yin(3) + yc00*yin(6) + c01*yin(5)
                                      zin(9) = c10*zin(3) + zc00*zin(6) + c01*zin(5)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    4
                                      ! i4 = i5 =    7

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   10

                                      xin(12) = c10*xin(6) + xc00*xin(9) + c01*xin(8)
                                      yin(12) = c10*yin(6) + yc00*yin(9) + c01*yin(8)
                                      zin(12) = c10*zin(6) + zc00*zin(9) + c01*zin(8)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =   10

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(16) = xc00
                                      yin(16) = yc00
                                      zin(16) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   14

                                      xin(14) = xcp00
                                      yin(14) = ycp00
                                      zin(14) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   17
                                      ! i2 =   16

                                      xin(17) = xcp00*xin(16) + cp10
                                      yin(17) = ycp00*yin(16) + cp10
                                      zin(17) = zcp00*zin(16) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   13
                                      ! i4 = i2 =   16

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   19
                                      ! i3 =   13
                                      ! i4 =   16

                                      xin(19) = c10*xin(13) + xc00*xin(16)
                                      yin(19) = c10*yin(13) + yc00*yin(16)
                                      zin(19) = c10*zin(13) + zc00*zin(16)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   20
                                      ! i5 =   19
                                      ! i4 =   16

                                      xin(20) = xcp00*xin(19) + cp10*xin(16)
                                      yin(20) = ycp00*yin(19) + cp10*yin(16)
                                      zin(20) = zcp00*zin(19) + cp10*zin(16)

                                      ! ------------------

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

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   23
                                      ! i5 =   22
                                      ! i4 =   19

                                      xin(23) = xcp00*xin(22) + cp10*xin(19)
                                      yin(23) = ycp00*yin(22) + cp10*yin(19)
                                      zin(23) = zcp00*zin(22) + cp10*zin(19)

                                      ! ------------------

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   22

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   13
                                      ! i4 = i1+k2 =   14

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   15
                                      ! i3 =   13
                                      ! i4 =   14

                                      xin(15) = cp01*xin(13) + xcp00*xin(14)
                                      yin(15) = cp01*yin(13) + ycp00*yin(14)
                                      zin(15) = cp01*zin(13) + zcp00*zin(14)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   18

                                      xin(18) = xc00*xin(15) + c01*xin(14)
                                      yin(18) = yc00*yin(15) + c01*yin(14)
                                      zin(18) = zc00*zin(15) + c01*zin(14)

                                      ! ------------------

                                      ! i3 = i4 =   14
                                      ! i4 = i5 =   15

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   13
                                      ! i4 = i2 =   16

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   19

                                      xin(21) = c10*xin(15) + xc00*xin(18) + c01*xin(17)
                                      yin(21) = c10*yin(15) + yc00*yin(18) + c01*yin(17)
                                      zin(21) = c10*zin(15) + zc00*zin(18) + c01*zin(17)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   16
                                      ! i4 = i5 =   19

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   22

                                      xin(24) = c10*xin(18) + xc00*xin(21) + c01*xin(20)
                                      yin(24) = c10*yin(18) + yc00*yin(21) + c01*yin(20)
                                      zin(24) = c10*zin(18) + zc00*zin(21) + c01*zin(20)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   22

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

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
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(28) = xc00
                                      yin(28) = yc00
                                      zin(28) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   26

                                      xin(26) = xcp00
                                      yin(26) = ycp00
                                      zin(26) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   29
                                      ! i2 =   28

                                      xin(29) = xcp00*xin(28) + cp10
                                      yin(29) = ycp00*yin(28) + cp10
                                      zin(29) = zcp00*zin(28) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   28

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   31
                                      ! i3 =   25
                                      ! i4 =   28

                                      xin(31) = c10*xin(25) + xc00*xin(28)
                                      yin(31) = c10*yin(25) + yc00*yin(28)
                                      zin(31) = c10*zin(25) + zc00*zin(28)

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

                                      ! n =    3

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

                                      ! n =    4

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   25
                                      ! i4 = i1+k2 =   26

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   27
                                      ! i3 =   25
                                      ! i4 =   26

                                      xin(27) = cp01*xin(25) + xcp00*xin(26)
                                      yin(27) = cp01*yin(25) + ycp00*yin(26)
                                      zin(27) = cp01*zin(25) + zcp00*zin(26)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   30

                                      xin(30) = xc00*xin(27) + c01*xin(26)
                                      yin(30) = yc00*yin(27) + c01*yin(26)
                                      zin(30) = zc00*zin(27) + c01*zin(26)

                                      ! ------------------

                                      ! i3 = i4 =   26
                                      ! i4 = i5 =   27

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   25
                                      ! i4 = i2 =   28

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    3

                                      ! i5 = in(nn+1) =   31

                                      xin(33) = c10*xin(27) + xc00*xin(30) + c01*xin(29)
                                      yin(33) = c10*yin(27) + yc00*yin(30) + c01*yin(29)
                                      zin(33) = c10*zin(27) + zc00*zin(30) + c01*zin(29)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   31

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   34

                                      xin(36) = c10*xin(30) + xc00*xin(33) + c01*xin(32)
                                      yin(36) = c10*yin(30) + yc00*yin(33) + c01*yin(32)
                                      zin(36) = c10*zin(30) + zc00*zin(33) + c01*zin(32)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   31
                                      ! i4 = i5 =   34

                                      ! nn =    4

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   36

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

       eri_value(1) = eri_value(1) + d03bra(1)*d02ket(1)*(xin(12)*yin(1)*zin(1) + xin(24)*yin(13)*zin(13) + xin(36)*yin(25)*zin(25))
       eri_value(2) = eri_value(2) + d03bra(1)*d02ket(2)*(xin(10)*yin(3)*zin(1) + xin(22)*yin(15)*zin(13) + xin(34)*yin(27)*zin(25))
       eri_value(3) = eri_value(3) + d03bra(1)*d02ket(3)*(xin(10)*yin(1)*zin(3) + xin(22)*yin(13)*zin(15) + xin(34)*yin(25)*zin(27))
       eri_value(4) = eri_value(4) + d03bra(1)*d02ket(4)*(xin(11)*yin(2)*zin(1) + xin(23)*yin(14)*zin(13) + xin(35)*yin(26)*zin(25))
       eri_value(5) = eri_value(5) + d03bra(1)*d02ket(5)*(xin(11)*yin(1)*zin(2) + xin(23)*yin(13)*zin(14) + xin(35)*yin(25)*zin(26))
       eri_value(6) = eri_value(6) + d03bra(1)*d02ket(6)*(xin(10)*yin(2)*zin(2) + xin(22)*yin(14)*zin(14) + xin(34)*yin(26)*zin(26))
       eri_value(7) = eri_value(7) + d03bra(2)*d02ket(1)*(xin(3)*yin(10)*zin(1) + xin(15)*yin(22)*zin(13) + xin(27)*yin(34)*zin(25))
       eri_value(8) = eri_value(8) + d03bra(2)*d02ket(2)*(xin(1)*yin(12)*zin(1) + xin(13)*yin(24)*zin(13) + xin(25)*yin(36)*zin(25))
       eri_value(9) = eri_value(9) + d03bra(2)*d02ket(3)*(xin(1)*yin(10)*zin(3) + xin(13)*yin(22)*zin(15) + xin(25)*yin(34)*zin(27))
     eri_value(10) = eri_value(10) + d03bra(2)*d02ket(4)*(xin(2)*yin(11)*zin(1) + xin(14)*yin(23)*zin(13) + xin(26)*yin(35)*zin(25))
     eri_value(11) = eri_value(11) + d03bra(2)*d02ket(5)*(xin(2)*yin(10)*zin(2) + xin(14)*yin(22)*zin(14) + xin(26)*yin(34)*zin(26))
     eri_value(12) = eri_value(12) + d03bra(2)*d02ket(6)*(xin(1)*yin(11)*zin(2) + xin(13)*yin(23)*zin(14) + xin(25)*yin(35)*zin(26))
     eri_value(13) = eri_value(13) + d03bra(3)*d02ket(1)*(xin(3)*yin(1)*zin(10) + xin(15)*yin(13)*zin(22) + xin(27)*yin(25)*zin(34))
     eri_value(14) = eri_value(14) + d03bra(3)*d02ket(2)*(xin(1)*yin(3)*zin(10) + xin(13)*yin(15)*zin(22) + xin(25)*yin(27)*zin(34))
     eri_value(15) = eri_value(15) + d03bra(3)*d02ket(3)*(xin(1)*yin(1)*zin(12) + xin(13)*yin(13)*zin(24) + xin(25)*yin(25)*zin(36))
     eri_value(16) = eri_value(16) + d03bra(3)*d02ket(4)*(xin(2)*yin(2)*zin(10) + xin(14)*yin(14)*zin(22) + xin(26)*yin(26)*zin(34))
     eri_value(17) = eri_value(17) + d03bra(3)*d02ket(5)*(xin(2)*yin(1)*zin(11) + xin(14)*yin(13)*zin(23) + xin(26)*yin(25)*zin(35))
     eri_value(18) = eri_value(18) + d03bra(3)*d02ket(6)*(xin(1)*yin(2)*zin(11) + xin(13)*yin(14)*zin(23) + xin(25)*yin(26)*zin(35))
      eri_value(19) = eri_value(19) + d03bra(4)*d02ket(1)*(xin(9)*yin(4)*zin(1) + xin(21)*yin(16)*zin(13) + xin(33)*yin(28)*zin(25))
      eri_value(20) = eri_value(20) + d03bra(4)*d02ket(2)*(xin(7)*yin(6)*zin(1) + xin(19)*yin(18)*zin(13) + xin(31)*yin(30)*zin(25))
      eri_value(21) = eri_value(21) + d03bra(4)*d02ket(3)*(xin(7)*yin(4)*zin(3) + xin(19)*yin(16)*zin(15) + xin(31)*yin(28)*zin(27))
      eri_value(22) = eri_value(22) + d03bra(4)*d02ket(4)*(xin(8)*yin(5)*zin(1) + xin(20)*yin(17)*zin(13) + xin(32)*yin(29)*zin(25))
      eri_value(23) = eri_value(23) + d03bra(4)*d02ket(5)*(xin(8)*yin(4)*zin(2) + xin(20)*yin(16)*zin(14) + xin(32)*yin(28)*zin(26))
      eri_value(24) = eri_value(24) + d03bra(4)*d02ket(6)*(xin(7)*yin(5)*zin(2) + xin(19)*yin(17)*zin(14) + xin(31)*yin(29)*zin(26))
      eri_value(25) = eri_value(25) + d03bra(5)*d02ket(1)*(xin(9)*yin(1)*zin(4) + xin(21)*yin(13)*zin(16) + xin(33)*yin(25)*zin(28))
      eri_value(26) = eri_value(26) + d03bra(5)*d02ket(2)*(xin(7)*yin(3)*zin(4) + xin(19)*yin(15)*zin(16) + xin(31)*yin(27)*zin(28))
      eri_value(27) = eri_value(27) + d03bra(5)*d02ket(3)*(xin(7)*yin(1)*zin(6) + xin(19)*yin(13)*zin(18) + xin(31)*yin(25)*zin(30))
      eri_value(28) = eri_value(28) + d03bra(5)*d02ket(4)*(xin(8)*yin(2)*zin(4) + xin(20)*yin(14)*zin(16) + xin(32)*yin(26)*zin(28))
      eri_value(29) = eri_value(29) + d03bra(5)*d02ket(5)*(xin(8)*yin(1)*zin(5) + xin(20)*yin(13)*zin(17) + xin(32)*yin(25)*zin(29))
      eri_value(30) = eri_value(30) + d03bra(5)*d02ket(6)*(xin(7)*yin(2)*zin(5) + xin(19)*yin(14)*zin(17) + xin(31)*yin(26)*zin(29))
      eri_value(31) = eri_value(31) + d03bra(6)*d02ket(1)*(xin(6)*yin(7)*zin(1) + xin(18)*yin(19)*zin(13) + xin(30)*yin(31)*zin(25))
      eri_value(32) = eri_value(32) + d03bra(6)*d02ket(2)*(xin(4)*yin(9)*zin(1) + xin(16)*yin(21)*zin(13) + xin(28)*yin(33)*zin(25))
      eri_value(33) = eri_value(33) + d03bra(6)*d02ket(3)*(xin(4)*yin(7)*zin(3) + xin(16)*yin(19)*zin(15) + xin(28)*yin(31)*zin(27))
      eri_value(34) = eri_value(34) + d03bra(6)*d02ket(4)*(xin(5)*yin(8)*zin(1) + xin(17)*yin(20)*zin(13) + xin(29)*yin(32)*zin(25))
      eri_value(35) = eri_value(35) + d03bra(6)*d02ket(5)*(xin(5)*yin(7)*zin(2) + xin(17)*yin(19)*zin(14) + xin(29)*yin(31)*zin(26))
      eri_value(36) = eri_value(36) + d03bra(6)*d02ket(6)*(xin(4)*yin(8)*zin(2) + xin(16)*yin(20)*zin(14) + xin(28)*yin(32)*zin(26))
      eri_value(37) = eri_value(37) + d03bra(7)*d02ket(1)*(xin(3)*yin(7)*zin(4) + xin(15)*yin(19)*zin(16) + xin(27)*yin(31)*zin(28))
      eri_value(38) = eri_value(38) + d03bra(7)*d02ket(2)*(xin(1)*yin(9)*zin(4) + xin(13)*yin(21)*zin(16) + xin(25)*yin(33)*zin(28))
      eri_value(39) = eri_value(39) + d03bra(7)*d02ket(3)*(xin(1)*yin(7)*zin(6) + xin(13)*yin(19)*zin(18) + xin(25)*yin(31)*zin(30))
      eri_value(40) = eri_value(40) + d03bra(7)*d02ket(4)*(xin(2)*yin(8)*zin(4) + xin(14)*yin(20)*zin(16) + xin(26)*yin(32)*zin(28))
      eri_value(41) = eri_value(41) + d03bra(7)*d02ket(5)*(xin(2)*yin(7)*zin(5) + xin(14)*yin(19)*zin(17) + xin(26)*yin(31)*zin(29))
      eri_value(42) = eri_value(42) + d03bra(7)*d02ket(6)*(xin(1)*yin(8)*zin(5) + xin(13)*yin(20)*zin(17) + xin(25)*yin(32)*zin(29))
      eri_value(43) = eri_value(43) + d03bra(8)*d02ket(1)*(xin(6)*yin(1)*zin(7) + xin(18)*yin(13)*zin(19) + xin(30)*yin(25)*zin(31))
      eri_value(44) = eri_value(44) + d03bra(8)*d02ket(2)*(xin(4)*yin(3)*zin(7) + xin(16)*yin(15)*zin(19) + xin(28)*yin(27)*zin(31))
      eri_value(45) = eri_value(45) + d03bra(8)*d02ket(3)*(xin(4)*yin(1)*zin(9) + xin(16)*yin(13)*zin(21) + xin(28)*yin(25)*zin(33))
      eri_value(46) = eri_value(46) + d03bra(8)*d02ket(4)*(xin(5)*yin(2)*zin(7) + xin(17)*yin(14)*zin(19) + xin(29)*yin(26)*zin(31))
      eri_value(47) = eri_value(47) + d03bra(8)*d02ket(5)*(xin(5)*yin(1)*zin(8) + xin(17)*yin(13)*zin(20) + xin(29)*yin(25)*zin(32))
      eri_value(48) = eri_value(48) + d03bra(8)*d02ket(6)*(xin(4)*yin(2)*zin(8) + xin(16)*yin(14)*zin(20) + xin(28)*yin(26)*zin(32))
      eri_value(49) = eri_value(49) + d03bra(9)*d02ket(1)*(xin(3)*yin(4)*zin(7) + xin(15)*yin(16)*zin(19) + xin(27)*yin(28)*zin(31))
      eri_value(50) = eri_value(50) + d03bra(9)*d02ket(2)*(xin(1)*yin(6)*zin(7) + xin(13)*yin(18)*zin(19) + xin(25)*yin(30)*zin(31))
      eri_value(51) = eri_value(51) + d03bra(9)*d02ket(3)*(xin(1)*yin(4)*zin(9) + xin(13)*yin(16)*zin(21) + xin(25)*yin(28)*zin(33))
      eri_value(52) = eri_value(52) + d03bra(9)*d02ket(4)*(xin(2)*yin(5)*zin(7) + xin(14)*yin(17)*zin(19) + xin(26)*yin(29)*zin(31))
      eri_value(53) = eri_value(53) + d03bra(9)*d02ket(5)*(xin(2)*yin(4)*zin(8) + xin(14)*yin(16)*zin(20) + xin(26)*yin(28)*zin(32))
      eri_value(54) = eri_value(54) + d03bra(9)*d02ket(6)*(xin(1)*yin(5)*zin(8) + xin(13)*yin(17)*zin(20) + xin(25)*yin(29)*zin(32))
     eri_value(55) = eri_value(55) + d03bra(10)*d02ket(1)*(xin(6)*yin(4)*zin(4) + xin(18)*yin(16)*zin(16) + xin(30)*yin(28)*zin(28))
     eri_value(56) = eri_value(56) + d03bra(10)*d02ket(2)*(xin(4)*yin(6)*zin(4) + xin(16)*yin(18)*zin(16) + xin(28)*yin(30)*zin(28))
     eri_value(57) = eri_value(57) + d03bra(10)*d02ket(3)*(xin(4)*yin(4)*zin(6) + xin(16)*yin(16)*zin(18) + xin(28)*yin(28)*zin(30))
     eri_value(58) = eri_value(58) + d03bra(10)*d02ket(4)*(xin(5)*yin(5)*zin(4) + xin(17)*yin(17)*zin(16) + xin(29)*yin(29)*zin(28))
     eri_value(59) = eri_value(59) + d03bra(10)*d02ket(5)*(xin(5)*yin(4)*zin(5) + xin(17)*yin(16)*zin(17) + xin(29)*yin(28)*zin(29))
     eri_value(60) = eri_value(60) + d03bra(10)*d02ket(6)*(xin(4)*yin(5)*zin(5) + xin(16)*yin(17)*zin(17) + xin(28)*yin(29)*zin(29))

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
                                    ip = (i - 1)*6 ! Stride between functions in i

                                    do j = 1, 1 ! # of cartesians in j

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

                              deallocate (n03bra)
                              deallocate (xint03bra)
                              deallocate (n02ket)
                              deallocate (xint02ket)

                              end subroutine int3020
                              end submodule
