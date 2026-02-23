! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3030_impl
contains
  module subroutine int3030(sf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n03bra(:)
    real(dp), allocatable :: xint03bra(:)
    integer(kind=int64) :: nsfbra
    real(dp) :: scutsfbra, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iijj, kkll
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2, nij, nkl, itmp
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(64), yin(64), zin(64)
    real(dp) :: eri_value(100)
    real(dp) :: d03bra(10), d03ket(10)
    integer(kind=int64) :: ix(10), jx(1), kx(10), lx(1)
    integer(kind=int64) :: iy(10), jy(1), ky(10), ly(1)
    integer(kind=int64) :: iz(10), jz(1), kz(10), lz(1)
    integer(kind=int64) :: in(4), in1(4), kn(4)
    integer(kind=int64) :: ijx(10), ijy(10), ijz(10)
    integer(kind=int64) :: klx(10), kly(10), klz(10)
    logical :: same

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 5
    in1(3) = 9
    in1(4) = 13

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

    jx(1) = 0

    ix(1) = 13
    ix(2) = 1
    ix(3) = 1
    ix(4) = 9
    ix(5) = 9
    ix(6) = 5
    ix(7) = 1
    ix(8) = 5
    ix(9) = 1
    ix(10) = 5

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

    iy(1) = 1
    iy(2) = 13
    iy(3) = 1
    iy(4) = 5
    iy(5) = 1
    iy(6) = 9
    iy(7) = 9
    iy(8) = 1
    iy(9) = 5
    iy(10) = 5

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

    iz(1) = 1
    iz(2) = 1
    iz(3) = 13
    iz(4) = 1
    iz(5) = 5
    iz(6) = 1
    iz(7) = 5
    iz(8) = 9
    iz(9) = 9
    iz(10) = 5

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 13
    ijx(2) = 1
    ijx(3) = 1
    ijx(4) = 9
    ijx(5) = 9
    ijx(6) = 5
    ijx(7) = 1
    ijx(8) = 5
    ijx(9) = 1
    ijx(10) = 5

    ijy(1) = 1
    ijy(2) = 13
    ijy(3) = 1
    ijy(4) = 5
    ijy(5) = 1
    ijy(6) = 9
    ijy(7) = 9
    ijy(8) = 1
    ijy(9) = 5
    ijy(10) = 5

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 13
    ijz(4) = 1
    ijz(5) = 5
    ijz(6) = 1
    ijz(7) = 5
    ijz(8) = 9
    ijz(9) = 9
    ijz(10) = 5

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

    allocate (n03bra(res%n_s_shl*res%n_f_shl))
    allocate (xint03bra(res%n_s_shl*res%n_f_shl))

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

    ! --multi-gpu--work
    nchunk = nsfbra/res%n_size
    nquart_start = nchunk*res%n_rank + 1
    nquart_end = nquart_start + nchunk - 1
    if (res%n_rank .EQ. res%n_size - 1) nquart_end = nsfbra

    ! Mappings to GPU

 !$omp target teams distribute parallel do collapse(2) default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nsfbra, xint03bra, n03bra, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d03ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d03bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,nij,nkl,itmp,same)
              do iijj = nquart_start, nquart_end
                do kkll = 1, nsfbra

                  if (kkll .gt. iijj) cycle

                  test = xint03bra(iijj)*xint03bra(kkll)
                  if (test .gt. cutoff_schwarz) then

                    ij = n03bra(iijj)
                    kl = n03bra(kkll)

                    ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                    jsh_tmp = (ij - 1)/res%n_f_shl + 1
                    ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                    lsh_tmp = (kl - 1)/res%n_f_shl + 1

                    ish = res%i_f_shl(ish_tmp)
                    jsh = res%i_s_shl(jsh_tmp)
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

30                            continue

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

100                           continue

                              do 240 l = 1, 4

                                jj = 0

105                             do 110 m = l, 4
                                  if (m .eq. 4) go to 120
                                  if (abs(wrk(m)) .le. (1.0D-14)*(abs(rts(m)) + abs(rts(m + 1)))) go to 120
110                               continue

120                               dpp = rts(l)
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
150                                 ds = df/dg
                                    dr = sqrt(ds*ds + 1.0D+00)
                                    wrk(mmii + 1) = dg*dr
                                    dc = 1.0D+00/dr
                                    ds = ds*dc
160                                 dg = rts(mmii + 1) - dpp
                                    dr = (rts(mmii) - dg)*ds + 2.0D+00*dc*db
                                    dpp = ds*dr
                                    rts(mmii + 1) = dg + dpp
                                    dg = dc*dr - db
                                    df = wts(mmii + 1)
                                    wts(mmii + 1) = ds*wts(mmii) + dc*df
                                    wts(mmii) = dc*wts(mmii) - ds*df

200                                 continue

                                    rts(l) = rts(l) - dpp
                                    wrk(l) = dg
                                    wrk(m) = 0.0D+00
                                    go to 105

240                                 continue

                                    do 300 ii = 2, 4

                                      iim1 = ii - 1
                                      kk = iim1
                                      dpp = rts(iim1)

                                      do 260 jj = ii, 4
                                        if (rts(jj) .ge. dpp) go to 260
                                        kk = jj
                                        dpp = rts(jj)
260                                     continue

                                        if (kk .eq. iim1) go to 300

                                        rts(kk) = rts(iim1)
                                        rts(iim1) = dpp
                                        dpp = wts(iim1)
                                        wts(iim1) = wts(kk)
                                        wts(kk) = dpp

300                                     continue

                                        do 310 kk = 1, 4
                                          wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                       continue

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

                                          ! i2 = in(2) =    5
                                          ! k2 = kn(2) =    1
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(5) = xc00
                                          yin(5) = yc00
                                          zin(5) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =    2

                                          xin(2) = xcp00
                                          yin(2) = ycp00
                                          zin(2) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =    6
                                          ! i2 =    5

                                          xin(6) = xcp00*xin(5) + cp10
                                          yin(6) = ycp00*yin(5) + cp10
                                          zin(6) = zcp00*zin(5) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =    1
                                          ! i4 = i2 =    5

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =    9
                                          ! i3 =    1
                                          ! i4 =    5

                                          xin(9) = c10*xin(1) + xc00*xin(5)
                                          yin(9) = c10*yin(1) + yc00*yin(5)
                                          zin(9) = c10*zin(1) + zc00*zin(5)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   10
                                          ! i5 =    9
                                          ! i4 =    5

                                          xin(10) = xcp00*xin(9) + cp10*xin(5)
                                          yin(10) = ycp00*yin(9) + cp10*yin(5)
                                          zin(10) = zcp00*zin(9) + cp10*zin(5)

                                          ! ------------------

                                          ! i3 = i4 =    5
                                          ! i4 = i5 =    9

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   13
                                          ! i3 =    5
                                          ! i4 =    9

                                          xin(13) = c10*xin(5) + xc00*xin(9)
                                          yin(13) = c10*yin(5) + yc00*yin(9)
                                          zin(13) = c10*zin(5) + zc00*zin(9)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   14
                                          ! i5 =   13
                                          ! i4 =    9

                                          xin(14) = xcp00*xin(13) + cp10*xin(9)
                                          yin(14) = ycp00*yin(13) + cp10*yin(9)
                                          zin(14) = zcp00*zin(13) + cp10*zin(9)

                                          ! ------------------

                                          ! i3 = i4 =    9
                                          ! i4 = i5 =   13

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

                                          ! i3 = i2+kn(n+1) =    7

                                          xin(7) = xc00*xin(3) + c01*xin(2)
                                          yin(7) = yc00*yin(3) + c01*yin(2)
                                          zin(7) = zc00*zin(3) + c01*zin(2)

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

                                          ! i3 = i2+kn(n+1) =    8

                                          xin(8) = xc00*xin(4) + c01*xin(3)
                                          yin(8) = yc00*yin(4) + c01*yin(3)
                                          zin(8) = zc00*zin(4) + c01*zin(3)

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
                                          ! i4 = i2 =    5

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =    9

                                          xin(11) = c10*xin(3) + xc00*xin(7) + c01*xin(6)
                                          yin(11) = c10*yin(3) + yc00*yin(7) + c01*yin(6)
                                          zin(11) = c10*zin(3) + zc00*zin(7) + c01*zin(6)

                                          c10 = c10 + b10

                                          ! i3 = i4 =    5
                                          ! i4 = i5 =    9

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   13

                                          xin(15) = c10*xin(7) + xc00*xin(11) + c01*xin(10)
                                          yin(15) = c10*yin(7) + yc00*yin(11) + c01*yin(10)
                                          zin(15) = c10*zin(7) + zc00*zin(11) + c01*zin(10)

                                          c10 = c10 + b10

                                          ! i3 = i4 =    9
                                          ! i4 = i5 =   13

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   2

                                          ! n =    3

                                          ! k4 = kn(n+1) =    3
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =    5

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =    9

                                          xin(12) = c10*xin(4) + xc00*xin(8) + c01*xin(7)
                                          yin(12) = c10*yin(4) + yc00*yin(8) + c01*yin(7)
                                          zin(12) = c10*zin(4) + zc00*zin(8) + c01*zin(7)

                                          c10 = c10 + b10

                                          ! i3 = i4 =    5
                                          ! i4 = i5 =    9

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   13

                                          xin(16) = c10*xin(8) + xc00*xin(12) + c01*xin(11)
                                          yin(16) = c10*yin(8) + yc00*yin(12) + c01*yin(11)
                                          zin(16) = c10*zin(8) + zc00*zin(12) + c01*zin(11)

                                          c10 = c10 + b10

                                          ! i3 = i4 =    9
                                          ! i4 = i5 =   13

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   3

                                          ! n =    4

                                          ! end do

                                          ! *** Now root =    2

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   16

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

                                          ! i1 = in(1) =   17

                                          xin(17) = 1.0_dp
                                          yin(17) = 1.0_dp
                                          zin(17) = f00

                                          ! i2 = in(2) =   21
                                          ! k2 = kn(2) =    1
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(21) = xc00
                                          yin(21) = yc00
                                          zin(21) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =   18

                                          xin(18) = xcp00
                                          yin(18) = ycp00
                                          zin(18) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   22
                                          ! i2 =   21

                                          xin(22) = xcp00*xin(21) + cp10
                                          yin(22) = ycp00*yin(21) + cp10
                                          zin(22) = zcp00*zin(21) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =   17
                                          ! i4 = i2 =   21

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   25
                                          ! i3 =   17
                                          ! i4 =   21

                                          xin(25) = c10*xin(17) + xc00*xin(21)
                                          yin(25) = c10*yin(17) + yc00*yin(21)
                                          zin(25) = c10*zin(17) + zc00*zin(21)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   26
                                          ! i5 =   25
                                          ! i4 =   21

                                          xin(26) = xcp00*xin(25) + cp10*xin(21)
                                          yin(26) = ycp00*yin(25) + cp10*yin(21)
                                          zin(26) = zcp00*zin(25) + cp10*zin(21)

                                          ! ------------------

                                          ! i3 = i4 =   21
                                          ! i4 = i5 =   25

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   29
                                          ! i3 =   21
                                          ! i4 =   25

                                          xin(29) = c10*xin(21) + xc00*xin(25)
                                          yin(29) = c10*yin(21) + yc00*yin(25)
                                          zin(29) = c10*zin(21) + zc00*zin(25)

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

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =   17
                                          ! i4 = i1+k2 =   18

                                          ! do n = 2,    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   19
                                          ! i3 =   17
                                          ! i4 =   18

                                          xin(19) = cp01*xin(17) + xcp00*xin(18)
                                          yin(19) = cp01*yin(17) + ycp00*yin(18)
                                          zin(19) = cp01*zin(17) + zcp00*zin(18)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   23

                                          xin(23) = xc00*xin(19) + c01*xin(18)
                                          yin(23) = yc00*yin(19) + c01*yin(18)
                                          zin(23) = zc00*zin(19) + c01*zin(18)

                                          ! ------------------

                                          ! i3 = i4 =   18
                                          ! i4 = i5 =   19

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   20
                                          ! i3 =   18
                                          ! i4 =   19

                                          xin(20) = cp01*xin(18) + xcp00*xin(19)
                                          yin(20) = cp01*yin(18) + ycp00*yin(19)
                                          zin(20) = cp01*zin(18) + zcp00*zin(19)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   24

                                          xin(24) = xc00*xin(20) + c01*xin(19)
                                          yin(24) = yc00*yin(20) + c01*yin(19)
                                          zin(24) = zc00*zin(20) + c01*zin(19)

                                          ! ------------------

                                          ! i3 = i4 =   19
                                          ! i4 = i5 =   20

                                          ! n =    4

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    1

                                          ! do n = 2,    3

                                          ! k4 = kn(n+1) =    2
                                          ! i3 = i1 =   17
                                          ! i4 = i2 =   21

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   25

                                          xin(27) = c10*xin(19) + xc00*xin(23) + c01*xin(22)
                                          yin(27) = c10*yin(19) + yc00*yin(23) + c01*yin(22)
                                          zin(27) = c10*zin(19) + zc00*zin(23) + c01*zin(22)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   21
                                          ! i4 = i5 =   25

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   29

                                          xin(31) = c10*xin(23) + xc00*xin(27) + c01*xin(26)
                                          yin(31) = c10*yin(23) + yc00*yin(27) + c01*yin(26)
                                          zin(31) = c10*zin(23) + zc00*zin(27) + c01*zin(26)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   25
                                          ! i4 = i5 =   29

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   2

                                          ! n =    3

                                          ! k4 = kn(n+1) =    3
                                          ! i3 = i1 =   17
                                          ! i4 = i2 =   21

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   25

                                          xin(28) = c10*xin(20) + xc00*xin(24) + c01*xin(23)
                                          yin(28) = c10*yin(20) + yc00*yin(24) + c01*yin(23)
                                          zin(28) = c10*zin(20) + zc00*zin(24) + c01*zin(23)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   21
                                          ! i4 = i5 =   25

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   29

                                          xin(32) = c10*xin(24) + xc00*xin(28) + c01*xin(27)
                                          yin(32) = c10*yin(24) + yc00*yin(28) + c01*yin(27)
                                          zin(32) = c10*zin(24) + zc00*zin(28) + c01*zin(27)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   25
                                          ! i4 = i5 =   29

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   3

                                          ! n =    4

                                          ! end do

                                          ! *** Now root =    3

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   32

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

                                          ! i1 = in(1) =   33

                                          xin(33) = 1.0_dp
                                          yin(33) = 1.0_dp
                                          zin(33) = f00

                                          ! i2 = in(2) =   37
                                          ! k2 = kn(2) =    1
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(37) = xc00
                                          yin(37) = yc00
                                          zin(37) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =   34

                                          xin(34) = xcp00
                                          yin(34) = ycp00
                                          zin(34) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   38
                                          ! i2 =   37

                                          xin(38) = xcp00*xin(37) + cp10
                                          yin(38) = ycp00*yin(37) + cp10
                                          zin(38) = zcp00*zin(37) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =   33
                                          ! i4 = i2 =   37

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   41
                                          ! i3 =   33
                                          ! i4 =   37

                                          xin(41) = c10*xin(33) + xc00*xin(37)
                                          yin(41) = c10*yin(33) + yc00*yin(37)
                                          zin(41) = c10*zin(33) + zc00*zin(37)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   42
                                          ! i5 =   41
                                          ! i4 =   37

                                          xin(42) = xcp00*xin(41) + cp10*xin(37)
                                          yin(42) = ycp00*yin(41) + cp10*yin(37)
                                          zin(42) = zcp00*zin(41) + cp10*zin(37)

                                          ! ------------------

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   41

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   45
                                          ! i3 =   37
                                          ! i4 =   41

                                          xin(45) = c10*xin(37) + xc00*xin(41)
                                          yin(45) = c10*yin(37) + yc00*yin(41)
                                          zin(45) = c10*zin(37) + zc00*zin(41)

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

                                          ! i3 = i1 =   33
                                          ! i4 = i1+k2 =   34

                                          ! do n = 2,    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   35
                                          ! i3 =   33
                                          ! i4 =   34

                                          xin(35) = cp01*xin(33) + xcp00*xin(34)
                                          yin(35) = cp01*yin(33) + ycp00*yin(34)
                                          zin(35) = cp01*zin(33) + zcp00*zin(34)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   39

                                          xin(39) = xc00*xin(35) + c01*xin(34)
                                          yin(39) = yc00*yin(35) + c01*yin(34)
                                          zin(39) = zc00*zin(35) + c01*zin(34)

                                          ! ------------------

                                          ! i3 = i4 =   34
                                          ! i4 = i5 =   35

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   36
                                          ! i3 =   34
                                          ! i4 =   35

                                          xin(36) = cp01*xin(34) + xcp00*xin(35)
                                          yin(36) = cp01*yin(34) + ycp00*yin(35)
                                          zin(36) = cp01*zin(34) + zcp00*zin(35)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   40

                                          xin(40) = xc00*xin(36) + c01*xin(35)
                                          yin(40) = yc00*yin(36) + c01*yin(35)
                                          zin(40) = zc00*zin(36) + c01*zin(35)

                                          ! ------------------

                                          ! i3 = i4 =   35
                                          ! i4 = i5 =   36

                                          ! n =    4

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    1

                                          ! do n = 2,    3

                                          ! k4 = kn(n+1) =    2
                                          ! i3 = i1 =   33
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   41

                                          xin(43) = c10*xin(35) + xc00*xin(39) + c01*xin(38)
                                          yin(43) = c10*yin(35) + yc00*yin(39) + c01*yin(38)
                                          zin(43) = c10*zin(35) + zc00*zin(39) + c01*zin(38)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   41

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   45

                                          xin(47) = c10*xin(39) + xc00*xin(43) + c01*xin(42)
                                          yin(47) = c10*yin(39) + yc00*yin(43) + c01*yin(42)
                                          zin(47) = c10*zin(39) + zc00*zin(43) + c01*zin(42)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   41
                                          ! i4 = i5 =   45

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   2

                                          ! n =    3

                                          ! k4 = kn(n+1) =    3
                                          ! i3 = i1 =   33
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   41

                                          xin(44) = c10*xin(36) + xc00*xin(40) + c01*xin(39)
                                          yin(44) = c10*yin(36) + yc00*yin(40) + c01*yin(39)
                                          zin(44) = c10*zin(36) + zc00*zin(40) + c01*zin(39)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   41

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   45

                                          xin(48) = c10*xin(40) + xc00*xin(44) + c01*xin(43)
                                          yin(48) = c10*yin(40) + yc00*yin(44) + c01*yin(43)
                                          zin(48) = c10*zin(40) + zc00*zin(44) + c01*zin(43)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   41
                                          ! i4 = i5 =   45

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   3

                                          ! n =    4

                                          ! end do

                                          ! *** Now root =    4

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   48

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

                                          ! i1 = in(1) =   49

                                          xin(49) = 1.0_dp
                                          yin(49) = 1.0_dp
                                          zin(49) = f00

                                          ! i2 = in(2) =   53
                                          ! k2 = kn(2) =    1
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(53) = xc00
                                          yin(53) = yc00
                                          zin(53) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =   50

                                          xin(50) = xcp00
                                          yin(50) = ycp00
                                          zin(50) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   54
                                          ! i2 =   53

                                          xin(54) = xcp00*xin(53) + cp10
                                          yin(54) = ycp00*yin(53) + cp10
                                          zin(54) = zcp00*zin(53) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =   49
                                          ! i4 = i2 =   53

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   57
                                          ! i3 =   49
                                          ! i4 =   53

                                          xin(57) = c10*xin(49) + xc00*xin(53)
                                          yin(57) = c10*yin(49) + yc00*yin(53)
                                          zin(57) = c10*zin(49) + zc00*zin(53)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   58
                                          ! i5 =   57
                                          ! i4 =   53

                                          xin(58) = xcp00*xin(57) + cp10*xin(53)
                                          yin(58) = ycp00*yin(57) + cp10*yin(53)
                                          zin(58) = zcp00*zin(57) + cp10*zin(53)

                                          ! ------------------

                                          ! i3 = i4 =   53
                                          ! i4 = i5 =   57

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   61
                                          ! i3 =   53
                                          ! i4 =   57

                                          xin(61) = c10*xin(53) + xc00*xin(57)
                                          yin(61) = c10*yin(53) + yc00*yin(57)
                                          zin(61) = c10*zin(53) + zc00*zin(57)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   62
                                          ! i5 =   61
                                          ! i4 =   57

                                          xin(62) = xcp00*xin(61) + cp10*xin(57)
                                          yin(62) = ycp00*yin(61) + cp10*yin(57)
                                          zin(62) = zcp00*zin(61) + cp10*zin(57)

                                          ! ------------------

                                          ! i3 = i4 =   57
                                          ! i4 = i5 =   61

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

                                          ! i3 = i2+kn(n+1) =   55

                                          xin(55) = xc00*xin(51) + c01*xin(50)
                                          yin(55) = yc00*yin(51) + c01*yin(50)
                                          zin(55) = zc00*zin(51) + c01*zin(50)

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

                                          ! i3 = i2+kn(n+1) =   56

                                          xin(56) = xc00*xin(52) + c01*xin(51)
                                          yin(56) = yc00*yin(52) + c01*yin(51)
                                          zin(56) = zc00*zin(52) + c01*zin(51)

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
                                          ! i4 = i2 =   53

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   57

                                          xin(59) = c10*xin(51) + xc00*xin(55) + c01*xin(54)
                                          yin(59) = c10*yin(51) + yc00*yin(55) + c01*yin(54)
                                          zin(59) = c10*zin(51) + zc00*zin(55) + c01*zin(54)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   53
                                          ! i4 = i5 =   57

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   61

                                          xin(63) = c10*xin(55) + xc00*xin(59) + c01*xin(58)
                                          yin(63) = c10*yin(55) + yc00*yin(59) + c01*yin(58)
                                          zin(63) = c10*zin(55) + zc00*zin(59) + c01*zin(58)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   57
                                          ! i4 = i5 =   61

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   2

                                          ! n =    3

                                          ! k4 = kn(n+1) =    3
                                          ! i3 = i1 =   49
                                          ! i4 = i2 =   53

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   57

                                          xin(60) = c10*xin(52) + xc00*xin(56) + c01*xin(55)
                                          yin(60) = c10*yin(52) + yc00*yin(56) + c01*yin(55)
                                          zin(60) = c10*zin(52) + zc00*zin(56) + c01*zin(55)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   53
                                          ! i4 = i5 =   57

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   61

                                          xin(64) = c10*xin(56) + xc00*xin(60) + c01*xin(59)
                                          yin(64) = c10*yin(56) + yc00*yin(60) + c01*yin(59)
                                          zin(64) = c10*zin(56) + zc00*zin(60) + c01*zin(59)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   57
                                          ! i4 = i5 =   61

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   3

                                          ! n =    4

                                          ! end do

                                          ! *** Now root =    5

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   64

                                          !                     --- END XYZINT ---

                                          !                       --- FORMS ---
                                          ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d03bra(  1)*d03ket(  1)*(xin(  16)*yin(   1)*zin(   1)+xin(  32)*yin(  17)*zin(  17)+xin(  48)*yin(  33)*zin(  33)+xin(  64)*yin(  49)*zin(  49))
          eri_value(    2)=eri_value(    2)+d03bra(  1)*d03ket(  2)*(xin(  13)*yin(   4)*zin(   1)+xin(  29)*yin(  20)*zin(  17)+xin(  45)*yin(  36)*zin(  33)+xin(  61)*yin(  52)*zin(  49))
          eri_value(    3)=eri_value(    3)+d03bra(  1)*d03ket(  3)*(xin(  13)*yin(   1)*zin(   4)+xin(  29)*yin(  17)*zin(  20)+xin(  45)*yin(  33)*zin(  36)+xin(  61)*yin(  49)*zin(  52))
          eri_value(    4)=eri_value(    4)+d03bra(  1)*d03ket(  4)*(xin(  15)*yin(   2)*zin(   1)+xin(  31)*yin(  18)*zin(  17)+xin(  47)*yin(  34)*zin(  33)+xin(  63)*yin(  50)*zin(  49))
          eri_value(    5)=eri_value(    5)+d03bra(  1)*d03ket(  5)*(xin(  15)*yin(   1)*zin(   2)+xin(  31)*yin(  17)*zin(  18)+xin(  47)*yin(  33)*zin(  34)+xin(  63)*yin(  49)*zin(  50))
          eri_value(    6)=eri_value(    6)+d03bra(  1)*d03ket(  6)*(xin(  14)*yin(   3)*zin(   1)+xin(  30)*yin(  19)*zin(  17)+xin(  46)*yin(  35)*zin(  33)+xin(  62)*yin(  51)*zin(  49))
          eri_value(    7)=eri_value(    7)+d03bra(  1)*d03ket(  7)*(xin(  13)*yin(   3)*zin(   2)+xin(  29)*yin(  19)*zin(  18)+xin(  45)*yin(  35)*zin(  34)+xin(  61)*yin(  51)*zin(  50))
          eri_value(    8)=eri_value(    8)+d03bra(  1)*d03ket(  8)*(xin(  14)*yin(   1)*zin(   3)+xin(  30)*yin(  17)*zin(  19)+xin(  46)*yin(  33)*zin(  35)+xin(  62)*yin(  49)*zin(  51))
          eri_value(    9)=eri_value(    9)+d03bra(  1)*d03ket(  9)*(xin(  13)*yin(   2)*zin(   3)+xin(  29)*yin(  18)*zin(  19)+xin(  45)*yin(  34)*zin(  35)+xin(  61)*yin(  50)*zin(  51))
          eri_value(   10)=eri_value(   10)+d03bra(  1)*d03ket( 10)*(xin(  14)*yin(   2)*zin(   2)+xin(  30)*yin(  18)*zin(  18)+xin(  46)*yin(  34)*zin(  34)+xin(  62)*yin(  50)*zin(  50))
          eri_value(   11)=eri_value(   11)+d03bra(  2)*d03ket(  1)*(xin(   4)*yin(  13)*zin(   1)+xin(  20)*yin(  29)*zin(  17)+xin(  36)*yin(  45)*zin(  33)+xin(  52)*yin(  61)*zin(  49))
          eri_value(   12)=eri_value(   12)+d03bra(  2)*d03ket(  2)*(xin(   1)*yin(  16)*zin(   1)+xin(  17)*yin(  32)*zin(  17)+xin(  33)*yin(  48)*zin(  33)+xin(  49)*yin(  64)*zin(  49))
          eri_value(   13)=eri_value(   13)+d03bra(  2)*d03ket(  3)*(xin(   1)*yin(  13)*zin(   4)+xin(  17)*yin(  29)*zin(  20)+xin(  33)*yin(  45)*zin(  36)+xin(  49)*yin(  61)*zin(  52))
          eri_value(   14)=eri_value(   14)+d03bra(  2)*d03ket(  4)*(xin(   3)*yin(  14)*zin(   1)+xin(  19)*yin(  30)*zin(  17)+xin(  35)*yin(  46)*zin(  33)+xin(  51)*yin(  62)*zin(  49))
          eri_value(   15)=eri_value(   15)+d03bra(  2)*d03ket(  5)*(xin(   3)*yin(  13)*zin(   2)+xin(  19)*yin(  29)*zin(  18)+xin(  35)*yin(  45)*zin(  34)+xin(  51)*yin(  61)*zin(  50))
          eri_value(   16)=eri_value(   16)+d03bra(  2)*d03ket(  6)*(xin(   2)*yin(  15)*zin(   1)+xin(  18)*yin(  31)*zin(  17)+xin(  34)*yin(  47)*zin(  33)+xin(  50)*yin(  63)*zin(  49))
          eri_value(   17)=eri_value(   17)+d03bra(  2)*d03ket(  7)*(xin(   1)*yin(  15)*zin(   2)+xin(  17)*yin(  31)*zin(  18)+xin(  33)*yin(  47)*zin(  34)+xin(  49)*yin(  63)*zin(  50))
          eri_value(   18)=eri_value(   18)+d03bra(  2)*d03ket(  8)*(xin(   2)*yin(  13)*zin(   3)+xin(  18)*yin(  29)*zin(  19)+xin(  34)*yin(  45)*zin(  35)+xin(  50)*yin(  61)*zin(  51))
          eri_value(   19)=eri_value(   19)+d03bra(  2)*d03ket(  9)*(xin(   1)*yin(  14)*zin(   3)+xin(  17)*yin(  30)*zin(  19)+xin(  33)*yin(  46)*zin(  35)+xin(  49)*yin(  62)*zin(  51))
          eri_value(   20)=eri_value(   20)+d03bra(  2)*d03ket( 10)*(xin(   2)*yin(  14)*zin(   2)+xin(  18)*yin(  30)*zin(  18)+xin(  34)*yin(  46)*zin(  34)+xin(  50)*yin(  62)*zin(  50))
          eri_value(   21)=eri_value(   21)+d03bra(  3)*d03ket(  1)*(xin(   4)*yin(   1)*zin(  13)+xin(  20)*yin(  17)*zin(  29)+xin(  36)*yin(  33)*zin(  45)+xin(  52)*yin(  49)*zin(  61))
          eri_value(   22)=eri_value(   22)+d03bra(  3)*d03ket(  2)*(xin(   1)*yin(   4)*zin(  13)+xin(  17)*yin(  20)*zin(  29)+xin(  33)*yin(  36)*zin(  45)+xin(  49)*yin(  52)*zin(  61))
          eri_value(   23)=eri_value(   23)+d03bra(  3)*d03ket(  3)*(xin(   1)*yin(   1)*zin(  16)+xin(  17)*yin(  17)*zin(  32)+xin(  33)*yin(  33)*zin(  48)+xin(  49)*yin(  49)*zin(  64))
          eri_value(   24)=eri_value(   24)+d03bra(  3)*d03ket(  4)*(xin(   3)*yin(   2)*zin(  13)+xin(  19)*yin(  18)*zin(  29)+xin(  35)*yin(  34)*zin(  45)+xin(  51)*yin(  50)*zin(  61))
          eri_value(   25)=eri_value(   25)+d03bra(  3)*d03ket(  5)*(xin(   3)*yin(   1)*zin(  14)+xin(  19)*yin(  17)*zin(  30)+xin(  35)*yin(  33)*zin(  46)+xin(  51)*yin(  49)*zin(  62))
          eri_value(   26)=eri_value(   26)+d03bra(  3)*d03ket(  6)*(xin(   2)*yin(   3)*zin(  13)+xin(  18)*yin(  19)*zin(  29)+xin(  34)*yin(  35)*zin(  45)+xin(  50)*yin(  51)*zin(  61))
          eri_value(   27)=eri_value(   27)+d03bra(  3)*d03ket(  7)*(xin(   1)*yin(   3)*zin(  14)+xin(  17)*yin(  19)*zin(  30)+xin(  33)*yin(  35)*zin(  46)+xin(  49)*yin(  51)*zin(  62))
          eri_value(   28)=eri_value(   28)+d03bra(  3)*d03ket(  8)*(xin(   2)*yin(   1)*zin(  15)+xin(  18)*yin(  17)*zin(  31)+xin(  34)*yin(  33)*zin(  47)+xin(  50)*yin(  49)*zin(  63))
          eri_value(   29)=eri_value(   29)+d03bra(  3)*d03ket(  9)*(xin(   1)*yin(   2)*zin(  15)+xin(  17)*yin(  18)*zin(  31)+xin(  33)*yin(  34)*zin(  47)+xin(  49)*yin(  50)*zin(  63))
          eri_value(   30)=eri_value(   30)+d03bra(  3)*d03ket( 10)*(xin(   2)*yin(   2)*zin(  14)+xin(  18)*yin(  18)*zin(  30)+xin(  34)*yin(  34)*zin(  46)+xin(  50)*yin(  50)*zin(  62))
          eri_value(   31)=eri_value(   31)+d03bra(  4)*d03ket(  1)*(xin(  12)*yin(   5)*zin(   1)+xin(  28)*yin(  21)*zin(  17)+xin(  44)*yin(  37)*zin(  33)+xin(  60)*yin(  53)*zin(  49))
          eri_value(   32)=eri_value(   32)+d03bra(  4)*d03ket(  2)*(xin(   9)*yin(   8)*zin(   1)+xin(  25)*yin(  24)*zin(  17)+xin(  41)*yin(  40)*zin(  33)+xin(  57)*yin(  56)*zin(  49))
          eri_value(   33)=eri_value(   33)+d03bra(  4)*d03ket(  3)*(xin(   9)*yin(   5)*zin(   4)+xin(  25)*yin(  21)*zin(  20)+xin(  41)*yin(  37)*zin(  36)+xin(  57)*yin(  53)*zin(  52))
          eri_value(   34)=eri_value(   34)+d03bra(  4)*d03ket(  4)*(xin(  11)*yin(   6)*zin(   1)+xin(  27)*yin(  22)*zin(  17)+xin(  43)*yin(  38)*zin(  33)+xin(  59)*yin(  54)*zin(  49))
          eri_value(   35)=eri_value(   35)+d03bra(  4)*d03ket(  5)*(xin(  11)*yin(   5)*zin(   2)+xin(  27)*yin(  21)*zin(  18)+xin(  43)*yin(  37)*zin(  34)+xin(  59)*yin(  53)*zin(  50))
          eri_value(   36)=eri_value(   36)+d03bra(  4)*d03ket(  6)*(xin(  10)*yin(   7)*zin(   1)+xin(  26)*yin(  23)*zin(  17)+xin(  42)*yin(  39)*zin(  33)+xin(  58)*yin(  55)*zin(  49))
          eri_value(   37)=eri_value(   37)+d03bra(  4)*d03ket(  7)*(xin(   9)*yin(   7)*zin(   2)+xin(  25)*yin(  23)*zin(  18)+xin(  41)*yin(  39)*zin(  34)+xin(  57)*yin(  55)*zin(  50))
          eri_value(   38)=eri_value(   38)+d03bra(  4)*d03ket(  8)*(xin(  10)*yin(   5)*zin(   3)+xin(  26)*yin(  21)*zin(  19)+xin(  42)*yin(  37)*zin(  35)+xin(  58)*yin(  53)*zin(  51))
          eri_value(   39)=eri_value(   39)+d03bra(  4)*d03ket(  9)*(xin(   9)*yin(   6)*zin(   3)+xin(  25)*yin(  22)*zin(  19)+xin(  41)*yin(  38)*zin(  35)+xin(  57)*yin(  54)*zin(  51))
          eri_value(   40)=eri_value(   40)+d03bra(  4)*d03ket( 10)*(xin(  10)*yin(   6)*zin(   2)+xin(  26)*yin(  22)*zin(  18)+xin(  42)*yin(  38)*zin(  34)+xin(  58)*yin(  54)*zin(  50))
          eri_value(   41)=eri_value(   41)+d03bra(  5)*d03ket(  1)*(xin(  12)*yin(   1)*zin(   5)+xin(  28)*yin(  17)*zin(  21)+xin(  44)*yin(  33)*zin(  37)+xin(  60)*yin(  49)*zin(  53))
          eri_value(   42)=eri_value(   42)+d03bra(  5)*d03ket(  2)*(xin(   9)*yin(   4)*zin(   5)+xin(  25)*yin(  20)*zin(  21)+xin(  41)*yin(  36)*zin(  37)+xin(  57)*yin(  52)*zin(  53))
          eri_value(   43)=eri_value(   43)+d03bra(  5)*d03ket(  3)*(xin(   9)*yin(   1)*zin(   8)+xin(  25)*yin(  17)*zin(  24)+xin(  41)*yin(  33)*zin(  40)+xin(  57)*yin(  49)*zin(  56))
          eri_value(   44)=eri_value(   44)+d03bra(  5)*d03ket(  4)*(xin(  11)*yin(   2)*zin(   5)+xin(  27)*yin(  18)*zin(  21)+xin(  43)*yin(  34)*zin(  37)+xin(  59)*yin(  50)*zin(  53))
          eri_value(   45)=eri_value(   45)+d03bra(  5)*d03ket(  5)*(xin(  11)*yin(   1)*zin(   6)+xin(  27)*yin(  17)*zin(  22)+xin(  43)*yin(  33)*zin(  38)+xin(  59)*yin(  49)*zin(  54))
          eri_value(   46)=eri_value(   46)+d03bra(  5)*d03ket(  6)*(xin(  10)*yin(   3)*zin(   5)+xin(  26)*yin(  19)*zin(  21)+xin(  42)*yin(  35)*zin(  37)+xin(  58)*yin(  51)*zin(  53))
          eri_value(   47)=eri_value(   47)+d03bra(  5)*d03ket(  7)*(xin(   9)*yin(   3)*zin(   6)+xin(  25)*yin(  19)*zin(  22)+xin(  41)*yin(  35)*zin(  38)+xin(  57)*yin(  51)*zin(  54))
          eri_value(   48)=eri_value(   48)+d03bra(  5)*d03ket(  8)*(xin(  10)*yin(   1)*zin(   7)+xin(  26)*yin(  17)*zin(  23)+xin(  42)*yin(  33)*zin(  39)+xin(  58)*yin(  49)*zin(  55))
          eri_value(   49)=eri_value(   49)+d03bra(  5)*d03ket(  9)*(xin(   9)*yin(   2)*zin(   7)+xin(  25)*yin(  18)*zin(  23)+xin(  41)*yin(  34)*zin(  39)+xin(  57)*yin(  50)*zin(  55))
          eri_value(   50)=eri_value(   50)+d03bra(  5)*d03ket( 10)*(xin(  10)*yin(   2)*zin(   6)+xin(  26)*yin(  18)*zin(  22)+xin(  42)*yin(  34)*zin(  38)+xin(  58)*yin(  50)*zin(  54))
          eri_value(   51)=eri_value(   51)+d03bra(  6)*d03ket(  1)*(xin(   8)*yin(   9)*zin(   1)+xin(  24)*yin(  25)*zin(  17)+xin(  40)*yin(  41)*zin(  33)+xin(  56)*yin(  57)*zin(  49))
          eri_value(   52)=eri_value(   52)+d03bra(  6)*d03ket(  2)*(xin(   5)*yin(  12)*zin(   1)+xin(  21)*yin(  28)*zin(  17)+xin(  37)*yin(  44)*zin(  33)+xin(  53)*yin(  60)*zin(  49))
          eri_value(   53)=eri_value(   53)+d03bra(  6)*d03ket(  3)*(xin(   5)*yin(   9)*zin(   4)+xin(  21)*yin(  25)*zin(  20)+xin(  37)*yin(  41)*zin(  36)+xin(  53)*yin(  57)*zin(  52))
          eri_value(   54)=eri_value(   54)+d03bra(  6)*d03ket(  4)*(xin(   7)*yin(  10)*zin(   1)+xin(  23)*yin(  26)*zin(  17)+xin(  39)*yin(  42)*zin(  33)+xin(  55)*yin(  58)*zin(  49))
          eri_value(   55)=eri_value(   55)+d03bra(  6)*d03ket(  5)*(xin(   7)*yin(   9)*zin(   2)+xin(  23)*yin(  25)*zin(  18)+xin(  39)*yin(  41)*zin(  34)+xin(  55)*yin(  57)*zin(  50))
          eri_value(   56)=eri_value(   56)+d03bra(  6)*d03ket(  6)*(xin(   6)*yin(  11)*zin(   1)+xin(  22)*yin(  27)*zin(  17)+xin(  38)*yin(  43)*zin(  33)+xin(  54)*yin(  59)*zin(  49))
          eri_value(   57)=eri_value(   57)+d03bra(  6)*d03ket(  7)*(xin(   5)*yin(  11)*zin(   2)+xin(  21)*yin(  27)*zin(  18)+xin(  37)*yin(  43)*zin(  34)+xin(  53)*yin(  59)*zin(  50))
          eri_value(   58)=eri_value(   58)+d03bra(  6)*d03ket(  8)*(xin(   6)*yin(   9)*zin(   3)+xin(  22)*yin(  25)*zin(  19)+xin(  38)*yin(  41)*zin(  35)+xin(  54)*yin(  57)*zin(  51))
          eri_value(   59)=eri_value(   59)+d03bra(  6)*d03ket(  9)*(xin(   5)*yin(  10)*zin(   3)+xin(  21)*yin(  26)*zin(  19)+xin(  37)*yin(  42)*zin(  35)+xin(  53)*yin(  58)*zin(  51))
          eri_value(   60)=eri_value(   60)+d03bra(  6)*d03ket( 10)*(xin(   6)*yin(  10)*zin(   2)+xin(  22)*yin(  26)*zin(  18)+xin(  38)*yin(  42)*zin(  34)+xin(  54)*yin(  58)*zin(  50))
          eri_value(   61)=eri_value(   61)+d03bra(  7)*d03ket(  1)*(xin(   4)*yin(   9)*zin(   5)+xin(  20)*yin(  25)*zin(  21)+xin(  36)*yin(  41)*zin(  37)+xin(  52)*yin(  57)*zin(  53))
          eri_value(   62)=eri_value(   62)+d03bra(  7)*d03ket(  2)*(xin(   1)*yin(  12)*zin(   5)+xin(  17)*yin(  28)*zin(  21)+xin(  33)*yin(  44)*zin(  37)+xin(  49)*yin(  60)*zin(  53))
          eri_value(   63)=eri_value(   63)+d03bra(  7)*d03ket(  3)*(xin(   1)*yin(   9)*zin(   8)+xin(  17)*yin(  25)*zin(  24)+xin(  33)*yin(  41)*zin(  40)+xin(  49)*yin(  57)*zin(  56))
          eri_value(   64)=eri_value(   64)+d03bra(  7)*d03ket(  4)*(xin(   3)*yin(  10)*zin(   5)+xin(  19)*yin(  26)*zin(  21)+xin(  35)*yin(  42)*zin(  37)+xin(  51)*yin(  58)*zin(  53))
          eri_value(   65)=eri_value(   65)+d03bra(  7)*d03ket(  5)*(xin(   3)*yin(   9)*zin(   6)+xin(  19)*yin(  25)*zin(  22)+xin(  35)*yin(  41)*zin(  38)+xin(  51)*yin(  57)*zin(  54))
          eri_value(   66)=eri_value(   66)+d03bra(  7)*d03ket(  6)*(xin(   2)*yin(  11)*zin(   5)+xin(  18)*yin(  27)*zin(  21)+xin(  34)*yin(  43)*zin(  37)+xin(  50)*yin(  59)*zin(  53))
          eri_value(   67)=eri_value(   67)+d03bra(  7)*d03ket(  7)*(xin(   1)*yin(  11)*zin(   6)+xin(  17)*yin(  27)*zin(  22)+xin(  33)*yin(  43)*zin(  38)+xin(  49)*yin(  59)*zin(  54))
          eri_value(   68)=eri_value(   68)+d03bra(  7)*d03ket(  8)*(xin(   2)*yin(   9)*zin(   7)+xin(  18)*yin(  25)*zin(  23)+xin(  34)*yin(  41)*zin(  39)+xin(  50)*yin(  57)*zin(  55))
          eri_value(   69)=eri_value(   69)+d03bra(  7)*d03ket(  9)*(xin(   1)*yin(  10)*zin(   7)+xin(  17)*yin(  26)*zin(  23)+xin(  33)*yin(  42)*zin(  39)+xin(  49)*yin(  58)*zin(  55))
          eri_value(   70)=eri_value(   70)+d03bra(  7)*d03ket( 10)*(xin(   2)*yin(  10)*zin(   6)+xin(  18)*yin(  26)*zin(  22)+xin(  34)*yin(  42)*zin(  38)+xin(  50)*yin(  58)*zin(  54))
          eri_value(   71)=eri_value(   71)+d03bra(  8)*d03ket(  1)*(xin(   8)*yin(   1)*zin(   9)+xin(  24)*yin(  17)*zin(  25)+xin(  40)*yin(  33)*zin(  41)+xin(  56)*yin(  49)*zin(  57))
          eri_value(   72)=eri_value(   72)+d03bra(  8)*d03ket(  2)*(xin(   5)*yin(   4)*zin(   9)+xin(  21)*yin(  20)*zin(  25)+xin(  37)*yin(  36)*zin(  41)+xin(  53)*yin(  52)*zin(  57))
          eri_value(   73)=eri_value(   73)+d03bra(  8)*d03ket(  3)*(xin(   5)*yin(   1)*zin(  12)+xin(  21)*yin(  17)*zin(  28)+xin(  37)*yin(  33)*zin(  44)+xin(  53)*yin(  49)*zin(  60))
          eri_value(   74)=eri_value(   74)+d03bra(  8)*d03ket(  4)*(xin(   7)*yin(   2)*zin(   9)+xin(  23)*yin(  18)*zin(  25)+xin(  39)*yin(  34)*zin(  41)+xin(  55)*yin(  50)*zin(  57))
          eri_value(   75)=eri_value(   75)+d03bra(  8)*d03ket(  5)*(xin(   7)*yin(   1)*zin(  10)+xin(  23)*yin(  17)*zin(  26)+xin(  39)*yin(  33)*zin(  42)+xin(  55)*yin(  49)*zin(  58))
          eri_value(   76)=eri_value(   76)+d03bra(  8)*d03ket(  6)*(xin(   6)*yin(   3)*zin(   9)+xin(  22)*yin(  19)*zin(  25)+xin(  38)*yin(  35)*zin(  41)+xin(  54)*yin(  51)*zin(  57))
          eri_value(   77)=eri_value(   77)+d03bra(  8)*d03ket(  7)*(xin(   5)*yin(   3)*zin(  10)+xin(  21)*yin(  19)*zin(  26)+xin(  37)*yin(  35)*zin(  42)+xin(  53)*yin(  51)*zin(  58))
          eri_value(   78)=eri_value(   78)+d03bra(  8)*d03ket(  8)*(xin(   6)*yin(   1)*zin(  11)+xin(  22)*yin(  17)*zin(  27)+xin(  38)*yin(  33)*zin(  43)+xin(  54)*yin(  49)*zin(  59))
          eri_value(   79)=eri_value(   79)+d03bra(  8)*d03ket(  9)*(xin(   5)*yin(   2)*zin(  11)+xin(  21)*yin(  18)*zin(  27)+xin(  37)*yin(  34)*zin(  43)+xin(  53)*yin(  50)*zin(  59))
          eri_value(   80)=eri_value(   80)+d03bra(  8)*d03ket( 10)*(xin(   6)*yin(   2)*zin(  10)+xin(  22)*yin(  18)*zin(  26)+xin(  38)*yin(  34)*zin(  42)+xin(  54)*yin(  50)*zin(  58))
          eri_value(   81)=eri_value(   81)+d03bra(  9)*d03ket(  1)*(xin(   4)*yin(   5)*zin(   9)+xin(  20)*yin(  21)*zin(  25)+xin(  36)*yin(  37)*zin(  41)+xin(  52)*yin(  53)*zin(  57))
          eri_value(   82)=eri_value(   82)+d03bra(  9)*d03ket(  2)*(xin(   1)*yin(   8)*zin(   9)+xin(  17)*yin(  24)*zin(  25)+xin(  33)*yin(  40)*zin(  41)+xin(  49)*yin(  56)*zin(  57))
          eri_value(   83)=eri_value(   83)+d03bra(  9)*d03ket(  3)*(xin(   1)*yin(   5)*zin(  12)+xin(  17)*yin(  21)*zin(  28)+xin(  33)*yin(  37)*zin(  44)+xin(  49)*yin(  53)*zin(  60))
          eri_value(   84)=eri_value(   84)+d03bra(  9)*d03ket(  4)*(xin(   3)*yin(   6)*zin(   9)+xin(  19)*yin(  22)*zin(  25)+xin(  35)*yin(  38)*zin(  41)+xin(  51)*yin(  54)*zin(  57))
          eri_value(   85)=eri_value(   85)+d03bra(  9)*d03ket(  5)*(xin(   3)*yin(   5)*zin(  10)+xin(  19)*yin(  21)*zin(  26)+xin(  35)*yin(  37)*zin(  42)+xin(  51)*yin(  53)*zin(  58))
          eri_value(   86)=eri_value(   86)+d03bra(  9)*d03ket(  6)*(xin(   2)*yin(   7)*zin(   9)+xin(  18)*yin(  23)*zin(  25)+xin(  34)*yin(  39)*zin(  41)+xin(  50)*yin(  55)*zin(  57))
          eri_value(   87)=eri_value(   87)+d03bra(  9)*d03ket(  7)*(xin(   1)*yin(   7)*zin(  10)+xin(  17)*yin(  23)*zin(  26)+xin(  33)*yin(  39)*zin(  42)+xin(  49)*yin(  55)*zin(  58))
          eri_value(   88)=eri_value(   88)+d03bra(  9)*d03ket(  8)*(xin(   2)*yin(   5)*zin(  11)+xin(  18)*yin(  21)*zin(  27)+xin(  34)*yin(  37)*zin(  43)+xin(  50)*yin(  53)*zin(  59))
          eri_value(   89)=eri_value(   89)+d03bra(  9)*d03ket(  9)*(xin(   1)*yin(   6)*zin(  11)+xin(  17)*yin(  22)*zin(  27)+xin(  33)*yin(  38)*zin(  43)+xin(  49)*yin(  54)*zin(  59))
          eri_value(   90)=eri_value(   90)+d03bra(  9)*d03ket( 10)*(xin(   2)*yin(   6)*zin(  10)+xin(  18)*yin(  22)*zin(  26)+xin(  34)*yin(  38)*zin(  42)+xin(  50)*yin(  54)*zin(  58))
          eri_value(   91)=eri_value(   91)+d03bra( 10)*d03ket(  1)*(xin(   8)*yin(   5)*zin(   5)+xin(  24)*yin(  21)*zin(  21)+xin(  40)*yin(  37)*zin(  37)+xin(  56)*yin(  53)*zin(  53))
          eri_value(   92)=eri_value(   92)+d03bra( 10)*d03ket(  2)*(xin(   5)*yin(   8)*zin(   5)+xin(  21)*yin(  24)*zin(  21)+xin(  37)*yin(  40)*zin(  37)+xin(  53)*yin(  56)*zin(  53))
          eri_value(   93)=eri_value(   93)+d03bra( 10)*d03ket(  3)*(xin(   5)*yin(   5)*zin(   8)+xin(  21)*yin(  21)*zin(  24)+xin(  37)*yin(  37)*zin(  40)+xin(  53)*yin(  53)*zin(  56))
          eri_value(   94)=eri_value(   94)+d03bra( 10)*d03ket(  4)*(xin(   7)*yin(   6)*zin(   5)+xin(  23)*yin(  22)*zin(  21)+xin(  39)*yin(  38)*zin(  37)+xin(  55)*yin(  54)*zin(  53))
          eri_value(   95)=eri_value(   95)+d03bra( 10)*d03ket(  5)*(xin(   7)*yin(   5)*zin(   6)+xin(  23)*yin(  21)*zin(  22)+xin(  39)*yin(  37)*zin(  38)+xin(  55)*yin(  53)*zin(  54))
          eri_value(   96)=eri_value(   96)+d03bra( 10)*d03ket(  6)*(xin(   6)*yin(   7)*zin(   5)+xin(  22)*yin(  23)*zin(  21)+xin(  38)*yin(  39)*zin(  37)+xin(  54)*yin(  55)*zin(  53))
          eri_value(   97)=eri_value(   97)+d03bra( 10)*d03ket(  7)*(xin(   5)*yin(   7)*zin(   6)+xin(  21)*yin(  23)*zin(  22)+xin(  37)*yin(  39)*zin(  38)+xin(  53)*yin(  55)*zin(  54))
          eri_value(   98)=eri_value(   98)+d03bra( 10)*d03ket(  8)*(xin(   6)*yin(   5)*zin(   7)+xin(  22)*yin(  21)*zin(  23)+xin(  38)*yin(  37)*zin(  39)+xin(  54)*yin(  53)*zin(  55))
          eri_value(   99)=eri_value(   99)+d03bra( 10)*d03ket(  9)*(xin(   5)*yin(   6)*zin(   7)+xin(  21)*yin(  22)*zin(  23)+xin(  37)*yin(  38)*zin(  39)+xin(  53)*yin(  54)*zin(  55))
          eri_value(  100)=eri_value(  100)+d03bra( 10)*d03ket( 10)*(xin(   6)*yin(   6)*zin(   6)+xin(  22)*yin(  22)*zin(  22)+xin(  38)*yin(  38)*zin(  38)+xin(  54)*yin(  54)*zin(  54))

                                          !                     --- END FORMS ---

                                        end do ! ij primitve loop

                                      end do ! kl primitve loop

                                      !                     --- DIRFCK_RHF ---
                                      !          Compute Fock matrix elements from 2EIs

                                      maxl = 1
                                      same = (ish .eq. ksh) .and. (jsh .eq. lsh)

                                      loci = res%atom_loc(ish) - 1
                                      locj = res%atom_loc(jsh) - 1
                                      lock = res%atom_loc(ksh) - 1
                                      locl = res%atom_loc(lsh) - 1

                                      nij = 0

                                      do i = 1, 10 ! # of cartesians in i

                                        ii1 = i + loci
                                        ip = (i - 1)*10 ! Stride between functions in i

                                        do j = 1, 1 ! # of cartesians in j

                                          nij = nij + 1

                                          maxl2 = maxl

                                          jj1 = j + locj
                                          i2 = ii1
                                          j2 = jj1
                                          if (ii1 .lt. jj1) then ! Sort <ij|
                                            i2 = jj1
                                            j2 = ii1
                                          end if

                                          ijp = (j - 1)*10 + ip ! Add stride between functions in j

                                          nkl = nij

                                          do k = 1, 10 ! # of cartesians in k

                                            kk1 = k + lock

                                            ijkp = (k - 1)*1 + ijp ! Add stride between functions in k

                                            if (same) then ! Account for non-unique permutations
                                              itmp = min(maxl2 - 1 + 1, nkl)
                                              if (itmp .eq. 0) exit ! Move to the next j iteration
                                              maxl2 = 1 + itmp - 1
                                              nkl = nkl - itmp
                                            end if

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

                                    end do ! kkll - or single loop for Do Concurrent

                                end do ! iijj
                                !$omp end target teams distribute parallel do


                                deallocate (n03bra)
                                deallocate (xint03bra)

                                end subroutine int3030
                                end submodule
