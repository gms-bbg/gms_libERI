! The total angular momentum of this class is:           3
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3000_impl
contains
  module subroutine int3000(sf_pair, ss_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sf_pair, ss_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n03bra(:), n00ket(:)
    real(dp), allocatable :: xint03bra(:), xint00ket(:)
    integer(kind=int64) :: nsfbra, nssket
    real(dp) :: scutsfbra, scutssket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
    integer(kind=int64) :: n, i1, i3, i4, i5
    real(dp) :: cp10, c10
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
    real(dp) :: roots(2), wghts(2)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(25), wgrid(25), p0(25), p1(25), p2(25)
    real(dp) :: rts(2), wts(2), alpha(2), beta(2), wrk(2)
    real(dp) :: xin(8), yin(8), zin(8)
    real(dp) :: eri_value(10)
    real(dp) :: d03bra(10), d00ket(1)
    integer(kind=int64) :: ix(10), jx(1), kx(1), lx(1)
    integer(kind=int64) :: iy(10), jy(1), ky(1), ly(1)
    integer(kind=int64) :: iz(10), jz(1), kz(1), lz(1)
    integer(kind=int64) :: in(4), in1(4), kn(1)
    integer(kind=int64) :: ijx(10), ijy(10), ijz(10)
    integer(kind=int64) :: klx(1), kly(1), klz(1)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 2
    in1(3) = 3
    in1(4) = 4

    kn(1) = 0

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 0

    jx(1) = 0

    ix(1) = 4
    ix(2) = 1
    ix(3) = 1
    ix(4) = 3
    ix(5) = 3
    ix(6) = 2
    ix(7) = 1
    ix(8) = 2
    ix(9) = 1
    ix(10) = 2

    ! y-arrays

    ly(1) = 0

    ky(1) = 0

    jy(1) = 0

    iy(1) = 1
    iy(2) = 4
    iy(3) = 1
    iy(4) = 2
    iy(5) = 1
    iy(6) = 3
    iy(7) = 3
    iy(8) = 1
    iy(9) = 2
    iy(10) = 2

    ! z-arrays

    lz(1) = 0

    kz(1) = 0

    jz(1) = 0

    iz(1) = 1
    iz(2) = 1
    iz(3) = 4
    iz(4) = 1
    iz(5) = 2
    iz(6) = 1
    iz(7) = 2
    iz(8) = 3
    iz(9) = 3
    iz(10) = 2

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 4
    ijx(2) = 1
    ijx(3) = 1
    ijx(4) = 3
    ijx(5) = 3
    ijx(6) = 2
    ijx(7) = 1
    ijx(8) = 2
    ijx(9) = 1
    ijx(10) = 2

    ijy(1) = 1
    ijy(2) = 4
    ijy(3) = 1
    ijy(4) = 2
    ijy(5) = 1
    ijy(6) = 3
    ijy(7) = 3
    ijy(8) = 1
    ijy(9) = 2
    ijy(10) = 2

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 4
    ijz(4) = 1
    ijz(5) = 2
    ijz(6) = 1
    ijz(7) = 2
    ijz(8) = 3
    ijz(9) = 3
    ijz(10) = 2

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 0

    kly(1) = 0

    klz(1) = 0

    allocate (n03bra(res%n_s_shl*res%n_f_shl))
    allocate (xint03bra(res%n_s_shl*res%n_f_shl))
    allocate (n00ket(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00ket(res%n_s_shl*(res%n_s_shl + 1)/2))

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

    if ((nsfbra*nssket) .le. nchunksize_int64) nchunksize_int64 = nsfbra*nssket
    ntile = int(nsfbra*nssket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nsfbra*nssket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nsfbra, xint03bra, n03bra, xint00ket, n00ket, sf_pair, ss_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d00ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d03bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nsfbra) + 1
              kl_tmp = (iquart - 1)/nsfbra + 1

              test = xint03bra(ij_tmp)*xint00ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n03bra(ij_tmp)
                kl = n00ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_s_shl(jsh_tmp)
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

                    if (xx .ge. 37.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.2752551286084109D+00
                      rts(2) = factr*0.2724744871391588D+01

                      wts(1) = factw*0.8049140900055122D+00
                      wts(2) = factw*0.8131283544724527D-01

                    else ! "regular" evaluation

                      rgrid(1) = 0.4935129360637439D-05
                      rgrid(2) = 0.1361431404118938D-03
                      rgrid(3) = 0.8129748816298266D-03
                      rgrid(4) = 0.2756670127399030D-02
                      rgrid(5) = 0.6935339478525350D-02
                      rgrid(6) = 0.1448902560832900D-01
                      rgrid(7) = 0.2663972894789983D-01
                      rgrid(8) = 0.4459215012310767D-01
                      rgrid(9) = 0.6943153026591924D-01
                      rgrid(10) = 0.1020252057162236D+00
                      rgrid(11) = 0.1429343223834524D+00
                      rgrid(12) = 0.1923415868672259D+00
                      rgrid(13) = 0.2499999999999999D+00
                      rgrid(14) = 0.3152062794779362D+00
                      rgrid(15) = 0.3868012061044409D+00
                      rgrid(16) = 0.4631975115256109D+00
                      rgrid(17) = 0.5424342617116348D+00
                      rgrid(18) = 0.6222550803643302D+00
                      rgrid(19) = 0.7002060974213685D+00
                      rgrid(20) = 0.7737482886456865D+00
                      rgrid(21) = 0.8403779682393594D+00
                      rgrid(22) = 0.8977486680056744D+00
                      rgrid(23) = 0.9437875461106039D+00
                      rgrid(24) = 0.9768000645999300D+00
                      rgrid(25) = 0.9955619049198584D+00

                      wgrid(1) = 0.5696899250513485D-02*exp(-xx*0.4935129360637439D-05)
                      wgrid(2) = 0.1317749330751579D-01*exp(-xx*0.1361431404118938D-03)
                      wgrid(3) = 0.2046957835065303D-01*exp(-xx*0.8129748816298266D-03)
                      wgrid(4) = 0.2745234798791782D-01*exp(-xx*0.2756670127399030D-02)
                      wgrid(5) = 0.3401916690617832D-01*exp(-xx*0.6935339478525350D-02)
                      wgrid(6) = 0.4007035016750040D-01*exp(-xx*0.1448902560832900D-01)
                      wgrid(7) = 0.4551413099148229D-01*exp(-xx*0.2663972894789983D-01)
                      wgrid(8) = 0.5026797453352529D-01*exp(-xx*0.4459215012310767D-01)
                      wgrid(9) = 0.5425981223713152D-01*exp(-xx*0.6943153026591924D-01)
                      wgrid(10) = 0.5742912957285592D-01*exp(-xx*0.1020252057162236D+00)
                      wgrid(11) = 0.5972788176789209D-01*exp(-xx*0.1429343223834524D+00)
                      wgrid(12) = 0.6112122149515455D-01*exp(-xx*0.1923415868672259D+00)
                      wgrid(13) = 0.6158802686335779D-01*exp(-xx*0.2499999999999999D+00)
                      wgrid(14) = 0.6112122149515536D-01*exp(-xx*0.3152062794779362D+00)
                      wgrid(15) = 0.5972788176789246D-01*exp(-xx*0.3868012061044409D+00)
                      wgrid(16) = 0.5742912957285590D-01*exp(-xx*0.4631975115256109D+00)
                      wgrid(17) = 0.5425981223713182D-01*exp(-xx*0.5424342617116348D+00)
                      wgrid(18) = 0.5026797453352500D-01*exp(-xx*0.6222550803643302D+00)
                      wgrid(19) = 0.4551413099148207D-01*exp(-xx*0.7002060974213685D+00)
                      wgrid(20) = 0.4007035016750037D-01*exp(-xx*0.7737482886456865D+00)
                      wgrid(21) = 0.3401916690617897D-01*exp(-xx*0.8403779682393594D+00)
                      wgrid(22) = 0.2745234798791809D-01*exp(-xx*0.8977486680056744D+00)
                      wgrid(23) = 0.2046957835065275D-01*exp(-xx*0.9437875461106039D+00)
                      wgrid(24) = 0.1317749330751527D-01*exp(-xx*0.9768000645999300D+00)
                      wgrid(25) = 0.5696899250513537D-02*exp(-xx*0.9955619049198584D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 25
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 25
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 1

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 25

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
                        wrk(2) = 0.0D+00
                        do 100 kk = 2, 2

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 2

                            jj = 0

105                         do 110 m = l, 2
                              if (m .eq. 2) go to 120
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

                                do 300 ii = 2, 2

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 2
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

                                    do 310 kk = 1, 2
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 2
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

                                      ! i2 = in(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(2) = xc00
                                      yin(2) = yc00
                                      zin(2) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    2

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    3
                                      ! i3 =    1
                                      ! i4 =    2

                                      xin(3) = c10*xin(1) + xc00*xin(2)
                                      yin(3) = c10*yin(1) + yc00*yin(2)
                                      zin(3) = c10*zin(1) + zc00*zin(2)

                                      ! i3 = i4 =    2
                                      ! i4 = i5 =    3

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    4
                                      ! i3 =    2
                                      ! i4 =    3

                                      xin(4) = c10*xin(2) + xc00*xin(3)
                                      yin(4) = c10*yin(2) + yc00*yin(3)
                                      zin(4) = c10*zin(2) + zc00*zin(3)

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    4

                                      ! n =    4

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =    4

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

                                      ! i1 = in(1) =    5

                                      xin(5) = 1.0_dp
                                      yin(5) = 1.0_dp
                                      zin(5) = f00

                                      ! i2 = in(2) =    6
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(6) = xc00
                                      yin(6) = yc00
                                      zin(6) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    5
                                      ! i4 = i2 =    6

                                      ! do n = 2,   3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    7
                                      ! i3 =    5
                                      ! i4 =    6

                                      xin(7) = c10*xin(5) + xc00*xin(6)
                                      yin(7) = c10*yin(5) + yc00*yin(6)
                                      zin(7) = c10*zin(5) + zc00*zin(6)

                                      ! i3 = i4 =    6
                                      ! i4 = i5 =    7

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    8
                                      ! i3 =    6
                                      ! i4 =    7

                                      xin(8) = c10*xin(6) + xc00*xin(7)
                                      yin(8) = c10*yin(6) + yc00*yin(7)
                                      zin(8) = c10*zin(6) + zc00*zin(7)

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    4

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =    8

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                     eri_value(1) = eri_value(1) + d03bra(1)*d00ket(1)*(xin(4)*yin(1)*zin(1) + xin(8)*yin(5)*zin(5))
                                     eri_value(2) = eri_value(2) + d03bra(2)*d00ket(1)*(xin(1)*yin(4)*zin(1) + xin(5)*yin(8)*zin(5))
                                     eri_value(3) = eri_value(3) + d03bra(3)*d00ket(1)*(xin(1)*yin(1)*zin(4) + xin(5)*yin(5)*zin(8))
                                     eri_value(4) = eri_value(4) + d03bra(4)*d00ket(1)*(xin(3)*yin(2)*zin(1) + xin(7)*yin(6)*zin(5))
                                     eri_value(5) = eri_value(5) + d03bra(5)*d00ket(1)*(xin(3)*yin(1)*zin(2) + xin(7)*yin(5)*zin(6))
                                     eri_value(6) = eri_value(6) + d03bra(6)*d00ket(1)*(xin(2)*yin(3)*zin(1) + xin(6)*yin(7)*zin(5))
                                     eri_value(7) = eri_value(7) + d03bra(7)*d00ket(1)*(xin(1)*yin(3)*zin(2) + xin(5)*yin(7)*zin(6))
                                     eri_value(8) = eri_value(8) + d03bra(8)*d00ket(1)*(xin(2)*yin(1)*zin(3) + xin(6)*yin(5)*zin(7))
                                     eri_value(9) = eri_value(9) + d03bra(9)*d00ket(1)*(xin(1)*yin(2)*zin(3) + xin(5)*yin(6)*zin(7))
                                  eri_value(10) = eri_value(10) + d03bra(10)*d00ket(1)*(xin(2)*yin(2)*zin(2) + xin(6)*yin(6)*zin(6))

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
                                    ip = (i - 1)*1 ! Stride between functions in i

                                    do j = 1, 1 ! # of cartesians in j

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

                              deallocate (n03bra)
                              deallocate (xint03bra)
                              deallocate (n00ket)
                              deallocate (xint00ket)

                              end subroutine int3000
                              end submodule
