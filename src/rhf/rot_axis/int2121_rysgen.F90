! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rot_axis_kernels) int2121_impl
contains
  module subroutine int2121(pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n12bra(:)
    real(dp), allocatable :: xint12bra(:)
    integer(kind=int64) :: npdbra
    real(dp) :: scutpdbra, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iijj, kkll
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2, nij, nkl, itmp
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
    real(dp) :: xin(144), yin(144), zin(144)
    real(dp) :: eri_value(324)
    real(dp) :: d12bra(18), d12ket(18)
    integer(kind=int64) :: ix(6), jx(3), kx(6), lx(3)
    integer(kind=int64) :: iy(6), jy(3), ky(6), ly(3)
    integer(kind=int64) :: iz(6), jz(3), kz(6), lz(3)
    integer(kind=int64) :: in(4), in1(4), kn(4)
    integer(kind=int64) :: ijx(18), ijy(18), ijz(18)
    integer(kind=int64) :: klx(18), kly(18), klz(18)
    logical :: same

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 31

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

    jx(1) = 6
    jx(2) = 0
    jx(3) = 0

    ix(1) = 25
    ix(2) = 1
    ix(3) = 1
    ix(4) = 13
    ix(5) = 13
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
    jy(2) = 6
    jy(3) = 0

    iy(1) = 1
    iy(2) = 25
    iy(3) = 1
    iy(4) = 13
    iy(5) = 1
    iy(6) = 13

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
    jz(3) = 6

    iz(1) = 1
    iz(2) = 1
    iz(3) = 25
    iz(4) = 1
    iz(5) = 13
    iz(6) = 13

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 31
    ijx(2) = 25
    ijx(3) = 25
    ijx(4) = 7
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 7
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 19
    ijx(11) = 13
    ijx(12) = 13
    ijx(13) = 19
    ijx(14) = 13
    ijx(15) = 13
    ijx(16) = 7
    ijx(17) = 1
    ijx(18) = 1

    ijy(1) = 1
    ijy(2) = 7
    ijy(3) = 1
    ijy(4) = 25
    ijy(5) = 31
    ijy(6) = 25
    ijy(7) = 1
    ijy(8) = 7
    ijy(9) = 1
    ijy(10) = 13
    ijy(11) = 19
    ijy(12) = 13
    ijy(13) = 1
    ijy(14) = 7
    ijy(15) = 1
    ijy(16) = 13
    ijy(17) = 19
    ijy(18) = 13

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 7
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 7
    ijz(7) = 25
    ijz(8) = 25
    ijz(9) = 31
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 7
    ijz(13) = 13
    ijz(14) = 13
    ijz(15) = 19
    ijz(16) = 13
    ijz(17) = 13
    ijz(18) = 19

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

    allocate (n12bra(res%n_p_shl*res%n_d_shl))
    allocate (xint12bra(res%n_p_shl*res%n_d_shl))

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

    ! --multi-gpu--work
    nchunk = npdbra/res%n_size
    nquart_start = nchunk*res%n_rank + 1
    nquart_end = nquart_start + nchunk - 1
    if (res%n_rank .EQ. res%n_size - 1) nquart_end = npdbra

    ! Mappings to GPU

 !$omp target teams distribute parallel do collapse(2) default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npdbra, xint12bra, n12bra, pd_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d12ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d12bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,nij,nkl,itmp,same)
              do iijj = nquart_start, nquart_end
                do kkll = 1, npdbra

                  if (kkll .gt. iijj) cycle

                  test = xint12bra(iijj)*xint12bra(kkll)
                  if (test .gt. cutoff_schwarz) then

                    ij = n12bra(iijj)
                    kl = n12bra(kkll)

                    ish_tmp = mod(ij - 1, res%n_d_shl) + 1
                    jsh_tmp = (ij - 1)/res%n_d_shl + 1
                    ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                    lsh_tmp = (kl - 1)/res%n_d_shl + 1

                    ish = res%i_d_shl(ish_tmp)
                    jsh = res%i_p_shl(jsh_tmp)
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

                                          ! do n = 2,   3

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

                                          ! i5 = in(n+1) =   31
                                          ! i3 =   13
                                          ! i4 =   25

                                          xin(31) = c10*xin(13) + xc00*xin(25)
                                          yin(31) = c10*yin(13) + yc00*yin(25)
                                          zin(31) = c10*zin(13) + zc00*zin(25)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   33
                                          ! i5 =   31
                                          ! i4 =   25

                                          xin(33) = xcp00*xin(31) + cp10*xin(25)
                                          yin(33) = ycp00*yin(31) + cp10*yin(25)
                                          zin(33) = zcp00*zin(31) + cp10*zin(25)

                                          ! ------------------

                                          ! i3 = i4 =   25
                                          ! i4 = i5 =   31

                                          ! n =    4

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

                                          ! i3 = i2+kn(n+1) =   17

                                          xin(17) = xc00*xin(5) + c01*xin(3)
                                          yin(17) = yc00*yin(5) + c01*yin(3)
                                          zin(17) = zc00*zin(5) + c01*zin(3)

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

                                          ! i3 = i2+kn(n+1) =   18

                                          xin(18) = xc00*xin(6) + c01*xin(5)
                                          yin(18) = yc00*yin(6) + c01*yin(5)
                                          zin(18) = zc00*zin(6) + c01*zin(5)

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
                                          ! i4 = i2 =   13

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   25

                                          xin(29) = c10*xin(5) + xc00*xin(17) + c01*xin(15)
                                          yin(29) = c10*yin(5) + yc00*yin(17) + c01*yin(15)
                                          zin(29) = c10*zin(5) + zc00*zin(17) + c01*zin(15)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   13
                                          ! i4 = i5 =   25

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   31

                                          xin(35) = c10*xin(17) + xc00*xin(29) + c01*xin(27)
                                          yin(35) = c10*yin(17) + yc00*yin(29) + c01*yin(27)
                                          zin(35) = c10*zin(17) + zc00*zin(29) + c01*zin(27)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   25
                                          ! i4 = i5 =   31

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   4

                                          ! n =    3

                                          ! k4 = kn(n+1) =    5
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   13

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   25

                                          xin(30) = c10*xin(6) + xc00*xin(18) + c01*xin(17)
                                          yin(30) = c10*yin(6) + yc00*yin(18) + c01*yin(17)
                                          zin(30) = c10*zin(6) + zc00*zin(18) + c01*zin(17)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   13
                                          ! i4 = i5 =   25

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   31

                                          xin(36) = c10*xin(18) + xc00*xin(30) + c01*xin(29)
                                          yin(36) = c10*yin(18) + yc00*yin(30) + c01*yin(29)
                                          zin(36) = c10*zin(18) + zc00*zin(30) + c01*zin(29)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   25
                                          ! i4 = i5 =   31

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   5

                                          ! n =    4

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =   31

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   31

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   25

                                          xin(31) = xin(31) + dxij*xin(25)
                                          yin(31) = yin(31) + dyij*yin(25)
                                          zin(31) = zin(31) + dzij*zin(25)

                                          ! i3 = i4 =   25
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                          ! do nj = 1,    1

                                          ! i4 = i3 =    7

                                          ! do ni = 1,    2

                                          xin(7) = xin(13) + dxij*xin(1)
                                          yin(7) = yin(13) + dyij*yin(1)
                                          zin(7) = zin(13) + dzij*zin(1)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   19

                                          ! ni =    2

                                          xin(19) = xin(25) + dxij*xin(13)
                                          yin(19) = yin(25) + dyij*yin(13)
                                          zin(19) = zin(25) + dzij*zin(13)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   13

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    2

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   33

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   27

                                          xin(33) = xin(33) + dxij*xin(27)
                                          yin(33) = yin(33) + dyij*yin(27)
                                          zin(33) = zin(33) + dzij*zin(27)

                                          ! i3 = i4 =   27
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                          ! do nj = 1,    1

                                          ! i4 = i3 =    9

                                          ! do ni = 1,    2

                                          xin(9) = xin(15) + dxij*xin(3)
                                          yin(9) = yin(15) + dyij*yin(3)
                                          zin(9) = zin(15) + dzij*zin(3)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                          ! ni =    2

                                          xin(21) = xin(27) + dxij*xin(15)
                                          yin(21) = yin(27) + dyij*yin(15)
                                          zin(21) = zin(27) + dzij*zin(15)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   15

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    4

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   35

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   29

                                          xin(35) = xin(35) + dxij*xin(29)
                                          yin(35) = yin(35) + dyij*yin(29)
                                          zin(35) = zin(35) + dzij*zin(29)

                                          ! i3 = i4 =   29
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   11

                                          ! do ni = 1,    2

                                          xin(11) = xin(17) + dxij*xin(5)
                                          yin(11) = yin(17) + dyij*yin(5)
                                          zin(11) = zin(17) + dzij*zin(5)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                          ! ni =    2

                                          xin(23) = xin(29) + dxij*xin(17)
                                          yin(23) = yin(29) + dyij*yin(17)
                                          zin(23) = zin(29) + dzij*zin(17)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   17

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    5

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   36

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   30

                                          xin(36) = xin(36) + dxij*xin(30)
                                          yin(36) = yin(36) + dyij*yin(30)
                                          zin(36) = zin(36) + dzij*zin(30)

                                          ! i3 = i4 =   30
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   12

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   12

                                          ! do ni = 1,    2

                                          xin(12) = xin(18) + dxij*xin(6)
                                          yin(12) = yin(18) + dyij*yin(6)
                                          zin(12) = zin(18) + dzij*zin(6)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                          ! ni =    2

                                          xin(24) = xin(30) + dxij*xin(18)
                                          yin(24) = yin(30) + dyij*yin(18)
                                          zin(24) = zin(30) + dzij*zin(18)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   18

                                          ! nj =    2

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

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   13

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   25

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   37

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
                                          ! k2 = kn(2) =    2
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(49) = xc00
                                          yin(49) = yc00
                                          zin(49) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =   39

                                          xin(39) = xcp00
                                          yin(39) = ycp00
                                          zin(39) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   51
                                          ! i2 =   49

                                          xin(51) = xcp00*xin(49) + cp10
                                          yin(51) = ycp00*yin(49) + cp10
                                          zin(51) = zcp00*zin(49) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =   37
                                          ! i4 = i2 =   49

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   61
                                          ! i3 =   37
                                          ! i4 =   49

                                          xin(61) = c10*xin(37) + xc00*xin(49)
                                          yin(61) = c10*yin(37) + yc00*yin(49)
                                          zin(61) = c10*zin(37) + zc00*zin(49)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   63
                                          ! i5 =   61
                                          ! i4 =   49

                                          xin(63) = xcp00*xin(61) + cp10*xin(49)
                                          yin(63) = ycp00*yin(61) + cp10*yin(49)
                                          zin(63) = zcp00*zin(61) + cp10*zin(49)

                                          ! ------------------

                                          ! i3 = i4 =   49
                                          ! i4 = i5 =   61

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   67
                                          ! i3 =   49
                                          ! i4 =   61

                                          xin(67) = c10*xin(49) + xc00*xin(61)
                                          yin(67) = c10*yin(49) + yc00*yin(61)
                                          zin(67) = c10*zin(49) + zc00*zin(61)

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

                                          ! n =    4

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =   37
                                          ! i4 = i1+k2 =   39

                                          ! do n = 2,    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   41
                                          ! i3 =   37
                                          ! i4 =   39

                                          xin(41) = cp01*xin(37) + xcp00*xin(39)
                                          yin(41) = cp01*yin(37) + ycp00*yin(39)
                                          zin(41) = cp01*zin(37) + zcp00*zin(39)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   53

                                          xin(53) = xc00*xin(41) + c01*xin(39)
                                          yin(53) = yc00*yin(41) + c01*yin(39)
                                          zin(53) = zc00*zin(41) + c01*zin(39)

                                          ! ------------------

                                          ! i3 = i4 =   39
                                          ! i4 = i5 =   41

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   42
                                          ! i3 =   39
                                          ! i4 =   41

                                          xin(42) = cp01*xin(39) + xcp00*xin(41)
                                          yin(42) = cp01*yin(39) + ycp00*yin(41)
                                          zin(42) = cp01*zin(39) + zcp00*zin(41)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   54

                                          xin(54) = xc00*xin(42) + c01*xin(41)
                                          yin(54) = yc00*yin(42) + c01*yin(41)
                                          zin(54) = zc00*zin(42) + c01*zin(41)

                                          ! ------------------

                                          ! i3 = i4 =   41
                                          ! i4 = i5 =   42

                                          ! n =    4

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    2

                                          ! do n = 2,    3

                                          ! k4 = kn(n+1) =    4
                                          ! i3 = i1 =   37
                                          ! i4 = i2 =   49

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   61

                                          xin(65) = c10*xin(41) + xc00*xin(53) + c01*xin(51)
                                          yin(65) = c10*yin(41) + yc00*yin(53) + c01*yin(51)
                                          zin(65) = c10*zin(41) + zc00*zin(53) + c01*zin(51)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   49
                                          ! i4 = i5 =   61

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   67

                                          xin(71) = c10*xin(53) + xc00*xin(65) + c01*xin(63)
                                          yin(71) = c10*yin(53) + yc00*yin(65) + c01*yin(63)
                                          zin(71) = c10*zin(53) + zc00*zin(65) + c01*zin(63)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   61
                                          ! i4 = i5 =   67

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   4

                                          ! n =    3

                                          ! k4 = kn(n+1) =    5
                                          ! i3 = i1 =   37
                                          ! i4 = i2 =   49

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   61

                                          xin(66) = c10*xin(42) + xc00*xin(54) + c01*xin(53)
                                          yin(66) = c10*yin(42) + yc00*yin(54) + c01*yin(53)
                                          zin(66) = c10*zin(42) + zc00*zin(54) + c01*zin(53)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   49
                                          ! i4 = i5 =   61

                                          ! nn =    3

                                          ! i5 = in(nn+1) =   67

                                          xin(72) = c10*xin(54) + xc00*xin(66) + c01*xin(65)
                                          yin(72) = c10*yin(54) + yc00*yin(66) + c01*yin(65)
                                          zin(72) = c10*zin(54) + zc00*zin(66) + c01*zin(65)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   61
                                          ! i4 = i5 =   67

                                          ! nn =    4

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

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   67

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   61

                                          xin(67) = xin(67) + dxij*xin(61)
                                          yin(67) = yin(67) + dyij*yin(61)
                                          zin(67) = zin(67) + dzij*zin(61)

                                          ! i3 = i4 =   61
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   43

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   43

                                          ! do ni = 1,    2

                                          xin(43) = xin(49) + dxij*xin(37)
                                          yin(43) = yin(49) + dyij*yin(37)
                                          zin(43) = zin(49) + dzij*zin(37)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                          ! ni =    2

                                          xin(55) = xin(61) + dxij*xin(49)
                                          yin(55) = yin(61) + dyij*yin(49)
                                          zin(55) = zin(61) + dzij*zin(49)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   49

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    2

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   69

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   63

                                          xin(69) = xin(69) + dxij*xin(63)
                                          yin(69) = yin(69) + dyij*yin(63)
                                          zin(69) = zin(69) + dzij*zin(63)

                                          ! i3 = i4 =   63
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   45

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   45

                                          ! do ni = 1,    2

                                          xin(45) = xin(51) + dxij*xin(39)
                                          yin(45) = yin(51) + dyij*yin(39)
                                          zin(45) = zin(51) + dzij*zin(39)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                          ! ni =    2

                                          xin(57) = xin(63) + dxij*xin(51)
                                          yin(57) = yin(63) + dyij*yin(51)
                                          zin(57) = zin(63) + dzij*zin(51)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   51

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    4

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   71

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   65

                                          xin(71) = xin(71) + dxij*xin(65)
                                          yin(71) = yin(71) + dyij*yin(65)
                                          zin(71) = zin(71) + dzij*zin(65)

                                          ! i3 = i4 =   65
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   47

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   47

                                          ! do ni = 1,    2

                                          xin(47) = xin(53) + dxij*xin(41)
                                          yin(47) = yin(53) + dyij*yin(41)
                                          zin(47) = zin(53) + dzij*zin(41)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                          ! ni =    2

                                          xin(59) = xin(65) + dxij*xin(53)
                                          yin(59) = yin(65) + dyij*yin(53)
                                          zin(59) = zin(65) + dzij*zin(53)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   53

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    5

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =   72

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   66

                                          xin(72) = xin(72) + dxij*xin(66)
                                          yin(72) = yin(72) + dyij*yin(66)
                                          zin(72) = zin(72) + dzij*zin(66)

                                          ! i3 = i4 =   66
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   48

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   48

                                          ! do ni = 1,    2

                                          xin(48) = xin(54) + dxij*xin(42)
                                          yin(48) = yin(54) + dyij*yin(42)
                                          zin(48) = zin(54) + dzij*zin(42)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                          ! ni =    2

                                          xin(60) = xin(66) + dxij*xin(54)
                                          yin(60) = yin(66) + dyij*yin(54)
                                          zin(60) = zin(66) + dzij*zin(54)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   54

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =    5

                                          ! iaa = i1 =   37

                                          ! ni = 0

                                          ! do while ni.le.iang

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

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   61

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

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
                                          ! k2 = kn(2) =    2
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(85) = xc00
                                          yin(85) = yc00
                                          zin(85) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =   75

                                          xin(75) = xcp00
                                          yin(75) = ycp00
                                          zin(75) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   87
                                          ! i2 =   85

                                          xin(87) = xcp00*xin(85) + cp10
                                          yin(87) = ycp00*yin(85) + cp10
                                          zin(87) = zcp00*zin(85) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =   73
                                          ! i4 = i2 =   85

                                          ! do n = 2,   3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   97
                                          ! i3 =   73
                                          ! i4 =   85

                                          xin(97) = c10*xin(73) + xc00*xin(85)
                                          yin(97) = c10*yin(73) + yc00*yin(85)
                                          zin(97) = c10*zin(73) + zc00*zin(85)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   99
                                          ! i5 =   97
                                          ! i4 =   85

                                          xin(99) = xcp00*xin(97) + cp10*xin(85)
                                          yin(99) = ycp00*yin(97) + cp10*yin(85)
                                          zin(99) = zcp00*zin(97) + cp10*zin(85)

                                          ! ------------------

                                          ! i3 = i4 =   85
                                          ! i4 = i5 =   97

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  103
                                          ! i3 =   85
                                          ! i4 =   97

                                          xin(103) = c10*xin(85) + xc00*xin(97)
                                          yin(103) = c10*yin(85) + yc00*yin(97)
                                          zin(103) = c10*zin(85) + zc00*zin(97)

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

                                          ! n =    4

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

                                          ! i3 = i2+kn(n+1) =   89

                                          xin(89) = xc00*xin(77) + c01*xin(75)
                                          yin(89) = yc00*yin(77) + c01*yin(75)
                                          zin(89) = zc00*zin(77) + c01*zin(75)

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

                                          ! i3 = i2+kn(n+1) =   90

                                          xin(90) = xc00*xin(78) + c01*xin(77)
                                          yin(90) = yc00*yin(78) + c01*yin(77)
                                          zin(90) = zc00*zin(78) + c01*zin(77)

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
                                          ! i4 = i2 =   85

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   97

                                          xin(101) = c10*xin(77) + xc00*xin(89) + c01*xin(87)
                                          yin(101) = c10*yin(77) + yc00*yin(89) + c01*yin(87)
                                          zin(101) = c10*zin(77) + zc00*zin(89) + c01*zin(87)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   85
                                          ! i4 = i5 =   97

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  103

                                          xin(107) = c10*xin(89) + xc00*xin(101) + c01*xin(99)
                                          yin(107) = c10*yin(89) + yc00*yin(101) + c01*yin(99)
                                          zin(107) = c10*zin(89) + zc00*zin(101) + c01*zin(99)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   97
                                          ! i4 = i5 =  103

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   4

                                          ! n =    3

                                          ! k4 = kn(n+1) =    5
                                          ! i3 = i1 =   73
                                          ! i4 = i2 =   85

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =   97

                                          xin(102) = c10*xin(78) + xc00*xin(90) + c01*xin(89)
                                          yin(102) = c10*yin(78) + yc00*yin(90) + c01*yin(89)
                                          zin(102) = c10*zin(78) + zc00*zin(90) + c01*zin(89)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   85
                                          ! i4 = i5 =   97

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  103

                                          xin(108) = c10*xin(90) + xc00*xin(102) + c01*xin(101)
                                          yin(108) = c10*yin(90) + yc00*yin(102) + c01*yin(101)
                                          zin(108) = c10*zin(90) + zc00*zin(102) + c01*zin(101)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   97
                                          ! i4 = i5 =  103

                                          ! nn =    4

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

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  103

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   97

                                          xin(103) = xin(103) + dxij*xin(97)
                                          yin(103) = yin(103) + dyij*yin(97)
                                          zin(103) = zin(103) + dzij*zin(97)

                                          ! i3 = i4 =   97
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   79

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   79

                                          ! do ni = 1,    2

                                          xin(79) = xin(85) + dxij*xin(73)
                                          yin(79) = yin(85) + dyij*yin(73)
                                          zin(79) = zin(85) + dzij*zin(73)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                          ! ni =    2

                                          xin(91) = xin(97) + dxij*xin(85)
                                          yin(91) = yin(97) + dyij*yin(85)
                                          zin(91) = zin(97) + dzij*zin(85)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   85

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    2

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  105

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =   99

                                          xin(105) = xin(105) + dxij*xin(99)
                                          yin(105) = yin(105) + dyij*yin(99)
                                          zin(105) = zin(105) + dzij*zin(99)

                                          ! i3 = i4 =   99
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   81

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   81

                                          ! do ni = 1,    2

                                          xin(81) = xin(87) + dxij*xin(75)
                                          yin(81) = yin(87) + dyij*yin(75)
                                          zin(81) = zin(87) + dzij*zin(75)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                          ! ni =    2

                                          xin(93) = xin(99) + dxij*xin(87)
                                          yin(93) = yin(99) + dyij*yin(87)
                                          zin(93) = zin(99) + dzij*zin(87)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   87

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    4

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  107

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  101

                                          xin(107) = xin(107) + dxij*xin(101)
                                          yin(107) = yin(107) + dyij*yin(101)
                                          zin(107) = zin(107) + dzij*zin(101)

                                          ! i3 = i4 =  101
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   83

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   83

                                          ! do ni = 1,    2

                                          xin(83) = xin(89) + dxij*xin(77)
                                          yin(83) = yin(89) + dyij*yin(77)
                                          zin(83) = zin(89) + dzij*zin(77)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                          ! ni =    2

                                          xin(95) = xin(101) + dxij*xin(89)
                                          yin(95) = yin(101) + dyij*yin(89)
                                          zin(95) = zin(101) + dzij*zin(89)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   89

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    5

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  108

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  102

                                          xin(108) = xin(108) + dxij*xin(102)
                                          yin(108) = yin(108) + dyij*yin(102)
                                          zin(108) = zin(108) + dzij*zin(102)

                                          ! i3 = i4 =  102
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   84

                                          ! do nj = 1,    1

                                          ! i4 = i3 =   84

                                          ! do ni = 1,    2

                                          xin(84) = xin(90) + dxij*xin(78)
                                          yin(84) = yin(90) + dyij*yin(78)
                                          zin(84) = zin(90) + dzij*zin(78)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                          ! ni =    2

                                          xin(96) = xin(102) + dxij*xin(90)
                                          yin(96) = yin(102) + dyij*yin(90)
                                          zin(96) = zin(102) + dzij*zin(90)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   90

                                          ! nj =    2

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

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   85

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

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
                                          ! k2 = kn(2) =    2
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(121) = xc00
                                          yin(121) = yc00
                                          zin(121) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  111

                                          xin(111) = xcp00
                                          yin(111) = ycp00
                                          zin(111) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  123
                                          ! i2 =  121

                                          xin(123) = xcp00*xin(121) + cp10
                                          yin(123) = ycp00*yin(121) + cp10
                                          zin(123) = zcp00*zin(121) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  109
                                          ! i4 = i2 =  121

                                          ! do n = 2,   3

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

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  139
                                          ! i3 =  121
                                          ! i4 =  133

                                          xin(139) = c10*xin(121) + xc00*xin(133)
                                          yin(139) = c10*yin(121) + yc00*yin(133)
                                          zin(139) = c10*zin(121) + zc00*zin(133)

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

                                          ! n =    4

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

                                          ! i3 = i2+kn(n+1) =  125

                                          xin(125) = xc00*xin(113) + c01*xin(111)
                                          yin(125) = yc00*yin(113) + c01*yin(111)
                                          zin(125) = zc00*zin(113) + c01*zin(111)

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

                                          ! i3 = i2+kn(n+1) =  126

                                          xin(126) = xc00*xin(114) + c01*xin(113)
                                          yin(126) = yc00*yin(114) + c01*yin(113)
                                          zin(126) = zc00*zin(114) + c01*zin(113)

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
                                          ! i4 = i2 =  121

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =  133

                                          xin(137) = c10*xin(113) + xc00*xin(125) + c01*xin(123)
                                          yin(137) = c10*yin(113) + yc00*yin(125) + c01*yin(123)
                                          zin(137) = c10*zin(113) + zc00*zin(125) + c01*zin(123)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  139

                                          xin(143) = c10*xin(125) + xc00*xin(137) + c01*xin(135)
                                          yin(143) = c10*yin(125) + yc00*yin(137) + c01*yin(135)
                                          zin(143) = c10*zin(125) + zc00*zin(137) + c01*zin(135)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  133
                                          ! i4 = i5 =  139

                                          ! nn =    4

                                          ! end do

                                          ! k3 = k4   4

                                          ! n =    3

                                          ! k4 = kn(n+1) =    5
                                          ! i3 = i1 =  109
                                          ! i4 = i2 =  121

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    3

                                          ! i5 = in(nn+1) =  133

                                          xin(138) = c10*xin(114) + xc00*xin(126) + c01*xin(125)
                                          yin(138) = c10*yin(114) + yc00*yin(126) + c01*yin(125)
                                          zin(138) = c10*zin(114) + zc00*zin(126) + c01*zin(125)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  139

                                          xin(144) = c10*xin(126) + xc00*xin(138) + c01*xin(137)
                                          yin(144) = c10*yin(126) + yc00*yin(138) + c01*yin(137)
                                          zin(144) = c10*zin(126) + zc00*zin(138) + c01*zin(137)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  133
                                          ! i4 = i5 =  139

                                          ! nn =    4

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

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  139

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  133

                                          xin(139) = xin(139) + dxij*xin(133)
                                          yin(139) = yin(139) + dyij*yin(133)
                                          zin(139) = zin(139) + dzij*zin(133)

                                          ! i3 = i4 =  133
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  115

                                          ! do nj = 1,    1

                                          ! i4 = i3 =  115

                                          ! do ni = 1,    2

                                          xin(115) = xin(121) + dxij*xin(109)
                                          yin(115) = yin(121) + dyij*yin(109)
                                          zin(115) = zin(121) + dzij*zin(109)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                          ! ni =    2

                                          xin(127) = xin(133) + dxij*xin(121)
                                          yin(127) = yin(133) + dyij*yin(121)
                                          zin(127) = zin(133) + dzij*zin(121)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  121

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    2

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  141

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  135

                                          xin(141) = xin(141) + dxij*xin(135)
                                          yin(141) = yin(141) + dyij*yin(135)
                                          zin(141) = zin(141) + dzij*zin(135)

                                          ! i3 = i4 =  135
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  117

                                          ! do nj = 1,    1

                                          ! i4 = i3 =  117

                                          ! do ni = 1,    2

                                          xin(117) = xin(123) + dxij*xin(111)
                                          yin(117) = yin(123) + dyij*yin(111)
                                          zin(117) = zin(123) + dzij*zin(111)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                          ! ni =    2

                                          xin(129) = xin(135) + dxij*xin(123)
                                          yin(129) = yin(135) + dyij*yin(123)
                                          zin(129) = zin(135) + dzij*zin(123)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  123

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    4

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  143

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  137

                                          xin(143) = xin(143) + dxij*xin(137)
                                          yin(143) = yin(143) + dyij*yin(137)
                                          zin(143) = zin(143) + dzij*zin(137)

                                          ! i3 = i4 =  137
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  119

                                          ! do nj = 1,    1

                                          ! i4 = i3 =  119

                                          ! do ni = 1,    2

                                          xin(119) = xin(125) + dxij*xin(113)
                                          yin(119) = yin(125) + dyij*yin(113)
                                          zin(119) = zin(125) + dzij*zin(113)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                          ! ni =    2

                                          xin(131) = xin(137) + dxij*xin(125)
                                          yin(131) = yin(137) + dyij*yin(125)
                                          zin(131) = zin(137) + dzij*zin(125)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  125

                                          ! nj =    2

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    5

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    3

                                          ! i3 = i5 + km =  144

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  138

                                          xin(144) = xin(144) + dxij*xin(138)
                                          yin(144) = yin(144) + dyij*yin(138)
                                          zin(144) = zin(144) + dzij*zin(138)

                                          ! i3 = i4 =  138
                                          ! nn = nn-1 =    2

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  120

                                          ! do nj = 1,    1

                                          ! i4 = i3 =  120

                                          ! do ni = 1,    2

                                          xin(120) = xin(126) + dxij*xin(114)
                                          yin(120) = yin(126) + dyij*yin(114)
                                          zin(120) = zin(126) + dzij*zin(114)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                          ! ni =    2

                                          xin(132) = xin(138) + dxij*xin(126)
                                          yin(132) = yin(138) + dyij*yin(126)
                                          zin(132) = zin(138) + dzij*zin(126)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                          ! ni =    3

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  126

                                          ! nj =    2

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

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  121

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  133

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

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

                                          ! nj = nj + 1 =    1

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

                                          ! nj = nj + 1 =    2

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                          ! end do

                                          ! *** Now root =    5

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

                                          !                     --- END XYZINT ---

                                          !                       --- FORMS ---
                                          ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d12bra(  1)*d12ket(  1)*(xin(  36)*yin(   1)*zin(   1)+xin(  72)*yin(  37)*zin(  37)+xin( 108)*yin(  73)*zin(  73)+xin( 144)*yin( 109)*zin( 109))
          eri_value(    2)=eri_value(    2)+d12bra(  1)*d12ket(  2)*(xin(  35)*yin(   2)*zin(   1)+xin(  71)*yin(  38)*zin(  37)+xin( 107)*yin(  74)*zin(  73)+xin( 143)*yin( 110)*zin( 109))
          eri_value(    3)=eri_value(    3)+d12bra(  1)*d12ket(  3)*(xin(  35)*yin(   1)*zin(   2)+xin(  71)*yin(  37)*zin(  38)+xin( 107)*yin(  73)*zin(  74)+xin( 143)*yin( 109)*zin( 110))
          eri_value(    4)=eri_value(    4)+d12bra(  1)*d12ket(  4)*(xin(  32)*yin(   5)*zin(   1)+xin(  68)*yin(  41)*zin(  37)+xin( 104)*yin(  77)*zin(  73)+xin( 140)*yin( 113)*zin( 109))
          eri_value(    5)=eri_value(    5)+d12bra(  1)*d12ket(  5)*(xin(  31)*yin(   6)*zin(   1)+xin(  67)*yin(  42)*zin(  37)+xin( 103)*yin(  78)*zin(  73)+xin( 139)*yin( 114)*zin( 109))
          eri_value(    6)=eri_value(    6)+d12bra(  1)*d12ket(  6)*(xin(  31)*yin(   5)*zin(   2)+xin(  67)*yin(  41)*zin(  38)+xin( 103)*yin(  77)*zin(  74)+xin( 139)*yin( 113)*zin( 110))
          eri_value(    7)=eri_value(    7)+d12bra(  1)*d12ket(  7)*(xin(  32)*yin(   1)*zin(   5)+xin(  68)*yin(  37)*zin(  41)+xin( 104)*yin(  73)*zin(  77)+xin( 140)*yin( 109)*zin( 113))
          eri_value(    8)=eri_value(    8)+d12bra(  1)*d12ket(  8)*(xin(  31)*yin(   2)*zin(   5)+xin(  67)*yin(  38)*zin(  41)+xin( 103)*yin(  74)*zin(  77)+xin( 139)*yin( 110)*zin( 113))
          eri_value(    9)=eri_value(    9)+d12bra(  1)*d12ket(  9)*(xin(  31)*yin(   1)*zin(   6)+xin(  67)*yin(  37)*zin(  42)+xin( 103)*yin(  73)*zin(  78)+xin( 139)*yin( 109)*zin( 114))
          eri_value(   10)=eri_value(   10)+d12bra(  1)*d12ket( 10)*(xin(  34)*yin(   3)*zin(   1)+xin(  70)*yin(  39)*zin(  37)+xin( 106)*yin(  75)*zin(  73)+xin( 142)*yin( 111)*zin( 109))
          eri_value(   11)=eri_value(   11)+d12bra(  1)*d12ket( 11)*(xin(  33)*yin(   4)*zin(   1)+xin(  69)*yin(  40)*zin(  37)+xin( 105)*yin(  76)*zin(  73)+xin( 141)*yin( 112)*zin( 109))
          eri_value(   12)=eri_value(   12)+d12bra(  1)*d12ket( 12)*(xin(  33)*yin(   3)*zin(   2)+xin(  69)*yin(  39)*zin(  38)+xin( 105)*yin(  75)*zin(  74)+xin( 141)*yin( 111)*zin( 110))
          eri_value(   13)=eri_value(   13)+d12bra(  1)*d12ket( 13)*(xin(  34)*yin(   1)*zin(   3)+xin(  70)*yin(  37)*zin(  39)+xin( 106)*yin(  73)*zin(  75)+xin( 142)*yin( 109)*zin( 111))
          eri_value(   14)=eri_value(   14)+d12bra(  1)*d12ket( 14)*(xin(  33)*yin(   2)*zin(   3)+xin(  69)*yin(  38)*zin(  39)+xin( 105)*yin(  74)*zin(  75)+xin( 141)*yin( 110)*zin( 111))
          eri_value(   15)=eri_value(   15)+d12bra(  1)*d12ket( 15)*(xin(  33)*yin(   1)*zin(   4)+xin(  69)*yin(  37)*zin(  40)+xin( 105)*yin(  73)*zin(  76)+xin( 141)*yin( 109)*zin( 112))
          eri_value(   16)=eri_value(   16)+d12bra(  1)*d12ket( 16)*(xin(  32)*yin(   3)*zin(   3)+xin(  68)*yin(  39)*zin(  39)+xin( 104)*yin(  75)*zin(  75)+xin( 140)*yin( 111)*zin( 111))
          eri_value(   17)=eri_value(   17)+d12bra(  1)*d12ket( 17)*(xin(  31)*yin(   4)*zin(   3)+xin(  67)*yin(  40)*zin(  39)+xin( 103)*yin(  76)*zin(  75)+xin( 139)*yin( 112)*zin( 111))
          eri_value(   18)=eri_value(   18)+d12bra(  1)*d12ket( 18)*(xin(  31)*yin(   3)*zin(   4)+xin(  67)*yin(  39)*zin(  40)+xin( 103)*yin(  75)*zin(  76)+xin( 139)*yin( 111)*zin( 112))
          eri_value(   19)=eri_value(   19)+d12bra(  2)*d12ket(  1)*(xin(  30)*yin(   7)*zin(   1)+xin(  66)*yin(  43)*zin(  37)+xin( 102)*yin(  79)*zin(  73)+xin( 138)*yin( 115)*zin( 109))
          eri_value(   20)=eri_value(   20)+d12bra(  2)*d12ket(  2)*(xin(  29)*yin(   8)*zin(   1)+xin(  65)*yin(  44)*zin(  37)+xin( 101)*yin(  80)*zin(  73)+xin( 137)*yin( 116)*zin( 109))
          eri_value(   21)=eri_value(   21)+d12bra(  2)*d12ket(  3)*(xin(  29)*yin(   7)*zin(   2)+xin(  65)*yin(  43)*zin(  38)+xin( 101)*yin(  79)*zin(  74)+xin( 137)*yin( 115)*zin( 110))
          eri_value(   22)=eri_value(   22)+d12bra(  2)*d12ket(  4)*(xin(  26)*yin(  11)*zin(   1)+xin(  62)*yin(  47)*zin(  37)+xin(  98)*yin(  83)*zin(  73)+xin( 134)*yin( 119)*zin( 109))
          eri_value(   23)=eri_value(   23)+d12bra(  2)*d12ket(  5)*(xin(  25)*yin(  12)*zin(   1)+xin(  61)*yin(  48)*zin(  37)+xin(  97)*yin(  84)*zin(  73)+xin( 133)*yin( 120)*zin( 109))
          eri_value(   24)=eri_value(   24)+d12bra(  2)*d12ket(  6)*(xin(  25)*yin(  11)*zin(   2)+xin(  61)*yin(  47)*zin(  38)+xin(  97)*yin(  83)*zin(  74)+xin( 133)*yin( 119)*zin( 110))
          eri_value(   25)=eri_value(   25)+d12bra(  2)*d12ket(  7)*(xin(  26)*yin(   7)*zin(   5)+xin(  62)*yin(  43)*zin(  41)+xin(  98)*yin(  79)*zin(  77)+xin( 134)*yin( 115)*zin( 113))
          eri_value(   26)=eri_value(   26)+d12bra(  2)*d12ket(  8)*(xin(  25)*yin(   8)*zin(   5)+xin(  61)*yin(  44)*zin(  41)+xin(  97)*yin(  80)*zin(  77)+xin( 133)*yin( 116)*zin( 113))
          eri_value(   27)=eri_value(   27)+d12bra(  2)*d12ket(  9)*(xin(  25)*yin(   7)*zin(   6)+xin(  61)*yin(  43)*zin(  42)+xin(  97)*yin(  79)*zin(  78)+xin( 133)*yin( 115)*zin( 114))
          eri_value(   28)=eri_value(   28)+d12bra(  2)*d12ket( 10)*(xin(  28)*yin(   9)*zin(   1)+xin(  64)*yin(  45)*zin(  37)+xin( 100)*yin(  81)*zin(  73)+xin( 136)*yin( 117)*zin( 109))
          eri_value(   29)=eri_value(   29)+d12bra(  2)*d12ket( 11)*(xin(  27)*yin(  10)*zin(   1)+xin(  63)*yin(  46)*zin(  37)+xin(  99)*yin(  82)*zin(  73)+xin( 135)*yin( 118)*zin( 109))
          eri_value(   30)=eri_value(   30)+d12bra(  2)*d12ket( 12)*(xin(  27)*yin(   9)*zin(   2)+xin(  63)*yin(  45)*zin(  38)+xin(  99)*yin(  81)*zin(  74)+xin( 135)*yin( 117)*zin( 110))
          eri_value(   31)=eri_value(   31)+d12bra(  2)*d12ket( 13)*(xin(  28)*yin(   7)*zin(   3)+xin(  64)*yin(  43)*zin(  39)+xin( 100)*yin(  79)*zin(  75)+xin( 136)*yin( 115)*zin( 111))
          eri_value(   32)=eri_value(   32)+d12bra(  2)*d12ket( 14)*(xin(  27)*yin(   8)*zin(   3)+xin(  63)*yin(  44)*zin(  39)+xin(  99)*yin(  80)*zin(  75)+xin( 135)*yin( 116)*zin( 111))
          eri_value(   33)=eri_value(   33)+d12bra(  2)*d12ket( 15)*(xin(  27)*yin(   7)*zin(   4)+xin(  63)*yin(  43)*zin(  40)+xin(  99)*yin(  79)*zin(  76)+xin( 135)*yin( 115)*zin( 112))
          eri_value(   34)=eri_value(   34)+d12bra(  2)*d12ket( 16)*(xin(  26)*yin(   9)*zin(   3)+xin(  62)*yin(  45)*zin(  39)+xin(  98)*yin(  81)*zin(  75)+xin( 134)*yin( 117)*zin( 111))
          eri_value(   35)=eri_value(   35)+d12bra(  2)*d12ket( 17)*(xin(  25)*yin(  10)*zin(   3)+xin(  61)*yin(  46)*zin(  39)+xin(  97)*yin(  82)*zin(  75)+xin( 133)*yin( 118)*zin( 111))
          eri_value(   36)=eri_value(   36)+d12bra(  2)*d12ket( 18)*(xin(  25)*yin(   9)*zin(   4)+xin(  61)*yin(  45)*zin(  40)+xin(  97)*yin(  81)*zin(  76)+xin( 133)*yin( 117)*zin( 112))
          eri_value(   37)=eri_value(   37)+d12bra(  3)*d12ket(  1)*(xin(  30)*yin(   1)*zin(   7)+xin(  66)*yin(  37)*zin(  43)+xin( 102)*yin(  73)*zin(  79)+xin( 138)*yin( 109)*zin( 115))
          eri_value(   38)=eri_value(   38)+d12bra(  3)*d12ket(  2)*(xin(  29)*yin(   2)*zin(   7)+xin(  65)*yin(  38)*zin(  43)+xin( 101)*yin(  74)*zin(  79)+xin( 137)*yin( 110)*zin( 115))
          eri_value(   39)=eri_value(   39)+d12bra(  3)*d12ket(  3)*(xin(  29)*yin(   1)*zin(   8)+xin(  65)*yin(  37)*zin(  44)+xin( 101)*yin(  73)*zin(  80)+xin( 137)*yin( 109)*zin( 116))
          eri_value(   40)=eri_value(   40)+d12bra(  3)*d12ket(  4)*(xin(  26)*yin(   5)*zin(   7)+xin(  62)*yin(  41)*zin(  43)+xin(  98)*yin(  77)*zin(  79)+xin( 134)*yin( 113)*zin( 115))
          eri_value(   41)=eri_value(   41)+d12bra(  3)*d12ket(  5)*(xin(  25)*yin(   6)*zin(   7)+xin(  61)*yin(  42)*zin(  43)+xin(  97)*yin(  78)*zin(  79)+xin( 133)*yin( 114)*zin( 115))
          eri_value(   42)=eri_value(   42)+d12bra(  3)*d12ket(  6)*(xin(  25)*yin(   5)*zin(   8)+xin(  61)*yin(  41)*zin(  44)+xin(  97)*yin(  77)*zin(  80)+xin( 133)*yin( 113)*zin( 116))
          eri_value(   43)=eri_value(   43)+d12bra(  3)*d12ket(  7)*(xin(  26)*yin(   1)*zin(  11)+xin(  62)*yin(  37)*zin(  47)+xin(  98)*yin(  73)*zin(  83)+xin( 134)*yin( 109)*zin( 119))
          eri_value(   44)=eri_value(   44)+d12bra(  3)*d12ket(  8)*(xin(  25)*yin(   2)*zin(  11)+xin(  61)*yin(  38)*zin(  47)+xin(  97)*yin(  74)*zin(  83)+xin( 133)*yin( 110)*zin( 119))
          eri_value(   45)=eri_value(   45)+d12bra(  3)*d12ket(  9)*(xin(  25)*yin(   1)*zin(  12)+xin(  61)*yin(  37)*zin(  48)+xin(  97)*yin(  73)*zin(  84)+xin( 133)*yin( 109)*zin( 120))
          eri_value(   46)=eri_value(   46)+d12bra(  3)*d12ket( 10)*(xin(  28)*yin(   3)*zin(   7)+xin(  64)*yin(  39)*zin(  43)+xin( 100)*yin(  75)*zin(  79)+xin( 136)*yin( 111)*zin( 115))
          eri_value(   47)=eri_value(   47)+d12bra(  3)*d12ket( 11)*(xin(  27)*yin(   4)*zin(   7)+xin(  63)*yin(  40)*zin(  43)+xin(  99)*yin(  76)*zin(  79)+xin( 135)*yin( 112)*zin( 115))
          eri_value(   48)=eri_value(   48)+d12bra(  3)*d12ket( 12)*(xin(  27)*yin(   3)*zin(   8)+xin(  63)*yin(  39)*zin(  44)+xin(  99)*yin(  75)*zin(  80)+xin( 135)*yin( 111)*zin( 116))
          eri_value(   49)=eri_value(   49)+d12bra(  3)*d12ket( 13)*(xin(  28)*yin(   1)*zin(   9)+xin(  64)*yin(  37)*zin(  45)+xin( 100)*yin(  73)*zin(  81)+xin( 136)*yin( 109)*zin( 117))
          eri_value(   50)=eri_value(   50)+d12bra(  3)*d12ket( 14)*(xin(  27)*yin(   2)*zin(   9)+xin(  63)*yin(  38)*zin(  45)+xin(  99)*yin(  74)*zin(  81)+xin( 135)*yin( 110)*zin( 117))
          eri_value(   51)=eri_value(   51)+d12bra(  3)*d12ket( 15)*(xin(  27)*yin(   1)*zin(  10)+xin(  63)*yin(  37)*zin(  46)+xin(  99)*yin(  73)*zin(  82)+xin( 135)*yin( 109)*zin( 118))
          eri_value(   52)=eri_value(   52)+d12bra(  3)*d12ket( 16)*(xin(  26)*yin(   3)*zin(   9)+xin(  62)*yin(  39)*zin(  45)+xin(  98)*yin(  75)*zin(  81)+xin( 134)*yin( 111)*zin( 117))
          eri_value(   53)=eri_value(   53)+d12bra(  3)*d12ket( 17)*(xin(  25)*yin(   4)*zin(   9)+xin(  61)*yin(  40)*zin(  45)+xin(  97)*yin(  76)*zin(  81)+xin( 133)*yin( 112)*zin( 117))
          eri_value(   54)=eri_value(   54)+d12bra(  3)*d12ket( 18)*(xin(  25)*yin(   3)*zin(  10)+xin(  61)*yin(  39)*zin(  46)+xin(  97)*yin(  75)*zin(  82)+xin( 133)*yin( 111)*zin( 118))
          eri_value(   55)=eri_value(   55)+d12bra(  4)*d12ket(  1)*(xin(  12)*yin(  25)*zin(   1)+xin(  48)*yin(  61)*zin(  37)+xin(  84)*yin(  97)*zin(  73)+xin( 120)*yin( 133)*zin( 109))
          eri_value(   56)=eri_value(   56)+d12bra(  4)*d12ket(  2)*(xin(  11)*yin(  26)*zin(   1)+xin(  47)*yin(  62)*zin(  37)+xin(  83)*yin(  98)*zin(  73)+xin( 119)*yin( 134)*zin( 109))
          eri_value(   57)=eri_value(   57)+d12bra(  4)*d12ket(  3)*(xin(  11)*yin(  25)*zin(   2)+xin(  47)*yin(  61)*zin(  38)+xin(  83)*yin(  97)*zin(  74)+xin( 119)*yin( 133)*zin( 110))
          eri_value(   58)=eri_value(   58)+d12bra(  4)*d12ket(  4)*(xin(   8)*yin(  29)*zin(   1)+xin(  44)*yin(  65)*zin(  37)+xin(  80)*yin( 101)*zin(  73)+xin( 116)*yin( 137)*zin( 109))
          eri_value(   59)=eri_value(   59)+d12bra(  4)*d12ket(  5)*(xin(   7)*yin(  30)*zin(   1)+xin(  43)*yin(  66)*zin(  37)+xin(  79)*yin( 102)*zin(  73)+xin( 115)*yin( 138)*zin( 109))
          eri_value(   60)=eri_value(   60)+d12bra(  4)*d12ket(  6)*(xin(   7)*yin(  29)*zin(   2)+xin(  43)*yin(  65)*zin(  38)+xin(  79)*yin( 101)*zin(  74)+xin( 115)*yin( 137)*zin( 110))
          eri_value(   61)=eri_value(   61)+d12bra(  4)*d12ket(  7)*(xin(   8)*yin(  25)*zin(   5)+xin(  44)*yin(  61)*zin(  41)+xin(  80)*yin(  97)*zin(  77)+xin( 116)*yin( 133)*zin( 113))
          eri_value(   62)=eri_value(   62)+d12bra(  4)*d12ket(  8)*(xin(   7)*yin(  26)*zin(   5)+xin(  43)*yin(  62)*zin(  41)+xin(  79)*yin(  98)*zin(  77)+xin( 115)*yin( 134)*zin( 113))
          eri_value(   63)=eri_value(   63)+d12bra(  4)*d12ket(  9)*(xin(   7)*yin(  25)*zin(   6)+xin(  43)*yin(  61)*zin(  42)+xin(  79)*yin(  97)*zin(  78)+xin( 115)*yin( 133)*zin( 114))
          eri_value(   64)=eri_value(   64)+d12bra(  4)*d12ket( 10)*(xin(  10)*yin(  27)*zin(   1)+xin(  46)*yin(  63)*zin(  37)+xin(  82)*yin(  99)*zin(  73)+xin( 118)*yin( 135)*zin( 109))
          eri_value(   65)=eri_value(   65)+d12bra(  4)*d12ket( 11)*(xin(   9)*yin(  28)*zin(   1)+xin(  45)*yin(  64)*zin(  37)+xin(  81)*yin( 100)*zin(  73)+xin( 117)*yin( 136)*zin( 109))
          eri_value(   66)=eri_value(   66)+d12bra(  4)*d12ket( 12)*(xin(   9)*yin(  27)*zin(   2)+xin(  45)*yin(  63)*zin(  38)+xin(  81)*yin(  99)*zin(  74)+xin( 117)*yin( 135)*zin( 110))
          eri_value(   67)=eri_value(   67)+d12bra(  4)*d12ket( 13)*(xin(  10)*yin(  25)*zin(   3)+xin(  46)*yin(  61)*zin(  39)+xin(  82)*yin(  97)*zin(  75)+xin( 118)*yin( 133)*zin( 111))
          eri_value(   68)=eri_value(   68)+d12bra(  4)*d12ket( 14)*(xin(   9)*yin(  26)*zin(   3)+xin(  45)*yin(  62)*zin(  39)+xin(  81)*yin(  98)*zin(  75)+xin( 117)*yin( 134)*zin( 111))
          eri_value(   69)=eri_value(   69)+d12bra(  4)*d12ket( 15)*(xin(   9)*yin(  25)*zin(   4)+xin(  45)*yin(  61)*zin(  40)+xin(  81)*yin(  97)*zin(  76)+xin( 117)*yin( 133)*zin( 112))
          eri_value(   70)=eri_value(   70)+d12bra(  4)*d12ket( 16)*(xin(   8)*yin(  27)*zin(   3)+xin(  44)*yin(  63)*zin(  39)+xin(  80)*yin(  99)*zin(  75)+xin( 116)*yin( 135)*zin( 111))
          eri_value(   71)=eri_value(   71)+d12bra(  4)*d12ket( 17)*(xin(   7)*yin(  28)*zin(   3)+xin(  43)*yin(  64)*zin(  39)+xin(  79)*yin( 100)*zin(  75)+xin( 115)*yin( 136)*zin( 111))
          eri_value(   72)=eri_value(   72)+d12bra(  4)*d12ket( 18)*(xin(   7)*yin(  27)*zin(   4)+xin(  43)*yin(  63)*zin(  40)+xin(  79)*yin(  99)*zin(  76)+xin( 115)*yin( 135)*zin( 112))
          eri_value(   73)=eri_value(   73)+d12bra(  5)*d12ket(  1)*(xin(   6)*yin(  31)*zin(   1)+xin(  42)*yin(  67)*zin(  37)+xin(  78)*yin( 103)*zin(  73)+xin( 114)*yin( 139)*zin( 109))
          eri_value(   74)=eri_value(   74)+d12bra(  5)*d12ket(  2)*(xin(   5)*yin(  32)*zin(   1)+xin(  41)*yin(  68)*zin(  37)+xin(  77)*yin( 104)*zin(  73)+xin( 113)*yin( 140)*zin( 109))
          eri_value(   75)=eri_value(   75)+d12bra(  5)*d12ket(  3)*(xin(   5)*yin(  31)*zin(   2)+xin(  41)*yin(  67)*zin(  38)+xin(  77)*yin( 103)*zin(  74)+xin( 113)*yin( 139)*zin( 110))
          eri_value(   76)=eri_value(   76)+d12bra(  5)*d12ket(  4)*(xin(   2)*yin(  35)*zin(   1)+xin(  38)*yin(  71)*zin(  37)+xin(  74)*yin( 107)*zin(  73)+xin( 110)*yin( 143)*zin( 109))
          eri_value(   77)=eri_value(   77)+d12bra(  5)*d12ket(  5)*(xin(   1)*yin(  36)*zin(   1)+xin(  37)*yin(  72)*zin(  37)+xin(  73)*yin( 108)*zin(  73)+xin( 109)*yin( 144)*zin( 109))
          eri_value(   78)=eri_value(   78)+d12bra(  5)*d12ket(  6)*(xin(   1)*yin(  35)*zin(   2)+xin(  37)*yin(  71)*zin(  38)+xin(  73)*yin( 107)*zin(  74)+xin( 109)*yin( 143)*zin( 110))
          eri_value(   79)=eri_value(   79)+d12bra(  5)*d12ket(  7)*(xin(   2)*yin(  31)*zin(   5)+xin(  38)*yin(  67)*zin(  41)+xin(  74)*yin( 103)*zin(  77)+xin( 110)*yin( 139)*zin( 113))
          eri_value(   80)=eri_value(   80)+d12bra(  5)*d12ket(  8)*(xin(   1)*yin(  32)*zin(   5)+xin(  37)*yin(  68)*zin(  41)+xin(  73)*yin( 104)*zin(  77)+xin( 109)*yin( 140)*zin( 113))
          eri_value(   81)=eri_value(   81)+d12bra(  5)*d12ket(  9)*(xin(   1)*yin(  31)*zin(   6)+xin(  37)*yin(  67)*zin(  42)+xin(  73)*yin( 103)*zin(  78)+xin( 109)*yin( 139)*zin( 114))
          eri_value(   82)=eri_value(   82)+d12bra(  5)*d12ket( 10)*(xin(   4)*yin(  33)*zin(   1)+xin(  40)*yin(  69)*zin(  37)+xin(  76)*yin( 105)*zin(  73)+xin( 112)*yin( 141)*zin( 109))
          eri_value(   83)=eri_value(   83)+d12bra(  5)*d12ket( 11)*(xin(   3)*yin(  34)*zin(   1)+xin(  39)*yin(  70)*zin(  37)+xin(  75)*yin( 106)*zin(  73)+xin( 111)*yin( 142)*zin( 109))
          eri_value(   84)=eri_value(   84)+d12bra(  5)*d12ket( 12)*(xin(   3)*yin(  33)*zin(   2)+xin(  39)*yin(  69)*zin(  38)+xin(  75)*yin( 105)*zin(  74)+xin( 111)*yin( 141)*zin( 110))
          eri_value(   85)=eri_value(   85)+d12bra(  5)*d12ket( 13)*(xin(   4)*yin(  31)*zin(   3)+xin(  40)*yin(  67)*zin(  39)+xin(  76)*yin( 103)*zin(  75)+xin( 112)*yin( 139)*zin( 111))
          eri_value(   86)=eri_value(   86)+d12bra(  5)*d12ket( 14)*(xin(   3)*yin(  32)*zin(   3)+xin(  39)*yin(  68)*zin(  39)+xin(  75)*yin( 104)*zin(  75)+xin( 111)*yin( 140)*zin( 111))
          eri_value(   87)=eri_value(   87)+d12bra(  5)*d12ket( 15)*(xin(   3)*yin(  31)*zin(   4)+xin(  39)*yin(  67)*zin(  40)+xin(  75)*yin( 103)*zin(  76)+xin( 111)*yin( 139)*zin( 112))
          eri_value(   88)=eri_value(   88)+d12bra(  5)*d12ket( 16)*(xin(   2)*yin(  33)*zin(   3)+xin(  38)*yin(  69)*zin(  39)+xin(  74)*yin( 105)*zin(  75)+xin( 110)*yin( 141)*zin( 111))
          eri_value(   89)=eri_value(   89)+d12bra(  5)*d12ket( 17)*(xin(   1)*yin(  34)*zin(   3)+xin(  37)*yin(  70)*zin(  39)+xin(  73)*yin( 106)*zin(  75)+xin( 109)*yin( 142)*zin( 111))
          eri_value(   90)=eri_value(   90)+d12bra(  5)*d12ket( 18)*(xin(   1)*yin(  33)*zin(   4)+xin(  37)*yin(  69)*zin(  40)+xin(  73)*yin( 105)*zin(  76)+xin( 109)*yin( 141)*zin( 112))
          eri_value(   91)=eri_value(   91)+d12bra(  6)*d12ket(  1)*(xin(   6)*yin(  25)*zin(   7)+xin(  42)*yin(  61)*zin(  43)+xin(  78)*yin(  97)*zin(  79)+xin( 114)*yin( 133)*zin( 115))
          eri_value(   92)=eri_value(   92)+d12bra(  6)*d12ket(  2)*(xin(   5)*yin(  26)*zin(   7)+xin(  41)*yin(  62)*zin(  43)+xin(  77)*yin(  98)*zin(  79)+xin( 113)*yin( 134)*zin( 115))
          eri_value(   93)=eri_value(   93)+d12bra(  6)*d12ket(  3)*(xin(   5)*yin(  25)*zin(   8)+xin(  41)*yin(  61)*zin(  44)+xin(  77)*yin(  97)*zin(  80)+xin( 113)*yin( 133)*zin( 116))
          eri_value(   94)=eri_value(   94)+d12bra(  6)*d12ket(  4)*(xin(   2)*yin(  29)*zin(   7)+xin(  38)*yin(  65)*zin(  43)+xin(  74)*yin( 101)*zin(  79)+xin( 110)*yin( 137)*zin( 115))
          eri_value(   95)=eri_value(   95)+d12bra(  6)*d12ket(  5)*(xin(   1)*yin(  30)*zin(   7)+xin(  37)*yin(  66)*zin(  43)+xin(  73)*yin( 102)*zin(  79)+xin( 109)*yin( 138)*zin( 115))
          eri_value(   96)=eri_value(   96)+d12bra(  6)*d12ket(  6)*(xin(   1)*yin(  29)*zin(   8)+xin(  37)*yin(  65)*zin(  44)+xin(  73)*yin( 101)*zin(  80)+xin( 109)*yin( 137)*zin( 116))
          eri_value(   97)=eri_value(   97)+d12bra(  6)*d12ket(  7)*(xin(   2)*yin(  25)*zin(  11)+xin(  38)*yin(  61)*zin(  47)+xin(  74)*yin(  97)*zin(  83)+xin( 110)*yin( 133)*zin( 119))
          eri_value(   98)=eri_value(   98)+d12bra(  6)*d12ket(  8)*(xin(   1)*yin(  26)*zin(  11)+xin(  37)*yin(  62)*zin(  47)+xin(  73)*yin(  98)*zin(  83)+xin( 109)*yin( 134)*zin( 119))
          eri_value(   99)=eri_value(   99)+d12bra(  6)*d12ket(  9)*(xin(   1)*yin(  25)*zin(  12)+xin(  37)*yin(  61)*zin(  48)+xin(  73)*yin(  97)*zin(  84)+xin( 109)*yin( 133)*zin( 120))
          eri_value(  100)=eri_value(  100)+d12bra(  6)*d12ket( 10)*(xin(   4)*yin(  27)*zin(   7)+xin(  40)*yin(  63)*zin(  43)+xin(  76)*yin(  99)*zin(  79)+xin( 112)*yin( 135)*zin( 115))
          eri_value(  101)=eri_value(  101)+d12bra(  6)*d12ket( 11)*(xin(   3)*yin(  28)*zin(   7)+xin(  39)*yin(  64)*zin(  43)+xin(  75)*yin( 100)*zin(  79)+xin( 111)*yin( 136)*zin( 115))
          eri_value(  102)=eri_value(  102)+d12bra(  6)*d12ket( 12)*(xin(   3)*yin(  27)*zin(   8)+xin(  39)*yin(  63)*zin(  44)+xin(  75)*yin(  99)*zin(  80)+xin( 111)*yin( 135)*zin( 116))
          eri_value(  103)=eri_value(  103)+d12bra(  6)*d12ket( 13)*(xin(   4)*yin(  25)*zin(   9)+xin(  40)*yin(  61)*zin(  45)+xin(  76)*yin(  97)*zin(  81)+xin( 112)*yin( 133)*zin( 117))
          eri_value(  104)=eri_value(  104)+d12bra(  6)*d12ket( 14)*(xin(   3)*yin(  26)*zin(   9)+xin(  39)*yin(  62)*zin(  45)+xin(  75)*yin(  98)*zin(  81)+xin( 111)*yin( 134)*zin( 117))
          eri_value(  105)=eri_value(  105)+d12bra(  6)*d12ket( 15)*(xin(   3)*yin(  25)*zin(  10)+xin(  39)*yin(  61)*zin(  46)+xin(  75)*yin(  97)*zin(  82)+xin( 111)*yin( 133)*zin( 118))
          eri_value(  106)=eri_value(  106)+d12bra(  6)*d12ket( 16)*(xin(   2)*yin(  27)*zin(   9)+xin(  38)*yin(  63)*zin(  45)+xin(  74)*yin(  99)*zin(  81)+xin( 110)*yin( 135)*zin( 117))
          eri_value(  107)=eri_value(  107)+d12bra(  6)*d12ket( 17)*(xin(   1)*yin(  28)*zin(   9)+xin(  37)*yin(  64)*zin(  45)+xin(  73)*yin( 100)*zin(  81)+xin( 109)*yin( 136)*zin( 117))
          eri_value(  108)=eri_value(  108)+d12bra(  6)*d12ket( 18)*(xin(   1)*yin(  27)*zin(  10)+xin(  37)*yin(  63)*zin(  46)+xin(  73)*yin(  99)*zin(  82)+xin( 109)*yin( 135)*zin( 118))
          eri_value(  109)=eri_value(  109)+d12bra(  7)*d12ket(  1)*(xin(  12)*yin(   1)*zin(  25)+xin(  48)*yin(  37)*zin(  61)+xin(  84)*yin(  73)*zin(  97)+xin( 120)*yin( 109)*zin( 133))
          eri_value(  110)=eri_value(  110)+d12bra(  7)*d12ket(  2)*(xin(  11)*yin(   2)*zin(  25)+xin(  47)*yin(  38)*zin(  61)+xin(  83)*yin(  74)*zin(  97)+xin( 119)*yin( 110)*zin( 133))
          eri_value(  111)=eri_value(  111)+d12bra(  7)*d12ket(  3)*(xin(  11)*yin(   1)*zin(  26)+xin(  47)*yin(  37)*zin(  62)+xin(  83)*yin(  73)*zin(  98)+xin( 119)*yin( 109)*zin( 134))
          eri_value(  112)=eri_value(  112)+d12bra(  7)*d12ket(  4)*(xin(   8)*yin(   5)*zin(  25)+xin(  44)*yin(  41)*zin(  61)+xin(  80)*yin(  77)*zin(  97)+xin( 116)*yin( 113)*zin( 133))
          eri_value(  113)=eri_value(  113)+d12bra(  7)*d12ket(  5)*(xin(   7)*yin(   6)*zin(  25)+xin(  43)*yin(  42)*zin(  61)+xin(  79)*yin(  78)*zin(  97)+xin( 115)*yin( 114)*zin( 133))
          eri_value(  114)=eri_value(  114)+d12bra(  7)*d12ket(  6)*(xin(   7)*yin(   5)*zin(  26)+xin(  43)*yin(  41)*zin(  62)+xin(  79)*yin(  77)*zin(  98)+xin( 115)*yin( 113)*zin( 134))
          eri_value(  115)=eri_value(  115)+d12bra(  7)*d12ket(  7)*(xin(   8)*yin(   1)*zin(  29)+xin(  44)*yin(  37)*zin(  65)+xin(  80)*yin(  73)*zin( 101)+xin( 116)*yin( 109)*zin( 137))
          eri_value(  116)=eri_value(  116)+d12bra(  7)*d12ket(  8)*(xin(   7)*yin(   2)*zin(  29)+xin(  43)*yin(  38)*zin(  65)+xin(  79)*yin(  74)*zin( 101)+xin( 115)*yin( 110)*zin( 137))
          eri_value(  117)=eri_value(  117)+d12bra(  7)*d12ket(  9)*(xin(   7)*yin(   1)*zin(  30)+xin(  43)*yin(  37)*zin(  66)+xin(  79)*yin(  73)*zin( 102)+xin( 115)*yin( 109)*zin( 138))
          eri_value(  118)=eri_value(  118)+d12bra(  7)*d12ket( 10)*(xin(  10)*yin(   3)*zin(  25)+xin(  46)*yin(  39)*zin(  61)+xin(  82)*yin(  75)*zin(  97)+xin( 118)*yin( 111)*zin( 133))
          eri_value(  119)=eri_value(  119)+d12bra(  7)*d12ket( 11)*(xin(   9)*yin(   4)*zin(  25)+xin(  45)*yin(  40)*zin(  61)+xin(  81)*yin(  76)*zin(  97)+xin( 117)*yin( 112)*zin( 133))
          eri_value(  120)=eri_value(  120)+d12bra(  7)*d12ket( 12)*(xin(   9)*yin(   3)*zin(  26)+xin(  45)*yin(  39)*zin(  62)+xin(  81)*yin(  75)*zin(  98)+xin( 117)*yin( 111)*zin( 134))
          eri_value(  121)=eri_value(  121)+d12bra(  7)*d12ket( 13)*(xin(  10)*yin(   1)*zin(  27)+xin(  46)*yin(  37)*zin(  63)+xin(  82)*yin(  73)*zin(  99)+xin( 118)*yin( 109)*zin( 135))
          eri_value(  122)=eri_value(  122)+d12bra(  7)*d12ket( 14)*(xin(   9)*yin(   2)*zin(  27)+xin(  45)*yin(  38)*zin(  63)+xin(  81)*yin(  74)*zin(  99)+xin( 117)*yin( 110)*zin( 135))
          eri_value(  123)=eri_value(  123)+d12bra(  7)*d12ket( 15)*(xin(   9)*yin(   1)*zin(  28)+xin(  45)*yin(  37)*zin(  64)+xin(  81)*yin(  73)*zin( 100)+xin( 117)*yin( 109)*zin( 136))
          eri_value(  124)=eri_value(  124)+d12bra(  7)*d12ket( 16)*(xin(   8)*yin(   3)*zin(  27)+xin(  44)*yin(  39)*zin(  63)+xin(  80)*yin(  75)*zin(  99)+xin( 116)*yin( 111)*zin( 135))
          eri_value(  125)=eri_value(  125)+d12bra(  7)*d12ket( 17)*(xin(   7)*yin(   4)*zin(  27)+xin(  43)*yin(  40)*zin(  63)+xin(  79)*yin(  76)*zin(  99)+xin( 115)*yin( 112)*zin( 135))
          eri_value(  126)=eri_value(  126)+d12bra(  7)*d12ket( 18)*(xin(   7)*yin(   3)*zin(  28)+xin(  43)*yin(  39)*zin(  64)+xin(  79)*yin(  75)*zin( 100)+xin( 115)*yin( 111)*zin( 136))
          eri_value(  127)=eri_value(  127)+d12bra(  8)*d12ket(  1)*(xin(   6)*yin(   7)*zin(  25)+xin(  42)*yin(  43)*zin(  61)+xin(  78)*yin(  79)*zin(  97)+xin( 114)*yin( 115)*zin( 133))
          eri_value(  128)=eri_value(  128)+d12bra(  8)*d12ket(  2)*(xin(   5)*yin(   8)*zin(  25)+xin(  41)*yin(  44)*zin(  61)+xin(  77)*yin(  80)*zin(  97)+xin( 113)*yin( 116)*zin( 133))
          eri_value(  129)=eri_value(  129)+d12bra(  8)*d12ket(  3)*(xin(   5)*yin(   7)*zin(  26)+xin(  41)*yin(  43)*zin(  62)+xin(  77)*yin(  79)*zin(  98)+xin( 113)*yin( 115)*zin( 134))
          eri_value(  130)=eri_value(  130)+d12bra(  8)*d12ket(  4)*(xin(   2)*yin(  11)*zin(  25)+xin(  38)*yin(  47)*zin(  61)+xin(  74)*yin(  83)*zin(  97)+xin( 110)*yin( 119)*zin( 133))
          eri_value(  131)=eri_value(  131)+d12bra(  8)*d12ket(  5)*(xin(   1)*yin(  12)*zin(  25)+xin(  37)*yin(  48)*zin(  61)+xin(  73)*yin(  84)*zin(  97)+xin( 109)*yin( 120)*zin( 133))
          eri_value(  132)=eri_value(  132)+d12bra(  8)*d12ket(  6)*(xin(   1)*yin(  11)*zin(  26)+xin(  37)*yin(  47)*zin(  62)+xin(  73)*yin(  83)*zin(  98)+xin( 109)*yin( 119)*zin( 134))
          eri_value(  133)=eri_value(  133)+d12bra(  8)*d12ket(  7)*(xin(   2)*yin(   7)*zin(  29)+xin(  38)*yin(  43)*zin(  65)+xin(  74)*yin(  79)*zin( 101)+xin( 110)*yin( 115)*zin( 137))
          eri_value(  134)=eri_value(  134)+d12bra(  8)*d12ket(  8)*(xin(   1)*yin(   8)*zin(  29)+xin(  37)*yin(  44)*zin(  65)+xin(  73)*yin(  80)*zin( 101)+xin( 109)*yin( 116)*zin( 137))
          eri_value(  135)=eri_value(  135)+d12bra(  8)*d12ket(  9)*(xin(   1)*yin(   7)*zin(  30)+xin(  37)*yin(  43)*zin(  66)+xin(  73)*yin(  79)*zin( 102)+xin( 109)*yin( 115)*zin( 138))
          eri_value(  136)=eri_value(  136)+d12bra(  8)*d12ket( 10)*(xin(   4)*yin(   9)*zin(  25)+xin(  40)*yin(  45)*zin(  61)+xin(  76)*yin(  81)*zin(  97)+xin( 112)*yin( 117)*zin( 133))
          eri_value(  137)=eri_value(  137)+d12bra(  8)*d12ket( 11)*(xin(   3)*yin(  10)*zin(  25)+xin(  39)*yin(  46)*zin(  61)+xin(  75)*yin(  82)*zin(  97)+xin( 111)*yin( 118)*zin( 133))
          eri_value(  138)=eri_value(  138)+d12bra(  8)*d12ket( 12)*(xin(   3)*yin(   9)*zin(  26)+xin(  39)*yin(  45)*zin(  62)+xin(  75)*yin(  81)*zin(  98)+xin( 111)*yin( 117)*zin( 134))
          eri_value(  139)=eri_value(  139)+d12bra(  8)*d12ket( 13)*(xin(   4)*yin(   7)*zin(  27)+xin(  40)*yin(  43)*zin(  63)+xin(  76)*yin(  79)*zin(  99)+xin( 112)*yin( 115)*zin( 135))
          eri_value(  140)=eri_value(  140)+d12bra(  8)*d12ket( 14)*(xin(   3)*yin(   8)*zin(  27)+xin(  39)*yin(  44)*zin(  63)+xin(  75)*yin(  80)*zin(  99)+xin( 111)*yin( 116)*zin( 135))
          eri_value(  141)=eri_value(  141)+d12bra(  8)*d12ket( 15)*(xin(   3)*yin(   7)*zin(  28)+xin(  39)*yin(  43)*zin(  64)+xin(  75)*yin(  79)*zin( 100)+xin( 111)*yin( 115)*zin( 136))
          eri_value(  142)=eri_value(  142)+d12bra(  8)*d12ket( 16)*(xin(   2)*yin(   9)*zin(  27)+xin(  38)*yin(  45)*zin(  63)+xin(  74)*yin(  81)*zin(  99)+xin( 110)*yin( 117)*zin( 135))
          eri_value(  143)=eri_value(  143)+d12bra(  8)*d12ket( 17)*(xin(   1)*yin(  10)*zin(  27)+xin(  37)*yin(  46)*zin(  63)+xin(  73)*yin(  82)*zin(  99)+xin( 109)*yin( 118)*zin( 135))
          eri_value(  144)=eri_value(  144)+d12bra(  8)*d12ket( 18)*(xin(   1)*yin(   9)*zin(  28)+xin(  37)*yin(  45)*zin(  64)+xin(  73)*yin(  81)*zin( 100)+xin( 109)*yin( 117)*zin( 136))
          eri_value(  145)=eri_value(  145)+d12bra(  9)*d12ket(  1)*(xin(   6)*yin(   1)*zin(  31)+xin(  42)*yin(  37)*zin(  67)+xin(  78)*yin(  73)*zin( 103)+xin( 114)*yin( 109)*zin( 139))
          eri_value(  146)=eri_value(  146)+d12bra(  9)*d12ket(  2)*(xin(   5)*yin(   2)*zin(  31)+xin(  41)*yin(  38)*zin(  67)+xin(  77)*yin(  74)*zin( 103)+xin( 113)*yin( 110)*zin( 139))
          eri_value(  147)=eri_value(  147)+d12bra(  9)*d12ket(  3)*(xin(   5)*yin(   1)*zin(  32)+xin(  41)*yin(  37)*zin(  68)+xin(  77)*yin(  73)*zin( 104)+xin( 113)*yin( 109)*zin( 140))
          eri_value(  148)=eri_value(  148)+d12bra(  9)*d12ket(  4)*(xin(   2)*yin(   5)*zin(  31)+xin(  38)*yin(  41)*zin(  67)+xin(  74)*yin(  77)*zin( 103)+xin( 110)*yin( 113)*zin( 139))
          eri_value(  149)=eri_value(  149)+d12bra(  9)*d12ket(  5)*(xin(   1)*yin(   6)*zin(  31)+xin(  37)*yin(  42)*zin(  67)+xin(  73)*yin(  78)*zin( 103)+xin( 109)*yin( 114)*zin( 139))
          eri_value(  150)=eri_value(  150)+d12bra(  9)*d12ket(  6)*(xin(   1)*yin(   5)*zin(  32)+xin(  37)*yin(  41)*zin(  68)+xin(  73)*yin(  77)*zin( 104)+xin( 109)*yin( 113)*zin( 140))
          eri_value(  151)=eri_value(  151)+d12bra(  9)*d12ket(  7)*(xin(   2)*yin(   1)*zin(  35)+xin(  38)*yin(  37)*zin(  71)+xin(  74)*yin(  73)*zin( 107)+xin( 110)*yin( 109)*zin( 143))
          eri_value(  152)=eri_value(  152)+d12bra(  9)*d12ket(  8)*(xin(   1)*yin(   2)*zin(  35)+xin(  37)*yin(  38)*zin(  71)+xin(  73)*yin(  74)*zin( 107)+xin( 109)*yin( 110)*zin( 143))
          eri_value(  153)=eri_value(  153)+d12bra(  9)*d12ket(  9)*(xin(   1)*yin(   1)*zin(  36)+xin(  37)*yin(  37)*zin(  72)+xin(  73)*yin(  73)*zin( 108)+xin( 109)*yin( 109)*zin( 144))
          eri_value(  154)=eri_value(  154)+d12bra(  9)*d12ket( 10)*(xin(   4)*yin(   3)*zin(  31)+xin(  40)*yin(  39)*zin(  67)+xin(  76)*yin(  75)*zin( 103)+xin( 112)*yin( 111)*zin( 139))
          eri_value(  155)=eri_value(  155)+d12bra(  9)*d12ket( 11)*(xin(   3)*yin(   4)*zin(  31)+xin(  39)*yin(  40)*zin(  67)+xin(  75)*yin(  76)*zin( 103)+xin( 111)*yin( 112)*zin( 139))
          eri_value(  156)=eri_value(  156)+d12bra(  9)*d12ket( 12)*(xin(   3)*yin(   3)*zin(  32)+xin(  39)*yin(  39)*zin(  68)+xin(  75)*yin(  75)*zin( 104)+xin( 111)*yin( 111)*zin( 140))
          eri_value(  157)=eri_value(  157)+d12bra(  9)*d12ket( 13)*(xin(   4)*yin(   1)*zin(  33)+xin(  40)*yin(  37)*zin(  69)+xin(  76)*yin(  73)*zin( 105)+xin( 112)*yin( 109)*zin( 141))
          eri_value(  158)=eri_value(  158)+d12bra(  9)*d12ket( 14)*(xin(   3)*yin(   2)*zin(  33)+xin(  39)*yin(  38)*zin(  69)+xin(  75)*yin(  74)*zin( 105)+xin( 111)*yin( 110)*zin( 141))
          eri_value(  159)=eri_value(  159)+d12bra(  9)*d12ket( 15)*(xin(   3)*yin(   1)*zin(  34)+xin(  39)*yin(  37)*zin(  70)+xin(  75)*yin(  73)*zin( 106)+xin( 111)*yin( 109)*zin( 142))
          eri_value(  160)=eri_value(  160)+d12bra(  9)*d12ket( 16)*(xin(   2)*yin(   3)*zin(  33)+xin(  38)*yin(  39)*zin(  69)+xin(  74)*yin(  75)*zin( 105)+xin( 110)*yin( 111)*zin( 141))
          eri_value(  161)=eri_value(  161)+d12bra(  9)*d12ket( 17)*(xin(   1)*yin(   4)*zin(  33)+xin(  37)*yin(  40)*zin(  69)+xin(  73)*yin(  76)*zin( 105)+xin( 109)*yin( 112)*zin( 141))
          eri_value(  162)=eri_value(  162)+d12bra(  9)*d12ket( 18)*(xin(   1)*yin(   3)*zin(  34)+xin(  37)*yin(  39)*zin(  70)+xin(  73)*yin(  75)*zin( 106)+xin( 109)*yin( 111)*zin( 142))
          eri_value(  163)=eri_value(  163)+d12bra( 10)*d12ket(  1)*(xin(  24)*yin(  13)*zin(   1)+xin(  60)*yin(  49)*zin(  37)+xin(  96)*yin(  85)*zin(  73)+xin( 132)*yin( 121)*zin( 109))
          eri_value(  164)=eri_value(  164)+d12bra( 10)*d12ket(  2)*(xin(  23)*yin(  14)*zin(   1)+xin(  59)*yin(  50)*zin(  37)+xin(  95)*yin(  86)*zin(  73)+xin( 131)*yin( 122)*zin( 109))
          eri_value(  165)=eri_value(  165)+d12bra( 10)*d12ket(  3)*(xin(  23)*yin(  13)*zin(   2)+xin(  59)*yin(  49)*zin(  38)+xin(  95)*yin(  85)*zin(  74)+xin( 131)*yin( 121)*zin( 110))
          eri_value(  166)=eri_value(  166)+d12bra( 10)*d12ket(  4)*(xin(  20)*yin(  17)*zin(   1)+xin(  56)*yin(  53)*zin(  37)+xin(  92)*yin(  89)*zin(  73)+xin( 128)*yin( 125)*zin( 109))
          eri_value(  167)=eri_value(  167)+d12bra( 10)*d12ket(  5)*(xin(  19)*yin(  18)*zin(   1)+xin(  55)*yin(  54)*zin(  37)+xin(  91)*yin(  90)*zin(  73)+xin( 127)*yin( 126)*zin( 109))
          eri_value(  168)=eri_value(  168)+d12bra( 10)*d12ket(  6)*(xin(  19)*yin(  17)*zin(   2)+xin(  55)*yin(  53)*zin(  38)+xin(  91)*yin(  89)*zin(  74)+xin( 127)*yin( 125)*zin( 110))
          eri_value(  169)=eri_value(  169)+d12bra( 10)*d12ket(  7)*(xin(  20)*yin(  13)*zin(   5)+xin(  56)*yin(  49)*zin(  41)+xin(  92)*yin(  85)*zin(  77)+xin( 128)*yin( 121)*zin( 113))
          eri_value(  170)=eri_value(  170)+d12bra( 10)*d12ket(  8)*(xin(  19)*yin(  14)*zin(   5)+xin(  55)*yin(  50)*zin(  41)+xin(  91)*yin(  86)*zin(  77)+xin( 127)*yin( 122)*zin( 113))
          eri_value(  171)=eri_value(  171)+d12bra( 10)*d12ket(  9)*(xin(  19)*yin(  13)*zin(   6)+xin(  55)*yin(  49)*zin(  42)+xin(  91)*yin(  85)*zin(  78)+xin( 127)*yin( 121)*zin( 114))
          eri_value(  172)=eri_value(  172)+d12bra( 10)*d12ket( 10)*(xin(  22)*yin(  15)*zin(   1)+xin(  58)*yin(  51)*zin(  37)+xin(  94)*yin(  87)*zin(  73)+xin( 130)*yin( 123)*zin( 109))
          eri_value(  173)=eri_value(  173)+d12bra( 10)*d12ket( 11)*(xin(  21)*yin(  16)*zin(   1)+xin(  57)*yin(  52)*zin(  37)+xin(  93)*yin(  88)*zin(  73)+xin( 129)*yin( 124)*zin( 109))
          eri_value(  174)=eri_value(  174)+d12bra( 10)*d12ket( 12)*(xin(  21)*yin(  15)*zin(   2)+xin(  57)*yin(  51)*zin(  38)+xin(  93)*yin(  87)*zin(  74)+xin( 129)*yin( 123)*zin( 110))
          eri_value(  175)=eri_value(  175)+d12bra( 10)*d12ket( 13)*(xin(  22)*yin(  13)*zin(   3)+xin(  58)*yin(  49)*zin(  39)+xin(  94)*yin(  85)*zin(  75)+xin( 130)*yin( 121)*zin( 111))
          eri_value(  176)=eri_value(  176)+d12bra( 10)*d12ket( 14)*(xin(  21)*yin(  14)*zin(   3)+xin(  57)*yin(  50)*zin(  39)+xin(  93)*yin(  86)*zin(  75)+xin( 129)*yin( 122)*zin( 111))
          eri_value(  177)=eri_value(  177)+d12bra( 10)*d12ket( 15)*(xin(  21)*yin(  13)*zin(   4)+xin(  57)*yin(  49)*zin(  40)+xin(  93)*yin(  85)*zin(  76)+xin( 129)*yin( 121)*zin( 112))
          eri_value(  178)=eri_value(  178)+d12bra( 10)*d12ket( 16)*(xin(  20)*yin(  15)*zin(   3)+xin(  56)*yin(  51)*zin(  39)+xin(  92)*yin(  87)*zin(  75)+xin( 128)*yin( 123)*zin( 111))
          eri_value(  179)=eri_value(  179)+d12bra( 10)*d12ket( 17)*(xin(  19)*yin(  16)*zin(   3)+xin(  55)*yin(  52)*zin(  39)+xin(  91)*yin(  88)*zin(  75)+xin( 127)*yin( 124)*zin( 111))
          eri_value(  180)=eri_value(  180)+d12bra( 10)*d12ket( 18)*(xin(  19)*yin(  15)*zin(   4)+xin(  55)*yin(  51)*zin(  40)+xin(  91)*yin(  87)*zin(  76)+xin( 127)*yin( 123)*zin( 112))
          eri_value(  181)=eri_value(  181)+d12bra( 11)*d12ket(  1)*(xin(  18)*yin(  19)*zin(   1)+xin(  54)*yin(  55)*zin(  37)+xin(  90)*yin(  91)*zin(  73)+xin( 126)*yin( 127)*zin( 109))
          eri_value(  182)=eri_value(  182)+d12bra( 11)*d12ket(  2)*(xin(  17)*yin(  20)*zin(   1)+xin(  53)*yin(  56)*zin(  37)+xin(  89)*yin(  92)*zin(  73)+xin( 125)*yin( 128)*zin( 109))
          eri_value(  183)=eri_value(  183)+d12bra( 11)*d12ket(  3)*(xin(  17)*yin(  19)*zin(   2)+xin(  53)*yin(  55)*zin(  38)+xin(  89)*yin(  91)*zin(  74)+xin( 125)*yin( 127)*zin( 110))
          eri_value(  184)=eri_value(  184)+d12bra( 11)*d12ket(  4)*(xin(  14)*yin(  23)*zin(   1)+xin(  50)*yin(  59)*zin(  37)+xin(  86)*yin(  95)*zin(  73)+xin( 122)*yin( 131)*zin( 109))
          eri_value(  185)=eri_value(  185)+d12bra( 11)*d12ket(  5)*(xin(  13)*yin(  24)*zin(   1)+xin(  49)*yin(  60)*zin(  37)+xin(  85)*yin(  96)*zin(  73)+xin( 121)*yin( 132)*zin( 109))
          eri_value(  186)=eri_value(  186)+d12bra( 11)*d12ket(  6)*(xin(  13)*yin(  23)*zin(   2)+xin(  49)*yin(  59)*zin(  38)+xin(  85)*yin(  95)*zin(  74)+xin( 121)*yin( 131)*zin( 110))
          eri_value(  187)=eri_value(  187)+d12bra( 11)*d12ket(  7)*(xin(  14)*yin(  19)*zin(   5)+xin(  50)*yin(  55)*zin(  41)+xin(  86)*yin(  91)*zin(  77)+xin( 122)*yin( 127)*zin( 113))
          eri_value(  188)=eri_value(  188)+d12bra( 11)*d12ket(  8)*(xin(  13)*yin(  20)*zin(   5)+xin(  49)*yin(  56)*zin(  41)+xin(  85)*yin(  92)*zin(  77)+xin( 121)*yin( 128)*zin( 113))
          eri_value(  189)=eri_value(  189)+d12bra( 11)*d12ket(  9)*(xin(  13)*yin(  19)*zin(   6)+xin(  49)*yin(  55)*zin(  42)+xin(  85)*yin(  91)*zin(  78)+xin( 121)*yin( 127)*zin( 114))
          eri_value(  190)=eri_value(  190)+d12bra( 11)*d12ket( 10)*(xin(  16)*yin(  21)*zin(   1)+xin(  52)*yin(  57)*zin(  37)+xin(  88)*yin(  93)*zin(  73)+xin( 124)*yin( 129)*zin( 109))
          eri_value(  191)=eri_value(  191)+d12bra( 11)*d12ket( 11)*(xin(  15)*yin(  22)*zin(   1)+xin(  51)*yin(  58)*zin(  37)+xin(  87)*yin(  94)*zin(  73)+xin( 123)*yin( 130)*zin( 109))
          eri_value(  192)=eri_value(  192)+d12bra( 11)*d12ket( 12)*(xin(  15)*yin(  21)*zin(   2)+xin(  51)*yin(  57)*zin(  38)+xin(  87)*yin(  93)*zin(  74)+xin( 123)*yin( 129)*zin( 110))
          eri_value(  193)=eri_value(  193)+d12bra( 11)*d12ket( 13)*(xin(  16)*yin(  19)*zin(   3)+xin(  52)*yin(  55)*zin(  39)+xin(  88)*yin(  91)*zin(  75)+xin( 124)*yin( 127)*zin( 111))
          eri_value(  194)=eri_value(  194)+d12bra( 11)*d12ket( 14)*(xin(  15)*yin(  20)*zin(   3)+xin(  51)*yin(  56)*zin(  39)+xin(  87)*yin(  92)*zin(  75)+xin( 123)*yin( 128)*zin( 111))
          eri_value(  195)=eri_value(  195)+d12bra( 11)*d12ket( 15)*(xin(  15)*yin(  19)*zin(   4)+xin(  51)*yin(  55)*zin(  40)+xin(  87)*yin(  91)*zin(  76)+xin( 123)*yin( 127)*zin( 112))
          eri_value(  196)=eri_value(  196)+d12bra( 11)*d12ket( 16)*(xin(  14)*yin(  21)*zin(   3)+xin(  50)*yin(  57)*zin(  39)+xin(  86)*yin(  93)*zin(  75)+xin( 122)*yin( 129)*zin( 111))
          eri_value(  197)=eri_value(  197)+d12bra( 11)*d12ket( 17)*(xin(  13)*yin(  22)*zin(   3)+xin(  49)*yin(  58)*zin(  39)+xin(  85)*yin(  94)*zin(  75)+xin( 121)*yin( 130)*zin( 111))
          eri_value(  198)=eri_value(  198)+d12bra( 11)*d12ket( 18)*(xin(  13)*yin(  21)*zin(   4)+xin(  49)*yin(  57)*zin(  40)+xin(  85)*yin(  93)*zin(  76)+xin( 121)*yin( 129)*zin( 112))
          eri_value(  199)=eri_value(  199)+d12bra( 12)*d12ket(  1)*(xin(  18)*yin(  13)*zin(   7)+xin(  54)*yin(  49)*zin(  43)+xin(  90)*yin(  85)*zin(  79)+xin( 126)*yin( 121)*zin( 115))
          eri_value(  200)=eri_value(  200)+d12bra( 12)*d12ket(  2)*(xin(  17)*yin(  14)*zin(   7)+xin(  53)*yin(  50)*zin(  43)+xin(  89)*yin(  86)*zin(  79)+xin( 125)*yin( 122)*zin( 115))
          eri_value(  201)=eri_value(  201)+d12bra( 12)*d12ket(  3)*(xin(  17)*yin(  13)*zin(   8)+xin(  53)*yin(  49)*zin(  44)+xin(  89)*yin(  85)*zin(  80)+xin( 125)*yin( 121)*zin( 116))
          eri_value(  202)=eri_value(  202)+d12bra( 12)*d12ket(  4)*(xin(  14)*yin(  17)*zin(   7)+xin(  50)*yin(  53)*zin(  43)+xin(  86)*yin(  89)*zin(  79)+xin( 122)*yin( 125)*zin( 115))
          eri_value(  203)=eri_value(  203)+d12bra( 12)*d12ket(  5)*(xin(  13)*yin(  18)*zin(   7)+xin(  49)*yin(  54)*zin(  43)+xin(  85)*yin(  90)*zin(  79)+xin( 121)*yin( 126)*zin( 115))
          eri_value(  204)=eri_value(  204)+d12bra( 12)*d12ket(  6)*(xin(  13)*yin(  17)*zin(   8)+xin(  49)*yin(  53)*zin(  44)+xin(  85)*yin(  89)*zin(  80)+xin( 121)*yin( 125)*zin( 116))
          eri_value(  205)=eri_value(  205)+d12bra( 12)*d12ket(  7)*(xin(  14)*yin(  13)*zin(  11)+xin(  50)*yin(  49)*zin(  47)+xin(  86)*yin(  85)*zin(  83)+xin( 122)*yin( 121)*zin( 119))
          eri_value(  206)=eri_value(  206)+d12bra( 12)*d12ket(  8)*(xin(  13)*yin(  14)*zin(  11)+xin(  49)*yin(  50)*zin(  47)+xin(  85)*yin(  86)*zin(  83)+xin( 121)*yin( 122)*zin( 119))
          eri_value(  207)=eri_value(  207)+d12bra( 12)*d12ket(  9)*(xin(  13)*yin(  13)*zin(  12)+xin(  49)*yin(  49)*zin(  48)+xin(  85)*yin(  85)*zin(  84)+xin( 121)*yin( 121)*zin( 120))
          eri_value(  208)=eri_value(  208)+d12bra( 12)*d12ket( 10)*(xin(  16)*yin(  15)*zin(   7)+xin(  52)*yin(  51)*zin(  43)+xin(  88)*yin(  87)*zin(  79)+xin( 124)*yin( 123)*zin( 115))
          eri_value(  209)=eri_value(  209)+d12bra( 12)*d12ket( 11)*(xin(  15)*yin(  16)*zin(   7)+xin(  51)*yin(  52)*zin(  43)+xin(  87)*yin(  88)*zin(  79)+xin( 123)*yin( 124)*zin( 115))
          eri_value(  210)=eri_value(  210)+d12bra( 12)*d12ket( 12)*(xin(  15)*yin(  15)*zin(   8)+xin(  51)*yin(  51)*zin(  44)+xin(  87)*yin(  87)*zin(  80)+xin( 123)*yin( 123)*zin( 116))
          eri_value(  211)=eri_value(  211)+d12bra( 12)*d12ket( 13)*(xin(  16)*yin(  13)*zin(   9)+xin(  52)*yin(  49)*zin(  45)+xin(  88)*yin(  85)*zin(  81)+xin( 124)*yin( 121)*zin( 117))
          eri_value(  212)=eri_value(  212)+d12bra( 12)*d12ket( 14)*(xin(  15)*yin(  14)*zin(   9)+xin(  51)*yin(  50)*zin(  45)+xin(  87)*yin(  86)*zin(  81)+xin( 123)*yin( 122)*zin( 117))
          eri_value(  213)=eri_value(  213)+d12bra( 12)*d12ket( 15)*(xin(  15)*yin(  13)*zin(  10)+xin(  51)*yin(  49)*zin(  46)+xin(  87)*yin(  85)*zin(  82)+xin( 123)*yin( 121)*zin( 118))
          eri_value(  214)=eri_value(  214)+d12bra( 12)*d12ket( 16)*(xin(  14)*yin(  15)*zin(   9)+xin(  50)*yin(  51)*zin(  45)+xin(  86)*yin(  87)*zin(  81)+xin( 122)*yin( 123)*zin( 117))
          eri_value(  215)=eri_value(  215)+d12bra( 12)*d12ket( 17)*(xin(  13)*yin(  16)*zin(   9)+xin(  49)*yin(  52)*zin(  45)+xin(  85)*yin(  88)*zin(  81)+xin( 121)*yin( 124)*zin( 117))
          eri_value(  216)=eri_value(  216)+d12bra( 12)*d12ket( 18)*(xin(  13)*yin(  15)*zin(  10)+xin(  49)*yin(  51)*zin(  46)+xin(  85)*yin(  87)*zin(  82)+xin( 121)*yin( 123)*zin( 118))
          eri_value(  217)=eri_value(  217)+d12bra( 13)*d12ket(  1)*(xin(  24)*yin(   1)*zin(  13)+xin(  60)*yin(  37)*zin(  49)+xin(  96)*yin(  73)*zin(  85)+xin( 132)*yin( 109)*zin( 121))
          eri_value(  218)=eri_value(  218)+d12bra( 13)*d12ket(  2)*(xin(  23)*yin(   2)*zin(  13)+xin(  59)*yin(  38)*zin(  49)+xin(  95)*yin(  74)*zin(  85)+xin( 131)*yin( 110)*zin( 121))
          eri_value(  219)=eri_value(  219)+d12bra( 13)*d12ket(  3)*(xin(  23)*yin(   1)*zin(  14)+xin(  59)*yin(  37)*zin(  50)+xin(  95)*yin(  73)*zin(  86)+xin( 131)*yin( 109)*zin( 122))
          eri_value(  220)=eri_value(  220)+d12bra( 13)*d12ket(  4)*(xin(  20)*yin(   5)*zin(  13)+xin(  56)*yin(  41)*zin(  49)+xin(  92)*yin(  77)*zin(  85)+xin( 128)*yin( 113)*zin( 121))
          eri_value(  221)=eri_value(  221)+d12bra( 13)*d12ket(  5)*(xin(  19)*yin(   6)*zin(  13)+xin(  55)*yin(  42)*zin(  49)+xin(  91)*yin(  78)*zin(  85)+xin( 127)*yin( 114)*zin( 121))
          eri_value(  222)=eri_value(  222)+d12bra( 13)*d12ket(  6)*(xin(  19)*yin(   5)*zin(  14)+xin(  55)*yin(  41)*zin(  50)+xin(  91)*yin(  77)*zin(  86)+xin( 127)*yin( 113)*zin( 122))
          eri_value(  223)=eri_value(  223)+d12bra( 13)*d12ket(  7)*(xin(  20)*yin(   1)*zin(  17)+xin(  56)*yin(  37)*zin(  53)+xin(  92)*yin(  73)*zin(  89)+xin( 128)*yin( 109)*zin( 125))
          eri_value(  224)=eri_value(  224)+d12bra( 13)*d12ket(  8)*(xin(  19)*yin(   2)*zin(  17)+xin(  55)*yin(  38)*zin(  53)+xin(  91)*yin(  74)*zin(  89)+xin( 127)*yin( 110)*zin( 125))
          eri_value(  225)=eri_value(  225)+d12bra( 13)*d12ket(  9)*(xin(  19)*yin(   1)*zin(  18)+xin(  55)*yin(  37)*zin(  54)+xin(  91)*yin(  73)*zin(  90)+xin( 127)*yin( 109)*zin( 126))
          eri_value(  226)=eri_value(  226)+d12bra( 13)*d12ket( 10)*(xin(  22)*yin(   3)*zin(  13)+xin(  58)*yin(  39)*zin(  49)+xin(  94)*yin(  75)*zin(  85)+xin( 130)*yin( 111)*zin( 121))
          eri_value(  227)=eri_value(  227)+d12bra( 13)*d12ket( 11)*(xin(  21)*yin(   4)*zin(  13)+xin(  57)*yin(  40)*zin(  49)+xin(  93)*yin(  76)*zin(  85)+xin( 129)*yin( 112)*zin( 121))
          eri_value(  228)=eri_value(  228)+d12bra( 13)*d12ket( 12)*(xin(  21)*yin(   3)*zin(  14)+xin(  57)*yin(  39)*zin(  50)+xin(  93)*yin(  75)*zin(  86)+xin( 129)*yin( 111)*zin( 122))
          eri_value(  229)=eri_value(  229)+d12bra( 13)*d12ket( 13)*(xin(  22)*yin(   1)*zin(  15)+xin(  58)*yin(  37)*zin(  51)+xin(  94)*yin(  73)*zin(  87)+xin( 130)*yin( 109)*zin( 123))
          eri_value(  230)=eri_value(  230)+d12bra( 13)*d12ket( 14)*(xin(  21)*yin(   2)*zin(  15)+xin(  57)*yin(  38)*zin(  51)+xin(  93)*yin(  74)*zin(  87)+xin( 129)*yin( 110)*zin( 123))
          eri_value(  231)=eri_value(  231)+d12bra( 13)*d12ket( 15)*(xin(  21)*yin(   1)*zin(  16)+xin(  57)*yin(  37)*zin(  52)+xin(  93)*yin(  73)*zin(  88)+xin( 129)*yin( 109)*zin( 124))
          eri_value(  232)=eri_value(  232)+d12bra( 13)*d12ket( 16)*(xin(  20)*yin(   3)*zin(  15)+xin(  56)*yin(  39)*zin(  51)+xin(  92)*yin(  75)*zin(  87)+xin( 128)*yin( 111)*zin( 123))
          eri_value(  233)=eri_value(  233)+d12bra( 13)*d12ket( 17)*(xin(  19)*yin(   4)*zin(  15)+xin(  55)*yin(  40)*zin(  51)+xin(  91)*yin(  76)*zin(  87)+xin( 127)*yin( 112)*zin( 123))
          eri_value(  234)=eri_value(  234)+d12bra( 13)*d12ket( 18)*(xin(  19)*yin(   3)*zin(  16)+xin(  55)*yin(  39)*zin(  52)+xin(  91)*yin(  75)*zin(  88)+xin( 127)*yin( 111)*zin( 124))
          eri_value(  235)=eri_value(  235)+d12bra( 14)*d12ket(  1)*(xin(  18)*yin(   7)*zin(  13)+xin(  54)*yin(  43)*zin(  49)+xin(  90)*yin(  79)*zin(  85)+xin( 126)*yin( 115)*zin( 121))
          eri_value(  236)=eri_value(  236)+d12bra( 14)*d12ket(  2)*(xin(  17)*yin(   8)*zin(  13)+xin(  53)*yin(  44)*zin(  49)+xin(  89)*yin(  80)*zin(  85)+xin( 125)*yin( 116)*zin( 121))
          eri_value(  237)=eri_value(  237)+d12bra( 14)*d12ket(  3)*(xin(  17)*yin(   7)*zin(  14)+xin(  53)*yin(  43)*zin(  50)+xin(  89)*yin(  79)*zin(  86)+xin( 125)*yin( 115)*zin( 122))
          eri_value(  238)=eri_value(  238)+d12bra( 14)*d12ket(  4)*(xin(  14)*yin(  11)*zin(  13)+xin(  50)*yin(  47)*zin(  49)+xin(  86)*yin(  83)*zin(  85)+xin( 122)*yin( 119)*zin( 121))
          eri_value(  239)=eri_value(  239)+d12bra( 14)*d12ket(  5)*(xin(  13)*yin(  12)*zin(  13)+xin(  49)*yin(  48)*zin(  49)+xin(  85)*yin(  84)*zin(  85)+xin( 121)*yin( 120)*zin( 121))
          eri_value(  240)=eri_value(  240)+d12bra( 14)*d12ket(  6)*(xin(  13)*yin(  11)*zin(  14)+xin(  49)*yin(  47)*zin(  50)+xin(  85)*yin(  83)*zin(  86)+xin( 121)*yin( 119)*zin( 122))
          eri_value(  241)=eri_value(  241)+d12bra( 14)*d12ket(  7)*(xin(  14)*yin(   7)*zin(  17)+xin(  50)*yin(  43)*zin(  53)+xin(  86)*yin(  79)*zin(  89)+xin( 122)*yin( 115)*zin( 125))
          eri_value(  242)=eri_value(  242)+d12bra( 14)*d12ket(  8)*(xin(  13)*yin(   8)*zin(  17)+xin(  49)*yin(  44)*zin(  53)+xin(  85)*yin(  80)*zin(  89)+xin( 121)*yin( 116)*zin( 125))
          eri_value(  243)=eri_value(  243)+d12bra( 14)*d12ket(  9)*(xin(  13)*yin(   7)*zin(  18)+xin(  49)*yin(  43)*zin(  54)+xin(  85)*yin(  79)*zin(  90)+xin( 121)*yin( 115)*zin( 126))
          eri_value(  244)=eri_value(  244)+d12bra( 14)*d12ket( 10)*(xin(  16)*yin(   9)*zin(  13)+xin(  52)*yin(  45)*zin(  49)+xin(  88)*yin(  81)*zin(  85)+xin( 124)*yin( 117)*zin( 121))
          eri_value(  245)=eri_value(  245)+d12bra( 14)*d12ket( 11)*(xin(  15)*yin(  10)*zin(  13)+xin(  51)*yin(  46)*zin(  49)+xin(  87)*yin(  82)*zin(  85)+xin( 123)*yin( 118)*zin( 121))
          eri_value(  246)=eri_value(  246)+d12bra( 14)*d12ket( 12)*(xin(  15)*yin(   9)*zin(  14)+xin(  51)*yin(  45)*zin(  50)+xin(  87)*yin(  81)*zin(  86)+xin( 123)*yin( 117)*zin( 122))
          eri_value(  247)=eri_value(  247)+d12bra( 14)*d12ket( 13)*(xin(  16)*yin(   7)*zin(  15)+xin(  52)*yin(  43)*zin(  51)+xin(  88)*yin(  79)*zin(  87)+xin( 124)*yin( 115)*zin( 123))
          eri_value(  248)=eri_value(  248)+d12bra( 14)*d12ket( 14)*(xin(  15)*yin(   8)*zin(  15)+xin(  51)*yin(  44)*zin(  51)+xin(  87)*yin(  80)*zin(  87)+xin( 123)*yin( 116)*zin( 123))
          eri_value(  249)=eri_value(  249)+d12bra( 14)*d12ket( 15)*(xin(  15)*yin(   7)*zin(  16)+xin(  51)*yin(  43)*zin(  52)+xin(  87)*yin(  79)*zin(  88)+xin( 123)*yin( 115)*zin( 124))
          eri_value(  250)=eri_value(  250)+d12bra( 14)*d12ket( 16)*(xin(  14)*yin(   9)*zin(  15)+xin(  50)*yin(  45)*zin(  51)+xin(  86)*yin(  81)*zin(  87)+xin( 122)*yin( 117)*zin( 123))
          eri_value(  251)=eri_value(  251)+d12bra( 14)*d12ket( 17)*(xin(  13)*yin(  10)*zin(  15)+xin(  49)*yin(  46)*zin(  51)+xin(  85)*yin(  82)*zin(  87)+xin( 121)*yin( 118)*zin( 123))
          eri_value(  252)=eri_value(  252)+d12bra( 14)*d12ket( 18)*(xin(  13)*yin(   9)*zin(  16)+xin(  49)*yin(  45)*zin(  52)+xin(  85)*yin(  81)*zin(  88)+xin( 121)*yin( 117)*zin( 124))
          eri_value(  253)=eri_value(  253)+d12bra( 15)*d12ket(  1)*(xin(  18)*yin(   1)*zin(  19)+xin(  54)*yin(  37)*zin(  55)+xin(  90)*yin(  73)*zin(  91)+xin( 126)*yin( 109)*zin( 127))
          eri_value(  254)=eri_value(  254)+d12bra( 15)*d12ket(  2)*(xin(  17)*yin(   2)*zin(  19)+xin(  53)*yin(  38)*zin(  55)+xin(  89)*yin(  74)*zin(  91)+xin( 125)*yin( 110)*zin( 127))
          eri_value(  255)=eri_value(  255)+d12bra( 15)*d12ket(  3)*(xin(  17)*yin(   1)*zin(  20)+xin(  53)*yin(  37)*zin(  56)+xin(  89)*yin(  73)*zin(  92)+xin( 125)*yin( 109)*zin( 128))
          eri_value(  256)=eri_value(  256)+d12bra( 15)*d12ket(  4)*(xin(  14)*yin(   5)*zin(  19)+xin(  50)*yin(  41)*zin(  55)+xin(  86)*yin(  77)*zin(  91)+xin( 122)*yin( 113)*zin( 127))
          eri_value(  257)=eri_value(  257)+d12bra( 15)*d12ket(  5)*(xin(  13)*yin(   6)*zin(  19)+xin(  49)*yin(  42)*zin(  55)+xin(  85)*yin(  78)*zin(  91)+xin( 121)*yin( 114)*zin( 127))
          eri_value(  258)=eri_value(  258)+d12bra( 15)*d12ket(  6)*(xin(  13)*yin(   5)*zin(  20)+xin(  49)*yin(  41)*zin(  56)+xin(  85)*yin(  77)*zin(  92)+xin( 121)*yin( 113)*zin( 128))
          eri_value(  259)=eri_value(  259)+d12bra( 15)*d12ket(  7)*(xin(  14)*yin(   1)*zin(  23)+xin(  50)*yin(  37)*zin(  59)+xin(  86)*yin(  73)*zin(  95)+xin( 122)*yin( 109)*zin( 131))
          eri_value(  260)=eri_value(  260)+d12bra( 15)*d12ket(  8)*(xin(  13)*yin(   2)*zin(  23)+xin(  49)*yin(  38)*zin(  59)+xin(  85)*yin(  74)*zin(  95)+xin( 121)*yin( 110)*zin( 131))
          eri_value(  261)=eri_value(  261)+d12bra( 15)*d12ket(  9)*(xin(  13)*yin(   1)*zin(  24)+xin(  49)*yin(  37)*zin(  60)+xin(  85)*yin(  73)*zin(  96)+xin( 121)*yin( 109)*zin( 132))
          eri_value(  262)=eri_value(  262)+d12bra( 15)*d12ket( 10)*(xin(  16)*yin(   3)*zin(  19)+xin(  52)*yin(  39)*zin(  55)+xin(  88)*yin(  75)*zin(  91)+xin( 124)*yin( 111)*zin( 127))
          eri_value(  263)=eri_value(  263)+d12bra( 15)*d12ket( 11)*(xin(  15)*yin(   4)*zin(  19)+xin(  51)*yin(  40)*zin(  55)+xin(  87)*yin(  76)*zin(  91)+xin( 123)*yin( 112)*zin( 127))
          eri_value(  264)=eri_value(  264)+d12bra( 15)*d12ket( 12)*(xin(  15)*yin(   3)*zin(  20)+xin(  51)*yin(  39)*zin(  56)+xin(  87)*yin(  75)*zin(  92)+xin( 123)*yin( 111)*zin( 128))
          eri_value(  265)=eri_value(  265)+d12bra( 15)*d12ket( 13)*(xin(  16)*yin(   1)*zin(  21)+xin(  52)*yin(  37)*zin(  57)+xin(  88)*yin(  73)*zin(  93)+xin( 124)*yin( 109)*zin( 129))
          eri_value(  266)=eri_value(  266)+d12bra( 15)*d12ket( 14)*(xin(  15)*yin(   2)*zin(  21)+xin(  51)*yin(  38)*zin(  57)+xin(  87)*yin(  74)*zin(  93)+xin( 123)*yin( 110)*zin( 129))
          eri_value(  267)=eri_value(  267)+d12bra( 15)*d12ket( 15)*(xin(  15)*yin(   1)*zin(  22)+xin(  51)*yin(  37)*zin(  58)+xin(  87)*yin(  73)*zin(  94)+xin( 123)*yin( 109)*zin( 130))
          eri_value(  268)=eri_value(  268)+d12bra( 15)*d12ket( 16)*(xin(  14)*yin(   3)*zin(  21)+xin(  50)*yin(  39)*zin(  57)+xin(  86)*yin(  75)*zin(  93)+xin( 122)*yin( 111)*zin( 129))
          eri_value(  269)=eri_value(  269)+d12bra( 15)*d12ket( 17)*(xin(  13)*yin(   4)*zin(  21)+xin(  49)*yin(  40)*zin(  57)+xin(  85)*yin(  76)*zin(  93)+xin( 121)*yin( 112)*zin( 129))
          eri_value(  270)=eri_value(  270)+d12bra( 15)*d12ket( 18)*(xin(  13)*yin(   3)*zin(  22)+xin(  49)*yin(  39)*zin(  58)+xin(  85)*yin(  75)*zin(  94)+xin( 121)*yin( 111)*zin( 130))
          eri_value(  271)=eri_value(  271)+d12bra( 16)*d12ket(  1)*(xin(  12)*yin(  13)*zin(  13)+xin(  48)*yin(  49)*zin(  49)+xin(  84)*yin(  85)*zin(  85)+xin( 120)*yin( 121)*zin( 121))
          eri_value(  272)=eri_value(  272)+d12bra( 16)*d12ket(  2)*(xin(  11)*yin(  14)*zin(  13)+xin(  47)*yin(  50)*zin(  49)+xin(  83)*yin(  86)*zin(  85)+xin( 119)*yin( 122)*zin( 121))
          eri_value(  273)=eri_value(  273)+d12bra( 16)*d12ket(  3)*(xin(  11)*yin(  13)*zin(  14)+xin(  47)*yin(  49)*zin(  50)+xin(  83)*yin(  85)*zin(  86)+xin( 119)*yin( 121)*zin( 122))
          eri_value(  274)=eri_value(  274)+d12bra( 16)*d12ket(  4)*(xin(   8)*yin(  17)*zin(  13)+xin(  44)*yin(  53)*zin(  49)+xin(  80)*yin(  89)*zin(  85)+xin( 116)*yin( 125)*zin( 121))
          eri_value(  275)=eri_value(  275)+d12bra( 16)*d12ket(  5)*(xin(   7)*yin(  18)*zin(  13)+xin(  43)*yin(  54)*zin(  49)+xin(  79)*yin(  90)*zin(  85)+xin( 115)*yin( 126)*zin( 121))
          eri_value(  276)=eri_value(  276)+d12bra( 16)*d12ket(  6)*(xin(   7)*yin(  17)*zin(  14)+xin(  43)*yin(  53)*zin(  50)+xin(  79)*yin(  89)*zin(  86)+xin( 115)*yin( 125)*zin( 122))
          eri_value(  277)=eri_value(  277)+d12bra( 16)*d12ket(  7)*(xin(   8)*yin(  13)*zin(  17)+xin(  44)*yin(  49)*zin(  53)+xin(  80)*yin(  85)*zin(  89)+xin( 116)*yin( 121)*zin( 125))
          eri_value(  278)=eri_value(  278)+d12bra( 16)*d12ket(  8)*(xin(   7)*yin(  14)*zin(  17)+xin(  43)*yin(  50)*zin(  53)+xin(  79)*yin(  86)*zin(  89)+xin( 115)*yin( 122)*zin( 125))
          eri_value(  279)=eri_value(  279)+d12bra( 16)*d12ket(  9)*(xin(   7)*yin(  13)*zin(  18)+xin(  43)*yin(  49)*zin(  54)+xin(  79)*yin(  85)*zin(  90)+xin( 115)*yin( 121)*zin( 126))
          eri_value(  280)=eri_value(  280)+d12bra( 16)*d12ket( 10)*(xin(  10)*yin(  15)*zin(  13)+xin(  46)*yin(  51)*zin(  49)+xin(  82)*yin(  87)*zin(  85)+xin( 118)*yin( 123)*zin( 121))
          eri_value(  281)=eri_value(  281)+d12bra( 16)*d12ket( 11)*(xin(   9)*yin(  16)*zin(  13)+xin(  45)*yin(  52)*zin(  49)+xin(  81)*yin(  88)*zin(  85)+xin( 117)*yin( 124)*zin( 121))
          eri_value(  282)=eri_value(  282)+d12bra( 16)*d12ket( 12)*(xin(   9)*yin(  15)*zin(  14)+xin(  45)*yin(  51)*zin(  50)+xin(  81)*yin(  87)*zin(  86)+xin( 117)*yin( 123)*zin( 122))
          eri_value(  283)=eri_value(  283)+d12bra( 16)*d12ket( 13)*(xin(  10)*yin(  13)*zin(  15)+xin(  46)*yin(  49)*zin(  51)+xin(  82)*yin(  85)*zin(  87)+xin( 118)*yin( 121)*zin( 123))
          eri_value(  284)=eri_value(  284)+d12bra( 16)*d12ket( 14)*(xin(   9)*yin(  14)*zin(  15)+xin(  45)*yin(  50)*zin(  51)+xin(  81)*yin(  86)*zin(  87)+xin( 117)*yin( 122)*zin( 123))
          eri_value(  285)=eri_value(  285)+d12bra( 16)*d12ket( 15)*(xin(   9)*yin(  13)*zin(  16)+xin(  45)*yin(  49)*zin(  52)+xin(  81)*yin(  85)*zin(  88)+xin( 117)*yin( 121)*zin( 124))
          eri_value(  286)=eri_value(  286)+d12bra( 16)*d12ket( 16)*(xin(   8)*yin(  15)*zin(  15)+xin(  44)*yin(  51)*zin(  51)+xin(  80)*yin(  87)*zin(  87)+xin( 116)*yin( 123)*zin( 123))
          eri_value(  287)=eri_value(  287)+d12bra( 16)*d12ket( 17)*(xin(   7)*yin(  16)*zin(  15)+xin(  43)*yin(  52)*zin(  51)+xin(  79)*yin(  88)*zin(  87)+xin( 115)*yin( 124)*zin( 123))
          eri_value(  288)=eri_value(  288)+d12bra( 16)*d12ket( 18)*(xin(   7)*yin(  15)*zin(  16)+xin(  43)*yin(  51)*zin(  52)+xin(  79)*yin(  87)*zin(  88)+xin( 115)*yin( 123)*zin( 124))
          eri_value(  289)=eri_value(  289)+d12bra( 17)*d12ket(  1)*(xin(   6)*yin(  19)*zin(  13)+xin(  42)*yin(  55)*zin(  49)+xin(  78)*yin(  91)*zin(  85)+xin( 114)*yin( 127)*zin( 121))
          eri_value(  290)=eri_value(  290)+d12bra( 17)*d12ket(  2)*(xin(   5)*yin(  20)*zin(  13)+xin(  41)*yin(  56)*zin(  49)+xin(  77)*yin(  92)*zin(  85)+xin( 113)*yin( 128)*zin( 121))
          eri_value(  291)=eri_value(  291)+d12bra( 17)*d12ket(  3)*(xin(   5)*yin(  19)*zin(  14)+xin(  41)*yin(  55)*zin(  50)+xin(  77)*yin(  91)*zin(  86)+xin( 113)*yin( 127)*zin( 122))
          eri_value(  292)=eri_value(  292)+d12bra( 17)*d12ket(  4)*(xin(   2)*yin(  23)*zin(  13)+xin(  38)*yin(  59)*zin(  49)+xin(  74)*yin(  95)*zin(  85)+xin( 110)*yin( 131)*zin( 121))
          eri_value(  293)=eri_value(  293)+d12bra( 17)*d12ket(  5)*(xin(   1)*yin(  24)*zin(  13)+xin(  37)*yin(  60)*zin(  49)+xin(  73)*yin(  96)*zin(  85)+xin( 109)*yin( 132)*zin( 121))
          eri_value(  294)=eri_value(  294)+d12bra( 17)*d12ket(  6)*(xin(   1)*yin(  23)*zin(  14)+xin(  37)*yin(  59)*zin(  50)+xin(  73)*yin(  95)*zin(  86)+xin( 109)*yin( 131)*zin( 122))
          eri_value(  295)=eri_value(  295)+d12bra( 17)*d12ket(  7)*(xin(   2)*yin(  19)*zin(  17)+xin(  38)*yin(  55)*zin(  53)+xin(  74)*yin(  91)*zin(  89)+xin( 110)*yin( 127)*zin( 125))
          eri_value(  296)=eri_value(  296)+d12bra( 17)*d12ket(  8)*(xin(   1)*yin(  20)*zin(  17)+xin(  37)*yin(  56)*zin(  53)+xin(  73)*yin(  92)*zin(  89)+xin( 109)*yin( 128)*zin( 125))
          eri_value(  297)=eri_value(  297)+d12bra( 17)*d12ket(  9)*(xin(   1)*yin(  19)*zin(  18)+xin(  37)*yin(  55)*zin(  54)+xin(  73)*yin(  91)*zin(  90)+xin( 109)*yin( 127)*zin( 126))
          eri_value(  298)=eri_value(  298)+d12bra( 17)*d12ket( 10)*(xin(   4)*yin(  21)*zin(  13)+xin(  40)*yin(  57)*zin(  49)+xin(  76)*yin(  93)*zin(  85)+xin( 112)*yin( 129)*zin( 121))
          eri_value(  299)=eri_value(  299)+d12bra( 17)*d12ket( 11)*(xin(   3)*yin(  22)*zin(  13)+xin(  39)*yin(  58)*zin(  49)+xin(  75)*yin(  94)*zin(  85)+xin( 111)*yin( 130)*zin( 121))
          eri_value(  300)=eri_value(  300)+d12bra( 17)*d12ket( 12)*(xin(   3)*yin(  21)*zin(  14)+xin(  39)*yin(  57)*zin(  50)+xin(  75)*yin(  93)*zin(  86)+xin( 111)*yin( 129)*zin( 122))
          eri_value(  301)=eri_value(  301)+d12bra( 17)*d12ket( 13)*(xin(   4)*yin(  19)*zin(  15)+xin(  40)*yin(  55)*zin(  51)+xin(  76)*yin(  91)*zin(  87)+xin( 112)*yin( 127)*zin( 123))
          eri_value(  302)=eri_value(  302)+d12bra( 17)*d12ket( 14)*(xin(   3)*yin(  20)*zin(  15)+xin(  39)*yin(  56)*zin(  51)+xin(  75)*yin(  92)*zin(  87)+xin( 111)*yin( 128)*zin( 123))
          eri_value(  303)=eri_value(  303)+d12bra( 17)*d12ket( 15)*(xin(   3)*yin(  19)*zin(  16)+xin(  39)*yin(  55)*zin(  52)+xin(  75)*yin(  91)*zin(  88)+xin( 111)*yin( 127)*zin( 124))
          eri_value(  304)=eri_value(  304)+d12bra( 17)*d12ket( 16)*(xin(   2)*yin(  21)*zin(  15)+xin(  38)*yin(  57)*zin(  51)+xin(  74)*yin(  93)*zin(  87)+xin( 110)*yin( 129)*zin( 123))
          eri_value(  305)=eri_value(  305)+d12bra( 17)*d12ket( 17)*(xin(   1)*yin(  22)*zin(  15)+xin(  37)*yin(  58)*zin(  51)+xin(  73)*yin(  94)*zin(  87)+xin( 109)*yin( 130)*zin( 123))
          eri_value(  306)=eri_value(  306)+d12bra( 17)*d12ket( 18)*(xin(   1)*yin(  21)*zin(  16)+xin(  37)*yin(  57)*zin(  52)+xin(  73)*yin(  93)*zin(  88)+xin( 109)*yin( 129)*zin( 124))
          eri_value(  307)=eri_value(  307)+d12bra( 18)*d12ket(  1)*(xin(   6)*yin(  13)*zin(  19)+xin(  42)*yin(  49)*zin(  55)+xin(  78)*yin(  85)*zin(  91)+xin( 114)*yin( 121)*zin( 127))
          eri_value(  308)=eri_value(  308)+d12bra( 18)*d12ket(  2)*(xin(   5)*yin(  14)*zin(  19)+xin(  41)*yin(  50)*zin(  55)+xin(  77)*yin(  86)*zin(  91)+xin( 113)*yin( 122)*zin( 127))
          eri_value(  309)=eri_value(  309)+d12bra( 18)*d12ket(  3)*(xin(   5)*yin(  13)*zin(  20)+xin(  41)*yin(  49)*zin(  56)+xin(  77)*yin(  85)*zin(  92)+xin( 113)*yin( 121)*zin( 128))
          eri_value(  310)=eri_value(  310)+d12bra( 18)*d12ket(  4)*(xin(   2)*yin(  17)*zin(  19)+xin(  38)*yin(  53)*zin(  55)+xin(  74)*yin(  89)*zin(  91)+xin( 110)*yin( 125)*zin( 127))
          eri_value(  311)=eri_value(  311)+d12bra( 18)*d12ket(  5)*(xin(   1)*yin(  18)*zin(  19)+xin(  37)*yin(  54)*zin(  55)+xin(  73)*yin(  90)*zin(  91)+xin( 109)*yin( 126)*zin( 127))
          eri_value(  312)=eri_value(  312)+d12bra( 18)*d12ket(  6)*(xin(   1)*yin(  17)*zin(  20)+xin(  37)*yin(  53)*zin(  56)+xin(  73)*yin(  89)*zin(  92)+xin( 109)*yin( 125)*zin( 128))
          eri_value(  313)=eri_value(  313)+d12bra( 18)*d12ket(  7)*(xin(   2)*yin(  13)*zin(  23)+xin(  38)*yin(  49)*zin(  59)+xin(  74)*yin(  85)*zin(  95)+xin( 110)*yin( 121)*zin( 131))
          eri_value(  314)=eri_value(  314)+d12bra( 18)*d12ket(  8)*(xin(   1)*yin(  14)*zin(  23)+xin(  37)*yin(  50)*zin(  59)+xin(  73)*yin(  86)*zin(  95)+xin( 109)*yin( 122)*zin( 131))
          eri_value(  315)=eri_value(  315)+d12bra( 18)*d12ket(  9)*(xin(   1)*yin(  13)*zin(  24)+xin(  37)*yin(  49)*zin(  60)+xin(  73)*yin(  85)*zin(  96)+xin( 109)*yin( 121)*zin( 132))
          eri_value(  316)=eri_value(  316)+d12bra( 18)*d12ket( 10)*(xin(   4)*yin(  15)*zin(  19)+xin(  40)*yin(  51)*zin(  55)+xin(  76)*yin(  87)*zin(  91)+xin( 112)*yin( 123)*zin( 127))
          eri_value(  317)=eri_value(  317)+d12bra( 18)*d12ket( 11)*(xin(   3)*yin(  16)*zin(  19)+xin(  39)*yin(  52)*zin(  55)+xin(  75)*yin(  88)*zin(  91)+xin( 111)*yin( 124)*zin( 127))
          eri_value(  318)=eri_value(  318)+d12bra( 18)*d12ket( 12)*(xin(   3)*yin(  15)*zin(  20)+xin(  39)*yin(  51)*zin(  56)+xin(  75)*yin(  87)*zin(  92)+xin( 111)*yin( 123)*zin( 128))
          eri_value(  319)=eri_value(  319)+d12bra( 18)*d12ket( 13)*(xin(   4)*yin(  13)*zin(  21)+xin(  40)*yin(  49)*zin(  57)+xin(  76)*yin(  85)*zin(  93)+xin( 112)*yin( 121)*zin( 129))
          eri_value(  320)=eri_value(  320)+d12bra( 18)*d12ket( 14)*(xin(   3)*yin(  14)*zin(  21)+xin(  39)*yin(  50)*zin(  57)+xin(  75)*yin(  86)*zin(  93)+xin( 111)*yin( 122)*zin( 129))
          eri_value(  321)=eri_value(  321)+d12bra( 18)*d12ket( 15)*(xin(   3)*yin(  13)*zin(  22)+xin(  39)*yin(  49)*zin(  58)+xin(  75)*yin(  85)*zin(  94)+xin( 111)*yin( 121)*zin( 130))
          eri_value(  322)=eri_value(  322)+d12bra( 18)*d12ket( 16)*(xin(   2)*yin(  15)*zin(  21)+xin(  38)*yin(  51)*zin(  57)+xin(  74)*yin(  87)*zin(  93)+xin( 110)*yin( 123)*zin( 129))
          eri_value(  323)=eri_value(  323)+d12bra( 18)*d12ket( 17)*(xin(   1)*yin(  16)*zin(  21)+xin(  37)*yin(  52)*zin(  57)+xin(  73)*yin(  88)*zin(  93)+xin( 109)*yin( 124)*zin( 129))
          eri_value(  324)=eri_value(  324)+d12bra( 18)*d12ket( 18)*(xin(   1)*yin(  15)*zin(  22)+xin(  37)*yin(  51)*zin(  58)+xin(  73)*yin(  87)*zin(  94)+xin( 109)*yin( 123)*zin( 130))

                                          !                     --- END FORMS ---

                                        end do ! ij primitve loop

                                      end do ! kl primitve loop

                                      !                     --- DIRFCK_RHF ---
                                      !          Compute Fock matrix elements from 2EIs

                                      maxl = 3
                                      same = (ish .eq. ksh) .and. (jsh .eq. lsh)

                                      loci = res%atom_loc(ish) - 1
                                      locj = res%atom_loc(jsh) - 1
                                      lock = res%atom_loc(ksh) - 1
                                      locl = res%atom_loc(lsh) - 1

                                      nij = 0

                                      do i = 1, 6 ! # of cartesians in i

                                        ii1 = i + loci
                                        ip = (i - 1)*54 ! Stride between functions in i

                                        do j = 1, 3 ! # of cartesians in j

                                          nij = nij + 1

                                          maxl2 = maxl

                                          jj1 = j + locj
                                          i2 = ii1
                                          j2 = jj1
                                          if (ii1 .lt. jj1) then ! Sort <ij|
                                            i2 = jj1
                                            j2 = ii1
                                          end if

                                          ijp = (j - 1)*18 + ip ! Add stride between functions in j

                                          nkl = nij

                                          do k = 1, 6 ! # of cartesians in k

                                            kk1 = k + lock

                                            ijkp = (k - 1)*3 + ijp ! Add stride between functions in k

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


                                deallocate (n12bra)
                                deallocate (xint12bra)

                                end subroutine int2121
                                end submodule
