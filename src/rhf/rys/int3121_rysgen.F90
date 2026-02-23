! The total angular momentum of this class is:           7
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3121_impl
contains
  module subroutine int3121(pf_pair, pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pf_pair, pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n13bra(:), n12ket(:)
    real(dp), allocatable :: xint13bra(:), xint12ket(:)
    integer(kind=int64) :: npfbra, npdket
    real(dp) :: scutpfbra, scutpdket, test
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
    real(dp) :: roots(4), wghts(4)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(30), wgrid(30), p0(30), p1(30), p2(30)
    real(dp) :: rts(4), wts(4), alpha(4), beta(4), wrk(4)
    real(dp) :: xin(192), yin(192), zin(192)
    real(dp) :: eri_value(540)
    real(dp) :: d13bra(30), d12ket(18)
    integer(kind=int64) :: ix(10), jx(3), kx(6), lx(3)
    integer(kind=int64) :: iy(10), jy(3), ky(6), ly(3)
    integer(kind=int64) :: iz(10), jz(3), kz(6), lz(3)
    integer(kind=int64) :: in(5), in1(5), kn(4)
    integer(kind=int64) :: ijx(30), ijy(30), ijz(30)
    integer(kind=int64) :: klx(18), kly(18), klz(18)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 37
    in1(5) = 43

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
    ky(2) = 4
    ky(3) = 0
    ky(4) = 2
    ky(5) = 0
    ky(6) = 2

    jy(1) = 0
    jy(2) = 6
    jy(3) = 0

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
    kz(3) = 4
    kz(4) = 0
    kz(5) = 2
    kz(6) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 6

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

    ijx(1) = 43
    ijx(2) = 37
    ijx(3) = 37
    ijx(4) = 7
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 7
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 31
    ijx(11) = 25
    ijx(12) = 25
    ijx(13) = 31
    ijx(14) = 25
    ijx(15) = 25
    ijx(16) = 19
    ijx(17) = 13
    ijx(18) = 13
    ijx(19) = 7
    ijx(20) = 1
    ijx(21) = 1
    ijx(22) = 19
    ijx(23) = 13
    ijx(24) = 13
    ijx(25) = 7
    ijx(26) = 1
    ijx(27) = 1
    ijx(28) = 19
    ijx(29) = 13
    ijx(30) = 13

    ijy(1) = 1
    ijy(2) = 7
    ijy(3) = 1
    ijy(4) = 37
    ijy(5) = 43
    ijy(6) = 37
    ijy(7) = 1
    ijy(8) = 7
    ijy(9) = 1
    ijy(10) = 13
    ijy(11) = 19
    ijy(12) = 13
    ijy(13) = 1
    ijy(14) = 7
    ijy(15) = 1
    ijy(16) = 25
    ijy(17) = 31
    ijy(18) = 25
    ijy(19) = 25
    ijy(20) = 31
    ijy(21) = 25
    ijy(22) = 1
    ijy(23) = 7
    ijy(24) = 1
    ijy(25) = 13
    ijy(26) = 19
    ijy(27) = 13
    ijy(28) = 13
    ijy(29) = 19
    ijy(30) = 13

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 7
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 7
    ijz(7) = 37
    ijz(8) = 37
    ijz(9) = 43
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 7
    ijz(13) = 13
    ijz(14) = 13
    ijz(15) = 19
    ijz(16) = 1
    ijz(17) = 1
    ijz(18) = 7
    ijz(19) = 13
    ijz(20) = 13
    ijz(21) = 19
    ijz(22) = 25
    ijz(23) = 25
    ijz(24) = 31
    ijz(25) = 25
    ijz(26) = 25
    ijz(27) = 31
    ijz(28) = 13
    ijz(29) = 13
    ijz(30) = 19

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

    allocate (n13bra(res%n_p_shl*res%n_f_shl))
    allocate (xint13bra(res%n_p_shl*res%n_f_shl))
    allocate (n12ket(res%n_p_shl*res%n_d_shl))
    allocate (xint12ket(res%n_p_shl*res%n_d_shl))

    ! Start screening

    scutpfbra = cutoff_schwarz/maxval(pf_pair%xints)
    npfbra = 0
    do ij = 1, res%n_p_shl*res%n_f_shl
      if (pf_pair%xints(ij) .ge. scutpfbra) then
        npfbra = npfbra + 1
        xint13bra(npfbra) = pf_pair%xints(ij)
        n13bra(npfbra) = ij
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

    if ((npfbra*npdket) .le. nchunksize_int64) nchunksize_int64 = npfbra*npdket
    ntile = int(npfbra*npdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = npfbra*npdket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npfbra, xint13bra, n13bra, xint12ket, n12ket, pd_pair, pf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d12ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d13bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, npfbra) + 1
              kl_tmp = (iquart - 1)/npfbra + 1

              test = xint13bra(ij_tmp)*xint12ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n13bra(ij_tmp)
                kl = n12ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_d_shl) + 1
                lsh_tmp = (kl - 1)/res%n_d_shl + 1

                ish = res%i_f_shl(ish_tmp)
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

                    t_expon_ab = pf_pair%t_expon_ab(pf_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = pf_pair%expon_b(pf_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = pf_pair%expon_a(pf_pair%pair_loc(ij) + bra_loop)
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

                    d13bra(1) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(2) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(3) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(4) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(5) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(6) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(7) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(8) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(9) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)
                    d13bra(10) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(11) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(12) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(13) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(14) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(15) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(16) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(17) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(18) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(19) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(20) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(21) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(22) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(23) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(24) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(25) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(26) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(27) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d13bra(28) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d13bra(29) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d13bra(30) = pf_pair%d_coeff_alt(pf_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3

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

                                      ! i5 = in(n+1) =   43
                                      ! i3 =   25
                                      ! i4 =   37

                                      xin(43) = c10*xin(25) + xc00*xin(37)
                                      yin(43) = c10*yin(25) + yc00*yin(37)
                                      zin(43) = c10*zin(25) + zc00*zin(37)

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

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   25

                                      xin(29) = c10*xin(5) + xc00*xin(17) + c01*xin(15)
                                      yin(29) = c10*yin(5) + yc00*yin(17) + c01*yin(15)
                                      zin(29) = c10*zin(5) + zc00*zin(17) + c01*zin(15)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(41) = c10*xin(17) + xc00*xin(29) + c01*xin(27)
                                      yin(41) = c10*yin(17) + yc00*yin(29) + c01*yin(27)
                                      zin(41) = c10*zin(17) + zc00*zin(29) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   43

                                      xin(47) = c10*xin(29) + xc00*xin(41) + c01*xin(39)
                                      yin(47) = c10*yin(29) + yc00*yin(41) + c01*yin(39)
                                      zin(47) = c10*zin(29) + zc00*zin(41) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   25

                                      xin(30) = c10*xin(6) + xc00*xin(18) + c01*xin(17)
                                      yin(30) = c10*yin(6) + yc00*yin(18) + c01*yin(17)
                                      zin(30) = c10*zin(6) + zc00*zin(18) + c01*zin(17)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(42) = c10*xin(18) + xc00*xin(30) + c01*xin(29)
                                      yin(42) = c10*yin(18) + yc00*yin(30) + c01*yin(29)
                                      zin(42) = c10*zin(18) + zc00*zin(30) + c01*zin(29)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   43

                                      xin(48) = c10*xin(30) + xc00*xin(42) + c01*xin(41)
                                      yin(48) = c10*yin(30) + yc00*yin(42) + c01*yin(41)
                                      zin(48) = c10*zin(30) + zc00*zin(42) + c01*zin(41)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   43

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   43

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   43

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   37

                                      xin(43) = xin(43) + dxij*xin(37)
                                      yin(43) = yin(43) + dyij*yin(37)
                                      zin(43) = zin(43) + dzij*zin(37)

                                      ! i3 = i4 =   37
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

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

                                      xin(31) = xin(37) + dxij*xin(25)
                                      yin(31) = yin(37) + dyij*yin(25)
                                      zin(31) = zin(37) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   39

                                      xin(45) = xin(45) + dxij*xin(39)
                                      yin(45) = yin(45) + dyij*yin(39)
                                      zin(45) = zin(45) + dzij*zin(39)

                                      ! i3 = i4 =   39
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

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

                                      xin(33) = xin(39) + dxij*xin(27)
                                      yin(33) = yin(39) + dyij*yin(27)
                                      zin(33) = zin(39) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(47) = xin(47) + dxij*xin(41)
                                      yin(47) = yin(47) + dyij*yin(41)
                                      zin(47) = zin(47) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

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

                                      xin(35) = xin(41) + dxij*xin(29)
                                      yin(35) = yin(41) + dyij*yin(29)
                                      zin(35) = zin(41) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   42

                                      xin(48) = xin(48) + dxij*xin(42)
                                      yin(48) = yin(48) + dyij*yin(42)
                                      zin(48) = zin(48) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   12

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

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

                                      xin(36) = xin(42) + dxij*xin(30)
                                      yin(36) = yin(42) + dyij*yin(30)
                                      zin(36) = zin(42) + dzij*zin(30)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

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

                                      ! do n = 2,   4

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

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   73
                                      ! i4 =   85

                                      xin(91) = c10*xin(73) + xc00*xin(85)
                                      yin(91) = c10*yin(73) + yc00*yin(85)
                                      zin(91) = c10*zin(73) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   93
                                      ! i5 =   91
                                      ! i4 =   85

                                      xin(93) = xcp00*xin(91) + cp10*xin(85)
                                      yin(93) = ycp00*yin(91) + cp10*yin(85)
                                      zin(93) = zcp00*zin(91) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   49
                                      ! i4 = i1+k2 =   51

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   53
                                      ! i3 =   49
                                      ! i4 =   51

                                      xin(53) = cp01*xin(49) + xcp00*xin(51)
                                      yin(53) = cp01*yin(49) + ycp00*yin(51)
                                      zin(53) = cp01*zin(49) + zcp00*zin(51)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   65

                                      xin(65) = xc00*xin(53) + c01*xin(51)
                                      yin(65) = yc00*yin(53) + c01*yin(51)
                                      zin(65) = zc00*zin(53) + c01*zin(51)

                                      ! ------------------

                                      ! i3 = i4 =   51
                                      ! i4 = i5 =   53

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   54
                                      ! i3 =   51
                                      ! i4 =   53

                                      xin(54) = cp01*xin(51) + xcp00*xin(53)
                                      yin(54) = cp01*yin(51) + ycp00*yin(53)
                                      zin(54) = cp01*zin(51) + zcp00*zin(53)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   66

                                      xin(66) = xc00*xin(54) + c01*xin(53)
                                      yin(66) = yc00*yin(54) + c01*yin(53)
                                      zin(66) = zc00*zin(54) + c01*zin(53)

                                      ! ------------------

                                      ! i3 = i4 =   53
                                      ! i4 = i5 =   54

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   73

                                      xin(77) = c10*xin(53) + xc00*xin(65) + c01*xin(63)
                                      yin(77) = c10*yin(53) + yc00*yin(65) + c01*yin(63)
                                      zin(77) = c10*zin(53) + zc00*zin(65) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(89) = c10*xin(65) + xc00*xin(77) + c01*xin(75)
                                      yin(89) = c10*yin(65) + yc00*yin(77) + c01*yin(75)
                                      zin(89) = c10*zin(65) + zc00*zin(77) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   91

                                      xin(95) = c10*xin(77) + xc00*xin(89) + c01*xin(87)
                                      yin(95) = c10*yin(77) + yc00*yin(89) + c01*yin(87)
                                      zin(95) = c10*zin(77) + zc00*zin(89) + c01*zin(87)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   73

                                      xin(78) = c10*xin(54) + xc00*xin(66) + c01*xin(65)
                                      yin(78) = c10*yin(54) + yc00*yin(66) + c01*yin(65)
                                      zin(78) = c10*zin(54) + zc00*zin(66) + c01*zin(65)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(90) = c10*xin(66) + xc00*xin(78) + c01*xin(77)
                                      yin(90) = c10*yin(66) + yc00*yin(78) + c01*yin(77)
                                      zin(90) = c10*zin(66) + zc00*zin(78) + c01*zin(77)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   91

                                      xin(96) = c10*xin(78) + xc00*xin(90) + c01*xin(89)
                                      yin(96) = c10*yin(78) + yc00*yin(90) + c01*yin(89)
                                      zin(96) = c10*zin(78) + zc00*zin(90) + c01*zin(89)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   91

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   91

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   91

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   85

                                      xin(91) = xin(91) + dxij*xin(85)
                                      yin(91) = yin(91) + dyij*yin(85)
                                      zin(91) = zin(91) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   55

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   55

                                      ! do ni = 1,    3

                                      xin(55) = xin(61) + dxij*xin(49)
                                      yin(55) = yin(61) + dyij*yin(49)
                                      zin(55) = zin(61) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    2

                                      xin(67) = xin(73) + dxij*xin(61)
                                      yin(67) = yin(73) + dyij*yin(61)
                                      zin(67) = zin(73) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(85) + dxij*xin(73)
                                      yin(79) = yin(85) + dyij*yin(73)
                                      zin(79) = zin(85) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   61

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   87

                                      xin(93) = xin(93) + dxij*xin(87)
                                      yin(93) = yin(93) + dyij*yin(87)
                                      zin(93) = zin(93) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   57

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   57

                                      ! do ni = 1,    3

                                      xin(57) = xin(63) + dxij*xin(51)
                                      yin(57) = yin(63) + dyij*yin(51)
                                      zin(57) = zin(63) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    2

                                      xin(69) = xin(75) + dxij*xin(63)
                                      yin(69) = yin(75) + dyij*yin(63)
                                      zin(69) = zin(75) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(87) + dxij*xin(75)
                                      yin(81) = yin(87) + dyij*yin(75)
                                      zin(81) = zin(87) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   63

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(95) = xin(95) + dxij*xin(89)
                                      yin(95) = yin(95) + dyij*yin(89)
                                      zin(95) = zin(95) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   59

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   59

                                      ! do ni = 1,    3

                                      xin(59) = xin(65) + dxij*xin(53)
                                      yin(59) = yin(65) + dyij*yin(53)
                                      zin(59) = zin(65) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    2

                                      xin(71) = xin(77) + dxij*xin(65)
                                      yin(71) = yin(77) + dyij*yin(65)
                                      zin(71) = zin(77) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    3

                                      xin(83) = xin(89) + dxij*xin(77)
                                      yin(83) = yin(89) + dyij*yin(77)
                                      zin(83) = zin(89) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   65

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(96) = xin(96) + dxij*xin(90)
                                      yin(96) = yin(96) + dyij*yin(90)
                                      zin(96) = zin(96) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   60

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   60

                                      ! do ni = 1,    3

                                      xin(60) = xin(66) + dxij*xin(54)
                                      yin(60) = yin(66) + dyij*yin(54)
                                      zin(60) = zin(66) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    2

                                      xin(72) = xin(78) + dxij*xin(66)
                                      yin(72) = yin(78) + dyij*yin(66)
                                      zin(72) = zin(78) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    3

                                      xin(84) = xin(90) + dxij*xin(78)
                                      yin(84) = yin(90) + dyij*yin(78)
                                      zin(84) = zin(90) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   66

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =   49

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! ni = ni + 1 =    1

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

                                      ! ni = ni + 1 =    2

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

                                      ! end do

                                      ! ni = ni + 1 =    3

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

                                      ! do n = 2,   4

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

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   99

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  101
                                      ! i3 =   97
                                      ! i4 =   99

                                      xin(101) = cp01*xin(97) + xcp00*xin(99)
                                      yin(101) = cp01*yin(97) + ycp00*yin(99)
                                      zin(101) = cp01*zin(97) + zcp00*zin(99)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  113

                                      xin(113) = xc00*xin(101) + c01*xin(99)
                                      yin(113) = yc00*yin(101) + c01*yin(99)
                                      zin(113) = zc00*zin(101) + c01*zin(99)

                                      ! ------------------

                                      ! i3 = i4 =   99
                                      ! i4 = i5 =  101

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  102
                                      ! i3 =   99
                                      ! i4 =  101

                                      xin(102) = cp01*xin(99) + xcp00*xin(101)
                                      yin(102) = cp01*yin(99) + ycp00*yin(101)
                                      zin(102) = cp01*zin(99) + zcp00*zin(101)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  114

                                      xin(114) = xc00*xin(102) + c01*xin(101)
                                      yin(114) = yc00*yin(102) + c01*yin(101)
                                      zin(114) = zc00*zin(102) + c01*zin(101)

                                      ! ------------------

                                      ! i3 = i4 =  101
                                      ! i4 = i5 =  102

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  121

                                      xin(125) = c10*xin(101) + xc00*xin(113) + c01*xin(111)
                                      yin(125) = c10*yin(101) + yc00*yin(113) + c01*yin(111)
                                      zin(125) = c10*zin(101) + zc00*zin(113) + c01*zin(111)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(137) = c10*xin(113) + xc00*xin(125) + c01*xin(123)
                                      yin(137) = c10*yin(113) + yc00*yin(125) + c01*yin(123)
                                      zin(137) = c10*zin(113) + zc00*zin(125) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  139

                                      xin(143) = c10*xin(125) + xc00*xin(137) + c01*xin(135)
                                      yin(143) = c10*yin(125) + yc00*yin(137) + c01*yin(135)
                                      zin(143) = c10*zin(125) + zc00*zin(137) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  139

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  121

                                      xin(126) = c10*xin(102) + xc00*xin(114) + c01*xin(113)
                                      yin(126) = c10*yin(102) + yc00*yin(114) + c01*yin(113)
                                      zin(126) = c10*zin(102) + zc00*zin(114) + c01*zin(113)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(138) = c10*xin(114) + xc00*xin(126) + c01*xin(125)
                                      yin(138) = c10*yin(114) + yc00*yin(126) + c01*yin(125)
                                      zin(138) = c10*zin(114) + zc00*zin(126) + c01*zin(125)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  139

                                      xin(144) = c10*xin(126) + xc00*xin(138) + c01*xin(137)
                                      yin(144) = c10*yin(126) + yc00*yin(138) + c01*yin(137)
                                      zin(144) = c10*zin(126) + zc00*zin(138) + c01*zin(137)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  139

                                      ! nn =    5

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

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(139) = xin(139) + dxij*xin(133)
                                      yin(139) = yin(139) + dyij*yin(133)
                                      zin(139) = zin(139) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  103

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(109) + dxij*xin(97)
                                      yin(103) = yin(109) + dyij*yin(97)
                                      zin(103) = zin(109) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    2

                                      xin(115) = xin(121) + dxij*xin(109)
                                      yin(115) = yin(121) + dyij*yin(109)
                                      zin(115) = zin(121) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    3

                                      xin(127) = xin(133) + dxij*xin(121)
                                      yin(127) = yin(133) + dyij*yin(121)
                                      zin(127) = zin(133) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  109

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(141) = xin(141) + dxij*xin(135)
                                      yin(141) = yin(141) + dyij*yin(135)
                                      zin(141) = zin(141) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  105

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  105

                                      ! do ni = 1,    3

                                      xin(105) = xin(111) + dxij*xin(99)
                                      yin(105) = yin(111) + dyij*yin(99)
                                      zin(105) = zin(111) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    2

                                      xin(117) = xin(123) + dxij*xin(111)
                                      yin(117) = yin(123) + dyij*yin(111)
                                      zin(117) = zin(123) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    3

                                      xin(129) = xin(135) + dxij*xin(123)
                                      yin(129) = yin(135) + dyij*yin(123)
                                      zin(129) = zin(135) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  111

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(143) = xin(143) + dxij*xin(137)
                                      yin(143) = yin(143) + dyij*yin(137)
                                      zin(143) = zin(143) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  107

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  107

                                      ! do ni = 1,    3

                                      xin(107) = xin(113) + dxij*xin(101)
                                      yin(107) = yin(113) + dyij*yin(101)
                                      zin(107) = zin(113) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    2

                                      xin(119) = xin(125) + dxij*xin(113)
                                      yin(119) = yin(125) + dyij*yin(113)
                                      zin(119) = zin(125) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      xin(131) = xin(137) + dxij*xin(125)
                                      yin(131) = yin(137) + dyij*yin(125)
                                      zin(131) = zin(137) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  113

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(144) = xin(144) + dxij*xin(138)
                                      yin(144) = yin(144) + dyij*yin(138)
                                      zin(144) = zin(144) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  108

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  108

                                      ! do ni = 1,    3

                                      xin(108) = xin(114) + dxij*xin(102)
                                      yin(108) = yin(114) + dyij*yin(102)
                                      zin(108) = zin(114) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    2

                                      xin(120) = xin(126) + dxij*xin(114)
                                      yin(120) = yin(126) + dyij*yin(114)
                                      zin(120) = zin(126) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    3

                                      xin(132) = xin(138) + dxij*xin(126)
                                      yin(132) = yin(138) + dyij*yin(126)
                                      zin(132) = zin(138) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  114

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    5

                                      ! iaa = i1 =   97

                                      ! ni = 0

                                      ! do while ni.le.iang

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

                                      ! ni = ni + 1 =    1

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

                                      ! end do

                                      ! ni = ni + 1 =    2

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

                                      ! ni = ni + 1 =    3

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

                                      ! do n = 2,   4

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

                                      ! i5 = in(n+1) =  187
                                      ! i3 =  169
                                      ! i4 =  181

                                      xin(187) = c10*xin(169) + xc00*xin(181)
                                      yin(187) = c10*yin(169) + yc00*yin(181)
                                      zin(187) = c10*zin(169) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  189
                                      ! i5 =  187
                                      ! i4 =  181

                                      xin(189) = xcp00*xin(187) + cp10*xin(181)
                                      yin(189) = ycp00*yin(187) + cp10*yin(181)
                                      zin(189) = zcp00*zin(187) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! n =    5

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

                                      ! i3 = i2+kn(n+1) =  161

                                      xin(161) = xc00*xin(149) + c01*xin(147)
                                      yin(161) = yc00*yin(149) + c01*yin(147)
                                      zin(161) = zc00*zin(149) + c01*zin(147)

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

                                      ! i3 = i2+kn(n+1) =  162

                                      xin(162) = xc00*xin(150) + c01*xin(149)
                                      yin(162) = yc00*yin(150) + c01*yin(149)
                                      zin(162) = zc00*zin(150) + c01*zin(149)

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
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  169

                                      xin(173) = c10*xin(149) + xc00*xin(161) + c01*xin(159)
                                      yin(173) = c10*yin(149) + yc00*yin(161) + c01*yin(159)
                                      zin(173) = c10*zin(149) + zc00*zin(161) + c01*zin(159)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(185) = c10*xin(161) + xc00*xin(173) + c01*xin(171)
                                      yin(185) = c10*yin(161) + yc00*yin(173) + c01*yin(171)
                                      zin(185) = c10*zin(161) + zc00*zin(173) + c01*zin(171)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  187

                                      xin(191) = c10*xin(173) + xc00*xin(185) + c01*xin(183)
                                      yin(191) = c10*yin(173) + yc00*yin(185) + c01*yin(183)
                                      zin(191) = c10*zin(173) + zc00*zin(185) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    5
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  169

                                      xin(174) = c10*xin(150) + xc00*xin(162) + c01*xin(161)
                                      yin(174) = c10*yin(150) + yc00*yin(162) + c01*yin(161)
                                      zin(174) = c10*zin(150) + zc00*zin(162) + c01*zin(161)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(186) = c10*xin(162) + xc00*xin(174) + c01*xin(173)
                                      yin(186) = c10*yin(162) + yc00*yin(174) + c01*yin(173)
                                      zin(186) = c10*zin(162) + zc00*zin(174) + c01*zin(173)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  187

                                      xin(192) = c10*xin(174) + xc00*xin(186) + c01*xin(185)
                                      yin(192) = c10*yin(174) + yc00*yin(186) + c01*yin(185)
                                      zin(192) = c10*zin(174) + zc00*zin(186) + c01*zin(185)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  187

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   5

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  187

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  187

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  181

                                      xin(187) = xin(187) + dxij*xin(181)
                                      yin(187) = yin(187) + dyij*yin(181)
                                      zin(187) = zin(187) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  151

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(157) + dxij*xin(145)
                                      yin(151) = yin(157) + dyij*yin(145)
                                      zin(151) = zin(157) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! ni =    2

                                      xin(163) = xin(169) + dxij*xin(157)
                                      yin(163) = yin(169) + dyij*yin(157)
                                      zin(163) = zin(169) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(181) + dxij*xin(169)
                                      yin(175) = yin(181) + dyij*yin(169)
                                      zin(175) = zin(181) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  183

                                      xin(189) = xin(189) + dxij*xin(183)
                                      yin(189) = yin(189) + dyij*yin(183)
                                      zin(189) = zin(189) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  153

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(159) + dxij*xin(147)
                                      yin(153) = yin(159) + dyij*yin(147)
                                      zin(153) = zin(159) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    2

                                      xin(165) = xin(171) + dxij*xin(159)
                                      yin(165) = yin(171) + dyij*yin(159)
                                      zin(165) = zin(171) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    3

                                      xin(177) = xin(183) + dxij*xin(171)
                                      yin(177) = yin(183) + dyij*yin(171)
                                      zin(177) = zin(183) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(191) = xin(191) + dxij*xin(185)
                                      yin(191) = yin(191) + dyij*yin(185)
                                      zin(191) = zin(191) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  155

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(161) + dxij*xin(149)
                                      yin(155) = yin(161) + dyij*yin(149)
                                      zin(155) = zin(161) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    2

                                      xin(167) = xin(173) + dxij*xin(161)
                                      yin(167) = yin(173) + dyij*yin(161)
                                      zin(167) = zin(173) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    3

                                      xin(179) = xin(185) + dxij*xin(173)
                                      yin(179) = yin(185) + dyij*yin(173)
                                      zin(179) = zin(185) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  161

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    5

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(192) = xin(192) + dxij*xin(186)
                                      yin(192) = yin(192) + dyij*yin(186)
                                      zin(192) = zin(192) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  156

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  156

                                      ! do ni = 1,    3

                                      xin(156) = xin(162) + dxij*xin(150)
                                      yin(156) = yin(162) + dyij*yin(150)
                                      zin(156) = zin(162) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    2

                                      xin(168) = xin(174) + dxij*xin(162)
                                      yin(168) = yin(174) + dyij*yin(162)
                                      zin(168) = zin(174) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    3

                                      xin(180) = xin(186) + dxij*xin(174)
                                      yin(180) = yin(186) + dyij*yin(174)
                                      zin(180) = zin(186) + dzij*zin(174)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  162

                                      ! nj =    2

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

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

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

                                      ! nj = nj + 1 =    1

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

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    3

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

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d13bra(  1)*d12ket(  1)*(xin(  48)*yin(   1)*zin(   1)+xin(  96)*yin(  49)*zin(  49)+xin( 144)*yin(  97)*zin(  97)+xin( 192)*yin( 145)*zin( 145))
          eri_value(    2)=eri_value(    2)+d13bra(  1)*d12ket(  2)*(xin(  47)*yin(   2)*zin(   1)+xin(  95)*yin(  50)*zin(  49)+xin( 143)*yin(  98)*zin(  97)+xin( 191)*yin( 146)*zin( 145))
          eri_value(    3)=eri_value(    3)+d13bra(  1)*d12ket(  3)*(xin(  47)*yin(   1)*zin(   2)+xin(  95)*yin(  49)*zin(  50)+xin( 143)*yin(  97)*zin(  98)+xin( 191)*yin( 145)*zin( 146))
          eri_value(    4)=eri_value(    4)+d13bra(  1)*d12ket(  4)*(xin(  44)*yin(   5)*zin(   1)+xin(  92)*yin(  53)*zin(  49)+xin( 140)*yin( 101)*zin(  97)+xin( 188)*yin( 149)*zin( 145))
          eri_value(    5)=eri_value(    5)+d13bra(  1)*d12ket(  5)*(xin(  43)*yin(   6)*zin(   1)+xin(  91)*yin(  54)*zin(  49)+xin( 139)*yin( 102)*zin(  97)+xin( 187)*yin( 150)*zin( 145))
          eri_value(    6)=eri_value(    6)+d13bra(  1)*d12ket(  6)*(xin(  43)*yin(   5)*zin(   2)+xin(  91)*yin(  53)*zin(  50)+xin( 139)*yin( 101)*zin(  98)+xin( 187)*yin( 149)*zin( 146))
          eri_value(    7)=eri_value(    7)+d13bra(  1)*d12ket(  7)*(xin(  44)*yin(   1)*zin(   5)+xin(  92)*yin(  49)*zin(  53)+xin( 140)*yin(  97)*zin( 101)+xin( 188)*yin( 145)*zin( 149))
          eri_value(    8)=eri_value(    8)+d13bra(  1)*d12ket(  8)*(xin(  43)*yin(   2)*zin(   5)+xin(  91)*yin(  50)*zin(  53)+xin( 139)*yin(  98)*zin( 101)+xin( 187)*yin( 146)*zin( 149))
          eri_value(    9)=eri_value(    9)+d13bra(  1)*d12ket(  9)*(xin(  43)*yin(   1)*zin(   6)+xin(  91)*yin(  49)*zin(  54)+xin( 139)*yin(  97)*zin( 102)+xin( 187)*yin( 145)*zin( 150))
          eri_value(   10)=eri_value(   10)+d13bra(  1)*d12ket( 10)*(xin(  46)*yin(   3)*zin(   1)+xin(  94)*yin(  51)*zin(  49)+xin( 142)*yin(  99)*zin(  97)+xin( 190)*yin( 147)*zin( 145))
          eri_value(   11)=eri_value(   11)+d13bra(  1)*d12ket( 11)*(xin(  45)*yin(   4)*zin(   1)+xin(  93)*yin(  52)*zin(  49)+xin( 141)*yin( 100)*zin(  97)+xin( 189)*yin( 148)*zin( 145))
          eri_value(   12)=eri_value(   12)+d13bra(  1)*d12ket( 12)*(xin(  45)*yin(   3)*zin(   2)+xin(  93)*yin(  51)*zin(  50)+xin( 141)*yin(  99)*zin(  98)+xin( 189)*yin( 147)*zin( 146))
          eri_value(   13)=eri_value(   13)+d13bra(  1)*d12ket( 13)*(xin(  46)*yin(   1)*zin(   3)+xin(  94)*yin(  49)*zin(  51)+xin( 142)*yin(  97)*zin(  99)+xin( 190)*yin( 145)*zin( 147))
          eri_value(   14)=eri_value(   14)+d13bra(  1)*d12ket( 14)*(xin(  45)*yin(   2)*zin(   3)+xin(  93)*yin(  50)*zin(  51)+xin( 141)*yin(  98)*zin(  99)+xin( 189)*yin( 146)*zin( 147))
          eri_value(   15)=eri_value(   15)+d13bra(  1)*d12ket( 15)*(xin(  45)*yin(   1)*zin(   4)+xin(  93)*yin(  49)*zin(  52)+xin( 141)*yin(  97)*zin( 100)+xin( 189)*yin( 145)*zin( 148))
          eri_value(   16)=eri_value(   16)+d13bra(  1)*d12ket( 16)*(xin(  44)*yin(   3)*zin(   3)+xin(  92)*yin(  51)*zin(  51)+xin( 140)*yin(  99)*zin(  99)+xin( 188)*yin( 147)*zin( 147))
          eri_value(   17)=eri_value(   17)+d13bra(  1)*d12ket( 17)*(xin(  43)*yin(   4)*zin(   3)+xin(  91)*yin(  52)*zin(  51)+xin( 139)*yin( 100)*zin(  99)+xin( 187)*yin( 148)*zin( 147))
          eri_value(   18)=eri_value(   18)+d13bra(  1)*d12ket( 18)*(xin(  43)*yin(   3)*zin(   4)+xin(  91)*yin(  51)*zin(  52)+xin( 139)*yin(  99)*zin( 100)+xin( 187)*yin( 147)*zin( 148))
          eri_value(   19)=eri_value(   19)+d13bra(  2)*d12ket(  1)*(xin(  42)*yin(   7)*zin(   1)+xin(  90)*yin(  55)*zin(  49)+xin( 138)*yin( 103)*zin(  97)+xin( 186)*yin( 151)*zin( 145))
          eri_value(   20)=eri_value(   20)+d13bra(  2)*d12ket(  2)*(xin(  41)*yin(   8)*zin(   1)+xin(  89)*yin(  56)*zin(  49)+xin( 137)*yin( 104)*zin(  97)+xin( 185)*yin( 152)*zin( 145))
          eri_value(   21)=eri_value(   21)+d13bra(  2)*d12ket(  3)*(xin(  41)*yin(   7)*zin(   2)+xin(  89)*yin(  55)*zin(  50)+xin( 137)*yin( 103)*zin(  98)+xin( 185)*yin( 151)*zin( 146))
          eri_value(   22)=eri_value(   22)+d13bra(  2)*d12ket(  4)*(xin(  38)*yin(  11)*zin(   1)+xin(  86)*yin(  59)*zin(  49)+xin( 134)*yin( 107)*zin(  97)+xin( 182)*yin( 155)*zin( 145))
          eri_value(   23)=eri_value(   23)+d13bra(  2)*d12ket(  5)*(xin(  37)*yin(  12)*zin(   1)+xin(  85)*yin(  60)*zin(  49)+xin( 133)*yin( 108)*zin(  97)+xin( 181)*yin( 156)*zin( 145))
          eri_value(   24)=eri_value(   24)+d13bra(  2)*d12ket(  6)*(xin(  37)*yin(  11)*zin(   2)+xin(  85)*yin(  59)*zin(  50)+xin( 133)*yin( 107)*zin(  98)+xin( 181)*yin( 155)*zin( 146))
          eri_value(   25)=eri_value(   25)+d13bra(  2)*d12ket(  7)*(xin(  38)*yin(   7)*zin(   5)+xin(  86)*yin(  55)*zin(  53)+xin( 134)*yin( 103)*zin( 101)+xin( 182)*yin( 151)*zin( 149))
          eri_value(   26)=eri_value(   26)+d13bra(  2)*d12ket(  8)*(xin(  37)*yin(   8)*zin(   5)+xin(  85)*yin(  56)*zin(  53)+xin( 133)*yin( 104)*zin( 101)+xin( 181)*yin( 152)*zin( 149))
          eri_value(   27)=eri_value(   27)+d13bra(  2)*d12ket(  9)*(xin(  37)*yin(   7)*zin(   6)+xin(  85)*yin(  55)*zin(  54)+xin( 133)*yin( 103)*zin( 102)+xin( 181)*yin( 151)*zin( 150))
          eri_value(   28)=eri_value(   28)+d13bra(  2)*d12ket( 10)*(xin(  40)*yin(   9)*zin(   1)+xin(  88)*yin(  57)*zin(  49)+xin( 136)*yin( 105)*zin(  97)+xin( 184)*yin( 153)*zin( 145))
          eri_value(   29)=eri_value(   29)+d13bra(  2)*d12ket( 11)*(xin(  39)*yin(  10)*zin(   1)+xin(  87)*yin(  58)*zin(  49)+xin( 135)*yin( 106)*zin(  97)+xin( 183)*yin( 154)*zin( 145))
          eri_value(   30)=eri_value(   30)+d13bra(  2)*d12ket( 12)*(xin(  39)*yin(   9)*zin(   2)+xin(  87)*yin(  57)*zin(  50)+xin( 135)*yin( 105)*zin(  98)+xin( 183)*yin( 153)*zin( 146))
          eri_value(   31)=eri_value(   31)+d13bra(  2)*d12ket( 13)*(xin(  40)*yin(   7)*zin(   3)+xin(  88)*yin(  55)*zin(  51)+xin( 136)*yin( 103)*zin(  99)+xin( 184)*yin( 151)*zin( 147))
          eri_value(   32)=eri_value(   32)+d13bra(  2)*d12ket( 14)*(xin(  39)*yin(   8)*zin(   3)+xin(  87)*yin(  56)*zin(  51)+xin( 135)*yin( 104)*zin(  99)+xin( 183)*yin( 152)*zin( 147))
          eri_value(   33)=eri_value(   33)+d13bra(  2)*d12ket( 15)*(xin(  39)*yin(   7)*zin(   4)+xin(  87)*yin(  55)*zin(  52)+xin( 135)*yin( 103)*zin( 100)+xin( 183)*yin( 151)*zin( 148))
          eri_value(   34)=eri_value(   34)+d13bra(  2)*d12ket( 16)*(xin(  38)*yin(   9)*zin(   3)+xin(  86)*yin(  57)*zin(  51)+xin( 134)*yin( 105)*zin(  99)+xin( 182)*yin( 153)*zin( 147))
          eri_value(   35)=eri_value(   35)+d13bra(  2)*d12ket( 17)*(xin(  37)*yin(  10)*zin(   3)+xin(  85)*yin(  58)*zin(  51)+xin( 133)*yin( 106)*zin(  99)+xin( 181)*yin( 154)*zin( 147))
          eri_value(   36)=eri_value(   36)+d13bra(  2)*d12ket( 18)*(xin(  37)*yin(   9)*zin(   4)+xin(  85)*yin(  57)*zin(  52)+xin( 133)*yin( 105)*zin( 100)+xin( 181)*yin( 153)*zin( 148))
          eri_value(   37)=eri_value(   37)+d13bra(  3)*d12ket(  1)*(xin(  42)*yin(   1)*zin(   7)+xin(  90)*yin(  49)*zin(  55)+xin( 138)*yin(  97)*zin( 103)+xin( 186)*yin( 145)*zin( 151))
          eri_value(   38)=eri_value(   38)+d13bra(  3)*d12ket(  2)*(xin(  41)*yin(   2)*zin(   7)+xin(  89)*yin(  50)*zin(  55)+xin( 137)*yin(  98)*zin( 103)+xin( 185)*yin( 146)*zin( 151))
          eri_value(   39)=eri_value(   39)+d13bra(  3)*d12ket(  3)*(xin(  41)*yin(   1)*zin(   8)+xin(  89)*yin(  49)*zin(  56)+xin( 137)*yin(  97)*zin( 104)+xin( 185)*yin( 145)*zin( 152))
          eri_value(   40)=eri_value(   40)+d13bra(  3)*d12ket(  4)*(xin(  38)*yin(   5)*zin(   7)+xin(  86)*yin(  53)*zin(  55)+xin( 134)*yin( 101)*zin( 103)+xin( 182)*yin( 149)*zin( 151))
          eri_value(   41)=eri_value(   41)+d13bra(  3)*d12ket(  5)*(xin(  37)*yin(   6)*zin(   7)+xin(  85)*yin(  54)*zin(  55)+xin( 133)*yin( 102)*zin( 103)+xin( 181)*yin( 150)*zin( 151))
          eri_value(   42)=eri_value(   42)+d13bra(  3)*d12ket(  6)*(xin(  37)*yin(   5)*zin(   8)+xin(  85)*yin(  53)*zin(  56)+xin( 133)*yin( 101)*zin( 104)+xin( 181)*yin( 149)*zin( 152))
          eri_value(   43)=eri_value(   43)+d13bra(  3)*d12ket(  7)*(xin(  38)*yin(   1)*zin(  11)+xin(  86)*yin(  49)*zin(  59)+xin( 134)*yin(  97)*zin( 107)+xin( 182)*yin( 145)*zin( 155))
          eri_value(   44)=eri_value(   44)+d13bra(  3)*d12ket(  8)*(xin(  37)*yin(   2)*zin(  11)+xin(  85)*yin(  50)*zin(  59)+xin( 133)*yin(  98)*zin( 107)+xin( 181)*yin( 146)*zin( 155))
          eri_value(   45)=eri_value(   45)+d13bra(  3)*d12ket(  9)*(xin(  37)*yin(   1)*zin(  12)+xin(  85)*yin(  49)*zin(  60)+xin( 133)*yin(  97)*zin( 108)+xin( 181)*yin( 145)*zin( 156))
          eri_value(   46)=eri_value(   46)+d13bra(  3)*d12ket( 10)*(xin(  40)*yin(   3)*zin(   7)+xin(  88)*yin(  51)*zin(  55)+xin( 136)*yin(  99)*zin( 103)+xin( 184)*yin( 147)*zin( 151))
          eri_value(   47)=eri_value(   47)+d13bra(  3)*d12ket( 11)*(xin(  39)*yin(   4)*zin(   7)+xin(  87)*yin(  52)*zin(  55)+xin( 135)*yin( 100)*zin( 103)+xin( 183)*yin( 148)*zin( 151))
          eri_value(   48)=eri_value(   48)+d13bra(  3)*d12ket( 12)*(xin(  39)*yin(   3)*zin(   8)+xin(  87)*yin(  51)*zin(  56)+xin( 135)*yin(  99)*zin( 104)+xin( 183)*yin( 147)*zin( 152))
          eri_value(   49)=eri_value(   49)+d13bra(  3)*d12ket( 13)*(xin(  40)*yin(   1)*zin(   9)+xin(  88)*yin(  49)*zin(  57)+xin( 136)*yin(  97)*zin( 105)+xin( 184)*yin( 145)*zin( 153))
          eri_value(   50)=eri_value(   50)+d13bra(  3)*d12ket( 14)*(xin(  39)*yin(   2)*zin(   9)+xin(  87)*yin(  50)*zin(  57)+xin( 135)*yin(  98)*zin( 105)+xin( 183)*yin( 146)*zin( 153))
          eri_value(   51)=eri_value(   51)+d13bra(  3)*d12ket( 15)*(xin(  39)*yin(   1)*zin(  10)+xin(  87)*yin(  49)*zin(  58)+xin( 135)*yin(  97)*zin( 106)+xin( 183)*yin( 145)*zin( 154))
          eri_value(   52)=eri_value(   52)+d13bra(  3)*d12ket( 16)*(xin(  38)*yin(   3)*zin(   9)+xin(  86)*yin(  51)*zin(  57)+xin( 134)*yin(  99)*zin( 105)+xin( 182)*yin( 147)*zin( 153))
          eri_value(   53)=eri_value(   53)+d13bra(  3)*d12ket( 17)*(xin(  37)*yin(   4)*zin(   9)+xin(  85)*yin(  52)*zin(  57)+xin( 133)*yin( 100)*zin( 105)+xin( 181)*yin( 148)*zin( 153))
          eri_value(   54)=eri_value(   54)+d13bra(  3)*d12ket( 18)*(xin(  37)*yin(   3)*zin(  10)+xin(  85)*yin(  51)*zin(  58)+xin( 133)*yin(  99)*zin( 106)+xin( 181)*yin( 147)*zin( 154))
          eri_value(   55)=eri_value(   55)+d13bra(  4)*d12ket(  1)*(xin(  12)*yin(  37)*zin(   1)+xin(  60)*yin(  85)*zin(  49)+xin( 108)*yin( 133)*zin(  97)+xin( 156)*yin( 181)*zin( 145))
          eri_value(   56)=eri_value(   56)+d13bra(  4)*d12ket(  2)*(xin(  11)*yin(  38)*zin(   1)+xin(  59)*yin(  86)*zin(  49)+xin( 107)*yin( 134)*zin(  97)+xin( 155)*yin( 182)*zin( 145))
          eri_value(   57)=eri_value(   57)+d13bra(  4)*d12ket(  3)*(xin(  11)*yin(  37)*zin(   2)+xin(  59)*yin(  85)*zin(  50)+xin( 107)*yin( 133)*zin(  98)+xin( 155)*yin( 181)*zin( 146))
          eri_value(   58)=eri_value(   58)+d13bra(  4)*d12ket(  4)*(xin(   8)*yin(  41)*zin(   1)+xin(  56)*yin(  89)*zin(  49)+xin( 104)*yin( 137)*zin(  97)+xin( 152)*yin( 185)*zin( 145))
          eri_value(   59)=eri_value(   59)+d13bra(  4)*d12ket(  5)*(xin(   7)*yin(  42)*zin(   1)+xin(  55)*yin(  90)*zin(  49)+xin( 103)*yin( 138)*zin(  97)+xin( 151)*yin( 186)*zin( 145))
          eri_value(   60)=eri_value(   60)+d13bra(  4)*d12ket(  6)*(xin(   7)*yin(  41)*zin(   2)+xin(  55)*yin(  89)*zin(  50)+xin( 103)*yin( 137)*zin(  98)+xin( 151)*yin( 185)*zin( 146))
          eri_value(   61)=eri_value(   61)+d13bra(  4)*d12ket(  7)*(xin(   8)*yin(  37)*zin(   5)+xin(  56)*yin(  85)*zin(  53)+xin( 104)*yin( 133)*zin( 101)+xin( 152)*yin( 181)*zin( 149))
          eri_value(   62)=eri_value(   62)+d13bra(  4)*d12ket(  8)*(xin(   7)*yin(  38)*zin(   5)+xin(  55)*yin(  86)*zin(  53)+xin( 103)*yin( 134)*zin( 101)+xin( 151)*yin( 182)*zin( 149))
          eri_value(   63)=eri_value(   63)+d13bra(  4)*d12ket(  9)*(xin(   7)*yin(  37)*zin(   6)+xin(  55)*yin(  85)*zin(  54)+xin( 103)*yin( 133)*zin( 102)+xin( 151)*yin( 181)*zin( 150))
          eri_value(   64)=eri_value(   64)+d13bra(  4)*d12ket( 10)*(xin(  10)*yin(  39)*zin(   1)+xin(  58)*yin(  87)*zin(  49)+xin( 106)*yin( 135)*zin(  97)+xin( 154)*yin( 183)*zin( 145))
          eri_value(   65)=eri_value(   65)+d13bra(  4)*d12ket( 11)*(xin(   9)*yin(  40)*zin(   1)+xin(  57)*yin(  88)*zin(  49)+xin( 105)*yin( 136)*zin(  97)+xin( 153)*yin( 184)*zin( 145))
          eri_value(   66)=eri_value(   66)+d13bra(  4)*d12ket( 12)*(xin(   9)*yin(  39)*zin(   2)+xin(  57)*yin(  87)*zin(  50)+xin( 105)*yin( 135)*zin(  98)+xin( 153)*yin( 183)*zin( 146))
          eri_value(   67)=eri_value(   67)+d13bra(  4)*d12ket( 13)*(xin(  10)*yin(  37)*zin(   3)+xin(  58)*yin(  85)*zin(  51)+xin( 106)*yin( 133)*zin(  99)+xin( 154)*yin( 181)*zin( 147))
          eri_value(   68)=eri_value(   68)+d13bra(  4)*d12ket( 14)*(xin(   9)*yin(  38)*zin(   3)+xin(  57)*yin(  86)*zin(  51)+xin( 105)*yin( 134)*zin(  99)+xin( 153)*yin( 182)*zin( 147))
          eri_value(   69)=eri_value(   69)+d13bra(  4)*d12ket( 15)*(xin(   9)*yin(  37)*zin(   4)+xin(  57)*yin(  85)*zin(  52)+xin( 105)*yin( 133)*zin( 100)+xin( 153)*yin( 181)*zin( 148))
          eri_value(   70)=eri_value(   70)+d13bra(  4)*d12ket( 16)*(xin(   8)*yin(  39)*zin(   3)+xin(  56)*yin(  87)*zin(  51)+xin( 104)*yin( 135)*zin(  99)+xin( 152)*yin( 183)*zin( 147))
          eri_value(   71)=eri_value(   71)+d13bra(  4)*d12ket( 17)*(xin(   7)*yin(  40)*zin(   3)+xin(  55)*yin(  88)*zin(  51)+xin( 103)*yin( 136)*zin(  99)+xin( 151)*yin( 184)*zin( 147))
          eri_value(   72)=eri_value(   72)+d13bra(  4)*d12ket( 18)*(xin(   7)*yin(  39)*zin(   4)+xin(  55)*yin(  87)*zin(  52)+xin( 103)*yin( 135)*zin( 100)+xin( 151)*yin( 183)*zin( 148))
          eri_value(   73)=eri_value(   73)+d13bra(  5)*d12ket(  1)*(xin(   6)*yin(  43)*zin(   1)+xin(  54)*yin(  91)*zin(  49)+xin( 102)*yin( 139)*zin(  97)+xin( 150)*yin( 187)*zin( 145))
          eri_value(   74)=eri_value(   74)+d13bra(  5)*d12ket(  2)*(xin(   5)*yin(  44)*zin(   1)+xin(  53)*yin(  92)*zin(  49)+xin( 101)*yin( 140)*zin(  97)+xin( 149)*yin( 188)*zin( 145))
          eri_value(   75)=eri_value(   75)+d13bra(  5)*d12ket(  3)*(xin(   5)*yin(  43)*zin(   2)+xin(  53)*yin(  91)*zin(  50)+xin( 101)*yin( 139)*zin(  98)+xin( 149)*yin( 187)*zin( 146))
          eri_value(   76)=eri_value(   76)+d13bra(  5)*d12ket(  4)*(xin(   2)*yin(  47)*zin(   1)+xin(  50)*yin(  95)*zin(  49)+xin(  98)*yin( 143)*zin(  97)+xin( 146)*yin( 191)*zin( 145))
          eri_value(   77)=eri_value(   77)+d13bra(  5)*d12ket(  5)*(xin(   1)*yin(  48)*zin(   1)+xin(  49)*yin(  96)*zin(  49)+xin(  97)*yin( 144)*zin(  97)+xin( 145)*yin( 192)*zin( 145))
          eri_value(   78)=eri_value(   78)+d13bra(  5)*d12ket(  6)*(xin(   1)*yin(  47)*zin(   2)+xin(  49)*yin(  95)*zin(  50)+xin(  97)*yin( 143)*zin(  98)+xin( 145)*yin( 191)*zin( 146))
          eri_value(   79)=eri_value(   79)+d13bra(  5)*d12ket(  7)*(xin(   2)*yin(  43)*zin(   5)+xin(  50)*yin(  91)*zin(  53)+xin(  98)*yin( 139)*zin( 101)+xin( 146)*yin( 187)*zin( 149))
          eri_value(   80)=eri_value(   80)+d13bra(  5)*d12ket(  8)*(xin(   1)*yin(  44)*zin(   5)+xin(  49)*yin(  92)*zin(  53)+xin(  97)*yin( 140)*zin( 101)+xin( 145)*yin( 188)*zin( 149))
          eri_value(   81)=eri_value(   81)+d13bra(  5)*d12ket(  9)*(xin(   1)*yin(  43)*zin(   6)+xin(  49)*yin(  91)*zin(  54)+xin(  97)*yin( 139)*zin( 102)+xin( 145)*yin( 187)*zin( 150))
          eri_value(   82)=eri_value(   82)+d13bra(  5)*d12ket( 10)*(xin(   4)*yin(  45)*zin(   1)+xin(  52)*yin(  93)*zin(  49)+xin( 100)*yin( 141)*zin(  97)+xin( 148)*yin( 189)*zin( 145))
          eri_value(   83)=eri_value(   83)+d13bra(  5)*d12ket( 11)*(xin(   3)*yin(  46)*zin(   1)+xin(  51)*yin(  94)*zin(  49)+xin(  99)*yin( 142)*zin(  97)+xin( 147)*yin( 190)*zin( 145))
          eri_value(   84)=eri_value(   84)+d13bra(  5)*d12ket( 12)*(xin(   3)*yin(  45)*zin(   2)+xin(  51)*yin(  93)*zin(  50)+xin(  99)*yin( 141)*zin(  98)+xin( 147)*yin( 189)*zin( 146))
          eri_value(   85)=eri_value(   85)+d13bra(  5)*d12ket( 13)*(xin(   4)*yin(  43)*zin(   3)+xin(  52)*yin(  91)*zin(  51)+xin( 100)*yin( 139)*zin(  99)+xin( 148)*yin( 187)*zin( 147))
          eri_value(   86)=eri_value(   86)+d13bra(  5)*d12ket( 14)*(xin(   3)*yin(  44)*zin(   3)+xin(  51)*yin(  92)*zin(  51)+xin(  99)*yin( 140)*zin(  99)+xin( 147)*yin( 188)*zin( 147))
          eri_value(   87)=eri_value(   87)+d13bra(  5)*d12ket( 15)*(xin(   3)*yin(  43)*zin(   4)+xin(  51)*yin(  91)*zin(  52)+xin(  99)*yin( 139)*zin( 100)+xin( 147)*yin( 187)*zin( 148))
          eri_value(   88)=eri_value(   88)+d13bra(  5)*d12ket( 16)*(xin(   2)*yin(  45)*zin(   3)+xin(  50)*yin(  93)*zin(  51)+xin(  98)*yin( 141)*zin(  99)+xin( 146)*yin( 189)*zin( 147))
          eri_value(   89)=eri_value(   89)+d13bra(  5)*d12ket( 17)*(xin(   1)*yin(  46)*zin(   3)+xin(  49)*yin(  94)*zin(  51)+xin(  97)*yin( 142)*zin(  99)+xin( 145)*yin( 190)*zin( 147))
          eri_value(   90)=eri_value(   90)+d13bra(  5)*d12ket( 18)*(xin(   1)*yin(  45)*zin(   4)+xin(  49)*yin(  93)*zin(  52)+xin(  97)*yin( 141)*zin( 100)+xin( 145)*yin( 189)*zin( 148))
          eri_value(   91)=eri_value(   91)+d13bra(  6)*d12ket(  1)*(xin(   6)*yin(  37)*zin(   7)+xin(  54)*yin(  85)*zin(  55)+xin( 102)*yin( 133)*zin( 103)+xin( 150)*yin( 181)*zin( 151))
          eri_value(   92)=eri_value(   92)+d13bra(  6)*d12ket(  2)*(xin(   5)*yin(  38)*zin(   7)+xin(  53)*yin(  86)*zin(  55)+xin( 101)*yin( 134)*zin( 103)+xin( 149)*yin( 182)*zin( 151))
          eri_value(   93)=eri_value(   93)+d13bra(  6)*d12ket(  3)*(xin(   5)*yin(  37)*zin(   8)+xin(  53)*yin(  85)*zin(  56)+xin( 101)*yin( 133)*zin( 104)+xin( 149)*yin( 181)*zin( 152))
          eri_value(   94)=eri_value(   94)+d13bra(  6)*d12ket(  4)*(xin(   2)*yin(  41)*zin(   7)+xin(  50)*yin(  89)*zin(  55)+xin(  98)*yin( 137)*zin( 103)+xin( 146)*yin( 185)*zin( 151))
          eri_value(   95)=eri_value(   95)+d13bra(  6)*d12ket(  5)*(xin(   1)*yin(  42)*zin(   7)+xin(  49)*yin(  90)*zin(  55)+xin(  97)*yin( 138)*zin( 103)+xin( 145)*yin( 186)*zin( 151))
          eri_value(   96)=eri_value(   96)+d13bra(  6)*d12ket(  6)*(xin(   1)*yin(  41)*zin(   8)+xin(  49)*yin(  89)*zin(  56)+xin(  97)*yin( 137)*zin( 104)+xin( 145)*yin( 185)*zin( 152))
          eri_value(   97)=eri_value(   97)+d13bra(  6)*d12ket(  7)*(xin(   2)*yin(  37)*zin(  11)+xin(  50)*yin(  85)*zin(  59)+xin(  98)*yin( 133)*zin( 107)+xin( 146)*yin( 181)*zin( 155))
          eri_value(   98)=eri_value(   98)+d13bra(  6)*d12ket(  8)*(xin(   1)*yin(  38)*zin(  11)+xin(  49)*yin(  86)*zin(  59)+xin(  97)*yin( 134)*zin( 107)+xin( 145)*yin( 182)*zin( 155))
          eri_value(   99)=eri_value(   99)+d13bra(  6)*d12ket(  9)*(xin(   1)*yin(  37)*zin(  12)+xin(  49)*yin(  85)*zin(  60)+xin(  97)*yin( 133)*zin( 108)+xin( 145)*yin( 181)*zin( 156))
          eri_value(  100)=eri_value(  100)+d13bra(  6)*d12ket( 10)*(xin(   4)*yin(  39)*zin(   7)+xin(  52)*yin(  87)*zin(  55)+xin( 100)*yin( 135)*zin( 103)+xin( 148)*yin( 183)*zin( 151))
          eri_value(  101)=eri_value(  101)+d13bra(  6)*d12ket( 11)*(xin(   3)*yin(  40)*zin(   7)+xin(  51)*yin(  88)*zin(  55)+xin(  99)*yin( 136)*zin( 103)+xin( 147)*yin( 184)*zin( 151))
          eri_value(  102)=eri_value(  102)+d13bra(  6)*d12ket( 12)*(xin(   3)*yin(  39)*zin(   8)+xin(  51)*yin(  87)*zin(  56)+xin(  99)*yin( 135)*zin( 104)+xin( 147)*yin( 183)*zin( 152))
          eri_value(  103)=eri_value(  103)+d13bra(  6)*d12ket( 13)*(xin(   4)*yin(  37)*zin(   9)+xin(  52)*yin(  85)*zin(  57)+xin( 100)*yin( 133)*zin( 105)+xin( 148)*yin( 181)*zin( 153))
          eri_value(  104)=eri_value(  104)+d13bra(  6)*d12ket( 14)*(xin(   3)*yin(  38)*zin(   9)+xin(  51)*yin(  86)*zin(  57)+xin(  99)*yin( 134)*zin( 105)+xin( 147)*yin( 182)*zin( 153))
          eri_value(  105)=eri_value(  105)+d13bra(  6)*d12ket( 15)*(xin(   3)*yin(  37)*zin(  10)+xin(  51)*yin(  85)*zin(  58)+xin(  99)*yin( 133)*zin( 106)+xin( 147)*yin( 181)*zin( 154))
          eri_value(  106)=eri_value(  106)+d13bra(  6)*d12ket( 16)*(xin(   2)*yin(  39)*zin(   9)+xin(  50)*yin(  87)*zin(  57)+xin(  98)*yin( 135)*zin( 105)+xin( 146)*yin( 183)*zin( 153))
          eri_value(  107)=eri_value(  107)+d13bra(  6)*d12ket( 17)*(xin(   1)*yin(  40)*zin(   9)+xin(  49)*yin(  88)*zin(  57)+xin(  97)*yin( 136)*zin( 105)+xin( 145)*yin( 184)*zin( 153))
          eri_value(  108)=eri_value(  108)+d13bra(  6)*d12ket( 18)*(xin(   1)*yin(  39)*zin(  10)+xin(  49)*yin(  87)*zin(  58)+xin(  97)*yin( 135)*zin( 106)+xin( 145)*yin( 183)*zin( 154))
          eri_value(  109)=eri_value(  109)+d13bra(  7)*d12ket(  1)*(xin(  12)*yin(   1)*zin(  37)+xin(  60)*yin(  49)*zin(  85)+xin( 108)*yin(  97)*zin( 133)+xin( 156)*yin( 145)*zin( 181))
          eri_value(  110)=eri_value(  110)+d13bra(  7)*d12ket(  2)*(xin(  11)*yin(   2)*zin(  37)+xin(  59)*yin(  50)*zin(  85)+xin( 107)*yin(  98)*zin( 133)+xin( 155)*yin( 146)*zin( 181))
          eri_value(  111)=eri_value(  111)+d13bra(  7)*d12ket(  3)*(xin(  11)*yin(   1)*zin(  38)+xin(  59)*yin(  49)*zin(  86)+xin( 107)*yin(  97)*zin( 134)+xin( 155)*yin( 145)*zin( 182))
          eri_value(  112)=eri_value(  112)+d13bra(  7)*d12ket(  4)*(xin(   8)*yin(   5)*zin(  37)+xin(  56)*yin(  53)*zin(  85)+xin( 104)*yin( 101)*zin( 133)+xin( 152)*yin( 149)*zin( 181))
          eri_value(  113)=eri_value(  113)+d13bra(  7)*d12ket(  5)*(xin(   7)*yin(   6)*zin(  37)+xin(  55)*yin(  54)*zin(  85)+xin( 103)*yin( 102)*zin( 133)+xin( 151)*yin( 150)*zin( 181))
          eri_value(  114)=eri_value(  114)+d13bra(  7)*d12ket(  6)*(xin(   7)*yin(   5)*zin(  38)+xin(  55)*yin(  53)*zin(  86)+xin( 103)*yin( 101)*zin( 134)+xin( 151)*yin( 149)*zin( 182))
          eri_value(  115)=eri_value(  115)+d13bra(  7)*d12ket(  7)*(xin(   8)*yin(   1)*zin(  41)+xin(  56)*yin(  49)*zin(  89)+xin( 104)*yin(  97)*zin( 137)+xin( 152)*yin( 145)*zin( 185))
          eri_value(  116)=eri_value(  116)+d13bra(  7)*d12ket(  8)*(xin(   7)*yin(   2)*zin(  41)+xin(  55)*yin(  50)*zin(  89)+xin( 103)*yin(  98)*zin( 137)+xin( 151)*yin( 146)*zin( 185))
          eri_value(  117)=eri_value(  117)+d13bra(  7)*d12ket(  9)*(xin(   7)*yin(   1)*zin(  42)+xin(  55)*yin(  49)*zin(  90)+xin( 103)*yin(  97)*zin( 138)+xin( 151)*yin( 145)*zin( 186))
          eri_value(  118)=eri_value(  118)+d13bra(  7)*d12ket( 10)*(xin(  10)*yin(   3)*zin(  37)+xin(  58)*yin(  51)*zin(  85)+xin( 106)*yin(  99)*zin( 133)+xin( 154)*yin( 147)*zin( 181))
          eri_value(  119)=eri_value(  119)+d13bra(  7)*d12ket( 11)*(xin(   9)*yin(   4)*zin(  37)+xin(  57)*yin(  52)*zin(  85)+xin( 105)*yin( 100)*zin( 133)+xin( 153)*yin( 148)*zin( 181))
          eri_value(  120)=eri_value(  120)+d13bra(  7)*d12ket( 12)*(xin(   9)*yin(   3)*zin(  38)+xin(  57)*yin(  51)*zin(  86)+xin( 105)*yin(  99)*zin( 134)+xin( 153)*yin( 147)*zin( 182))
          eri_value(  121)=eri_value(  121)+d13bra(  7)*d12ket( 13)*(xin(  10)*yin(   1)*zin(  39)+xin(  58)*yin(  49)*zin(  87)+xin( 106)*yin(  97)*zin( 135)+xin( 154)*yin( 145)*zin( 183))
          eri_value(  122)=eri_value(  122)+d13bra(  7)*d12ket( 14)*(xin(   9)*yin(   2)*zin(  39)+xin(  57)*yin(  50)*zin(  87)+xin( 105)*yin(  98)*zin( 135)+xin( 153)*yin( 146)*zin( 183))
          eri_value(  123)=eri_value(  123)+d13bra(  7)*d12ket( 15)*(xin(   9)*yin(   1)*zin(  40)+xin(  57)*yin(  49)*zin(  88)+xin( 105)*yin(  97)*zin( 136)+xin( 153)*yin( 145)*zin( 184))
          eri_value(  124)=eri_value(  124)+d13bra(  7)*d12ket( 16)*(xin(   8)*yin(   3)*zin(  39)+xin(  56)*yin(  51)*zin(  87)+xin( 104)*yin(  99)*zin( 135)+xin( 152)*yin( 147)*zin( 183))
          eri_value(  125)=eri_value(  125)+d13bra(  7)*d12ket( 17)*(xin(   7)*yin(   4)*zin(  39)+xin(  55)*yin(  52)*zin(  87)+xin( 103)*yin( 100)*zin( 135)+xin( 151)*yin( 148)*zin( 183))
          eri_value(  126)=eri_value(  126)+d13bra(  7)*d12ket( 18)*(xin(   7)*yin(   3)*zin(  40)+xin(  55)*yin(  51)*zin(  88)+xin( 103)*yin(  99)*zin( 136)+xin( 151)*yin( 147)*zin( 184))
          eri_value(  127)=eri_value(  127)+d13bra(  8)*d12ket(  1)*(xin(   6)*yin(   7)*zin(  37)+xin(  54)*yin(  55)*zin(  85)+xin( 102)*yin( 103)*zin( 133)+xin( 150)*yin( 151)*zin( 181))
          eri_value(  128)=eri_value(  128)+d13bra(  8)*d12ket(  2)*(xin(   5)*yin(   8)*zin(  37)+xin(  53)*yin(  56)*zin(  85)+xin( 101)*yin( 104)*zin( 133)+xin( 149)*yin( 152)*zin( 181))
          eri_value(  129)=eri_value(  129)+d13bra(  8)*d12ket(  3)*(xin(   5)*yin(   7)*zin(  38)+xin(  53)*yin(  55)*zin(  86)+xin( 101)*yin( 103)*zin( 134)+xin( 149)*yin( 151)*zin( 182))
          eri_value(  130)=eri_value(  130)+d13bra(  8)*d12ket(  4)*(xin(   2)*yin(  11)*zin(  37)+xin(  50)*yin(  59)*zin(  85)+xin(  98)*yin( 107)*zin( 133)+xin( 146)*yin( 155)*zin( 181))
          eri_value(  131)=eri_value(  131)+d13bra(  8)*d12ket(  5)*(xin(   1)*yin(  12)*zin(  37)+xin(  49)*yin(  60)*zin(  85)+xin(  97)*yin( 108)*zin( 133)+xin( 145)*yin( 156)*zin( 181))
          eri_value(  132)=eri_value(  132)+d13bra(  8)*d12ket(  6)*(xin(   1)*yin(  11)*zin(  38)+xin(  49)*yin(  59)*zin(  86)+xin(  97)*yin( 107)*zin( 134)+xin( 145)*yin( 155)*zin( 182))
          eri_value(  133)=eri_value(  133)+d13bra(  8)*d12ket(  7)*(xin(   2)*yin(   7)*zin(  41)+xin(  50)*yin(  55)*zin(  89)+xin(  98)*yin( 103)*zin( 137)+xin( 146)*yin( 151)*zin( 185))
          eri_value(  134)=eri_value(  134)+d13bra(  8)*d12ket(  8)*(xin(   1)*yin(   8)*zin(  41)+xin(  49)*yin(  56)*zin(  89)+xin(  97)*yin( 104)*zin( 137)+xin( 145)*yin( 152)*zin( 185))
          eri_value(  135)=eri_value(  135)+d13bra(  8)*d12ket(  9)*(xin(   1)*yin(   7)*zin(  42)+xin(  49)*yin(  55)*zin(  90)+xin(  97)*yin( 103)*zin( 138)+xin( 145)*yin( 151)*zin( 186))
          eri_value(  136)=eri_value(  136)+d13bra(  8)*d12ket( 10)*(xin(   4)*yin(   9)*zin(  37)+xin(  52)*yin(  57)*zin(  85)+xin( 100)*yin( 105)*zin( 133)+xin( 148)*yin( 153)*zin( 181))
          eri_value(  137)=eri_value(  137)+d13bra(  8)*d12ket( 11)*(xin(   3)*yin(  10)*zin(  37)+xin(  51)*yin(  58)*zin(  85)+xin(  99)*yin( 106)*zin( 133)+xin( 147)*yin( 154)*zin( 181))
          eri_value(  138)=eri_value(  138)+d13bra(  8)*d12ket( 12)*(xin(   3)*yin(   9)*zin(  38)+xin(  51)*yin(  57)*zin(  86)+xin(  99)*yin( 105)*zin( 134)+xin( 147)*yin( 153)*zin( 182))
          eri_value(  139)=eri_value(  139)+d13bra(  8)*d12ket( 13)*(xin(   4)*yin(   7)*zin(  39)+xin(  52)*yin(  55)*zin(  87)+xin( 100)*yin( 103)*zin( 135)+xin( 148)*yin( 151)*zin( 183))
          eri_value(  140)=eri_value(  140)+d13bra(  8)*d12ket( 14)*(xin(   3)*yin(   8)*zin(  39)+xin(  51)*yin(  56)*zin(  87)+xin(  99)*yin( 104)*zin( 135)+xin( 147)*yin( 152)*zin( 183))
          eri_value(  141)=eri_value(  141)+d13bra(  8)*d12ket( 15)*(xin(   3)*yin(   7)*zin(  40)+xin(  51)*yin(  55)*zin(  88)+xin(  99)*yin( 103)*zin( 136)+xin( 147)*yin( 151)*zin( 184))
          eri_value(  142)=eri_value(  142)+d13bra(  8)*d12ket( 16)*(xin(   2)*yin(   9)*zin(  39)+xin(  50)*yin(  57)*zin(  87)+xin(  98)*yin( 105)*zin( 135)+xin( 146)*yin( 153)*zin( 183))
          eri_value(  143)=eri_value(  143)+d13bra(  8)*d12ket( 17)*(xin(   1)*yin(  10)*zin(  39)+xin(  49)*yin(  58)*zin(  87)+xin(  97)*yin( 106)*zin( 135)+xin( 145)*yin( 154)*zin( 183))
          eri_value(  144)=eri_value(  144)+d13bra(  8)*d12ket( 18)*(xin(   1)*yin(   9)*zin(  40)+xin(  49)*yin(  57)*zin(  88)+xin(  97)*yin( 105)*zin( 136)+xin( 145)*yin( 153)*zin( 184))
          eri_value(  145)=eri_value(  145)+d13bra(  9)*d12ket(  1)*(xin(   6)*yin(   1)*zin(  43)+xin(  54)*yin(  49)*zin(  91)+xin( 102)*yin(  97)*zin( 139)+xin( 150)*yin( 145)*zin( 187))
          eri_value(  146)=eri_value(  146)+d13bra(  9)*d12ket(  2)*(xin(   5)*yin(   2)*zin(  43)+xin(  53)*yin(  50)*zin(  91)+xin( 101)*yin(  98)*zin( 139)+xin( 149)*yin( 146)*zin( 187))
          eri_value(  147)=eri_value(  147)+d13bra(  9)*d12ket(  3)*(xin(   5)*yin(   1)*zin(  44)+xin(  53)*yin(  49)*zin(  92)+xin( 101)*yin(  97)*zin( 140)+xin( 149)*yin( 145)*zin( 188))
          eri_value(  148)=eri_value(  148)+d13bra(  9)*d12ket(  4)*(xin(   2)*yin(   5)*zin(  43)+xin(  50)*yin(  53)*zin(  91)+xin(  98)*yin( 101)*zin( 139)+xin( 146)*yin( 149)*zin( 187))
          eri_value(  149)=eri_value(  149)+d13bra(  9)*d12ket(  5)*(xin(   1)*yin(   6)*zin(  43)+xin(  49)*yin(  54)*zin(  91)+xin(  97)*yin( 102)*zin( 139)+xin( 145)*yin( 150)*zin( 187))
          eri_value(  150)=eri_value(  150)+d13bra(  9)*d12ket(  6)*(xin(   1)*yin(   5)*zin(  44)+xin(  49)*yin(  53)*zin(  92)+xin(  97)*yin( 101)*zin( 140)+xin( 145)*yin( 149)*zin( 188))
          eri_value(  151)=eri_value(  151)+d13bra(  9)*d12ket(  7)*(xin(   2)*yin(   1)*zin(  47)+xin(  50)*yin(  49)*zin(  95)+xin(  98)*yin(  97)*zin( 143)+xin( 146)*yin( 145)*zin( 191))
          eri_value(  152)=eri_value(  152)+d13bra(  9)*d12ket(  8)*(xin(   1)*yin(   2)*zin(  47)+xin(  49)*yin(  50)*zin(  95)+xin(  97)*yin(  98)*zin( 143)+xin( 145)*yin( 146)*zin( 191))
          eri_value(  153)=eri_value(  153)+d13bra(  9)*d12ket(  9)*(xin(   1)*yin(   1)*zin(  48)+xin(  49)*yin(  49)*zin(  96)+xin(  97)*yin(  97)*zin( 144)+xin( 145)*yin( 145)*zin( 192))
          eri_value(  154)=eri_value(  154)+d13bra(  9)*d12ket( 10)*(xin(   4)*yin(   3)*zin(  43)+xin(  52)*yin(  51)*zin(  91)+xin( 100)*yin(  99)*zin( 139)+xin( 148)*yin( 147)*zin( 187))
          eri_value(  155)=eri_value(  155)+d13bra(  9)*d12ket( 11)*(xin(   3)*yin(   4)*zin(  43)+xin(  51)*yin(  52)*zin(  91)+xin(  99)*yin( 100)*zin( 139)+xin( 147)*yin( 148)*zin( 187))
          eri_value(  156)=eri_value(  156)+d13bra(  9)*d12ket( 12)*(xin(   3)*yin(   3)*zin(  44)+xin(  51)*yin(  51)*zin(  92)+xin(  99)*yin(  99)*zin( 140)+xin( 147)*yin( 147)*zin( 188))
          eri_value(  157)=eri_value(  157)+d13bra(  9)*d12ket( 13)*(xin(   4)*yin(   1)*zin(  45)+xin(  52)*yin(  49)*zin(  93)+xin( 100)*yin(  97)*zin( 141)+xin( 148)*yin( 145)*zin( 189))
          eri_value(  158)=eri_value(  158)+d13bra(  9)*d12ket( 14)*(xin(   3)*yin(   2)*zin(  45)+xin(  51)*yin(  50)*zin(  93)+xin(  99)*yin(  98)*zin( 141)+xin( 147)*yin( 146)*zin( 189))
          eri_value(  159)=eri_value(  159)+d13bra(  9)*d12ket( 15)*(xin(   3)*yin(   1)*zin(  46)+xin(  51)*yin(  49)*zin(  94)+xin(  99)*yin(  97)*zin( 142)+xin( 147)*yin( 145)*zin( 190))
          eri_value(  160)=eri_value(  160)+d13bra(  9)*d12ket( 16)*(xin(   2)*yin(   3)*zin(  45)+xin(  50)*yin(  51)*zin(  93)+xin(  98)*yin(  99)*zin( 141)+xin( 146)*yin( 147)*zin( 189))
          eri_value(  161)=eri_value(  161)+d13bra(  9)*d12ket( 17)*(xin(   1)*yin(   4)*zin(  45)+xin(  49)*yin(  52)*zin(  93)+xin(  97)*yin( 100)*zin( 141)+xin( 145)*yin( 148)*zin( 189))
          eri_value(  162)=eri_value(  162)+d13bra(  9)*d12ket( 18)*(xin(   1)*yin(   3)*zin(  46)+xin(  49)*yin(  51)*zin(  94)+xin(  97)*yin(  99)*zin( 142)+xin( 145)*yin( 147)*zin( 190))
          eri_value(  163)=eri_value(  163)+d13bra( 10)*d12ket(  1)*(xin(  36)*yin(  13)*zin(   1)+xin(  84)*yin(  61)*zin(  49)+xin( 132)*yin( 109)*zin(  97)+xin( 180)*yin( 157)*zin( 145))
          eri_value(  164)=eri_value(  164)+d13bra( 10)*d12ket(  2)*(xin(  35)*yin(  14)*zin(   1)+xin(  83)*yin(  62)*zin(  49)+xin( 131)*yin( 110)*zin(  97)+xin( 179)*yin( 158)*zin( 145))
          eri_value(  165)=eri_value(  165)+d13bra( 10)*d12ket(  3)*(xin(  35)*yin(  13)*zin(   2)+xin(  83)*yin(  61)*zin(  50)+xin( 131)*yin( 109)*zin(  98)+xin( 179)*yin( 157)*zin( 146))
          eri_value(  166)=eri_value(  166)+d13bra( 10)*d12ket(  4)*(xin(  32)*yin(  17)*zin(   1)+xin(  80)*yin(  65)*zin(  49)+xin( 128)*yin( 113)*zin(  97)+xin( 176)*yin( 161)*zin( 145))
          eri_value(  167)=eri_value(  167)+d13bra( 10)*d12ket(  5)*(xin(  31)*yin(  18)*zin(   1)+xin(  79)*yin(  66)*zin(  49)+xin( 127)*yin( 114)*zin(  97)+xin( 175)*yin( 162)*zin( 145))
          eri_value(  168)=eri_value(  168)+d13bra( 10)*d12ket(  6)*(xin(  31)*yin(  17)*zin(   2)+xin(  79)*yin(  65)*zin(  50)+xin( 127)*yin( 113)*zin(  98)+xin( 175)*yin( 161)*zin( 146))
          eri_value(  169)=eri_value(  169)+d13bra( 10)*d12ket(  7)*(xin(  32)*yin(  13)*zin(   5)+xin(  80)*yin(  61)*zin(  53)+xin( 128)*yin( 109)*zin( 101)+xin( 176)*yin( 157)*zin( 149))
          eri_value(  170)=eri_value(  170)+d13bra( 10)*d12ket(  8)*(xin(  31)*yin(  14)*zin(   5)+xin(  79)*yin(  62)*zin(  53)+xin( 127)*yin( 110)*zin( 101)+xin( 175)*yin( 158)*zin( 149))
          eri_value(  171)=eri_value(  171)+d13bra( 10)*d12ket(  9)*(xin(  31)*yin(  13)*zin(   6)+xin(  79)*yin(  61)*zin(  54)+xin( 127)*yin( 109)*zin( 102)+xin( 175)*yin( 157)*zin( 150))
          eri_value(  172)=eri_value(  172)+d13bra( 10)*d12ket( 10)*(xin(  34)*yin(  15)*zin(   1)+xin(  82)*yin(  63)*zin(  49)+xin( 130)*yin( 111)*zin(  97)+xin( 178)*yin( 159)*zin( 145))
          eri_value(  173)=eri_value(  173)+d13bra( 10)*d12ket( 11)*(xin(  33)*yin(  16)*zin(   1)+xin(  81)*yin(  64)*zin(  49)+xin( 129)*yin( 112)*zin(  97)+xin( 177)*yin( 160)*zin( 145))
          eri_value(  174)=eri_value(  174)+d13bra( 10)*d12ket( 12)*(xin(  33)*yin(  15)*zin(   2)+xin(  81)*yin(  63)*zin(  50)+xin( 129)*yin( 111)*zin(  98)+xin( 177)*yin( 159)*zin( 146))
          eri_value(  175)=eri_value(  175)+d13bra( 10)*d12ket( 13)*(xin(  34)*yin(  13)*zin(   3)+xin(  82)*yin(  61)*zin(  51)+xin( 130)*yin( 109)*zin(  99)+xin( 178)*yin( 157)*zin( 147))
          eri_value(  176)=eri_value(  176)+d13bra( 10)*d12ket( 14)*(xin(  33)*yin(  14)*zin(   3)+xin(  81)*yin(  62)*zin(  51)+xin( 129)*yin( 110)*zin(  99)+xin( 177)*yin( 158)*zin( 147))
          eri_value(  177)=eri_value(  177)+d13bra( 10)*d12ket( 15)*(xin(  33)*yin(  13)*zin(   4)+xin(  81)*yin(  61)*zin(  52)+xin( 129)*yin( 109)*zin( 100)+xin( 177)*yin( 157)*zin( 148))
          eri_value(  178)=eri_value(  178)+d13bra( 10)*d12ket( 16)*(xin(  32)*yin(  15)*zin(   3)+xin(  80)*yin(  63)*zin(  51)+xin( 128)*yin( 111)*zin(  99)+xin( 176)*yin( 159)*zin( 147))
          eri_value(  179)=eri_value(  179)+d13bra( 10)*d12ket( 17)*(xin(  31)*yin(  16)*zin(   3)+xin(  79)*yin(  64)*zin(  51)+xin( 127)*yin( 112)*zin(  99)+xin( 175)*yin( 160)*zin( 147))
          eri_value(  180)=eri_value(  180)+d13bra( 10)*d12ket( 18)*(xin(  31)*yin(  15)*zin(   4)+xin(  79)*yin(  63)*zin(  52)+xin( 127)*yin( 111)*zin( 100)+xin( 175)*yin( 159)*zin( 148))
          eri_value(  181)=eri_value(  181)+d13bra( 11)*d12ket(  1)*(xin(  30)*yin(  19)*zin(   1)+xin(  78)*yin(  67)*zin(  49)+xin( 126)*yin( 115)*zin(  97)+xin( 174)*yin( 163)*zin( 145))
          eri_value(  182)=eri_value(  182)+d13bra( 11)*d12ket(  2)*(xin(  29)*yin(  20)*zin(   1)+xin(  77)*yin(  68)*zin(  49)+xin( 125)*yin( 116)*zin(  97)+xin( 173)*yin( 164)*zin( 145))
          eri_value(  183)=eri_value(  183)+d13bra( 11)*d12ket(  3)*(xin(  29)*yin(  19)*zin(   2)+xin(  77)*yin(  67)*zin(  50)+xin( 125)*yin( 115)*zin(  98)+xin( 173)*yin( 163)*zin( 146))
          eri_value(  184)=eri_value(  184)+d13bra( 11)*d12ket(  4)*(xin(  26)*yin(  23)*zin(   1)+xin(  74)*yin(  71)*zin(  49)+xin( 122)*yin( 119)*zin(  97)+xin( 170)*yin( 167)*zin( 145))
          eri_value(  185)=eri_value(  185)+d13bra( 11)*d12ket(  5)*(xin(  25)*yin(  24)*zin(   1)+xin(  73)*yin(  72)*zin(  49)+xin( 121)*yin( 120)*zin(  97)+xin( 169)*yin( 168)*zin( 145))
          eri_value(  186)=eri_value(  186)+d13bra( 11)*d12ket(  6)*(xin(  25)*yin(  23)*zin(   2)+xin(  73)*yin(  71)*zin(  50)+xin( 121)*yin( 119)*zin(  98)+xin( 169)*yin( 167)*zin( 146))
          eri_value(  187)=eri_value(  187)+d13bra( 11)*d12ket(  7)*(xin(  26)*yin(  19)*zin(   5)+xin(  74)*yin(  67)*zin(  53)+xin( 122)*yin( 115)*zin( 101)+xin( 170)*yin( 163)*zin( 149))
          eri_value(  188)=eri_value(  188)+d13bra( 11)*d12ket(  8)*(xin(  25)*yin(  20)*zin(   5)+xin(  73)*yin(  68)*zin(  53)+xin( 121)*yin( 116)*zin( 101)+xin( 169)*yin( 164)*zin( 149))
          eri_value(  189)=eri_value(  189)+d13bra( 11)*d12ket(  9)*(xin(  25)*yin(  19)*zin(   6)+xin(  73)*yin(  67)*zin(  54)+xin( 121)*yin( 115)*zin( 102)+xin( 169)*yin( 163)*zin( 150))
          eri_value(  190)=eri_value(  190)+d13bra( 11)*d12ket( 10)*(xin(  28)*yin(  21)*zin(   1)+xin(  76)*yin(  69)*zin(  49)+xin( 124)*yin( 117)*zin(  97)+xin( 172)*yin( 165)*zin( 145))
          eri_value(  191)=eri_value(  191)+d13bra( 11)*d12ket( 11)*(xin(  27)*yin(  22)*zin(   1)+xin(  75)*yin(  70)*zin(  49)+xin( 123)*yin( 118)*zin(  97)+xin( 171)*yin( 166)*zin( 145))
          eri_value(  192)=eri_value(  192)+d13bra( 11)*d12ket( 12)*(xin(  27)*yin(  21)*zin(   2)+xin(  75)*yin(  69)*zin(  50)+xin( 123)*yin( 117)*zin(  98)+xin( 171)*yin( 165)*zin( 146))
          eri_value(  193)=eri_value(  193)+d13bra( 11)*d12ket( 13)*(xin(  28)*yin(  19)*zin(   3)+xin(  76)*yin(  67)*zin(  51)+xin( 124)*yin( 115)*zin(  99)+xin( 172)*yin( 163)*zin( 147))
          eri_value(  194)=eri_value(  194)+d13bra( 11)*d12ket( 14)*(xin(  27)*yin(  20)*zin(   3)+xin(  75)*yin(  68)*zin(  51)+xin( 123)*yin( 116)*zin(  99)+xin( 171)*yin( 164)*zin( 147))
          eri_value(  195)=eri_value(  195)+d13bra( 11)*d12ket( 15)*(xin(  27)*yin(  19)*zin(   4)+xin(  75)*yin(  67)*zin(  52)+xin( 123)*yin( 115)*zin( 100)+xin( 171)*yin( 163)*zin( 148))
          eri_value(  196)=eri_value(  196)+d13bra( 11)*d12ket( 16)*(xin(  26)*yin(  21)*zin(   3)+xin(  74)*yin(  69)*zin(  51)+xin( 122)*yin( 117)*zin(  99)+xin( 170)*yin( 165)*zin( 147))
          eri_value(  197)=eri_value(  197)+d13bra( 11)*d12ket( 17)*(xin(  25)*yin(  22)*zin(   3)+xin(  73)*yin(  70)*zin(  51)+xin( 121)*yin( 118)*zin(  99)+xin( 169)*yin( 166)*zin( 147))
          eri_value(  198)=eri_value(  198)+d13bra( 11)*d12ket( 18)*(xin(  25)*yin(  21)*zin(   4)+xin(  73)*yin(  69)*zin(  52)+xin( 121)*yin( 117)*zin( 100)+xin( 169)*yin( 165)*zin( 148))
          eri_value(  199)=eri_value(  199)+d13bra( 12)*d12ket(  1)*(xin(  30)*yin(  13)*zin(   7)+xin(  78)*yin(  61)*zin(  55)+xin( 126)*yin( 109)*zin( 103)+xin( 174)*yin( 157)*zin( 151))
          eri_value(  200)=eri_value(  200)+d13bra( 12)*d12ket(  2)*(xin(  29)*yin(  14)*zin(   7)+xin(  77)*yin(  62)*zin(  55)+xin( 125)*yin( 110)*zin( 103)+xin( 173)*yin( 158)*zin( 151))
          eri_value(  201)=eri_value(  201)+d13bra( 12)*d12ket(  3)*(xin(  29)*yin(  13)*zin(   8)+xin(  77)*yin(  61)*zin(  56)+xin( 125)*yin( 109)*zin( 104)+xin( 173)*yin( 157)*zin( 152))
          eri_value(  202)=eri_value(  202)+d13bra( 12)*d12ket(  4)*(xin(  26)*yin(  17)*zin(   7)+xin(  74)*yin(  65)*zin(  55)+xin( 122)*yin( 113)*zin( 103)+xin( 170)*yin( 161)*zin( 151))
          eri_value(  203)=eri_value(  203)+d13bra( 12)*d12ket(  5)*(xin(  25)*yin(  18)*zin(   7)+xin(  73)*yin(  66)*zin(  55)+xin( 121)*yin( 114)*zin( 103)+xin( 169)*yin( 162)*zin( 151))
          eri_value(  204)=eri_value(  204)+d13bra( 12)*d12ket(  6)*(xin(  25)*yin(  17)*zin(   8)+xin(  73)*yin(  65)*zin(  56)+xin( 121)*yin( 113)*zin( 104)+xin( 169)*yin( 161)*zin( 152))
          eri_value(  205)=eri_value(  205)+d13bra( 12)*d12ket(  7)*(xin(  26)*yin(  13)*zin(  11)+xin(  74)*yin(  61)*zin(  59)+xin( 122)*yin( 109)*zin( 107)+xin( 170)*yin( 157)*zin( 155))
          eri_value(  206)=eri_value(  206)+d13bra( 12)*d12ket(  8)*(xin(  25)*yin(  14)*zin(  11)+xin(  73)*yin(  62)*zin(  59)+xin( 121)*yin( 110)*zin( 107)+xin( 169)*yin( 158)*zin( 155))
          eri_value(  207)=eri_value(  207)+d13bra( 12)*d12ket(  9)*(xin(  25)*yin(  13)*zin(  12)+xin(  73)*yin(  61)*zin(  60)+xin( 121)*yin( 109)*zin( 108)+xin( 169)*yin( 157)*zin( 156))
          eri_value(  208)=eri_value(  208)+d13bra( 12)*d12ket( 10)*(xin(  28)*yin(  15)*zin(   7)+xin(  76)*yin(  63)*zin(  55)+xin( 124)*yin( 111)*zin( 103)+xin( 172)*yin( 159)*zin( 151))
          eri_value(  209)=eri_value(  209)+d13bra( 12)*d12ket( 11)*(xin(  27)*yin(  16)*zin(   7)+xin(  75)*yin(  64)*zin(  55)+xin( 123)*yin( 112)*zin( 103)+xin( 171)*yin( 160)*zin( 151))
          eri_value(  210)=eri_value(  210)+d13bra( 12)*d12ket( 12)*(xin(  27)*yin(  15)*zin(   8)+xin(  75)*yin(  63)*zin(  56)+xin( 123)*yin( 111)*zin( 104)+xin( 171)*yin( 159)*zin( 152))
          eri_value(  211)=eri_value(  211)+d13bra( 12)*d12ket( 13)*(xin(  28)*yin(  13)*zin(   9)+xin(  76)*yin(  61)*zin(  57)+xin( 124)*yin( 109)*zin( 105)+xin( 172)*yin( 157)*zin( 153))
          eri_value(  212)=eri_value(  212)+d13bra( 12)*d12ket( 14)*(xin(  27)*yin(  14)*zin(   9)+xin(  75)*yin(  62)*zin(  57)+xin( 123)*yin( 110)*zin( 105)+xin( 171)*yin( 158)*zin( 153))
          eri_value(  213)=eri_value(  213)+d13bra( 12)*d12ket( 15)*(xin(  27)*yin(  13)*zin(  10)+xin(  75)*yin(  61)*zin(  58)+xin( 123)*yin( 109)*zin( 106)+xin( 171)*yin( 157)*zin( 154))
          eri_value(  214)=eri_value(  214)+d13bra( 12)*d12ket( 16)*(xin(  26)*yin(  15)*zin(   9)+xin(  74)*yin(  63)*zin(  57)+xin( 122)*yin( 111)*zin( 105)+xin( 170)*yin( 159)*zin( 153))
          eri_value(  215)=eri_value(  215)+d13bra( 12)*d12ket( 17)*(xin(  25)*yin(  16)*zin(   9)+xin(  73)*yin(  64)*zin(  57)+xin( 121)*yin( 112)*zin( 105)+xin( 169)*yin( 160)*zin( 153))
          eri_value(  216)=eri_value(  216)+d13bra( 12)*d12ket( 18)*(xin(  25)*yin(  15)*zin(  10)+xin(  73)*yin(  63)*zin(  58)+xin( 121)*yin( 111)*zin( 106)+xin( 169)*yin( 159)*zin( 154))
          eri_value(  217)=eri_value(  217)+d13bra( 13)*d12ket(  1)*(xin(  36)*yin(   1)*zin(  13)+xin(  84)*yin(  49)*zin(  61)+xin( 132)*yin(  97)*zin( 109)+xin( 180)*yin( 145)*zin( 157))
          eri_value(  218)=eri_value(  218)+d13bra( 13)*d12ket(  2)*(xin(  35)*yin(   2)*zin(  13)+xin(  83)*yin(  50)*zin(  61)+xin( 131)*yin(  98)*zin( 109)+xin( 179)*yin( 146)*zin( 157))
          eri_value(  219)=eri_value(  219)+d13bra( 13)*d12ket(  3)*(xin(  35)*yin(   1)*zin(  14)+xin(  83)*yin(  49)*zin(  62)+xin( 131)*yin(  97)*zin( 110)+xin( 179)*yin( 145)*zin( 158))
          eri_value(  220)=eri_value(  220)+d13bra( 13)*d12ket(  4)*(xin(  32)*yin(   5)*zin(  13)+xin(  80)*yin(  53)*zin(  61)+xin( 128)*yin( 101)*zin( 109)+xin( 176)*yin( 149)*zin( 157))
          eri_value(  221)=eri_value(  221)+d13bra( 13)*d12ket(  5)*(xin(  31)*yin(   6)*zin(  13)+xin(  79)*yin(  54)*zin(  61)+xin( 127)*yin( 102)*zin( 109)+xin( 175)*yin( 150)*zin( 157))
          eri_value(  222)=eri_value(  222)+d13bra( 13)*d12ket(  6)*(xin(  31)*yin(   5)*zin(  14)+xin(  79)*yin(  53)*zin(  62)+xin( 127)*yin( 101)*zin( 110)+xin( 175)*yin( 149)*zin( 158))
          eri_value(  223)=eri_value(  223)+d13bra( 13)*d12ket(  7)*(xin(  32)*yin(   1)*zin(  17)+xin(  80)*yin(  49)*zin(  65)+xin( 128)*yin(  97)*zin( 113)+xin( 176)*yin( 145)*zin( 161))
          eri_value(  224)=eri_value(  224)+d13bra( 13)*d12ket(  8)*(xin(  31)*yin(   2)*zin(  17)+xin(  79)*yin(  50)*zin(  65)+xin( 127)*yin(  98)*zin( 113)+xin( 175)*yin( 146)*zin( 161))
          eri_value(  225)=eri_value(  225)+d13bra( 13)*d12ket(  9)*(xin(  31)*yin(   1)*zin(  18)+xin(  79)*yin(  49)*zin(  66)+xin( 127)*yin(  97)*zin( 114)+xin( 175)*yin( 145)*zin( 162))
          eri_value(  226)=eri_value(  226)+d13bra( 13)*d12ket( 10)*(xin(  34)*yin(   3)*zin(  13)+xin(  82)*yin(  51)*zin(  61)+xin( 130)*yin(  99)*zin( 109)+xin( 178)*yin( 147)*zin( 157))
          eri_value(  227)=eri_value(  227)+d13bra( 13)*d12ket( 11)*(xin(  33)*yin(   4)*zin(  13)+xin(  81)*yin(  52)*zin(  61)+xin( 129)*yin( 100)*zin( 109)+xin( 177)*yin( 148)*zin( 157))
          eri_value(  228)=eri_value(  228)+d13bra( 13)*d12ket( 12)*(xin(  33)*yin(   3)*zin(  14)+xin(  81)*yin(  51)*zin(  62)+xin( 129)*yin(  99)*zin( 110)+xin( 177)*yin( 147)*zin( 158))
          eri_value(  229)=eri_value(  229)+d13bra( 13)*d12ket( 13)*(xin(  34)*yin(   1)*zin(  15)+xin(  82)*yin(  49)*zin(  63)+xin( 130)*yin(  97)*zin( 111)+xin( 178)*yin( 145)*zin( 159))
          eri_value(  230)=eri_value(  230)+d13bra( 13)*d12ket( 14)*(xin(  33)*yin(   2)*zin(  15)+xin(  81)*yin(  50)*zin(  63)+xin( 129)*yin(  98)*zin( 111)+xin( 177)*yin( 146)*zin( 159))
          eri_value(  231)=eri_value(  231)+d13bra( 13)*d12ket( 15)*(xin(  33)*yin(   1)*zin(  16)+xin(  81)*yin(  49)*zin(  64)+xin( 129)*yin(  97)*zin( 112)+xin( 177)*yin( 145)*zin( 160))
          eri_value(  232)=eri_value(  232)+d13bra( 13)*d12ket( 16)*(xin(  32)*yin(   3)*zin(  15)+xin(  80)*yin(  51)*zin(  63)+xin( 128)*yin(  99)*zin( 111)+xin( 176)*yin( 147)*zin( 159))
          eri_value(  233)=eri_value(  233)+d13bra( 13)*d12ket( 17)*(xin(  31)*yin(   4)*zin(  15)+xin(  79)*yin(  52)*zin(  63)+xin( 127)*yin( 100)*zin( 111)+xin( 175)*yin( 148)*zin( 159))
          eri_value(  234)=eri_value(  234)+d13bra( 13)*d12ket( 18)*(xin(  31)*yin(   3)*zin(  16)+xin(  79)*yin(  51)*zin(  64)+xin( 127)*yin(  99)*zin( 112)+xin( 175)*yin( 147)*zin( 160))
          eri_value(  235)=eri_value(  235)+d13bra( 14)*d12ket(  1)*(xin(  30)*yin(   7)*zin(  13)+xin(  78)*yin(  55)*zin(  61)+xin( 126)*yin( 103)*zin( 109)+xin( 174)*yin( 151)*zin( 157))
          eri_value(  236)=eri_value(  236)+d13bra( 14)*d12ket(  2)*(xin(  29)*yin(   8)*zin(  13)+xin(  77)*yin(  56)*zin(  61)+xin( 125)*yin( 104)*zin( 109)+xin( 173)*yin( 152)*zin( 157))
          eri_value(  237)=eri_value(  237)+d13bra( 14)*d12ket(  3)*(xin(  29)*yin(   7)*zin(  14)+xin(  77)*yin(  55)*zin(  62)+xin( 125)*yin( 103)*zin( 110)+xin( 173)*yin( 151)*zin( 158))
          eri_value(  238)=eri_value(  238)+d13bra( 14)*d12ket(  4)*(xin(  26)*yin(  11)*zin(  13)+xin(  74)*yin(  59)*zin(  61)+xin( 122)*yin( 107)*zin( 109)+xin( 170)*yin( 155)*zin( 157))
          eri_value(  239)=eri_value(  239)+d13bra( 14)*d12ket(  5)*(xin(  25)*yin(  12)*zin(  13)+xin(  73)*yin(  60)*zin(  61)+xin( 121)*yin( 108)*zin( 109)+xin( 169)*yin( 156)*zin( 157))
          eri_value(  240)=eri_value(  240)+d13bra( 14)*d12ket(  6)*(xin(  25)*yin(  11)*zin(  14)+xin(  73)*yin(  59)*zin(  62)+xin( 121)*yin( 107)*zin( 110)+xin( 169)*yin( 155)*zin( 158))
          eri_value(  241)=eri_value(  241)+d13bra( 14)*d12ket(  7)*(xin(  26)*yin(   7)*zin(  17)+xin(  74)*yin(  55)*zin(  65)+xin( 122)*yin( 103)*zin( 113)+xin( 170)*yin( 151)*zin( 161))
          eri_value(  242)=eri_value(  242)+d13bra( 14)*d12ket(  8)*(xin(  25)*yin(   8)*zin(  17)+xin(  73)*yin(  56)*zin(  65)+xin( 121)*yin( 104)*zin( 113)+xin( 169)*yin( 152)*zin( 161))
          eri_value(  243)=eri_value(  243)+d13bra( 14)*d12ket(  9)*(xin(  25)*yin(   7)*zin(  18)+xin(  73)*yin(  55)*zin(  66)+xin( 121)*yin( 103)*zin( 114)+xin( 169)*yin( 151)*zin( 162))
          eri_value(  244)=eri_value(  244)+d13bra( 14)*d12ket( 10)*(xin(  28)*yin(   9)*zin(  13)+xin(  76)*yin(  57)*zin(  61)+xin( 124)*yin( 105)*zin( 109)+xin( 172)*yin( 153)*zin( 157))
          eri_value(  245)=eri_value(  245)+d13bra( 14)*d12ket( 11)*(xin(  27)*yin(  10)*zin(  13)+xin(  75)*yin(  58)*zin(  61)+xin( 123)*yin( 106)*zin( 109)+xin( 171)*yin( 154)*zin( 157))
          eri_value(  246)=eri_value(  246)+d13bra( 14)*d12ket( 12)*(xin(  27)*yin(   9)*zin(  14)+xin(  75)*yin(  57)*zin(  62)+xin( 123)*yin( 105)*zin( 110)+xin( 171)*yin( 153)*zin( 158))
          eri_value(  247)=eri_value(  247)+d13bra( 14)*d12ket( 13)*(xin(  28)*yin(   7)*zin(  15)+xin(  76)*yin(  55)*zin(  63)+xin( 124)*yin( 103)*zin( 111)+xin( 172)*yin( 151)*zin( 159))
          eri_value(  248)=eri_value(  248)+d13bra( 14)*d12ket( 14)*(xin(  27)*yin(   8)*zin(  15)+xin(  75)*yin(  56)*zin(  63)+xin( 123)*yin( 104)*zin( 111)+xin( 171)*yin( 152)*zin( 159))
          eri_value(  249)=eri_value(  249)+d13bra( 14)*d12ket( 15)*(xin(  27)*yin(   7)*zin(  16)+xin(  75)*yin(  55)*zin(  64)+xin( 123)*yin( 103)*zin( 112)+xin( 171)*yin( 151)*zin( 160))
          eri_value(  250)=eri_value(  250)+d13bra( 14)*d12ket( 16)*(xin(  26)*yin(   9)*zin(  15)+xin(  74)*yin(  57)*zin(  63)+xin( 122)*yin( 105)*zin( 111)+xin( 170)*yin( 153)*zin( 159))
          eri_value(  251)=eri_value(  251)+d13bra( 14)*d12ket( 17)*(xin(  25)*yin(  10)*zin(  15)+xin(  73)*yin(  58)*zin(  63)+xin( 121)*yin( 106)*zin( 111)+xin( 169)*yin( 154)*zin( 159))
          eri_value(  252)=eri_value(  252)+d13bra( 14)*d12ket( 18)*(xin(  25)*yin(   9)*zin(  16)+xin(  73)*yin(  57)*zin(  64)+xin( 121)*yin( 105)*zin( 112)+xin( 169)*yin( 153)*zin( 160))
          eri_value(  253)=eri_value(  253)+d13bra( 15)*d12ket(  1)*(xin(  30)*yin(   1)*zin(  19)+xin(  78)*yin(  49)*zin(  67)+xin( 126)*yin(  97)*zin( 115)+xin( 174)*yin( 145)*zin( 163))
          eri_value(  254)=eri_value(  254)+d13bra( 15)*d12ket(  2)*(xin(  29)*yin(   2)*zin(  19)+xin(  77)*yin(  50)*zin(  67)+xin( 125)*yin(  98)*zin( 115)+xin( 173)*yin( 146)*zin( 163))
          eri_value(  255)=eri_value(  255)+d13bra( 15)*d12ket(  3)*(xin(  29)*yin(   1)*zin(  20)+xin(  77)*yin(  49)*zin(  68)+xin( 125)*yin(  97)*zin( 116)+xin( 173)*yin( 145)*zin( 164))
          eri_value(  256)=eri_value(  256)+d13bra( 15)*d12ket(  4)*(xin(  26)*yin(   5)*zin(  19)+xin(  74)*yin(  53)*zin(  67)+xin( 122)*yin( 101)*zin( 115)+xin( 170)*yin( 149)*zin( 163))
          eri_value(  257)=eri_value(  257)+d13bra( 15)*d12ket(  5)*(xin(  25)*yin(   6)*zin(  19)+xin(  73)*yin(  54)*zin(  67)+xin( 121)*yin( 102)*zin( 115)+xin( 169)*yin( 150)*zin( 163))
          eri_value(  258)=eri_value(  258)+d13bra( 15)*d12ket(  6)*(xin(  25)*yin(   5)*zin(  20)+xin(  73)*yin(  53)*zin(  68)+xin( 121)*yin( 101)*zin( 116)+xin( 169)*yin( 149)*zin( 164))
          eri_value(  259)=eri_value(  259)+d13bra( 15)*d12ket(  7)*(xin(  26)*yin(   1)*zin(  23)+xin(  74)*yin(  49)*zin(  71)+xin( 122)*yin(  97)*zin( 119)+xin( 170)*yin( 145)*zin( 167))
          eri_value(  260)=eri_value(  260)+d13bra( 15)*d12ket(  8)*(xin(  25)*yin(   2)*zin(  23)+xin(  73)*yin(  50)*zin(  71)+xin( 121)*yin(  98)*zin( 119)+xin( 169)*yin( 146)*zin( 167))
          eri_value(  261)=eri_value(  261)+d13bra( 15)*d12ket(  9)*(xin(  25)*yin(   1)*zin(  24)+xin(  73)*yin(  49)*zin(  72)+xin( 121)*yin(  97)*zin( 120)+xin( 169)*yin( 145)*zin( 168))
          eri_value(  262)=eri_value(  262)+d13bra( 15)*d12ket( 10)*(xin(  28)*yin(   3)*zin(  19)+xin(  76)*yin(  51)*zin(  67)+xin( 124)*yin(  99)*zin( 115)+xin( 172)*yin( 147)*zin( 163))
          eri_value(  263)=eri_value(  263)+d13bra( 15)*d12ket( 11)*(xin(  27)*yin(   4)*zin(  19)+xin(  75)*yin(  52)*zin(  67)+xin( 123)*yin( 100)*zin( 115)+xin( 171)*yin( 148)*zin( 163))
          eri_value(  264)=eri_value(  264)+d13bra( 15)*d12ket( 12)*(xin(  27)*yin(   3)*zin(  20)+xin(  75)*yin(  51)*zin(  68)+xin( 123)*yin(  99)*zin( 116)+xin( 171)*yin( 147)*zin( 164))
          eri_value(  265)=eri_value(  265)+d13bra( 15)*d12ket( 13)*(xin(  28)*yin(   1)*zin(  21)+xin(  76)*yin(  49)*zin(  69)+xin( 124)*yin(  97)*zin( 117)+xin( 172)*yin( 145)*zin( 165))
          eri_value(  266)=eri_value(  266)+d13bra( 15)*d12ket( 14)*(xin(  27)*yin(   2)*zin(  21)+xin(  75)*yin(  50)*zin(  69)+xin( 123)*yin(  98)*zin( 117)+xin( 171)*yin( 146)*zin( 165))
          eri_value(  267)=eri_value(  267)+d13bra( 15)*d12ket( 15)*(xin(  27)*yin(   1)*zin(  22)+xin(  75)*yin(  49)*zin(  70)+xin( 123)*yin(  97)*zin( 118)+xin( 171)*yin( 145)*zin( 166))
          eri_value(  268)=eri_value(  268)+d13bra( 15)*d12ket( 16)*(xin(  26)*yin(   3)*zin(  21)+xin(  74)*yin(  51)*zin(  69)+xin( 122)*yin(  99)*zin( 117)+xin( 170)*yin( 147)*zin( 165))
          eri_value(  269)=eri_value(  269)+d13bra( 15)*d12ket( 17)*(xin(  25)*yin(   4)*zin(  21)+xin(  73)*yin(  52)*zin(  69)+xin( 121)*yin( 100)*zin( 117)+xin( 169)*yin( 148)*zin( 165))
          eri_value(  270)=eri_value(  270)+d13bra( 15)*d12ket( 18)*(xin(  25)*yin(   3)*zin(  22)+xin(  73)*yin(  51)*zin(  70)+xin( 121)*yin(  99)*zin( 118)+xin( 169)*yin( 147)*zin( 166))
          eri_value(  271)=eri_value(  271)+d13bra( 16)*d12ket(  1)*(xin(  24)*yin(  25)*zin(   1)+xin(  72)*yin(  73)*zin(  49)+xin( 120)*yin( 121)*zin(  97)+xin( 168)*yin( 169)*zin( 145))
          eri_value(  272)=eri_value(  272)+d13bra( 16)*d12ket(  2)*(xin(  23)*yin(  26)*zin(   1)+xin(  71)*yin(  74)*zin(  49)+xin( 119)*yin( 122)*zin(  97)+xin( 167)*yin( 170)*zin( 145))
          eri_value(  273)=eri_value(  273)+d13bra( 16)*d12ket(  3)*(xin(  23)*yin(  25)*zin(   2)+xin(  71)*yin(  73)*zin(  50)+xin( 119)*yin( 121)*zin(  98)+xin( 167)*yin( 169)*zin( 146))
          eri_value(  274)=eri_value(  274)+d13bra( 16)*d12ket(  4)*(xin(  20)*yin(  29)*zin(   1)+xin(  68)*yin(  77)*zin(  49)+xin( 116)*yin( 125)*zin(  97)+xin( 164)*yin( 173)*zin( 145))
          eri_value(  275)=eri_value(  275)+d13bra( 16)*d12ket(  5)*(xin(  19)*yin(  30)*zin(   1)+xin(  67)*yin(  78)*zin(  49)+xin( 115)*yin( 126)*zin(  97)+xin( 163)*yin( 174)*zin( 145))
          eri_value(  276)=eri_value(  276)+d13bra( 16)*d12ket(  6)*(xin(  19)*yin(  29)*zin(   2)+xin(  67)*yin(  77)*zin(  50)+xin( 115)*yin( 125)*zin(  98)+xin( 163)*yin( 173)*zin( 146))
          eri_value(  277)=eri_value(  277)+d13bra( 16)*d12ket(  7)*(xin(  20)*yin(  25)*zin(   5)+xin(  68)*yin(  73)*zin(  53)+xin( 116)*yin( 121)*zin( 101)+xin( 164)*yin( 169)*zin( 149))
          eri_value(  278)=eri_value(  278)+d13bra( 16)*d12ket(  8)*(xin(  19)*yin(  26)*zin(   5)+xin(  67)*yin(  74)*zin(  53)+xin( 115)*yin( 122)*zin( 101)+xin( 163)*yin( 170)*zin( 149))
          eri_value(  279)=eri_value(  279)+d13bra( 16)*d12ket(  9)*(xin(  19)*yin(  25)*zin(   6)+xin(  67)*yin(  73)*zin(  54)+xin( 115)*yin( 121)*zin( 102)+xin( 163)*yin( 169)*zin( 150))
          eri_value(  280)=eri_value(  280)+d13bra( 16)*d12ket( 10)*(xin(  22)*yin(  27)*zin(   1)+xin(  70)*yin(  75)*zin(  49)+xin( 118)*yin( 123)*zin(  97)+xin( 166)*yin( 171)*zin( 145))
          eri_value(  281)=eri_value(  281)+d13bra( 16)*d12ket( 11)*(xin(  21)*yin(  28)*zin(   1)+xin(  69)*yin(  76)*zin(  49)+xin( 117)*yin( 124)*zin(  97)+xin( 165)*yin( 172)*zin( 145))
          eri_value(  282)=eri_value(  282)+d13bra( 16)*d12ket( 12)*(xin(  21)*yin(  27)*zin(   2)+xin(  69)*yin(  75)*zin(  50)+xin( 117)*yin( 123)*zin(  98)+xin( 165)*yin( 171)*zin( 146))
          eri_value(  283)=eri_value(  283)+d13bra( 16)*d12ket( 13)*(xin(  22)*yin(  25)*zin(   3)+xin(  70)*yin(  73)*zin(  51)+xin( 118)*yin( 121)*zin(  99)+xin( 166)*yin( 169)*zin( 147))
          eri_value(  284)=eri_value(  284)+d13bra( 16)*d12ket( 14)*(xin(  21)*yin(  26)*zin(   3)+xin(  69)*yin(  74)*zin(  51)+xin( 117)*yin( 122)*zin(  99)+xin( 165)*yin( 170)*zin( 147))
          eri_value(  285)=eri_value(  285)+d13bra( 16)*d12ket( 15)*(xin(  21)*yin(  25)*zin(   4)+xin(  69)*yin(  73)*zin(  52)+xin( 117)*yin( 121)*zin( 100)+xin( 165)*yin( 169)*zin( 148))
          eri_value(  286)=eri_value(  286)+d13bra( 16)*d12ket( 16)*(xin(  20)*yin(  27)*zin(   3)+xin(  68)*yin(  75)*zin(  51)+xin( 116)*yin( 123)*zin(  99)+xin( 164)*yin( 171)*zin( 147))
          eri_value(  287)=eri_value(  287)+d13bra( 16)*d12ket( 17)*(xin(  19)*yin(  28)*zin(   3)+xin(  67)*yin(  76)*zin(  51)+xin( 115)*yin( 124)*zin(  99)+xin( 163)*yin( 172)*zin( 147))
          eri_value(  288)=eri_value(  288)+d13bra( 16)*d12ket( 18)*(xin(  19)*yin(  27)*zin(   4)+xin(  67)*yin(  75)*zin(  52)+xin( 115)*yin( 123)*zin( 100)+xin( 163)*yin( 171)*zin( 148))
          eri_value(  289)=eri_value(  289)+d13bra( 17)*d12ket(  1)*(xin(  18)*yin(  31)*zin(   1)+xin(  66)*yin(  79)*zin(  49)+xin( 114)*yin( 127)*zin(  97)+xin( 162)*yin( 175)*zin( 145))
          eri_value(  290)=eri_value(  290)+d13bra( 17)*d12ket(  2)*(xin(  17)*yin(  32)*zin(   1)+xin(  65)*yin(  80)*zin(  49)+xin( 113)*yin( 128)*zin(  97)+xin( 161)*yin( 176)*zin( 145))
          eri_value(  291)=eri_value(  291)+d13bra( 17)*d12ket(  3)*(xin(  17)*yin(  31)*zin(   2)+xin(  65)*yin(  79)*zin(  50)+xin( 113)*yin( 127)*zin(  98)+xin( 161)*yin( 175)*zin( 146))
          eri_value(  292)=eri_value(  292)+d13bra( 17)*d12ket(  4)*(xin(  14)*yin(  35)*zin(   1)+xin(  62)*yin(  83)*zin(  49)+xin( 110)*yin( 131)*zin(  97)+xin( 158)*yin( 179)*zin( 145))
          eri_value(  293)=eri_value(  293)+d13bra( 17)*d12ket(  5)*(xin(  13)*yin(  36)*zin(   1)+xin(  61)*yin(  84)*zin(  49)+xin( 109)*yin( 132)*zin(  97)+xin( 157)*yin( 180)*zin( 145))
          eri_value(  294)=eri_value(  294)+d13bra( 17)*d12ket(  6)*(xin(  13)*yin(  35)*zin(   2)+xin(  61)*yin(  83)*zin(  50)+xin( 109)*yin( 131)*zin(  98)+xin( 157)*yin( 179)*zin( 146))
          eri_value(  295)=eri_value(  295)+d13bra( 17)*d12ket(  7)*(xin(  14)*yin(  31)*zin(   5)+xin(  62)*yin(  79)*zin(  53)+xin( 110)*yin( 127)*zin( 101)+xin( 158)*yin( 175)*zin( 149))
          eri_value(  296)=eri_value(  296)+d13bra( 17)*d12ket(  8)*(xin(  13)*yin(  32)*zin(   5)+xin(  61)*yin(  80)*zin(  53)+xin( 109)*yin( 128)*zin( 101)+xin( 157)*yin( 176)*zin( 149))
          eri_value(  297)=eri_value(  297)+d13bra( 17)*d12ket(  9)*(xin(  13)*yin(  31)*zin(   6)+xin(  61)*yin(  79)*zin(  54)+xin( 109)*yin( 127)*zin( 102)+xin( 157)*yin( 175)*zin( 150))
          eri_value(  298)=eri_value(  298)+d13bra( 17)*d12ket( 10)*(xin(  16)*yin(  33)*zin(   1)+xin(  64)*yin(  81)*zin(  49)+xin( 112)*yin( 129)*zin(  97)+xin( 160)*yin( 177)*zin( 145))
          eri_value(  299)=eri_value(  299)+d13bra( 17)*d12ket( 11)*(xin(  15)*yin(  34)*zin(   1)+xin(  63)*yin(  82)*zin(  49)+xin( 111)*yin( 130)*zin(  97)+xin( 159)*yin( 178)*zin( 145))
          eri_value(  300)=eri_value(  300)+d13bra( 17)*d12ket( 12)*(xin(  15)*yin(  33)*zin(   2)+xin(  63)*yin(  81)*zin(  50)+xin( 111)*yin( 129)*zin(  98)+xin( 159)*yin( 177)*zin( 146))
          eri_value(  301)=eri_value(  301)+d13bra( 17)*d12ket( 13)*(xin(  16)*yin(  31)*zin(   3)+xin(  64)*yin(  79)*zin(  51)+xin( 112)*yin( 127)*zin(  99)+xin( 160)*yin( 175)*zin( 147))
          eri_value(  302)=eri_value(  302)+d13bra( 17)*d12ket( 14)*(xin(  15)*yin(  32)*zin(   3)+xin(  63)*yin(  80)*zin(  51)+xin( 111)*yin( 128)*zin(  99)+xin( 159)*yin( 176)*zin( 147))
          eri_value(  303)=eri_value(  303)+d13bra( 17)*d12ket( 15)*(xin(  15)*yin(  31)*zin(   4)+xin(  63)*yin(  79)*zin(  52)+xin( 111)*yin( 127)*zin( 100)+xin( 159)*yin( 175)*zin( 148))
          eri_value(  304)=eri_value(  304)+d13bra( 17)*d12ket( 16)*(xin(  14)*yin(  33)*zin(   3)+xin(  62)*yin(  81)*zin(  51)+xin( 110)*yin( 129)*zin(  99)+xin( 158)*yin( 177)*zin( 147))
          eri_value(  305)=eri_value(  305)+d13bra( 17)*d12ket( 17)*(xin(  13)*yin(  34)*zin(   3)+xin(  61)*yin(  82)*zin(  51)+xin( 109)*yin( 130)*zin(  99)+xin( 157)*yin( 178)*zin( 147))
          eri_value(  306)=eri_value(  306)+d13bra( 17)*d12ket( 18)*(xin(  13)*yin(  33)*zin(   4)+xin(  61)*yin(  81)*zin(  52)+xin( 109)*yin( 129)*zin( 100)+xin( 157)*yin( 177)*zin( 148))
          eri_value(  307)=eri_value(  307)+d13bra( 18)*d12ket(  1)*(xin(  18)*yin(  25)*zin(   7)+xin(  66)*yin(  73)*zin(  55)+xin( 114)*yin( 121)*zin( 103)+xin( 162)*yin( 169)*zin( 151))
          eri_value(  308)=eri_value(  308)+d13bra( 18)*d12ket(  2)*(xin(  17)*yin(  26)*zin(   7)+xin(  65)*yin(  74)*zin(  55)+xin( 113)*yin( 122)*zin( 103)+xin( 161)*yin( 170)*zin( 151))
          eri_value(  309)=eri_value(  309)+d13bra( 18)*d12ket(  3)*(xin(  17)*yin(  25)*zin(   8)+xin(  65)*yin(  73)*zin(  56)+xin( 113)*yin( 121)*zin( 104)+xin( 161)*yin( 169)*zin( 152))
          eri_value(  310)=eri_value(  310)+d13bra( 18)*d12ket(  4)*(xin(  14)*yin(  29)*zin(   7)+xin(  62)*yin(  77)*zin(  55)+xin( 110)*yin( 125)*zin( 103)+xin( 158)*yin( 173)*zin( 151))
          eri_value(  311)=eri_value(  311)+d13bra( 18)*d12ket(  5)*(xin(  13)*yin(  30)*zin(   7)+xin(  61)*yin(  78)*zin(  55)+xin( 109)*yin( 126)*zin( 103)+xin( 157)*yin( 174)*zin( 151))
          eri_value(  312)=eri_value(  312)+d13bra( 18)*d12ket(  6)*(xin(  13)*yin(  29)*zin(   8)+xin(  61)*yin(  77)*zin(  56)+xin( 109)*yin( 125)*zin( 104)+xin( 157)*yin( 173)*zin( 152))
          eri_value(  313)=eri_value(  313)+d13bra( 18)*d12ket(  7)*(xin(  14)*yin(  25)*zin(  11)+xin(  62)*yin(  73)*zin(  59)+xin( 110)*yin( 121)*zin( 107)+xin( 158)*yin( 169)*zin( 155))
          eri_value(  314)=eri_value(  314)+d13bra( 18)*d12ket(  8)*(xin(  13)*yin(  26)*zin(  11)+xin(  61)*yin(  74)*zin(  59)+xin( 109)*yin( 122)*zin( 107)+xin( 157)*yin( 170)*zin( 155))
          eri_value(  315)=eri_value(  315)+d13bra( 18)*d12ket(  9)*(xin(  13)*yin(  25)*zin(  12)+xin(  61)*yin(  73)*zin(  60)+xin( 109)*yin( 121)*zin( 108)+xin( 157)*yin( 169)*zin( 156))
          eri_value(  316)=eri_value(  316)+d13bra( 18)*d12ket( 10)*(xin(  16)*yin(  27)*zin(   7)+xin(  64)*yin(  75)*zin(  55)+xin( 112)*yin( 123)*zin( 103)+xin( 160)*yin( 171)*zin( 151))
          eri_value(  317)=eri_value(  317)+d13bra( 18)*d12ket( 11)*(xin(  15)*yin(  28)*zin(   7)+xin(  63)*yin(  76)*zin(  55)+xin( 111)*yin( 124)*zin( 103)+xin( 159)*yin( 172)*zin( 151))
          eri_value(  318)=eri_value(  318)+d13bra( 18)*d12ket( 12)*(xin(  15)*yin(  27)*zin(   8)+xin(  63)*yin(  75)*zin(  56)+xin( 111)*yin( 123)*zin( 104)+xin( 159)*yin( 171)*zin( 152))
          eri_value(  319)=eri_value(  319)+d13bra( 18)*d12ket( 13)*(xin(  16)*yin(  25)*zin(   9)+xin(  64)*yin(  73)*zin(  57)+xin( 112)*yin( 121)*zin( 105)+xin( 160)*yin( 169)*zin( 153))
          eri_value(  320)=eri_value(  320)+d13bra( 18)*d12ket( 14)*(xin(  15)*yin(  26)*zin(   9)+xin(  63)*yin(  74)*zin(  57)+xin( 111)*yin( 122)*zin( 105)+xin( 159)*yin( 170)*zin( 153))
          eri_value(  321)=eri_value(  321)+d13bra( 18)*d12ket( 15)*(xin(  15)*yin(  25)*zin(  10)+xin(  63)*yin(  73)*zin(  58)+xin( 111)*yin( 121)*zin( 106)+xin( 159)*yin( 169)*zin( 154))
          eri_value(  322)=eri_value(  322)+d13bra( 18)*d12ket( 16)*(xin(  14)*yin(  27)*zin(   9)+xin(  62)*yin(  75)*zin(  57)+xin( 110)*yin( 123)*zin( 105)+xin( 158)*yin( 171)*zin( 153))
          eri_value(  323)=eri_value(  323)+d13bra( 18)*d12ket( 17)*(xin(  13)*yin(  28)*zin(   9)+xin(  61)*yin(  76)*zin(  57)+xin( 109)*yin( 124)*zin( 105)+xin( 157)*yin( 172)*zin( 153))
          eri_value(  324)=eri_value(  324)+d13bra( 18)*d12ket( 18)*(xin(  13)*yin(  27)*zin(  10)+xin(  61)*yin(  75)*zin(  58)+xin( 109)*yin( 123)*zin( 106)+xin( 157)*yin( 171)*zin( 154))
          eri_value(  325)=eri_value(  325)+d13bra( 19)*d12ket(  1)*(xin(  12)*yin(  25)*zin(  13)+xin(  60)*yin(  73)*zin(  61)+xin( 108)*yin( 121)*zin( 109)+xin( 156)*yin( 169)*zin( 157))
          eri_value(  326)=eri_value(  326)+d13bra( 19)*d12ket(  2)*(xin(  11)*yin(  26)*zin(  13)+xin(  59)*yin(  74)*zin(  61)+xin( 107)*yin( 122)*zin( 109)+xin( 155)*yin( 170)*zin( 157))
          eri_value(  327)=eri_value(  327)+d13bra( 19)*d12ket(  3)*(xin(  11)*yin(  25)*zin(  14)+xin(  59)*yin(  73)*zin(  62)+xin( 107)*yin( 121)*zin( 110)+xin( 155)*yin( 169)*zin( 158))
          eri_value(  328)=eri_value(  328)+d13bra( 19)*d12ket(  4)*(xin(   8)*yin(  29)*zin(  13)+xin(  56)*yin(  77)*zin(  61)+xin( 104)*yin( 125)*zin( 109)+xin( 152)*yin( 173)*zin( 157))
          eri_value(  329)=eri_value(  329)+d13bra( 19)*d12ket(  5)*(xin(   7)*yin(  30)*zin(  13)+xin(  55)*yin(  78)*zin(  61)+xin( 103)*yin( 126)*zin( 109)+xin( 151)*yin( 174)*zin( 157))
          eri_value(  330)=eri_value(  330)+d13bra( 19)*d12ket(  6)*(xin(   7)*yin(  29)*zin(  14)+xin(  55)*yin(  77)*zin(  62)+xin( 103)*yin( 125)*zin( 110)+xin( 151)*yin( 173)*zin( 158))
          eri_value(  331)=eri_value(  331)+d13bra( 19)*d12ket(  7)*(xin(   8)*yin(  25)*zin(  17)+xin(  56)*yin(  73)*zin(  65)+xin( 104)*yin( 121)*zin( 113)+xin( 152)*yin( 169)*zin( 161))
          eri_value(  332)=eri_value(  332)+d13bra( 19)*d12ket(  8)*(xin(   7)*yin(  26)*zin(  17)+xin(  55)*yin(  74)*zin(  65)+xin( 103)*yin( 122)*zin( 113)+xin( 151)*yin( 170)*zin( 161))
          eri_value(  333)=eri_value(  333)+d13bra( 19)*d12ket(  9)*(xin(   7)*yin(  25)*zin(  18)+xin(  55)*yin(  73)*zin(  66)+xin( 103)*yin( 121)*zin( 114)+xin( 151)*yin( 169)*zin( 162))
          eri_value(  334)=eri_value(  334)+d13bra( 19)*d12ket( 10)*(xin(  10)*yin(  27)*zin(  13)+xin(  58)*yin(  75)*zin(  61)+xin( 106)*yin( 123)*zin( 109)+xin( 154)*yin( 171)*zin( 157))
          eri_value(  335)=eri_value(  335)+d13bra( 19)*d12ket( 11)*(xin(   9)*yin(  28)*zin(  13)+xin(  57)*yin(  76)*zin(  61)+xin( 105)*yin( 124)*zin( 109)+xin( 153)*yin( 172)*zin( 157))
          eri_value(  336)=eri_value(  336)+d13bra( 19)*d12ket( 12)*(xin(   9)*yin(  27)*zin(  14)+xin(  57)*yin(  75)*zin(  62)+xin( 105)*yin( 123)*zin( 110)+xin( 153)*yin( 171)*zin( 158))
          eri_value(  337)=eri_value(  337)+d13bra( 19)*d12ket( 13)*(xin(  10)*yin(  25)*zin(  15)+xin(  58)*yin(  73)*zin(  63)+xin( 106)*yin( 121)*zin( 111)+xin( 154)*yin( 169)*zin( 159))
          eri_value(  338)=eri_value(  338)+d13bra( 19)*d12ket( 14)*(xin(   9)*yin(  26)*zin(  15)+xin(  57)*yin(  74)*zin(  63)+xin( 105)*yin( 122)*zin( 111)+xin( 153)*yin( 170)*zin( 159))
          eri_value(  339)=eri_value(  339)+d13bra( 19)*d12ket( 15)*(xin(   9)*yin(  25)*zin(  16)+xin(  57)*yin(  73)*zin(  64)+xin( 105)*yin( 121)*zin( 112)+xin( 153)*yin( 169)*zin( 160))
          eri_value(  340)=eri_value(  340)+d13bra( 19)*d12ket( 16)*(xin(   8)*yin(  27)*zin(  15)+xin(  56)*yin(  75)*zin(  63)+xin( 104)*yin( 123)*zin( 111)+xin( 152)*yin( 171)*zin( 159))
          eri_value(  341)=eri_value(  341)+d13bra( 19)*d12ket( 17)*(xin(   7)*yin(  28)*zin(  15)+xin(  55)*yin(  76)*zin(  63)+xin( 103)*yin( 124)*zin( 111)+xin( 151)*yin( 172)*zin( 159))
          eri_value(  342)=eri_value(  342)+d13bra( 19)*d12ket( 18)*(xin(   7)*yin(  27)*zin(  16)+xin(  55)*yin(  75)*zin(  64)+xin( 103)*yin( 123)*zin( 112)+xin( 151)*yin( 171)*zin( 160))
          eri_value(  343)=eri_value(  343)+d13bra( 20)*d12ket(  1)*(xin(   6)*yin(  31)*zin(  13)+xin(  54)*yin(  79)*zin(  61)+xin( 102)*yin( 127)*zin( 109)+xin( 150)*yin( 175)*zin( 157))
          eri_value(  344)=eri_value(  344)+d13bra( 20)*d12ket(  2)*(xin(   5)*yin(  32)*zin(  13)+xin(  53)*yin(  80)*zin(  61)+xin( 101)*yin( 128)*zin( 109)+xin( 149)*yin( 176)*zin( 157))
          eri_value(  345)=eri_value(  345)+d13bra( 20)*d12ket(  3)*(xin(   5)*yin(  31)*zin(  14)+xin(  53)*yin(  79)*zin(  62)+xin( 101)*yin( 127)*zin( 110)+xin( 149)*yin( 175)*zin( 158))
          eri_value(  346)=eri_value(  346)+d13bra( 20)*d12ket(  4)*(xin(   2)*yin(  35)*zin(  13)+xin(  50)*yin(  83)*zin(  61)+xin(  98)*yin( 131)*zin( 109)+xin( 146)*yin( 179)*zin( 157))
          eri_value(  347)=eri_value(  347)+d13bra( 20)*d12ket(  5)*(xin(   1)*yin(  36)*zin(  13)+xin(  49)*yin(  84)*zin(  61)+xin(  97)*yin( 132)*zin( 109)+xin( 145)*yin( 180)*zin( 157))
          eri_value(  348)=eri_value(  348)+d13bra( 20)*d12ket(  6)*(xin(   1)*yin(  35)*zin(  14)+xin(  49)*yin(  83)*zin(  62)+xin(  97)*yin( 131)*zin( 110)+xin( 145)*yin( 179)*zin( 158))
          eri_value(  349)=eri_value(  349)+d13bra( 20)*d12ket(  7)*(xin(   2)*yin(  31)*zin(  17)+xin(  50)*yin(  79)*zin(  65)+xin(  98)*yin( 127)*zin( 113)+xin( 146)*yin( 175)*zin( 161))
          eri_value(  350)=eri_value(  350)+d13bra( 20)*d12ket(  8)*(xin(   1)*yin(  32)*zin(  17)+xin(  49)*yin(  80)*zin(  65)+xin(  97)*yin( 128)*zin( 113)+xin( 145)*yin( 176)*zin( 161))
          eri_value(  351)=eri_value(  351)+d13bra( 20)*d12ket(  9)*(xin(   1)*yin(  31)*zin(  18)+xin(  49)*yin(  79)*zin(  66)+xin(  97)*yin( 127)*zin( 114)+xin( 145)*yin( 175)*zin( 162))
          eri_value(  352)=eri_value(  352)+d13bra( 20)*d12ket( 10)*(xin(   4)*yin(  33)*zin(  13)+xin(  52)*yin(  81)*zin(  61)+xin( 100)*yin( 129)*zin( 109)+xin( 148)*yin( 177)*zin( 157))
          eri_value(  353)=eri_value(  353)+d13bra( 20)*d12ket( 11)*(xin(   3)*yin(  34)*zin(  13)+xin(  51)*yin(  82)*zin(  61)+xin(  99)*yin( 130)*zin( 109)+xin( 147)*yin( 178)*zin( 157))
          eri_value(  354)=eri_value(  354)+d13bra( 20)*d12ket( 12)*(xin(   3)*yin(  33)*zin(  14)+xin(  51)*yin(  81)*zin(  62)+xin(  99)*yin( 129)*zin( 110)+xin( 147)*yin( 177)*zin( 158))
          eri_value(  355)=eri_value(  355)+d13bra( 20)*d12ket( 13)*(xin(   4)*yin(  31)*zin(  15)+xin(  52)*yin(  79)*zin(  63)+xin( 100)*yin( 127)*zin( 111)+xin( 148)*yin( 175)*zin( 159))
          eri_value(  356)=eri_value(  356)+d13bra( 20)*d12ket( 14)*(xin(   3)*yin(  32)*zin(  15)+xin(  51)*yin(  80)*zin(  63)+xin(  99)*yin( 128)*zin( 111)+xin( 147)*yin( 176)*zin( 159))
          eri_value(  357)=eri_value(  357)+d13bra( 20)*d12ket( 15)*(xin(   3)*yin(  31)*zin(  16)+xin(  51)*yin(  79)*zin(  64)+xin(  99)*yin( 127)*zin( 112)+xin( 147)*yin( 175)*zin( 160))
          eri_value(  358)=eri_value(  358)+d13bra( 20)*d12ket( 16)*(xin(   2)*yin(  33)*zin(  15)+xin(  50)*yin(  81)*zin(  63)+xin(  98)*yin( 129)*zin( 111)+xin( 146)*yin( 177)*zin( 159))
          eri_value(  359)=eri_value(  359)+d13bra( 20)*d12ket( 17)*(xin(   1)*yin(  34)*zin(  15)+xin(  49)*yin(  82)*zin(  63)+xin(  97)*yin( 130)*zin( 111)+xin( 145)*yin( 178)*zin( 159))
          eri_value(  360)=eri_value(  360)+d13bra( 20)*d12ket( 18)*(xin(   1)*yin(  33)*zin(  16)+xin(  49)*yin(  81)*zin(  64)+xin(  97)*yin( 129)*zin( 112)+xin( 145)*yin( 177)*zin( 160))
          eri_value(  361)=eri_value(  361)+d13bra( 21)*d12ket(  1)*(xin(   6)*yin(  25)*zin(  19)+xin(  54)*yin(  73)*zin(  67)+xin( 102)*yin( 121)*zin( 115)+xin( 150)*yin( 169)*zin( 163))
          eri_value(  362)=eri_value(  362)+d13bra( 21)*d12ket(  2)*(xin(   5)*yin(  26)*zin(  19)+xin(  53)*yin(  74)*zin(  67)+xin( 101)*yin( 122)*zin( 115)+xin( 149)*yin( 170)*zin( 163))
          eri_value(  363)=eri_value(  363)+d13bra( 21)*d12ket(  3)*(xin(   5)*yin(  25)*zin(  20)+xin(  53)*yin(  73)*zin(  68)+xin( 101)*yin( 121)*zin( 116)+xin( 149)*yin( 169)*zin( 164))
          eri_value(  364)=eri_value(  364)+d13bra( 21)*d12ket(  4)*(xin(   2)*yin(  29)*zin(  19)+xin(  50)*yin(  77)*zin(  67)+xin(  98)*yin( 125)*zin( 115)+xin( 146)*yin( 173)*zin( 163))
          eri_value(  365)=eri_value(  365)+d13bra( 21)*d12ket(  5)*(xin(   1)*yin(  30)*zin(  19)+xin(  49)*yin(  78)*zin(  67)+xin(  97)*yin( 126)*zin( 115)+xin( 145)*yin( 174)*zin( 163))
          eri_value(  366)=eri_value(  366)+d13bra( 21)*d12ket(  6)*(xin(   1)*yin(  29)*zin(  20)+xin(  49)*yin(  77)*zin(  68)+xin(  97)*yin( 125)*zin( 116)+xin( 145)*yin( 173)*zin( 164))
          eri_value(  367)=eri_value(  367)+d13bra( 21)*d12ket(  7)*(xin(   2)*yin(  25)*zin(  23)+xin(  50)*yin(  73)*zin(  71)+xin(  98)*yin( 121)*zin( 119)+xin( 146)*yin( 169)*zin( 167))
          eri_value(  368)=eri_value(  368)+d13bra( 21)*d12ket(  8)*(xin(   1)*yin(  26)*zin(  23)+xin(  49)*yin(  74)*zin(  71)+xin(  97)*yin( 122)*zin( 119)+xin( 145)*yin( 170)*zin( 167))
          eri_value(  369)=eri_value(  369)+d13bra( 21)*d12ket(  9)*(xin(   1)*yin(  25)*zin(  24)+xin(  49)*yin(  73)*zin(  72)+xin(  97)*yin( 121)*zin( 120)+xin( 145)*yin( 169)*zin( 168))
          eri_value(  370)=eri_value(  370)+d13bra( 21)*d12ket( 10)*(xin(   4)*yin(  27)*zin(  19)+xin(  52)*yin(  75)*zin(  67)+xin( 100)*yin( 123)*zin( 115)+xin( 148)*yin( 171)*zin( 163))
          eri_value(  371)=eri_value(  371)+d13bra( 21)*d12ket( 11)*(xin(   3)*yin(  28)*zin(  19)+xin(  51)*yin(  76)*zin(  67)+xin(  99)*yin( 124)*zin( 115)+xin( 147)*yin( 172)*zin( 163))
          eri_value(  372)=eri_value(  372)+d13bra( 21)*d12ket( 12)*(xin(   3)*yin(  27)*zin(  20)+xin(  51)*yin(  75)*zin(  68)+xin(  99)*yin( 123)*zin( 116)+xin( 147)*yin( 171)*zin( 164))
          eri_value(  373)=eri_value(  373)+d13bra( 21)*d12ket( 13)*(xin(   4)*yin(  25)*zin(  21)+xin(  52)*yin(  73)*zin(  69)+xin( 100)*yin( 121)*zin( 117)+xin( 148)*yin( 169)*zin( 165))
          eri_value(  374)=eri_value(  374)+d13bra( 21)*d12ket( 14)*(xin(   3)*yin(  26)*zin(  21)+xin(  51)*yin(  74)*zin(  69)+xin(  99)*yin( 122)*zin( 117)+xin( 147)*yin( 170)*zin( 165))
          eri_value(  375)=eri_value(  375)+d13bra( 21)*d12ket( 15)*(xin(   3)*yin(  25)*zin(  22)+xin(  51)*yin(  73)*zin(  70)+xin(  99)*yin( 121)*zin( 118)+xin( 147)*yin( 169)*zin( 166))
          eri_value(  376)=eri_value(  376)+d13bra( 21)*d12ket( 16)*(xin(   2)*yin(  27)*zin(  21)+xin(  50)*yin(  75)*zin(  69)+xin(  98)*yin( 123)*zin( 117)+xin( 146)*yin( 171)*zin( 165))
          eri_value(  377)=eri_value(  377)+d13bra( 21)*d12ket( 17)*(xin(   1)*yin(  28)*zin(  21)+xin(  49)*yin(  76)*zin(  69)+xin(  97)*yin( 124)*zin( 117)+xin( 145)*yin( 172)*zin( 165))
          eri_value(  378)=eri_value(  378)+d13bra( 21)*d12ket( 18)*(xin(   1)*yin(  27)*zin(  22)+xin(  49)*yin(  75)*zin(  70)+xin(  97)*yin( 123)*zin( 118)+xin( 145)*yin( 171)*zin( 166))
          eri_value(  379)=eri_value(  379)+d13bra( 22)*d12ket(  1)*(xin(  24)*yin(   1)*zin(  25)+xin(  72)*yin(  49)*zin(  73)+xin( 120)*yin(  97)*zin( 121)+xin( 168)*yin( 145)*zin( 169))
          eri_value(  380)=eri_value(  380)+d13bra( 22)*d12ket(  2)*(xin(  23)*yin(   2)*zin(  25)+xin(  71)*yin(  50)*zin(  73)+xin( 119)*yin(  98)*zin( 121)+xin( 167)*yin( 146)*zin( 169))
          eri_value(  381)=eri_value(  381)+d13bra( 22)*d12ket(  3)*(xin(  23)*yin(   1)*zin(  26)+xin(  71)*yin(  49)*zin(  74)+xin( 119)*yin(  97)*zin( 122)+xin( 167)*yin( 145)*zin( 170))
          eri_value(  382)=eri_value(  382)+d13bra( 22)*d12ket(  4)*(xin(  20)*yin(   5)*zin(  25)+xin(  68)*yin(  53)*zin(  73)+xin( 116)*yin( 101)*zin( 121)+xin( 164)*yin( 149)*zin( 169))
          eri_value(  383)=eri_value(  383)+d13bra( 22)*d12ket(  5)*(xin(  19)*yin(   6)*zin(  25)+xin(  67)*yin(  54)*zin(  73)+xin( 115)*yin( 102)*zin( 121)+xin( 163)*yin( 150)*zin( 169))
          eri_value(  384)=eri_value(  384)+d13bra( 22)*d12ket(  6)*(xin(  19)*yin(   5)*zin(  26)+xin(  67)*yin(  53)*zin(  74)+xin( 115)*yin( 101)*zin( 122)+xin( 163)*yin( 149)*zin( 170))
          eri_value(  385)=eri_value(  385)+d13bra( 22)*d12ket(  7)*(xin(  20)*yin(   1)*zin(  29)+xin(  68)*yin(  49)*zin(  77)+xin( 116)*yin(  97)*zin( 125)+xin( 164)*yin( 145)*zin( 173))
          eri_value(  386)=eri_value(  386)+d13bra( 22)*d12ket(  8)*(xin(  19)*yin(   2)*zin(  29)+xin(  67)*yin(  50)*zin(  77)+xin( 115)*yin(  98)*zin( 125)+xin( 163)*yin( 146)*zin( 173))
          eri_value(  387)=eri_value(  387)+d13bra( 22)*d12ket(  9)*(xin(  19)*yin(   1)*zin(  30)+xin(  67)*yin(  49)*zin(  78)+xin( 115)*yin(  97)*zin( 126)+xin( 163)*yin( 145)*zin( 174))
          eri_value(  388)=eri_value(  388)+d13bra( 22)*d12ket( 10)*(xin(  22)*yin(   3)*zin(  25)+xin(  70)*yin(  51)*zin(  73)+xin( 118)*yin(  99)*zin( 121)+xin( 166)*yin( 147)*zin( 169))
          eri_value(  389)=eri_value(  389)+d13bra( 22)*d12ket( 11)*(xin(  21)*yin(   4)*zin(  25)+xin(  69)*yin(  52)*zin(  73)+xin( 117)*yin( 100)*zin( 121)+xin( 165)*yin( 148)*zin( 169))
          eri_value(  390)=eri_value(  390)+d13bra( 22)*d12ket( 12)*(xin(  21)*yin(   3)*zin(  26)+xin(  69)*yin(  51)*zin(  74)+xin( 117)*yin(  99)*zin( 122)+xin( 165)*yin( 147)*zin( 170))
          eri_value(  391)=eri_value(  391)+d13bra( 22)*d12ket( 13)*(xin(  22)*yin(   1)*zin(  27)+xin(  70)*yin(  49)*zin(  75)+xin( 118)*yin(  97)*zin( 123)+xin( 166)*yin( 145)*zin( 171))
          eri_value(  392)=eri_value(  392)+d13bra( 22)*d12ket( 14)*(xin(  21)*yin(   2)*zin(  27)+xin(  69)*yin(  50)*zin(  75)+xin( 117)*yin(  98)*zin( 123)+xin( 165)*yin( 146)*zin( 171))
          eri_value(  393)=eri_value(  393)+d13bra( 22)*d12ket( 15)*(xin(  21)*yin(   1)*zin(  28)+xin(  69)*yin(  49)*zin(  76)+xin( 117)*yin(  97)*zin( 124)+xin( 165)*yin( 145)*zin( 172))
          eri_value(  394)=eri_value(  394)+d13bra( 22)*d12ket( 16)*(xin(  20)*yin(   3)*zin(  27)+xin(  68)*yin(  51)*zin(  75)+xin( 116)*yin(  99)*zin( 123)+xin( 164)*yin( 147)*zin( 171))
          eri_value(  395)=eri_value(  395)+d13bra( 22)*d12ket( 17)*(xin(  19)*yin(   4)*zin(  27)+xin(  67)*yin(  52)*zin(  75)+xin( 115)*yin( 100)*zin( 123)+xin( 163)*yin( 148)*zin( 171))
          eri_value(  396)=eri_value(  396)+d13bra( 22)*d12ket( 18)*(xin(  19)*yin(   3)*zin(  28)+xin(  67)*yin(  51)*zin(  76)+xin( 115)*yin(  99)*zin( 124)+xin( 163)*yin( 147)*zin( 172))
          eri_value(  397)=eri_value(  397)+d13bra( 23)*d12ket(  1)*(xin(  18)*yin(   7)*zin(  25)+xin(  66)*yin(  55)*zin(  73)+xin( 114)*yin( 103)*zin( 121)+xin( 162)*yin( 151)*zin( 169))
          eri_value(  398)=eri_value(  398)+d13bra( 23)*d12ket(  2)*(xin(  17)*yin(   8)*zin(  25)+xin(  65)*yin(  56)*zin(  73)+xin( 113)*yin( 104)*zin( 121)+xin( 161)*yin( 152)*zin( 169))
          eri_value(  399)=eri_value(  399)+d13bra( 23)*d12ket(  3)*(xin(  17)*yin(   7)*zin(  26)+xin(  65)*yin(  55)*zin(  74)+xin( 113)*yin( 103)*zin( 122)+xin( 161)*yin( 151)*zin( 170))
          eri_value(  400)=eri_value(  400)+d13bra( 23)*d12ket(  4)*(xin(  14)*yin(  11)*zin(  25)+xin(  62)*yin(  59)*zin(  73)+xin( 110)*yin( 107)*zin( 121)+xin( 158)*yin( 155)*zin( 169))
          eri_value(  401)=eri_value(  401)+d13bra( 23)*d12ket(  5)*(xin(  13)*yin(  12)*zin(  25)+xin(  61)*yin(  60)*zin(  73)+xin( 109)*yin( 108)*zin( 121)+xin( 157)*yin( 156)*zin( 169))
          eri_value(  402)=eri_value(  402)+d13bra( 23)*d12ket(  6)*(xin(  13)*yin(  11)*zin(  26)+xin(  61)*yin(  59)*zin(  74)+xin( 109)*yin( 107)*zin( 122)+xin( 157)*yin( 155)*zin( 170))
          eri_value(  403)=eri_value(  403)+d13bra( 23)*d12ket(  7)*(xin(  14)*yin(   7)*zin(  29)+xin(  62)*yin(  55)*zin(  77)+xin( 110)*yin( 103)*zin( 125)+xin( 158)*yin( 151)*zin( 173))
          eri_value(  404)=eri_value(  404)+d13bra( 23)*d12ket(  8)*(xin(  13)*yin(   8)*zin(  29)+xin(  61)*yin(  56)*zin(  77)+xin( 109)*yin( 104)*zin( 125)+xin( 157)*yin( 152)*zin( 173))
          eri_value(  405)=eri_value(  405)+d13bra( 23)*d12ket(  9)*(xin(  13)*yin(   7)*zin(  30)+xin(  61)*yin(  55)*zin(  78)+xin( 109)*yin( 103)*zin( 126)+xin( 157)*yin( 151)*zin( 174))
          eri_value(  406)=eri_value(  406)+d13bra( 23)*d12ket( 10)*(xin(  16)*yin(   9)*zin(  25)+xin(  64)*yin(  57)*zin(  73)+xin( 112)*yin( 105)*zin( 121)+xin( 160)*yin( 153)*zin( 169))
          eri_value(  407)=eri_value(  407)+d13bra( 23)*d12ket( 11)*(xin(  15)*yin(  10)*zin(  25)+xin(  63)*yin(  58)*zin(  73)+xin( 111)*yin( 106)*zin( 121)+xin( 159)*yin( 154)*zin( 169))
          eri_value(  408)=eri_value(  408)+d13bra( 23)*d12ket( 12)*(xin(  15)*yin(   9)*zin(  26)+xin(  63)*yin(  57)*zin(  74)+xin( 111)*yin( 105)*zin( 122)+xin( 159)*yin( 153)*zin( 170))
          eri_value(  409)=eri_value(  409)+d13bra( 23)*d12ket( 13)*(xin(  16)*yin(   7)*zin(  27)+xin(  64)*yin(  55)*zin(  75)+xin( 112)*yin( 103)*zin( 123)+xin( 160)*yin( 151)*zin( 171))
          eri_value(  410)=eri_value(  410)+d13bra( 23)*d12ket( 14)*(xin(  15)*yin(   8)*zin(  27)+xin(  63)*yin(  56)*zin(  75)+xin( 111)*yin( 104)*zin( 123)+xin( 159)*yin( 152)*zin( 171))
          eri_value(  411)=eri_value(  411)+d13bra( 23)*d12ket( 15)*(xin(  15)*yin(   7)*zin(  28)+xin(  63)*yin(  55)*zin(  76)+xin( 111)*yin( 103)*zin( 124)+xin( 159)*yin( 151)*zin( 172))
          eri_value(  412)=eri_value(  412)+d13bra( 23)*d12ket( 16)*(xin(  14)*yin(   9)*zin(  27)+xin(  62)*yin(  57)*zin(  75)+xin( 110)*yin( 105)*zin( 123)+xin( 158)*yin( 153)*zin( 171))
          eri_value(  413)=eri_value(  413)+d13bra( 23)*d12ket( 17)*(xin(  13)*yin(  10)*zin(  27)+xin(  61)*yin(  58)*zin(  75)+xin( 109)*yin( 106)*zin( 123)+xin( 157)*yin( 154)*zin( 171))
          eri_value(  414)=eri_value(  414)+d13bra( 23)*d12ket( 18)*(xin(  13)*yin(   9)*zin(  28)+xin(  61)*yin(  57)*zin(  76)+xin( 109)*yin( 105)*zin( 124)+xin( 157)*yin( 153)*zin( 172))
          eri_value(  415)=eri_value(  415)+d13bra( 24)*d12ket(  1)*(xin(  18)*yin(   1)*zin(  31)+xin(  66)*yin(  49)*zin(  79)+xin( 114)*yin(  97)*zin( 127)+xin( 162)*yin( 145)*zin( 175))
          eri_value(  416)=eri_value(  416)+d13bra( 24)*d12ket(  2)*(xin(  17)*yin(   2)*zin(  31)+xin(  65)*yin(  50)*zin(  79)+xin( 113)*yin(  98)*zin( 127)+xin( 161)*yin( 146)*zin( 175))
          eri_value(  417)=eri_value(  417)+d13bra( 24)*d12ket(  3)*(xin(  17)*yin(   1)*zin(  32)+xin(  65)*yin(  49)*zin(  80)+xin( 113)*yin(  97)*zin( 128)+xin( 161)*yin( 145)*zin( 176))
          eri_value(  418)=eri_value(  418)+d13bra( 24)*d12ket(  4)*(xin(  14)*yin(   5)*zin(  31)+xin(  62)*yin(  53)*zin(  79)+xin( 110)*yin( 101)*zin( 127)+xin( 158)*yin( 149)*zin( 175))
          eri_value(  419)=eri_value(  419)+d13bra( 24)*d12ket(  5)*(xin(  13)*yin(   6)*zin(  31)+xin(  61)*yin(  54)*zin(  79)+xin( 109)*yin( 102)*zin( 127)+xin( 157)*yin( 150)*zin( 175))
          eri_value(  420)=eri_value(  420)+d13bra( 24)*d12ket(  6)*(xin(  13)*yin(   5)*zin(  32)+xin(  61)*yin(  53)*zin(  80)+xin( 109)*yin( 101)*zin( 128)+xin( 157)*yin( 149)*zin( 176))
          eri_value(  421)=eri_value(  421)+d13bra( 24)*d12ket(  7)*(xin(  14)*yin(   1)*zin(  35)+xin(  62)*yin(  49)*zin(  83)+xin( 110)*yin(  97)*zin( 131)+xin( 158)*yin( 145)*zin( 179))
          eri_value(  422)=eri_value(  422)+d13bra( 24)*d12ket(  8)*(xin(  13)*yin(   2)*zin(  35)+xin(  61)*yin(  50)*zin(  83)+xin( 109)*yin(  98)*zin( 131)+xin( 157)*yin( 146)*zin( 179))
          eri_value(  423)=eri_value(  423)+d13bra( 24)*d12ket(  9)*(xin(  13)*yin(   1)*zin(  36)+xin(  61)*yin(  49)*zin(  84)+xin( 109)*yin(  97)*zin( 132)+xin( 157)*yin( 145)*zin( 180))
          eri_value(  424)=eri_value(  424)+d13bra( 24)*d12ket( 10)*(xin(  16)*yin(   3)*zin(  31)+xin(  64)*yin(  51)*zin(  79)+xin( 112)*yin(  99)*zin( 127)+xin( 160)*yin( 147)*zin( 175))
          eri_value(  425)=eri_value(  425)+d13bra( 24)*d12ket( 11)*(xin(  15)*yin(   4)*zin(  31)+xin(  63)*yin(  52)*zin(  79)+xin( 111)*yin( 100)*zin( 127)+xin( 159)*yin( 148)*zin( 175))
          eri_value(  426)=eri_value(  426)+d13bra( 24)*d12ket( 12)*(xin(  15)*yin(   3)*zin(  32)+xin(  63)*yin(  51)*zin(  80)+xin( 111)*yin(  99)*zin( 128)+xin( 159)*yin( 147)*zin( 176))
          eri_value(  427)=eri_value(  427)+d13bra( 24)*d12ket( 13)*(xin(  16)*yin(   1)*zin(  33)+xin(  64)*yin(  49)*zin(  81)+xin( 112)*yin(  97)*zin( 129)+xin( 160)*yin( 145)*zin( 177))
          eri_value(  428)=eri_value(  428)+d13bra( 24)*d12ket( 14)*(xin(  15)*yin(   2)*zin(  33)+xin(  63)*yin(  50)*zin(  81)+xin( 111)*yin(  98)*zin( 129)+xin( 159)*yin( 146)*zin( 177))
          eri_value(  429)=eri_value(  429)+d13bra( 24)*d12ket( 15)*(xin(  15)*yin(   1)*zin(  34)+xin(  63)*yin(  49)*zin(  82)+xin( 111)*yin(  97)*zin( 130)+xin( 159)*yin( 145)*zin( 178))
          eri_value(  430)=eri_value(  430)+d13bra( 24)*d12ket( 16)*(xin(  14)*yin(   3)*zin(  33)+xin(  62)*yin(  51)*zin(  81)+xin( 110)*yin(  99)*zin( 129)+xin( 158)*yin( 147)*zin( 177))
          eri_value(  431)=eri_value(  431)+d13bra( 24)*d12ket( 17)*(xin(  13)*yin(   4)*zin(  33)+xin(  61)*yin(  52)*zin(  81)+xin( 109)*yin( 100)*zin( 129)+xin( 157)*yin( 148)*zin( 177))
          eri_value(  432)=eri_value(  432)+d13bra( 24)*d12ket( 18)*(xin(  13)*yin(   3)*zin(  34)+xin(  61)*yin(  51)*zin(  82)+xin( 109)*yin(  99)*zin( 130)+xin( 157)*yin( 147)*zin( 178))
          eri_value(  433)=eri_value(  433)+d13bra( 25)*d12ket(  1)*(xin(  12)*yin(  13)*zin(  25)+xin(  60)*yin(  61)*zin(  73)+xin( 108)*yin( 109)*zin( 121)+xin( 156)*yin( 157)*zin( 169))
          eri_value(  434)=eri_value(  434)+d13bra( 25)*d12ket(  2)*(xin(  11)*yin(  14)*zin(  25)+xin(  59)*yin(  62)*zin(  73)+xin( 107)*yin( 110)*zin( 121)+xin( 155)*yin( 158)*zin( 169))
          eri_value(  435)=eri_value(  435)+d13bra( 25)*d12ket(  3)*(xin(  11)*yin(  13)*zin(  26)+xin(  59)*yin(  61)*zin(  74)+xin( 107)*yin( 109)*zin( 122)+xin( 155)*yin( 157)*zin( 170))
          eri_value(  436)=eri_value(  436)+d13bra( 25)*d12ket(  4)*(xin(   8)*yin(  17)*zin(  25)+xin(  56)*yin(  65)*zin(  73)+xin( 104)*yin( 113)*zin( 121)+xin( 152)*yin( 161)*zin( 169))
          eri_value(  437)=eri_value(  437)+d13bra( 25)*d12ket(  5)*(xin(   7)*yin(  18)*zin(  25)+xin(  55)*yin(  66)*zin(  73)+xin( 103)*yin( 114)*zin( 121)+xin( 151)*yin( 162)*zin( 169))
          eri_value(  438)=eri_value(  438)+d13bra( 25)*d12ket(  6)*(xin(   7)*yin(  17)*zin(  26)+xin(  55)*yin(  65)*zin(  74)+xin( 103)*yin( 113)*zin( 122)+xin( 151)*yin( 161)*zin( 170))
          eri_value(  439)=eri_value(  439)+d13bra( 25)*d12ket(  7)*(xin(   8)*yin(  13)*zin(  29)+xin(  56)*yin(  61)*zin(  77)+xin( 104)*yin( 109)*zin( 125)+xin( 152)*yin( 157)*zin( 173))
          eri_value(  440)=eri_value(  440)+d13bra( 25)*d12ket(  8)*(xin(   7)*yin(  14)*zin(  29)+xin(  55)*yin(  62)*zin(  77)+xin( 103)*yin( 110)*zin( 125)+xin( 151)*yin( 158)*zin( 173))
          eri_value(  441)=eri_value(  441)+d13bra( 25)*d12ket(  9)*(xin(   7)*yin(  13)*zin(  30)+xin(  55)*yin(  61)*zin(  78)+xin( 103)*yin( 109)*zin( 126)+xin( 151)*yin( 157)*zin( 174))
          eri_value(  442)=eri_value(  442)+d13bra( 25)*d12ket( 10)*(xin(  10)*yin(  15)*zin(  25)+xin(  58)*yin(  63)*zin(  73)+xin( 106)*yin( 111)*zin( 121)+xin( 154)*yin( 159)*zin( 169))
          eri_value(  443)=eri_value(  443)+d13bra( 25)*d12ket( 11)*(xin(   9)*yin(  16)*zin(  25)+xin(  57)*yin(  64)*zin(  73)+xin( 105)*yin( 112)*zin( 121)+xin( 153)*yin( 160)*zin( 169))
          eri_value(  444)=eri_value(  444)+d13bra( 25)*d12ket( 12)*(xin(   9)*yin(  15)*zin(  26)+xin(  57)*yin(  63)*zin(  74)+xin( 105)*yin( 111)*zin( 122)+xin( 153)*yin( 159)*zin( 170))
          eri_value(  445)=eri_value(  445)+d13bra( 25)*d12ket( 13)*(xin(  10)*yin(  13)*zin(  27)+xin(  58)*yin(  61)*zin(  75)+xin( 106)*yin( 109)*zin( 123)+xin( 154)*yin( 157)*zin( 171))
          eri_value(  446)=eri_value(  446)+d13bra( 25)*d12ket( 14)*(xin(   9)*yin(  14)*zin(  27)+xin(  57)*yin(  62)*zin(  75)+xin( 105)*yin( 110)*zin( 123)+xin( 153)*yin( 158)*zin( 171))
          eri_value(  447)=eri_value(  447)+d13bra( 25)*d12ket( 15)*(xin(   9)*yin(  13)*zin(  28)+xin(  57)*yin(  61)*zin(  76)+xin( 105)*yin( 109)*zin( 124)+xin( 153)*yin( 157)*zin( 172))
          eri_value(  448)=eri_value(  448)+d13bra( 25)*d12ket( 16)*(xin(   8)*yin(  15)*zin(  27)+xin(  56)*yin(  63)*zin(  75)+xin( 104)*yin( 111)*zin( 123)+xin( 152)*yin( 159)*zin( 171))
          eri_value(  449)=eri_value(  449)+d13bra( 25)*d12ket( 17)*(xin(   7)*yin(  16)*zin(  27)+xin(  55)*yin(  64)*zin(  75)+xin( 103)*yin( 112)*zin( 123)+xin( 151)*yin( 160)*zin( 171))
          eri_value(  450)=eri_value(  450)+d13bra( 25)*d12ket( 18)*(xin(   7)*yin(  15)*zin(  28)+xin(  55)*yin(  63)*zin(  76)+xin( 103)*yin( 111)*zin( 124)+xin( 151)*yin( 159)*zin( 172))
          eri_value(  451)=eri_value(  451)+d13bra( 26)*d12ket(  1)*(xin(   6)*yin(  19)*zin(  25)+xin(  54)*yin(  67)*zin(  73)+xin( 102)*yin( 115)*zin( 121)+xin( 150)*yin( 163)*zin( 169))
          eri_value(  452)=eri_value(  452)+d13bra( 26)*d12ket(  2)*(xin(   5)*yin(  20)*zin(  25)+xin(  53)*yin(  68)*zin(  73)+xin( 101)*yin( 116)*zin( 121)+xin( 149)*yin( 164)*zin( 169))
          eri_value(  453)=eri_value(  453)+d13bra( 26)*d12ket(  3)*(xin(   5)*yin(  19)*zin(  26)+xin(  53)*yin(  67)*zin(  74)+xin( 101)*yin( 115)*zin( 122)+xin( 149)*yin( 163)*zin( 170))
          eri_value(  454)=eri_value(  454)+d13bra( 26)*d12ket(  4)*(xin(   2)*yin(  23)*zin(  25)+xin(  50)*yin(  71)*zin(  73)+xin(  98)*yin( 119)*zin( 121)+xin( 146)*yin( 167)*zin( 169))
          eri_value(  455)=eri_value(  455)+d13bra( 26)*d12ket(  5)*(xin(   1)*yin(  24)*zin(  25)+xin(  49)*yin(  72)*zin(  73)+xin(  97)*yin( 120)*zin( 121)+xin( 145)*yin( 168)*zin( 169))
          eri_value(  456)=eri_value(  456)+d13bra( 26)*d12ket(  6)*(xin(   1)*yin(  23)*zin(  26)+xin(  49)*yin(  71)*zin(  74)+xin(  97)*yin( 119)*zin( 122)+xin( 145)*yin( 167)*zin( 170))
          eri_value(  457)=eri_value(  457)+d13bra( 26)*d12ket(  7)*(xin(   2)*yin(  19)*zin(  29)+xin(  50)*yin(  67)*zin(  77)+xin(  98)*yin( 115)*zin( 125)+xin( 146)*yin( 163)*zin( 173))
          eri_value(  458)=eri_value(  458)+d13bra( 26)*d12ket(  8)*(xin(   1)*yin(  20)*zin(  29)+xin(  49)*yin(  68)*zin(  77)+xin(  97)*yin( 116)*zin( 125)+xin( 145)*yin( 164)*zin( 173))
          eri_value(  459)=eri_value(  459)+d13bra( 26)*d12ket(  9)*(xin(   1)*yin(  19)*zin(  30)+xin(  49)*yin(  67)*zin(  78)+xin(  97)*yin( 115)*zin( 126)+xin( 145)*yin( 163)*zin( 174))
          eri_value(  460)=eri_value(  460)+d13bra( 26)*d12ket( 10)*(xin(   4)*yin(  21)*zin(  25)+xin(  52)*yin(  69)*zin(  73)+xin( 100)*yin( 117)*zin( 121)+xin( 148)*yin( 165)*zin( 169))
          eri_value(  461)=eri_value(  461)+d13bra( 26)*d12ket( 11)*(xin(   3)*yin(  22)*zin(  25)+xin(  51)*yin(  70)*zin(  73)+xin(  99)*yin( 118)*zin( 121)+xin( 147)*yin( 166)*zin( 169))
          eri_value(  462)=eri_value(  462)+d13bra( 26)*d12ket( 12)*(xin(   3)*yin(  21)*zin(  26)+xin(  51)*yin(  69)*zin(  74)+xin(  99)*yin( 117)*zin( 122)+xin( 147)*yin( 165)*zin( 170))
          eri_value(  463)=eri_value(  463)+d13bra( 26)*d12ket( 13)*(xin(   4)*yin(  19)*zin(  27)+xin(  52)*yin(  67)*zin(  75)+xin( 100)*yin( 115)*zin( 123)+xin( 148)*yin( 163)*zin( 171))
          eri_value(  464)=eri_value(  464)+d13bra( 26)*d12ket( 14)*(xin(   3)*yin(  20)*zin(  27)+xin(  51)*yin(  68)*zin(  75)+xin(  99)*yin( 116)*zin( 123)+xin( 147)*yin( 164)*zin( 171))
          eri_value(  465)=eri_value(  465)+d13bra( 26)*d12ket( 15)*(xin(   3)*yin(  19)*zin(  28)+xin(  51)*yin(  67)*zin(  76)+xin(  99)*yin( 115)*zin( 124)+xin( 147)*yin( 163)*zin( 172))
          eri_value(  466)=eri_value(  466)+d13bra( 26)*d12ket( 16)*(xin(   2)*yin(  21)*zin(  27)+xin(  50)*yin(  69)*zin(  75)+xin(  98)*yin( 117)*zin( 123)+xin( 146)*yin( 165)*zin( 171))
          eri_value(  467)=eri_value(  467)+d13bra( 26)*d12ket( 17)*(xin(   1)*yin(  22)*zin(  27)+xin(  49)*yin(  70)*zin(  75)+xin(  97)*yin( 118)*zin( 123)+xin( 145)*yin( 166)*zin( 171))
          eri_value(  468)=eri_value(  468)+d13bra( 26)*d12ket( 18)*(xin(   1)*yin(  21)*zin(  28)+xin(  49)*yin(  69)*zin(  76)+xin(  97)*yin( 117)*zin( 124)+xin( 145)*yin( 165)*zin( 172))
          eri_value(  469)=eri_value(  469)+d13bra( 27)*d12ket(  1)*(xin(   6)*yin(  13)*zin(  31)+xin(  54)*yin(  61)*zin(  79)+xin( 102)*yin( 109)*zin( 127)+xin( 150)*yin( 157)*zin( 175))
          eri_value(  470)=eri_value(  470)+d13bra( 27)*d12ket(  2)*(xin(   5)*yin(  14)*zin(  31)+xin(  53)*yin(  62)*zin(  79)+xin( 101)*yin( 110)*zin( 127)+xin( 149)*yin( 158)*zin( 175))
          eri_value(  471)=eri_value(  471)+d13bra( 27)*d12ket(  3)*(xin(   5)*yin(  13)*zin(  32)+xin(  53)*yin(  61)*zin(  80)+xin( 101)*yin( 109)*zin( 128)+xin( 149)*yin( 157)*zin( 176))
          eri_value(  472)=eri_value(  472)+d13bra( 27)*d12ket(  4)*(xin(   2)*yin(  17)*zin(  31)+xin(  50)*yin(  65)*zin(  79)+xin(  98)*yin( 113)*zin( 127)+xin( 146)*yin( 161)*zin( 175))
          eri_value(  473)=eri_value(  473)+d13bra( 27)*d12ket(  5)*(xin(   1)*yin(  18)*zin(  31)+xin(  49)*yin(  66)*zin(  79)+xin(  97)*yin( 114)*zin( 127)+xin( 145)*yin( 162)*zin( 175))
          eri_value(  474)=eri_value(  474)+d13bra( 27)*d12ket(  6)*(xin(   1)*yin(  17)*zin(  32)+xin(  49)*yin(  65)*zin(  80)+xin(  97)*yin( 113)*zin( 128)+xin( 145)*yin( 161)*zin( 176))
          eri_value(  475)=eri_value(  475)+d13bra( 27)*d12ket(  7)*(xin(   2)*yin(  13)*zin(  35)+xin(  50)*yin(  61)*zin(  83)+xin(  98)*yin( 109)*zin( 131)+xin( 146)*yin( 157)*zin( 179))
          eri_value(  476)=eri_value(  476)+d13bra( 27)*d12ket(  8)*(xin(   1)*yin(  14)*zin(  35)+xin(  49)*yin(  62)*zin(  83)+xin(  97)*yin( 110)*zin( 131)+xin( 145)*yin( 158)*zin( 179))
          eri_value(  477)=eri_value(  477)+d13bra( 27)*d12ket(  9)*(xin(   1)*yin(  13)*zin(  36)+xin(  49)*yin(  61)*zin(  84)+xin(  97)*yin( 109)*zin( 132)+xin( 145)*yin( 157)*zin( 180))
          eri_value(  478)=eri_value(  478)+d13bra( 27)*d12ket( 10)*(xin(   4)*yin(  15)*zin(  31)+xin(  52)*yin(  63)*zin(  79)+xin( 100)*yin( 111)*zin( 127)+xin( 148)*yin( 159)*zin( 175))
          eri_value(  479)=eri_value(  479)+d13bra( 27)*d12ket( 11)*(xin(   3)*yin(  16)*zin(  31)+xin(  51)*yin(  64)*zin(  79)+xin(  99)*yin( 112)*zin( 127)+xin( 147)*yin( 160)*zin( 175))
          eri_value(  480)=eri_value(  480)+d13bra( 27)*d12ket( 12)*(xin(   3)*yin(  15)*zin(  32)+xin(  51)*yin(  63)*zin(  80)+xin(  99)*yin( 111)*zin( 128)+xin( 147)*yin( 159)*zin( 176))
          eri_value(  481)=eri_value(  481)+d13bra( 27)*d12ket( 13)*(xin(   4)*yin(  13)*zin(  33)+xin(  52)*yin(  61)*zin(  81)+xin( 100)*yin( 109)*zin( 129)+xin( 148)*yin( 157)*zin( 177))
          eri_value(  482)=eri_value(  482)+d13bra( 27)*d12ket( 14)*(xin(   3)*yin(  14)*zin(  33)+xin(  51)*yin(  62)*zin(  81)+xin(  99)*yin( 110)*zin( 129)+xin( 147)*yin( 158)*zin( 177))
          eri_value(  483)=eri_value(  483)+d13bra( 27)*d12ket( 15)*(xin(   3)*yin(  13)*zin(  34)+xin(  51)*yin(  61)*zin(  82)+xin(  99)*yin( 109)*zin( 130)+xin( 147)*yin( 157)*zin( 178))
          eri_value(  484)=eri_value(  484)+d13bra( 27)*d12ket( 16)*(xin(   2)*yin(  15)*zin(  33)+xin(  50)*yin(  63)*zin(  81)+xin(  98)*yin( 111)*zin( 129)+xin( 146)*yin( 159)*zin( 177))
          eri_value(  485)=eri_value(  485)+d13bra( 27)*d12ket( 17)*(xin(   1)*yin(  16)*zin(  33)+xin(  49)*yin(  64)*zin(  81)+xin(  97)*yin( 112)*zin( 129)+xin( 145)*yin( 160)*zin( 177))
          eri_value(  486)=eri_value(  486)+d13bra( 27)*d12ket( 18)*(xin(   1)*yin(  15)*zin(  34)+xin(  49)*yin(  63)*zin(  82)+xin(  97)*yin( 111)*zin( 130)+xin( 145)*yin( 159)*zin( 178))
          eri_value(  487)=eri_value(  487)+d13bra( 28)*d12ket(  1)*(xin(  24)*yin(  13)*zin(  13)+xin(  72)*yin(  61)*zin(  61)+xin( 120)*yin( 109)*zin( 109)+xin( 168)*yin( 157)*zin( 157))
          eri_value(  488)=eri_value(  488)+d13bra( 28)*d12ket(  2)*(xin(  23)*yin(  14)*zin(  13)+xin(  71)*yin(  62)*zin(  61)+xin( 119)*yin( 110)*zin( 109)+xin( 167)*yin( 158)*zin( 157))
          eri_value(  489)=eri_value(  489)+d13bra( 28)*d12ket(  3)*(xin(  23)*yin(  13)*zin(  14)+xin(  71)*yin(  61)*zin(  62)+xin( 119)*yin( 109)*zin( 110)+xin( 167)*yin( 157)*zin( 158))
          eri_value(  490)=eri_value(  490)+d13bra( 28)*d12ket(  4)*(xin(  20)*yin(  17)*zin(  13)+xin(  68)*yin(  65)*zin(  61)+xin( 116)*yin( 113)*zin( 109)+xin( 164)*yin( 161)*zin( 157))
          eri_value(  491)=eri_value(  491)+d13bra( 28)*d12ket(  5)*(xin(  19)*yin(  18)*zin(  13)+xin(  67)*yin(  66)*zin(  61)+xin( 115)*yin( 114)*zin( 109)+xin( 163)*yin( 162)*zin( 157))
          eri_value(  492)=eri_value(  492)+d13bra( 28)*d12ket(  6)*(xin(  19)*yin(  17)*zin(  14)+xin(  67)*yin(  65)*zin(  62)+xin( 115)*yin( 113)*zin( 110)+xin( 163)*yin( 161)*zin( 158))
          eri_value(  493)=eri_value(  493)+d13bra( 28)*d12ket(  7)*(xin(  20)*yin(  13)*zin(  17)+xin(  68)*yin(  61)*zin(  65)+xin( 116)*yin( 109)*zin( 113)+xin( 164)*yin( 157)*zin( 161))
          eri_value(  494)=eri_value(  494)+d13bra( 28)*d12ket(  8)*(xin(  19)*yin(  14)*zin(  17)+xin(  67)*yin(  62)*zin(  65)+xin( 115)*yin( 110)*zin( 113)+xin( 163)*yin( 158)*zin( 161))
          eri_value(  495)=eri_value(  495)+d13bra( 28)*d12ket(  9)*(xin(  19)*yin(  13)*zin(  18)+xin(  67)*yin(  61)*zin(  66)+xin( 115)*yin( 109)*zin( 114)+xin( 163)*yin( 157)*zin( 162))
          eri_value(  496)=eri_value(  496)+d13bra( 28)*d12ket( 10)*(xin(  22)*yin(  15)*zin(  13)+xin(  70)*yin(  63)*zin(  61)+xin( 118)*yin( 111)*zin( 109)+xin( 166)*yin( 159)*zin( 157))
          eri_value(  497)=eri_value(  497)+d13bra( 28)*d12ket( 11)*(xin(  21)*yin(  16)*zin(  13)+xin(  69)*yin(  64)*zin(  61)+xin( 117)*yin( 112)*zin( 109)+xin( 165)*yin( 160)*zin( 157))
          eri_value(  498)=eri_value(  498)+d13bra( 28)*d12ket( 12)*(xin(  21)*yin(  15)*zin(  14)+xin(  69)*yin(  63)*zin(  62)+xin( 117)*yin( 111)*zin( 110)+xin( 165)*yin( 159)*zin( 158))
          eri_value(  499)=eri_value(  499)+d13bra( 28)*d12ket( 13)*(xin(  22)*yin(  13)*zin(  15)+xin(  70)*yin(  61)*zin(  63)+xin( 118)*yin( 109)*zin( 111)+xin( 166)*yin( 157)*zin( 159))
          eri_value(  500)=eri_value(  500)+d13bra( 28)*d12ket( 14)*(xin(  21)*yin(  14)*zin(  15)+xin(  69)*yin(  62)*zin(  63)+xin( 117)*yin( 110)*zin( 111)+xin( 165)*yin( 158)*zin( 159))
          eri_value(  501)=eri_value(  501)+d13bra( 28)*d12ket( 15)*(xin(  21)*yin(  13)*zin(  16)+xin(  69)*yin(  61)*zin(  64)+xin( 117)*yin( 109)*zin( 112)+xin( 165)*yin( 157)*zin( 160))
          eri_value(  502)=eri_value(  502)+d13bra( 28)*d12ket( 16)*(xin(  20)*yin(  15)*zin(  15)+xin(  68)*yin(  63)*zin(  63)+xin( 116)*yin( 111)*zin( 111)+xin( 164)*yin( 159)*zin( 159))
          eri_value(  503)=eri_value(  503)+d13bra( 28)*d12ket( 17)*(xin(  19)*yin(  16)*zin(  15)+xin(  67)*yin(  64)*zin(  63)+xin( 115)*yin( 112)*zin( 111)+xin( 163)*yin( 160)*zin( 159))
          eri_value(  504)=eri_value(  504)+d13bra( 28)*d12ket( 18)*(xin(  19)*yin(  15)*zin(  16)+xin(  67)*yin(  63)*zin(  64)+xin( 115)*yin( 111)*zin( 112)+xin( 163)*yin( 159)*zin( 160))
          eri_value(  505)=eri_value(  505)+d13bra( 29)*d12ket(  1)*(xin(  18)*yin(  19)*zin(  13)+xin(  66)*yin(  67)*zin(  61)+xin( 114)*yin( 115)*zin( 109)+xin( 162)*yin( 163)*zin( 157))
          eri_value(  506)=eri_value(  506)+d13bra( 29)*d12ket(  2)*(xin(  17)*yin(  20)*zin(  13)+xin(  65)*yin(  68)*zin(  61)+xin( 113)*yin( 116)*zin( 109)+xin( 161)*yin( 164)*zin( 157))
          eri_value(  507)=eri_value(  507)+d13bra( 29)*d12ket(  3)*(xin(  17)*yin(  19)*zin(  14)+xin(  65)*yin(  67)*zin(  62)+xin( 113)*yin( 115)*zin( 110)+xin( 161)*yin( 163)*zin( 158))
          eri_value(  508)=eri_value(  508)+d13bra( 29)*d12ket(  4)*(xin(  14)*yin(  23)*zin(  13)+xin(  62)*yin(  71)*zin(  61)+xin( 110)*yin( 119)*zin( 109)+xin( 158)*yin( 167)*zin( 157))
          eri_value(  509)=eri_value(  509)+d13bra( 29)*d12ket(  5)*(xin(  13)*yin(  24)*zin(  13)+xin(  61)*yin(  72)*zin(  61)+xin( 109)*yin( 120)*zin( 109)+xin( 157)*yin( 168)*zin( 157))
          eri_value(  510)=eri_value(  510)+d13bra( 29)*d12ket(  6)*(xin(  13)*yin(  23)*zin(  14)+xin(  61)*yin(  71)*zin(  62)+xin( 109)*yin( 119)*zin( 110)+xin( 157)*yin( 167)*zin( 158))
          eri_value(  511)=eri_value(  511)+d13bra( 29)*d12ket(  7)*(xin(  14)*yin(  19)*zin(  17)+xin(  62)*yin(  67)*zin(  65)+xin( 110)*yin( 115)*zin( 113)+xin( 158)*yin( 163)*zin( 161))
          eri_value(  512)=eri_value(  512)+d13bra( 29)*d12ket(  8)*(xin(  13)*yin(  20)*zin(  17)+xin(  61)*yin(  68)*zin(  65)+xin( 109)*yin( 116)*zin( 113)+xin( 157)*yin( 164)*zin( 161))
          eri_value(  513)=eri_value(  513)+d13bra( 29)*d12ket(  9)*(xin(  13)*yin(  19)*zin(  18)+xin(  61)*yin(  67)*zin(  66)+xin( 109)*yin( 115)*zin( 114)+xin( 157)*yin( 163)*zin( 162))
          eri_value(  514)=eri_value(  514)+d13bra( 29)*d12ket( 10)*(xin(  16)*yin(  21)*zin(  13)+xin(  64)*yin(  69)*zin(  61)+xin( 112)*yin( 117)*zin( 109)+xin( 160)*yin( 165)*zin( 157))
          eri_value(  515)=eri_value(  515)+d13bra( 29)*d12ket( 11)*(xin(  15)*yin(  22)*zin(  13)+xin(  63)*yin(  70)*zin(  61)+xin( 111)*yin( 118)*zin( 109)+xin( 159)*yin( 166)*zin( 157))
          eri_value(  516)=eri_value(  516)+d13bra( 29)*d12ket( 12)*(xin(  15)*yin(  21)*zin(  14)+xin(  63)*yin(  69)*zin(  62)+xin( 111)*yin( 117)*zin( 110)+xin( 159)*yin( 165)*zin( 158))
          eri_value(  517)=eri_value(  517)+d13bra( 29)*d12ket( 13)*(xin(  16)*yin(  19)*zin(  15)+xin(  64)*yin(  67)*zin(  63)+xin( 112)*yin( 115)*zin( 111)+xin( 160)*yin( 163)*zin( 159))
          eri_value(  518)=eri_value(  518)+d13bra( 29)*d12ket( 14)*(xin(  15)*yin(  20)*zin(  15)+xin(  63)*yin(  68)*zin(  63)+xin( 111)*yin( 116)*zin( 111)+xin( 159)*yin( 164)*zin( 159))
          eri_value(  519)=eri_value(  519)+d13bra( 29)*d12ket( 15)*(xin(  15)*yin(  19)*zin(  16)+xin(  63)*yin(  67)*zin(  64)+xin( 111)*yin( 115)*zin( 112)+xin( 159)*yin( 163)*zin( 160))
          eri_value(  520)=eri_value(  520)+d13bra( 29)*d12ket( 16)*(xin(  14)*yin(  21)*zin(  15)+xin(  62)*yin(  69)*zin(  63)+xin( 110)*yin( 117)*zin( 111)+xin( 158)*yin( 165)*zin( 159))
          eri_value(  521)=eri_value(  521)+d13bra( 29)*d12ket( 17)*(xin(  13)*yin(  22)*zin(  15)+xin(  61)*yin(  70)*zin(  63)+xin( 109)*yin( 118)*zin( 111)+xin( 157)*yin( 166)*zin( 159))
          eri_value(  522)=eri_value(  522)+d13bra( 29)*d12ket( 18)*(xin(  13)*yin(  21)*zin(  16)+xin(  61)*yin(  69)*zin(  64)+xin( 109)*yin( 117)*zin( 112)+xin( 157)*yin( 165)*zin( 160))
          eri_value(  523)=eri_value(  523)+d13bra( 30)*d12ket(  1)*(xin(  18)*yin(  13)*zin(  19)+xin(  66)*yin(  61)*zin(  67)+xin( 114)*yin( 109)*zin( 115)+xin( 162)*yin( 157)*zin( 163))
          eri_value(  524)=eri_value(  524)+d13bra( 30)*d12ket(  2)*(xin(  17)*yin(  14)*zin(  19)+xin(  65)*yin(  62)*zin(  67)+xin( 113)*yin( 110)*zin( 115)+xin( 161)*yin( 158)*zin( 163))
          eri_value(  525)=eri_value(  525)+d13bra( 30)*d12ket(  3)*(xin(  17)*yin(  13)*zin(  20)+xin(  65)*yin(  61)*zin(  68)+xin( 113)*yin( 109)*zin( 116)+xin( 161)*yin( 157)*zin( 164))
          eri_value(  526)=eri_value(  526)+d13bra( 30)*d12ket(  4)*(xin(  14)*yin(  17)*zin(  19)+xin(  62)*yin(  65)*zin(  67)+xin( 110)*yin( 113)*zin( 115)+xin( 158)*yin( 161)*zin( 163))
          eri_value(  527)=eri_value(  527)+d13bra( 30)*d12ket(  5)*(xin(  13)*yin(  18)*zin(  19)+xin(  61)*yin(  66)*zin(  67)+xin( 109)*yin( 114)*zin( 115)+xin( 157)*yin( 162)*zin( 163))
          eri_value(  528)=eri_value(  528)+d13bra( 30)*d12ket(  6)*(xin(  13)*yin(  17)*zin(  20)+xin(  61)*yin(  65)*zin(  68)+xin( 109)*yin( 113)*zin( 116)+xin( 157)*yin( 161)*zin( 164))
          eri_value(  529)=eri_value(  529)+d13bra( 30)*d12ket(  7)*(xin(  14)*yin(  13)*zin(  23)+xin(  62)*yin(  61)*zin(  71)+xin( 110)*yin( 109)*zin( 119)+xin( 158)*yin( 157)*zin( 167))
          eri_value(  530)=eri_value(  530)+d13bra( 30)*d12ket(  8)*(xin(  13)*yin(  14)*zin(  23)+xin(  61)*yin(  62)*zin(  71)+xin( 109)*yin( 110)*zin( 119)+xin( 157)*yin( 158)*zin( 167))
          eri_value(  531)=eri_value(  531)+d13bra( 30)*d12ket(  9)*(xin(  13)*yin(  13)*zin(  24)+xin(  61)*yin(  61)*zin(  72)+xin( 109)*yin( 109)*zin( 120)+xin( 157)*yin( 157)*zin( 168))
          eri_value(  532)=eri_value(  532)+d13bra( 30)*d12ket( 10)*(xin(  16)*yin(  15)*zin(  19)+xin(  64)*yin(  63)*zin(  67)+xin( 112)*yin( 111)*zin( 115)+xin( 160)*yin( 159)*zin( 163))
          eri_value(  533)=eri_value(  533)+d13bra( 30)*d12ket( 11)*(xin(  15)*yin(  16)*zin(  19)+xin(  63)*yin(  64)*zin(  67)+xin( 111)*yin( 112)*zin( 115)+xin( 159)*yin( 160)*zin( 163))
          eri_value(  534)=eri_value(  534)+d13bra( 30)*d12ket( 12)*(xin(  15)*yin(  15)*zin(  20)+xin(  63)*yin(  63)*zin(  68)+xin( 111)*yin( 111)*zin( 116)+xin( 159)*yin( 159)*zin( 164))
          eri_value(  535)=eri_value(  535)+d13bra( 30)*d12ket( 13)*(xin(  16)*yin(  13)*zin(  21)+xin(  64)*yin(  61)*zin(  69)+xin( 112)*yin( 109)*zin( 117)+xin( 160)*yin( 157)*zin( 165))
          eri_value(  536)=eri_value(  536)+d13bra( 30)*d12ket( 14)*(xin(  15)*yin(  14)*zin(  21)+xin(  63)*yin(  62)*zin(  69)+xin( 111)*yin( 110)*zin( 117)+xin( 159)*yin( 158)*zin( 165))
          eri_value(  537)=eri_value(  537)+d13bra( 30)*d12ket( 15)*(xin(  15)*yin(  13)*zin(  22)+xin(  63)*yin(  61)*zin(  70)+xin( 111)*yin( 109)*zin( 118)+xin( 159)*yin( 157)*zin( 166))
          eri_value(  538)=eri_value(  538)+d13bra( 30)*d12ket( 16)*(xin(  14)*yin(  15)*zin(  21)+xin(  62)*yin(  63)*zin(  69)+xin( 110)*yin( 111)*zin( 117)+xin( 158)*yin( 159)*zin( 165))
          eri_value(  539)=eri_value(  539)+d13bra( 30)*d12ket( 17)*(xin(  13)*yin(  16)*zin(  21)+xin(  61)*yin(  64)*zin(  69)+xin( 109)*yin( 112)*zin( 117)+xin( 157)*yin( 160)*zin( 165))
          eri_value(  540)=eri_value(  540)+d13bra( 30)*d12ket( 18)*(xin(  13)*yin(  15)*zin(  22)+xin(  61)*yin(  63)*zin(  70)+xin( 109)*yin( 111)*zin( 118)+xin( 157)*yin( 159)*zin( 166))

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
                                    ip = (i - 1)*54 ! Stride between functions in i

                                    do j = 1, 3 ! # of cartesians in j

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

                              deallocate (n13bra)
                              deallocate (xint13bra)
                              deallocate (n12ket)
                              deallocate (xint12ket)

                              end subroutine int3121
                              end submodule
