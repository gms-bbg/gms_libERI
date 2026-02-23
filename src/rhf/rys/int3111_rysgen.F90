! The total angular momentum of this class is:           6
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3111_impl
contains
  module subroutine int3111(pf_pair, pp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pf_pair, pp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n13bra(:), n11ket(:)
    real(dp), allocatable :: xint13bra(:), xint11ket(:)
    integer(kind=int64) :: npfbra, nppket
    real(dp) :: scutpfbra, scutppket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
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
    real(dp) :: xin(128), yin(128), zin(128)
    real(dp) :: eri_value(270)
    real(dp) :: d13bra(30), d11ket(9)
    integer(kind=int64) :: ix(10), jx(3), kx(3), lx(3)
    integer(kind=int64) :: iy(10), jy(3), ky(3), ly(3)
    integer(kind=int64) :: iz(10), jz(3), kz(3), lz(3)
    integer(kind=int64) :: in(5), in1(5), kn(3)
    integer(kind=int64) :: ijx(30), ijy(30), ijz(30)
    integer(kind=int64) :: klx(9), kly(9), klz(9)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 9
    in1(3) = 17
    in1(4) = 25
    in1(5) = 29

    kn(1) = 0
    kn(2) = 2
    kn(3) = 3

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 1
    lx(2) = 0
    lx(3) = 0

    kx(1) = 2
    kx(2) = 0
    kx(3) = 0

    jx(1) = 4
    jx(2) = 0
    jx(3) = 0

    ix(1) = 25
    ix(2) = 1
    ix(3) = 1
    ix(4) = 17
    ix(5) = 17
    ix(6) = 9
    ix(7) = 1
    ix(8) = 9
    ix(9) = 1
    ix(10) = 9

    ! y-arrays

    ly(1) = 0
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 2
    ky(3) = 0

    jy(1) = 0
    jy(2) = 4
    jy(3) = 0

    iy(1) = 1
    iy(2) = 25
    iy(3) = 1
    iy(4) = 9
    iy(5) = 1
    iy(6) = 17
    iy(7) = 17
    iy(8) = 1
    iy(9) = 9
    iy(10) = 9

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 4

    iz(1) = 1
    iz(2) = 1
    iz(3) = 25
    iz(4) = 1
    iz(5) = 9
    iz(6) = 1
    iz(7) = 9
    iz(8) = 17
    iz(9) = 17
    iz(10) = 9

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 29
    ijx(2) = 25
    ijx(3) = 25
    ijx(4) = 5
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 5
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 21
    ijx(11) = 17
    ijx(12) = 17
    ijx(13) = 21
    ijx(14) = 17
    ijx(15) = 17
    ijx(16) = 13
    ijx(17) = 9
    ijx(18) = 9
    ijx(19) = 5
    ijx(20) = 1
    ijx(21) = 1
    ijx(22) = 13
    ijx(23) = 9
    ijx(24) = 9
    ijx(25) = 5
    ijx(26) = 1
    ijx(27) = 1
    ijx(28) = 13
    ijx(29) = 9
    ijx(30) = 9

    ijy(1) = 1
    ijy(2) = 5
    ijy(3) = 1
    ijy(4) = 25
    ijy(5) = 29
    ijy(6) = 25
    ijy(7) = 1
    ijy(8) = 5
    ijy(9) = 1
    ijy(10) = 9
    ijy(11) = 13
    ijy(12) = 9
    ijy(13) = 1
    ijy(14) = 5
    ijy(15) = 1
    ijy(16) = 17
    ijy(17) = 21
    ijy(18) = 17
    ijy(19) = 17
    ijy(20) = 21
    ijy(21) = 17
    ijy(22) = 1
    ijy(23) = 5
    ijy(24) = 1
    ijy(25) = 9
    ijy(26) = 13
    ijy(27) = 9
    ijy(28) = 9
    ijy(29) = 13
    ijy(30) = 9

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 5
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 5
    ijz(7) = 25
    ijz(8) = 25
    ijz(9) = 29
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 5
    ijz(13) = 9
    ijz(14) = 9
    ijz(15) = 13
    ijz(16) = 1
    ijz(17) = 1
    ijz(18) = 5
    ijz(19) = 9
    ijz(20) = 9
    ijz(21) = 13
    ijz(22) = 17
    ijz(23) = 17
    ijz(24) = 21
    ijz(25) = 17
    ijz(26) = 17
    ijz(27) = 21
    ijz(28) = 9
    ijz(29) = 9
    ijz(30) = 13

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 3
    klx(2) = 2
    klx(3) = 2
    klx(4) = 1
    klx(5) = 0
    klx(6) = 0
    klx(7) = 1
    klx(8) = 0
    klx(9) = 0

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0
    kly(4) = 2
    kly(5) = 3
    kly(6) = 2
    kly(7) = 0
    kly(8) = 1
    kly(9) = 0

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1
    klz(4) = 0
    klz(5) = 0
    klz(6) = 1
    klz(7) = 2
    klz(8) = 2
    klz(9) = 3

    allocate (n13bra(res%n_p_shl*res%n_f_shl))
    allocate (xint13bra(res%n_p_shl*res%n_f_shl))
    allocate (n11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11ket(res%n_p_shl*(res%n_p_shl + 1)/2))

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

    scutppket = cutoff_schwarz/maxval(pp_pair%xints)
    nppket = 0
    do ij = 1, res%n_p_shl*(res%n_p_shl + 1)/2
      if (pp_pair%xints(ij) .ge. scutppket) then
        nppket = nppket + 1
        xint11ket(nppket) = pp_pair%xints(ij)
        n11ket(nppket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((npfbra*nppket) .le. nchunksize_int64) nchunksize_int64 = npfbra*nppket
    ntile = int(npfbra*nppket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = npfbra*nppket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npfbra, xint13bra, n13bra, xint11ket, n11ket, pf_pair, pp_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d11ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d13bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, npfbra) + 1
              kl_tmp = (iquart - 1)/npfbra + 1

              test = xint13bra(ij_tmp)*xint11ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n13bra(ij_tmp)
                kl = n11ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_p_shl(jsh_tmp)
                ksh = res%i_p_shl(ksh_tmp)
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

                  t_expon_cd = pp_pair%t_expon_ab(pp_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = pp_pair%expon_a(pp_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = pp_pair%expon_b(pp_pair%pair_loc(kl) + ket_loop)
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

                  d11ket(1) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(2) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(3) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(4) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(5) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(6) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(7) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(8) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d11ket(9) = pp_pair%d_coeff_alt(pp_pair%pair_loc(kl) + ket_loop)*twopi_5_2

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

                                      ! i2 = in(2) =    9
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(9) = xc00
                                      yin(9) = yc00
                                      zin(9) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   11
                                      ! i2 =    9

                                      xin(11) = xcp00*xin(9) + cp10
                                      yin(11) = ycp00*yin(9) + cp10
                                      zin(11) = zcp00*zin(9) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    9

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   17
                                      ! i3 =    1
                                      ! i4 =    9

                                      xin(17) = c10*xin(1) + xc00*xin(9)
                                      yin(17) = c10*yin(1) + yc00*yin(9)
                                      zin(17) = c10*zin(1) + zc00*zin(9)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   19
                                      ! i5 =   17
                                      ! i4 =    9

                                      xin(19) = xcp00*xin(17) + cp10*xin(9)
                                      yin(19) = ycp00*yin(17) + cp10*yin(9)
                                      zin(19) = zcp00*zin(17) + cp10*zin(9)

                                      ! ------------------

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   17

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   25
                                      ! i3 =    9
                                      ! i4 =   17

                                      xin(25) = c10*xin(9) + xc00*xin(17)
                                      yin(25) = c10*yin(9) + yc00*yin(17)
                                      zin(25) = c10*zin(9) + zc00*zin(17)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   27
                                      ! i5 =   25
                                      ! i4 =   17

                                      xin(27) = xcp00*xin(25) + cp10*xin(17)
                                      yin(27) = ycp00*yin(25) + cp10*yin(17)
                                      zin(27) = zcp00*zin(25) + cp10*zin(17)

                                      ! ------------------

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   25

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   29
                                      ! i3 =   17
                                      ! i4 =   25

                                      xin(29) = c10*xin(17) + xc00*xin(25)
                                      yin(29) = c10*yin(17) + yc00*yin(25)
                                      zin(29) = c10*zin(17) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   31
                                      ! i5 =   29
                                      ! i4 =   25

                                      xin(31) = xcp00*xin(29) + cp10*xin(25)
                                      yin(31) = ycp00*yin(29) + cp10*yin(25)
                                      zin(31) = zcp00*zin(29) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    3

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    4
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(4) = cp01*xin(1) + xcp00*xin(3)
                                      yin(4) = cp01*yin(1) + ycp00*yin(3)
                                      zin(4) = cp01*zin(1) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   12

                                      xin(12) = xc00*xin(4) + c01*xin(3)
                                      yin(12) = yc00*yin(4) + c01*yin(3)
                                      zin(12) = zc00*zin(4) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    4

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    9

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   17

                                      xin(20) = c10*xin(4) + xc00*xin(12) + c01*xin(11)
                                      yin(20) = c10*yin(4) + yc00*yin(12) + c01*yin(11)
                                      zin(20) = c10*zin(4) + zc00*zin(12) + c01*zin(11)

                                      c10 = c10 + b10

                                      ! i3 = i4 =    9
                                      ! i4 = i5 =   17

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   25

                                      xin(28) = c10*xin(12) + xc00*xin(20) + c01*xin(19)
                                      yin(28) = c10*yin(12) + yc00*yin(20) + c01*yin(19)
                                      zin(28) = c10*zin(12) + zc00*zin(20) + c01*zin(19)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   17
                                      ! i4 = i5 =   25

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   29

                                      xin(32) = c10*xin(20) + xc00*xin(28) + c01*xin(27)
                                      yin(32) = c10*yin(20) + yc00*yin(28) + c01*yin(27)
                                      zin(32) = c10*zin(20) + zc00*zin(28) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   29

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   29

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   29

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   25

                                      xin(29) = xin(29) + dxij*xin(25)
                                      yin(29) = yin(29) + dyij*yin(25)
                                      zin(29) = zin(29) + dzij*zin(25)

                                      ! i3 = i4 =   25
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

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

                                      xin(21) = xin(25) + dxij*xin(17)
                                      yin(21) = yin(25) + dyij*yin(17)
                                      zin(21) = zin(25) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   31

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   27

                                      xin(31) = xin(31) + dxij*xin(27)
                                      yin(31) = yin(31) + dyij*yin(27)
                                      zin(31) = zin(31) + dzij*zin(27)

                                      ! i3 = i4 =   27
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

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

                                      xin(23) = xin(27) + dxij*xin(19)
                                      yin(23) = yin(27) + dyij*yin(19)
                                      zin(23) = zin(27) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   32

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   28

                                      xin(32) = xin(32) + dxij*xin(28)
                                      yin(32) = yin(32) + dyij*yin(28)
                                      zin(32) = zin(32) + dzij*zin(28)

                                      ! i3 = i4 =   28
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

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

                                      xin(24) = xin(28) + dxij*xin(20)
                                      yin(24) = yin(28) + dyij*yin(20)
                                      zin(24) = zin(28) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =    4

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    3

                                      xin(4) = xin(4) + dxkl*xin(3)
                                      yin(4) = yin(4) + dykl*yin(3)
                                      zin(4) = zin(4) + dzkl*zin(3)

                                      ! i3 = i4 =    3
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    1

                                      xin(2) = xin(3) + dxkl*xin(1)
                                      yin(2) = yin(3) + dykl*yin(1)
                                      zin(2) = zin(3) + dzkl*zin(1)
                                      ! i4 = i4 + lang+1 =    4

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    5

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =    8

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =    6

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    6

                                      ! do nk = 1,    1

                                      xin(6) = xin(7) + dxkl*xin(5)
                                      yin(6) = yin(7) + dykl*yin(5)
                                      zin(6) = zin(7) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =    7

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    9

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =    9

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   12

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   11

                                      xin(12) = xin(12) + dxkl*xin(11)
                                      yin(12) = yin(12) + dykl*yin(11)
                                      zin(12) = zin(12) + dzkl*zin(11)

                                      ! i3 = i4 =   11
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   10

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   10

                                      ! do nk = 1,    1

                                      xin(10) = xin(11) + dxkl*xin(9)
                                      yin(10) = yin(11) + dykl*yin(9)
                                      zin(10) = zin(11) + dzkl*zin(9)
                                      ! i4 = i4 + lang+1 =   12

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   11

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   13

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   16

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   15

                                      xin(16) = xin(16) + dxkl*xin(15)
                                      yin(16) = yin(16) + dykl*yin(15)
                                      zin(16) = zin(16) + dzkl*zin(15)

                                      ! i3 = i4 =   15
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   14

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   14

                                      ! do nk = 1,    1

                                      xin(14) = xin(15) + dxkl*xin(13)
                                      yin(14) = yin(15) + dykl*yin(13)
                                      zin(14) = zin(15) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   16

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   15

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   17

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   20

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   19

                                      xin(20) = xin(20) + dxkl*xin(19)
                                      yin(20) = yin(20) + dykl*yin(19)
                                      zin(20) = zin(20) + dzkl*zin(19)

                                      ! i3 = i4 =   19
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   18

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   18

                                      ! do nk = 1,    1

                                      xin(18) = xin(19) + dxkl*xin(17)
                                      yin(18) = yin(19) + dykl*yin(17)
                                      zin(18) = zin(19) + dzkl*zin(17)
                                      ! i4 = i4 + lang+1 =   20

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   19

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   21

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   24

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   23

                                      xin(24) = xin(24) + dxkl*xin(23)
                                      yin(24) = yin(24) + dykl*yin(23)
                                      zin(24) = zin(24) + dzkl*zin(23)

                                      ! i3 = i4 =   23
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   22

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   22

                                      ! do nk = 1,    1

                                      xin(22) = xin(23) + dxkl*xin(21)
                                      yin(22) = yin(23) + dykl*yin(21)
                                      zin(22) = zin(23) + dzkl*zin(21)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   23

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   25

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   28

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   27

                                      xin(28) = xin(28) + dxkl*xin(27)
                                      yin(28) = yin(28) + dykl*yin(27)
                                      zin(28) = zin(28) + dzkl*zin(27)

                                      ! i3 = i4 =   27
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   26

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   26

                                      ! do nk = 1,    1

                                      xin(26) = xin(27) + dxkl*xin(25)
                                      yin(26) = yin(27) + dykl*yin(25)
                                      zin(26) = zin(27) + dzkl*zin(25)
                                      ! i4 = i4 + lang+1 =   28

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   27

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   29

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   32

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   31

                                      xin(32) = xin(32) + dxkl*xin(31)
                                      yin(32) = yin(32) + dykl*yin(31)
                                      zin(32) = zin(32) + dzkl*zin(31)

                                      ! i3 = i4 =   31
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   30

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   30

                                      ! do nk = 1,    1

                                      xin(30) = xin(31) + dxkl*xin(29)
                                      yin(30) = yin(31) + dykl*yin(29)
                                      zin(30) = zin(31) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   31

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   33

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   32

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

                                      ! i1 = in(1) =   33

                                      xin(33) = 1.0_dp
                                      yin(33) = 1.0_dp
                                      zin(33) = f00

                                      ! i2 = in(2) =   41
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(41) = xc00
                                      yin(41) = yc00
                                      zin(41) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   35

                                      xin(35) = xcp00
                                      yin(35) = ycp00
                                      zin(35) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   43
                                      ! i2 =   41

                                      xin(43) = xcp00*xin(41) + cp10
                                      yin(43) = ycp00*yin(41) + cp10
                                      zin(43) = zcp00*zin(41) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   33
                                      ! i4 = i2 =   41

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =   33
                                      ! i4 =   41

                                      xin(49) = c10*xin(33) + xc00*xin(41)
                                      yin(49) = c10*yin(33) + yc00*yin(41)
                                      zin(49) = c10*zin(33) + zc00*zin(41)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   51
                                      ! i5 =   49
                                      ! i4 =   41

                                      xin(51) = xcp00*xin(49) + cp10*xin(41)
                                      yin(51) = ycp00*yin(49) + cp10*yin(41)
                                      zin(51) = zcp00*zin(49) + cp10*zin(41)

                                      ! ------------------

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   49

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   57
                                      ! i3 =   41
                                      ! i4 =   49

                                      xin(57) = c10*xin(41) + xc00*xin(49)
                                      yin(57) = c10*yin(41) + yc00*yin(49)
                                      zin(57) = c10*zin(41) + zc00*zin(49)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   59
                                      ! i5 =   57
                                      ! i4 =   49

                                      xin(59) = xcp00*xin(57) + cp10*xin(49)
                                      yin(59) = ycp00*yin(57) + cp10*yin(49)
                                      zin(59) = zcp00*zin(57) + cp10*zin(49)

                                      ! ------------------

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   57

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   61
                                      ! i3 =   49
                                      ! i4 =   57

                                      xin(61) = c10*xin(49) + xc00*xin(57)
                                      yin(61) = c10*yin(49) + yc00*yin(57)
                                      zin(61) = c10*zin(49) + zc00*zin(57)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   63
                                      ! i5 =   61
                                      ! i4 =   57

                                      xin(63) = xcp00*xin(61) + cp10*xin(57)
                                      yin(63) = ycp00*yin(61) + cp10*yin(57)
                                      zin(63) = zcp00*zin(61) + cp10*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   33
                                      ! i4 = i1+k2 =   35

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   36
                                      ! i3 =   33
                                      ! i4 =   35

                                      xin(36) = cp01*xin(33) + xcp00*xin(35)
                                      yin(36) = cp01*yin(33) + ycp00*yin(35)
                                      zin(36) = cp01*zin(33) + zcp00*zin(35)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   44

                                      xin(44) = xc00*xin(36) + c01*xin(35)
                                      yin(44) = yc00*yin(36) + c01*yin(35)
                                      zin(44) = zc00*zin(36) + c01*zin(35)

                                      ! ------------------

                                      ! i3 = i4 =   35
                                      ! i4 = i5 =   36

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   33
                                      ! i4 = i2 =   41

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   49

                                      xin(52) = c10*xin(36) + xc00*xin(44) + c01*xin(43)
                                      yin(52) = c10*yin(36) + yc00*yin(44) + c01*yin(43)
                                      zin(52) = c10*zin(36) + zc00*zin(44) + c01*zin(43)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   57

                                      xin(60) = c10*xin(44) + xc00*xin(52) + c01*xin(51)
                                      yin(60) = c10*yin(44) + yc00*yin(52) + c01*yin(51)
                                      zin(60) = c10*zin(44) + zc00*zin(52) + c01*zin(51)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   57

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   61

                                      xin(64) = c10*xin(52) + xc00*xin(60) + c01*xin(59)
                                      yin(64) = c10*yin(52) + yc00*yin(60) + c01*yin(59)
                                      zin(64) = c10*zin(52) + zc00*zin(60) + c01*zin(59)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   61

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   61

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   61

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(61) = xin(61) + dxij*xin(57)
                                      yin(61) = yin(61) + dyij*yin(57)
                                      zin(61) = zin(61) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   37

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   37

                                      ! do ni = 1,    3

                                      xin(37) = xin(41) + dxij*xin(33)
                                      yin(37) = yin(41) + dyij*yin(33)
                                      zin(37) = zin(41) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    2

                                      xin(45) = xin(49) + dxij*xin(41)
                                      yin(45) = yin(49) + dyij*yin(41)
                                      zin(45) = zin(49) + dzij*zin(41)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    3

                                      xin(53) = xin(57) + dxij*xin(49)
                                      yin(53) = yin(57) + dyij*yin(49)
                                      zin(53) = zin(57) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   41

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   63

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(63) = xin(63) + dxij*xin(59)
                                      yin(63) = yin(63) + dyij*yin(59)
                                      zin(63) = zin(63) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   39

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   39

                                      ! do ni = 1,    3

                                      xin(39) = xin(43) + dxij*xin(35)
                                      yin(39) = yin(43) + dyij*yin(35)
                                      zin(39) = zin(43) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    2

                                      xin(47) = xin(51) + dxij*xin(43)
                                      yin(47) = yin(51) + dyij*yin(43)
                                      zin(47) = zin(51) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! ni =    3

                                      xin(55) = xin(59) + dxij*xin(51)
                                      yin(55) = yin(59) + dyij*yin(51)
                                      zin(55) = zin(59) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   43

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   64

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   60

                                      xin(64) = xin(64) + dxij*xin(60)
                                      yin(64) = yin(64) + dyij*yin(60)
                                      zin(64) = zin(64) + dzij*zin(60)

                                      ! i3 = i4 =   60
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   40

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   40

                                      ! do ni = 1,    3

                                      xin(40) = xin(44) + dxij*xin(36)
                                      yin(40) = yin(44) + dyij*yin(36)
                                      zin(40) = zin(44) + dzij*zin(36)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    2

                                      xin(48) = xin(52) + dxij*xin(44)
                                      yin(48) = yin(52) + dyij*yin(44)
                                      zin(48) = zin(52) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   56

                                      ! ni =    3

                                      xin(56) = xin(60) + dxij*xin(52)
                                      yin(56) = yin(60) + dyij*yin(52)
                                      zin(56) = zin(60) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   44

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   33

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   34

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   34

                                      ! do nk = 1,    1

                                      xin(34) = xin(35) + dxkl*xin(33)
                                      yin(34) = yin(35) + dykl*yin(33)
                                      zin(34) = zin(35) + dzkl*zin(33)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   35

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   37

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   40

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   39

                                      xin(40) = xin(40) + dxkl*xin(39)
                                      yin(40) = yin(40) + dykl*yin(39)
                                      zin(40) = zin(40) + dzkl*zin(39)

                                      ! i3 = i4 =   39
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   38

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   38

                                      ! do nk = 1,    1

                                      xin(38) = xin(39) + dxkl*xin(37)
                                      yin(38) = yin(39) + dykl*yin(37)
                                      zin(38) = zin(39) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   40

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   39

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   41

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   44

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   43

                                      xin(44) = xin(44) + dxkl*xin(43)
                                      yin(44) = yin(44) + dykl*yin(43)
                                      zin(44) = zin(44) + dzkl*zin(43)

                                      ! i3 = i4 =   43
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   42

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   42

                                      ! do nk = 1,    1

                                      xin(42) = xin(43) + dxkl*xin(41)
                                      yin(42) = yin(43) + dykl*yin(41)
                                      zin(42) = zin(43) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   43

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   45

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   48

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   47

                                      xin(48) = xin(48) + dxkl*xin(47)
                                      yin(48) = yin(48) + dykl*yin(47)
                                      zin(48) = zin(48) + dzkl*zin(47)

                                      ! i3 = i4 =   47
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   46

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   46

                                      ! do nk = 1,    1

                                      xin(46) = xin(47) + dxkl*xin(45)
                                      yin(46) = yin(47) + dykl*yin(45)
                                      zin(46) = zin(47) + dzkl*zin(45)
                                      ! i4 = i4 + lang+1 =   48

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   47

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   49

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   52

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   51

                                      xin(52) = xin(52) + dxkl*xin(51)
                                      yin(52) = yin(52) + dykl*yin(51)
                                      zin(52) = zin(52) + dzkl*zin(51)

                                      ! i3 = i4 =   51
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   50

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   50

                                      ! do nk = 1,    1

                                      xin(50) = xin(51) + dxkl*xin(49)
                                      yin(50) = yin(51) + dykl*yin(49)
                                      zin(50) = zin(51) + dzkl*zin(49)
                                      ! i4 = i4 + lang+1 =   52

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   51

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   53

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   56

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   55

                                      xin(56) = xin(56) + dxkl*xin(55)
                                      yin(56) = yin(56) + dykl*yin(55)
                                      zin(56) = zin(56) + dzkl*zin(55)

                                      ! i3 = i4 =   55
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   54

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   54

                                      ! do nk = 1,    1

                                      xin(54) = xin(55) + dxkl*xin(53)
                                      yin(54) = yin(55) + dykl*yin(53)
                                      zin(54) = zin(55) + dzkl*zin(53)
                                      ! i4 = i4 + lang+1 =   56

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   55

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   57

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   60

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   59

                                      xin(60) = xin(60) + dxkl*xin(59)
                                      yin(60) = yin(60) + dykl*yin(59)
                                      zin(60) = zin(60) + dzkl*zin(59)

                                      ! i3 = i4 =   59
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   58

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   58

                                      ! do nk = 1,    1

                                      xin(58) = xin(59) + dxkl*xin(57)
                                      yin(58) = yin(59) + dykl*yin(57)
                                      zin(58) = zin(59) + dzkl*zin(57)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   59

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   61

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   64

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   63

                                      xin(64) = xin(64) + dxkl*xin(63)
                                      yin(64) = yin(64) + dykl*yin(63)
                                      zin(64) = zin(64) + dzkl*zin(63)

                                      ! i3 = i4 =   63
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   62

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   62

                                      ! do nk = 1,    1

                                      xin(62) = xin(63) + dxkl*xin(61)
                                      yin(62) = yin(63) + dykl*yin(61)
                                      zin(62) = zin(63) + dzkl*zin(61)
                                      ! i4 = i4 + lang+1 =   64

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   63

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   65

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   64

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

                                      ! i1 = in(1) =   65

                                      xin(65) = 1.0_dp
                                      yin(65) = 1.0_dp
                                      zin(65) = f00

                                      ! i2 = in(2) =   73
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(73) = xc00
                                      yin(73) = yc00
                                      zin(73) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   67

                                      xin(67) = xcp00
                                      yin(67) = ycp00
                                      zin(67) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   75
                                      ! i2 =   73

                                      xin(75) = xcp00*xin(73) + cp10
                                      yin(75) = ycp00*yin(73) + cp10
                                      zin(75) = zcp00*zin(73) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   65
                                      ! i4 = i2 =   73

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   81
                                      ! i3 =   65
                                      ! i4 =   73

                                      xin(81) = c10*xin(65) + xc00*xin(73)
                                      yin(81) = c10*yin(65) + yc00*yin(73)
                                      zin(81) = c10*zin(65) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   83
                                      ! i5 =   81
                                      ! i4 =   73

                                      xin(83) = xcp00*xin(81) + cp10*xin(73)
                                      yin(83) = ycp00*yin(81) + cp10*yin(73)
                                      zin(83) = zcp00*zin(81) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   81

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   89
                                      ! i3 =   73
                                      ! i4 =   81

                                      xin(89) = c10*xin(73) + xc00*xin(81)
                                      yin(89) = c10*yin(73) + yc00*yin(81)
                                      zin(89) = c10*zin(73) + zc00*zin(81)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   91
                                      ! i5 =   89
                                      ! i4 =   81

                                      xin(91) = xcp00*xin(89) + cp10*xin(81)
                                      yin(91) = ycp00*yin(89) + cp10*yin(81)
                                      zin(91) = zcp00*zin(89) + cp10*zin(81)

                                      ! ------------------

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   89

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   93
                                      ! i3 =   81
                                      ! i4 =   89

                                      xin(93) = c10*xin(81) + xc00*xin(89)
                                      yin(93) = c10*yin(81) + yc00*yin(89)
                                      zin(93) = c10*zin(81) + zc00*zin(89)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   95
                                      ! i5 =   93
                                      ! i4 =   89

                                      xin(95) = xcp00*xin(93) + cp10*xin(89)
                                      yin(95) = ycp00*yin(93) + cp10*yin(89)
                                      zin(95) = zcp00*zin(93) + cp10*zin(89)

                                      ! ------------------

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   65
                                      ! i4 = i1+k2 =   67

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   68
                                      ! i3 =   65
                                      ! i4 =   67

                                      xin(68) = cp01*xin(65) + xcp00*xin(67)
                                      yin(68) = cp01*yin(65) + ycp00*yin(67)
                                      zin(68) = cp01*zin(65) + zcp00*zin(67)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   76

                                      xin(76) = xc00*xin(68) + c01*xin(67)
                                      yin(76) = yc00*yin(68) + c01*yin(67)
                                      zin(76) = zc00*zin(68) + c01*zin(67)

                                      ! ------------------

                                      ! i3 = i4 =   67
                                      ! i4 = i5 =   68

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   65
                                      ! i4 = i2 =   73

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   81

                                      xin(84) = c10*xin(68) + xc00*xin(76) + c01*xin(75)
                                      yin(84) = c10*yin(68) + yc00*yin(76) + c01*yin(75)
                                      zin(84) = c10*zin(68) + zc00*zin(76) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   81

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   89

                                      xin(92) = c10*xin(76) + xc00*xin(84) + c01*xin(83)
                                      yin(92) = c10*yin(76) + yc00*yin(84) + c01*yin(83)
                                      zin(92) = c10*zin(76) + zc00*zin(84) + c01*zin(83)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   81
                                      ! i4 = i5 =   89

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   93

                                      xin(96) = c10*xin(84) + xc00*xin(92) + c01*xin(91)
                                      yin(96) = c10*yin(84) + yc00*yin(92) + c01*yin(91)
                                      zin(96) = c10*zin(84) + zc00*zin(92) + c01*zin(91)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   93

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   69

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   69

                                      ! do ni = 1,    3

                                      xin(69) = xin(73) + dxij*xin(65)
                                      yin(69) = yin(73) + dyij*yin(65)
                                      zin(69) = zin(73) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   77

                                      ! ni =    2

                                      xin(77) = xin(81) + dxij*xin(73)
                                      yin(77) = yin(81) + dyij*yin(73)
                                      zin(77) = zin(81) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                      ! ni =    3

                                      xin(85) = xin(89) + dxij*xin(81)
                                      yin(85) = yin(89) + dyij*yin(81)
                                      zin(85) = zin(89) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   73

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   71

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   71

                                      ! do ni = 1,    3

                                      xin(71) = xin(75) + dxij*xin(67)
                                      yin(71) = yin(75) + dyij*yin(67)
                                      zin(71) = zin(75) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    2

                                      xin(79) = xin(83) + dxij*xin(75)
                                      yin(79) = yin(83) + dyij*yin(75)
                                      zin(79) = zin(83) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   87

                                      ! ni =    3

                                      xin(87) = xin(91) + dxij*xin(83)
                                      yin(87) = yin(91) + dyij*yin(83)
                                      zin(87) = zin(91) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   75

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   72

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   72

                                      ! do ni = 1,    3

                                      xin(72) = xin(76) + dxij*xin(68)
                                      yin(72) = yin(76) + dyij*yin(68)
                                      zin(72) = zin(76) + dzij*zin(68)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    2

                                      xin(80) = xin(84) + dxij*xin(76)
                                      yin(80) = yin(84) + dyij*yin(76)
                                      zin(80) = zin(84) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                      ! ni =    3

                                      xin(88) = xin(92) + dxij*xin(84)
                                      yin(88) = yin(92) + dyij*yin(84)
                                      zin(88) = zin(92) + dzij*zin(84)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   76

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   65

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   68

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   67

                                      xin(68) = xin(68) + dxkl*xin(67)
                                      yin(68) = yin(68) + dykl*yin(67)
                                      zin(68) = zin(68) + dzkl*zin(67)

                                      ! i3 = i4 =   67
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   66

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    1

                                      xin(66) = xin(67) + dxkl*xin(65)
                                      yin(66) = yin(67) + dykl*yin(65)
                                      zin(66) = zin(67) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   69

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   70

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   70

                                      ! do nk = 1,    1

                                      xin(70) = xin(71) + dxkl*xin(69)
                                      yin(70) = yin(71) + dykl*yin(69)
                                      zin(70) = zin(71) + dzkl*zin(69)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   71

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   76

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   75

                                      xin(76) = xin(76) + dxkl*xin(75)
                                      yin(76) = yin(76) + dykl*yin(75)
                                      zin(76) = zin(76) + dzkl*zin(75)

                                      ! i3 = i4 =   75
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    1

                                      xin(74) = xin(75) + dxkl*xin(73)
                                      yin(74) = yin(75) + dykl*yin(73)
                                      zin(74) = zin(75) + dzkl*zin(73)
                                      ! i4 = i4 + lang+1 =   76

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   77

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   80

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   78

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   78

                                      ! do nk = 1,    1

                                      xin(78) = xin(79) + dxkl*xin(77)
                                      yin(78) = yin(79) + dykl*yin(77)
                                      zin(78) = zin(79) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   79

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   81

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   84

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   83

                                      xin(84) = xin(84) + dxkl*xin(83)
                                      yin(84) = yin(84) + dykl*yin(83)
                                      zin(84) = zin(84) + dzkl*zin(83)

                                      ! i3 = i4 =   83
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   82

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   82

                                      ! do nk = 1,    1

                                      xin(82) = xin(83) + dxkl*xin(81)
                                      yin(82) = yin(83) + dykl*yin(81)
                                      zin(82) = zin(83) + dzkl*zin(81)
                                      ! i4 = i4 + lang+1 =   84

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   83

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   85

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   88

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   87

                                      xin(88) = xin(88) + dxkl*xin(87)
                                      yin(88) = yin(88) + dykl*yin(87)
                                      zin(88) = zin(88) + dzkl*zin(87)

                                      ! i3 = i4 =   87
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   86

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   86

                                      ! do nk = 1,    1

                                      xin(86) = xin(87) + dxkl*xin(85)
                                      yin(86) = yin(87) + dykl*yin(85)
                                      zin(86) = zin(87) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   88

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   87

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   89

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   92

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   91

                                      xin(92) = xin(92) + dxkl*xin(91)
                                      yin(92) = yin(92) + dykl*yin(91)
                                      zin(92) = zin(92) + dzkl*zin(91)

                                      ! i3 = i4 =   91
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   90

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   90

                                      ! do nk = 1,    1

                                      xin(90) = xin(91) + dxkl*xin(89)
                                      yin(90) = yin(91) + dykl*yin(89)
                                      zin(90) = zin(91) + dzkl*zin(89)
                                      ! i4 = i4 + lang+1 =   92

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   91

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   93

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =   96

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   95

                                      xin(96) = xin(96) + dxkl*xin(95)
                                      yin(96) = yin(96) + dykl*yin(95)
                                      zin(96) = zin(96) + dzkl*zin(95)

                                      ! i3 = i4 =   95
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   94

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   94

                                      ! do nk = 1,    1

                                      xin(94) = xin(95) + dxkl*xin(93)
                                      yin(94) = yin(95) + dykl*yin(93)
                                      zin(94) = zin(95) + dzkl*zin(93)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   95

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   97

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

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

                                      ! i1 = in(1) =   97

                                      xin(97) = 1.0_dp
                                      yin(97) = 1.0_dp
                                      zin(97) = f00

                                      ! i2 = in(2) =  105
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(105) = xc00
                                      yin(105) = yc00
                                      zin(105) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   99

                                      xin(99) = xcp00
                                      yin(99) = ycp00
                                      zin(99) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  107
                                      ! i2 =  105

                                      xin(107) = xcp00*xin(105) + cp10
                                      yin(107) = ycp00*yin(105) + cp10
                                      zin(107) = zcp00*zin(105) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  105

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  113
                                      ! i3 =   97
                                      ! i4 =  105

                                      xin(113) = c10*xin(97) + xc00*xin(105)
                                      yin(113) = c10*yin(97) + yc00*yin(105)
                                      zin(113) = c10*zin(97) + zc00*zin(105)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  115
                                      ! i5 =  113
                                      ! i4 =  105

                                      xin(115) = xcp00*xin(113) + cp10*xin(105)
                                      yin(115) = ycp00*yin(113) + cp10*yin(105)
                                      zin(115) = zcp00*zin(113) + cp10*zin(105)

                                      ! ------------------

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =  105
                                      ! i4 =  113

                                      xin(121) = c10*xin(105) + xc00*xin(113)
                                      yin(121) = c10*yin(105) + yc00*yin(113)
                                      zin(121) = c10*zin(105) + zc00*zin(113)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  123
                                      ! i5 =  121
                                      ! i4 =  113

                                      xin(123) = xcp00*xin(121) + cp10*xin(113)
                                      yin(123) = ycp00*yin(121) + cp10*yin(113)
                                      zin(123) = zcp00*zin(121) + cp10*zin(113)

                                      ! ------------------

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  125
                                      ! i3 =  113
                                      ! i4 =  121

                                      xin(125) = c10*xin(113) + xc00*xin(121)
                                      yin(125) = c10*yin(113) + yc00*yin(121)
                                      zin(125) = c10*zin(113) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  127
                                      ! i5 =  125
                                      ! i4 =  121

                                      xin(127) = xcp00*xin(125) + cp10*xin(121)
                                      yin(127) = ycp00*yin(125) + cp10*yin(121)
                                      zin(127) = zcp00*zin(125) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   99

                                      ! do n = 2,    2

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  100
                                      ! i3 =   97
                                      ! i4 =   99

                                      xin(100) = cp01*xin(97) + xcp00*xin(99)
                                      yin(100) = cp01*yin(97) + ycp00*yin(99)
                                      zin(100) = cp01*zin(97) + zcp00*zin(99)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  108

                                      xin(108) = xc00*xin(100) + c01*xin(99)
                                      yin(108) = yc00*yin(100) + c01*yin(99)
                                      zin(108) = zc00*zin(100) + c01*zin(99)

                                      ! ------------------

                                      ! i3 = i4 =   99
                                      ! i4 = i5 =  100

                                      ! n =    3

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    2

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  105

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  113

                                      xin(116) = c10*xin(100) + xc00*xin(108) + c01*xin(107)
                                      yin(116) = c10*yin(100) + yc00*yin(108) + c01*yin(107)
                                      zin(116) = c10*zin(100) + zc00*zin(108) + c01*zin(107)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  105
                                      ! i4 = i5 =  113

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  121

                                      xin(124) = c10*xin(108) + xc00*xin(116) + c01*xin(115)
                                      yin(124) = c10*yin(108) + yc00*yin(116) + c01*yin(115)
                                      zin(124) = c10*zin(108) + zc00*zin(116) + c01*zin(115)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  113
                                      ! i4 = i5 =  121

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  125

                                      xin(128) = c10*xin(116) + xc00*xin(124) + c01*xin(123)
                                      yin(128) = c10*yin(116) + yc00*yin(124) + c01*yin(123)
                                      zin(128) = c10*zin(116) + zc00*zin(124) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  125

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    3

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  125

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  125

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  121

                                      xin(125) = xin(125) + dxij*xin(121)
                                      yin(125) = yin(125) + dyij*yin(121)
                                      zin(125) = zin(125) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  101

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  101

                                      ! do ni = 1,    3

                                      xin(101) = xin(105) + dxij*xin(97)
                                      yin(101) = yin(105) + dyij*yin(97)
                                      zin(101) = zin(105) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! ni =    2

                                      xin(109) = xin(113) + dxij*xin(105)
                                      yin(109) = yin(113) + dyij*yin(105)
                                      zin(109) = zin(113) + dzij*zin(105)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    3

                                      xin(117) = xin(121) + dxij*xin(113)
                                      yin(117) = yin(121) + dyij*yin(113)
                                      zin(117) = zin(121) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  105

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  127

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  123

                                      xin(127) = xin(127) + dxij*xin(123)
                                      yin(127) = yin(127) + dyij*yin(123)
                                      zin(127) = zin(127) + dzij*zin(123)

                                      ! i3 = i4 =  123
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  103

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(107) + dxij*xin(99)
                                      yin(103) = yin(107) + dyij*yin(99)
                                      zin(103) = zin(107) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  111

                                      ! ni =    2

                                      xin(111) = xin(115) + dxij*xin(107)
                                      yin(111) = yin(115) + dyij*yin(107)
                                      zin(111) = zin(115) + dzij*zin(107)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    3

                                      xin(119) = xin(123) + dxij*xin(115)
                                      yin(119) = yin(123) + dyij*yin(115)
                                      zin(119) = zin(123) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  107

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  128

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  124

                                      xin(128) = xin(128) + dxij*xin(124)
                                      yin(128) = yin(128) + dyij*yin(124)
                                      zin(128) = zin(128) + dzij*zin(124)

                                      ! i3 = i4 =  124
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  104

                                      ! do nj = 1,    1

                                      ! i4 = i3 =  104

                                      ! do ni = 1,    3

                                      xin(104) = xin(108) + dxij*xin(100)
                                      yin(104) = yin(108) + dyij*yin(100)
                                      zin(104) = zin(108) + dzij*zin(100)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  112

                                      ! ni =    2

                                      xin(112) = xin(116) + dxij*xin(108)
                                      yin(112) = yin(116) + dyij*yin(108)
                                      zin(112) = zin(116) + dzij*zin(108)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    3

                                      xin(120) = xin(124) + dxij*xin(116)
                                      yin(120) = yin(124) + dyij*yin(116)
                                      zin(120) = zin(124) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  108

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    3

                                      ! iaa = i1 =   97

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  100

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   99

                                      xin(100) = xin(100) + dxkl*xin(99)
                                      yin(100) = yin(100) + dykl*yin(99)
                                      zin(100) = zin(100) + dzkl*zin(99)

                                      ! i3 = i4 =   99
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =   98

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   98

                                      ! do nk = 1,    1

                                      xin(98) = xin(99) + dxkl*xin(97)
                                      yin(98) = yin(99) + dykl*yin(97)
                                      zin(98) = zin(99) + dzkl*zin(97)
                                      ! i4 = i4 + lang+1 =  100

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =   99

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  101

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  104

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  103

                                      xin(104) = xin(104) + dxkl*xin(103)
                                      yin(104) = yin(104) + dykl*yin(103)
                                      zin(104) = zin(104) + dzkl*zin(103)

                                      ! i3 = i4 =  103
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  102

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  102

                                      ! do nk = 1,    1

                                      xin(102) = xin(103) + dxkl*xin(101)
                                      yin(102) = yin(103) + dykl*yin(101)
                                      zin(102) = zin(103) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  103

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  105

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  106

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  106

                                      ! do nk = 1,    1

                                      xin(106) = xin(107) + dxkl*xin(105)
                                      yin(106) = yin(107) + dykl*yin(105)
                                      zin(106) = zin(107) + dzkl*zin(105)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  107

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  109

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  112

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  111

                                      xin(112) = xin(112) + dxkl*xin(111)
                                      yin(112) = yin(112) + dykl*yin(111)
                                      zin(112) = zin(112) + dzkl*zin(111)

                                      ! i3 = i4 =  111
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  110

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  110

                                      ! do nk = 1,    1

                                      xin(110) = xin(111) + dxkl*xin(109)
                                      yin(110) = yin(111) + dykl*yin(109)
                                      zin(110) = zin(111) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  112

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  111

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  113

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  116

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  115

                                      xin(116) = xin(116) + dxkl*xin(115)
                                      yin(116) = yin(116) + dykl*yin(115)
                                      zin(116) = zin(116) + dzkl*zin(115)

                                      ! i3 = i4 =  115
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  114

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  114

                                      ! do nk = 1,    1

                                      xin(114) = xin(115) + dxkl*xin(113)
                                      yin(114) = yin(115) + dykl*yin(113)
                                      zin(114) = zin(115) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  115

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  117

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  120

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  119

                                      xin(120) = xin(120) + dxkl*xin(119)
                                      yin(120) = yin(120) + dykl*yin(119)
                                      zin(120) = zin(120) + dzkl*zin(119)

                                      ! i3 = i4 =  119
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  118

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  118

                                      ! do nk = 1,    1

                                      xin(118) = xin(119) + dxkl*xin(117)
                                      yin(118) = yin(119) + dykl*yin(117)
                                      zin(118) = zin(119) + dzkl*zin(117)
                                      ! i4 = i4 + lang+1 =  120

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  119

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  121

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  124

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  123

                                      xin(124) = xin(124) + dxkl*xin(123)
                                      yin(124) = yin(124) + dykl*yin(123)
                                      zin(124) = zin(124) + dzkl*zin(123)

                                      ! i3 = i4 =  123
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  122

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  122

                                      ! do nk = 1,    1

                                      xin(122) = xin(123) + dxkl*xin(121)
                                      yin(122) = yin(123) + dykl*yin(121)
                                      zin(122) = zin(123) + dzkl*zin(121)
                                      ! i4 = i4 + lang+1 =  124

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  123

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  125

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    2

                                      ! i3 = ib+i5 =  128

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  127

                                      xin(128) = xin(128) + dxkl*xin(127)
                                      yin(128) = yin(128) + dykl*yin(127)
                                      zin(128) = zin(128) + dzkl*zin(127)

                                      ! i3 = i4 =  127
                                      ! nm = nm -1 =    1

                                      ! end do

                                      ! min = min + 1 =    2

                                      ! end do

                                      ! i3 = ib + 1 =  126

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  126

                                      ! do nk = 1,    1

                                      xin(126) = xin(127) + dxkl*xin(125)
                                      yin(126) = yin(127) + dykl*yin(125)
                                      zin(126) = zin(127) + dzkl*zin(125)
                                      ! i4 = i4 + lang+1 =  128

                                      ! nk =    2

                                      ! end do

                                      ! i3 = i3 + 1 =  127

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  129

                                      ! nj = nj + 1 =    2

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  128

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(    1)=eri_value(    1)+d13bra(  1)*d11ket(  1)*(xin(  32)*yin(   1)*zin(   1)+xin(  64)*yin(  33)*zin(  33)+xin(  96)*yin(  65)*zin(  65)+xin( 128)*yin(  97)*zin(  97))
          eri_value(    2)=eri_value(    2)+d13bra(  1)*d11ket(  2)*(xin(  31)*yin(   2)*zin(   1)+xin(  63)*yin(  34)*zin(  33)+xin(  95)*yin(  66)*zin(  65)+xin( 127)*yin(  98)*zin(  97))
          eri_value(    3)=eri_value(    3)+d13bra(  1)*d11ket(  3)*(xin(  31)*yin(   1)*zin(   2)+xin(  63)*yin(  33)*zin(  34)+xin(  95)*yin(  65)*zin(  66)+xin( 127)*yin(  97)*zin(  98))
          eri_value(    4)=eri_value(    4)+d13bra(  1)*d11ket(  4)*(xin(  30)*yin(   3)*zin(   1)+xin(  62)*yin(  35)*zin(  33)+xin(  94)*yin(  67)*zin(  65)+xin( 126)*yin(  99)*zin(  97))
          eri_value(    5)=eri_value(    5)+d13bra(  1)*d11ket(  5)*(xin(  29)*yin(   4)*zin(   1)+xin(  61)*yin(  36)*zin(  33)+xin(  93)*yin(  68)*zin(  65)+xin( 125)*yin( 100)*zin(  97))
          eri_value(    6)=eri_value(    6)+d13bra(  1)*d11ket(  6)*(xin(  29)*yin(   3)*zin(   2)+xin(  61)*yin(  35)*zin(  34)+xin(  93)*yin(  67)*zin(  66)+xin( 125)*yin(  99)*zin(  98))
          eri_value(    7)=eri_value(    7)+d13bra(  1)*d11ket(  7)*(xin(  30)*yin(   1)*zin(   3)+xin(  62)*yin(  33)*zin(  35)+xin(  94)*yin(  65)*zin(  67)+xin( 126)*yin(  97)*zin(  99))
          eri_value(    8)=eri_value(    8)+d13bra(  1)*d11ket(  8)*(xin(  29)*yin(   2)*zin(   3)+xin(  61)*yin(  34)*zin(  35)+xin(  93)*yin(  66)*zin(  67)+xin( 125)*yin(  98)*zin(  99))
          eri_value(    9)=eri_value(    9)+d13bra(  1)*d11ket(  9)*(xin(  29)*yin(   1)*zin(   4)+xin(  61)*yin(  33)*zin(  36)+xin(  93)*yin(  65)*zin(  68)+xin( 125)*yin(  97)*zin( 100))
          eri_value(   10)=eri_value(   10)+d13bra(  2)*d11ket(  1)*(xin(  28)*yin(   5)*zin(   1)+xin(  60)*yin(  37)*zin(  33)+xin(  92)*yin(  69)*zin(  65)+xin( 124)*yin( 101)*zin(  97))
          eri_value(   11)=eri_value(   11)+d13bra(  2)*d11ket(  2)*(xin(  27)*yin(   6)*zin(   1)+xin(  59)*yin(  38)*zin(  33)+xin(  91)*yin(  70)*zin(  65)+xin( 123)*yin( 102)*zin(  97))
          eri_value(   12)=eri_value(   12)+d13bra(  2)*d11ket(  3)*(xin(  27)*yin(   5)*zin(   2)+xin(  59)*yin(  37)*zin(  34)+xin(  91)*yin(  69)*zin(  66)+xin( 123)*yin( 101)*zin(  98))
          eri_value(   13)=eri_value(   13)+d13bra(  2)*d11ket(  4)*(xin(  26)*yin(   7)*zin(   1)+xin(  58)*yin(  39)*zin(  33)+xin(  90)*yin(  71)*zin(  65)+xin( 122)*yin( 103)*zin(  97))
          eri_value(   14)=eri_value(   14)+d13bra(  2)*d11ket(  5)*(xin(  25)*yin(   8)*zin(   1)+xin(  57)*yin(  40)*zin(  33)+xin(  89)*yin(  72)*zin(  65)+xin( 121)*yin( 104)*zin(  97))
          eri_value(   15)=eri_value(   15)+d13bra(  2)*d11ket(  6)*(xin(  25)*yin(   7)*zin(   2)+xin(  57)*yin(  39)*zin(  34)+xin(  89)*yin(  71)*zin(  66)+xin( 121)*yin( 103)*zin(  98))
          eri_value(   16)=eri_value(   16)+d13bra(  2)*d11ket(  7)*(xin(  26)*yin(   5)*zin(   3)+xin(  58)*yin(  37)*zin(  35)+xin(  90)*yin(  69)*zin(  67)+xin( 122)*yin( 101)*zin(  99))
          eri_value(   17)=eri_value(   17)+d13bra(  2)*d11ket(  8)*(xin(  25)*yin(   6)*zin(   3)+xin(  57)*yin(  38)*zin(  35)+xin(  89)*yin(  70)*zin(  67)+xin( 121)*yin( 102)*zin(  99))
          eri_value(   18)=eri_value(   18)+d13bra(  2)*d11ket(  9)*(xin(  25)*yin(   5)*zin(   4)+xin(  57)*yin(  37)*zin(  36)+xin(  89)*yin(  69)*zin(  68)+xin( 121)*yin( 101)*zin( 100))
          eri_value(   19)=eri_value(   19)+d13bra(  3)*d11ket(  1)*(xin(  28)*yin(   1)*zin(   5)+xin(  60)*yin(  33)*zin(  37)+xin(  92)*yin(  65)*zin(  69)+xin( 124)*yin(  97)*zin( 101))
          eri_value(   20)=eri_value(   20)+d13bra(  3)*d11ket(  2)*(xin(  27)*yin(   2)*zin(   5)+xin(  59)*yin(  34)*zin(  37)+xin(  91)*yin(  66)*zin(  69)+xin( 123)*yin(  98)*zin( 101))
          eri_value(   21)=eri_value(   21)+d13bra(  3)*d11ket(  3)*(xin(  27)*yin(   1)*zin(   6)+xin(  59)*yin(  33)*zin(  38)+xin(  91)*yin(  65)*zin(  70)+xin( 123)*yin(  97)*zin( 102))
          eri_value(   22)=eri_value(   22)+d13bra(  3)*d11ket(  4)*(xin(  26)*yin(   3)*zin(   5)+xin(  58)*yin(  35)*zin(  37)+xin(  90)*yin(  67)*zin(  69)+xin( 122)*yin(  99)*zin( 101))
          eri_value(   23)=eri_value(   23)+d13bra(  3)*d11ket(  5)*(xin(  25)*yin(   4)*zin(   5)+xin(  57)*yin(  36)*zin(  37)+xin(  89)*yin(  68)*zin(  69)+xin( 121)*yin( 100)*zin( 101))
          eri_value(   24)=eri_value(   24)+d13bra(  3)*d11ket(  6)*(xin(  25)*yin(   3)*zin(   6)+xin(  57)*yin(  35)*zin(  38)+xin(  89)*yin(  67)*zin(  70)+xin( 121)*yin(  99)*zin( 102))
          eri_value(   25)=eri_value(   25)+d13bra(  3)*d11ket(  7)*(xin(  26)*yin(   1)*zin(   7)+xin(  58)*yin(  33)*zin(  39)+xin(  90)*yin(  65)*zin(  71)+xin( 122)*yin(  97)*zin( 103))
          eri_value(   26)=eri_value(   26)+d13bra(  3)*d11ket(  8)*(xin(  25)*yin(   2)*zin(   7)+xin(  57)*yin(  34)*zin(  39)+xin(  89)*yin(  66)*zin(  71)+xin( 121)*yin(  98)*zin( 103))
          eri_value(   27)=eri_value(   27)+d13bra(  3)*d11ket(  9)*(xin(  25)*yin(   1)*zin(   8)+xin(  57)*yin(  33)*zin(  40)+xin(  89)*yin(  65)*zin(  72)+xin( 121)*yin(  97)*zin( 104))
          eri_value(   28)=eri_value(   28)+d13bra(  4)*d11ket(  1)*(xin(   8)*yin(  25)*zin(   1)+xin(  40)*yin(  57)*zin(  33)+xin(  72)*yin(  89)*zin(  65)+xin( 104)*yin( 121)*zin(  97))
          eri_value(   29)=eri_value(   29)+d13bra(  4)*d11ket(  2)*(xin(   7)*yin(  26)*zin(   1)+xin(  39)*yin(  58)*zin(  33)+xin(  71)*yin(  90)*zin(  65)+xin( 103)*yin( 122)*zin(  97))
          eri_value(   30)=eri_value(   30)+d13bra(  4)*d11ket(  3)*(xin(   7)*yin(  25)*zin(   2)+xin(  39)*yin(  57)*zin(  34)+xin(  71)*yin(  89)*zin(  66)+xin( 103)*yin( 121)*zin(  98))
          eri_value(   31)=eri_value(   31)+d13bra(  4)*d11ket(  4)*(xin(   6)*yin(  27)*zin(   1)+xin(  38)*yin(  59)*zin(  33)+xin(  70)*yin(  91)*zin(  65)+xin( 102)*yin( 123)*zin(  97))
          eri_value(   32)=eri_value(   32)+d13bra(  4)*d11ket(  5)*(xin(   5)*yin(  28)*zin(   1)+xin(  37)*yin(  60)*zin(  33)+xin(  69)*yin(  92)*zin(  65)+xin( 101)*yin( 124)*zin(  97))
          eri_value(   33)=eri_value(   33)+d13bra(  4)*d11ket(  6)*(xin(   5)*yin(  27)*zin(   2)+xin(  37)*yin(  59)*zin(  34)+xin(  69)*yin(  91)*zin(  66)+xin( 101)*yin( 123)*zin(  98))
          eri_value(   34)=eri_value(   34)+d13bra(  4)*d11ket(  7)*(xin(   6)*yin(  25)*zin(   3)+xin(  38)*yin(  57)*zin(  35)+xin(  70)*yin(  89)*zin(  67)+xin( 102)*yin( 121)*zin(  99))
          eri_value(   35)=eri_value(   35)+d13bra(  4)*d11ket(  8)*(xin(   5)*yin(  26)*zin(   3)+xin(  37)*yin(  58)*zin(  35)+xin(  69)*yin(  90)*zin(  67)+xin( 101)*yin( 122)*zin(  99))
          eri_value(   36)=eri_value(   36)+d13bra(  4)*d11ket(  9)*(xin(   5)*yin(  25)*zin(   4)+xin(  37)*yin(  57)*zin(  36)+xin(  69)*yin(  89)*zin(  68)+xin( 101)*yin( 121)*zin( 100))
          eri_value(   37)=eri_value(   37)+d13bra(  5)*d11ket(  1)*(xin(   4)*yin(  29)*zin(   1)+xin(  36)*yin(  61)*zin(  33)+xin(  68)*yin(  93)*zin(  65)+xin( 100)*yin( 125)*zin(  97))
          eri_value(   38)=eri_value(   38)+d13bra(  5)*d11ket(  2)*(xin(   3)*yin(  30)*zin(   1)+xin(  35)*yin(  62)*zin(  33)+xin(  67)*yin(  94)*zin(  65)+xin(  99)*yin( 126)*zin(  97))
          eri_value(   39)=eri_value(   39)+d13bra(  5)*d11ket(  3)*(xin(   3)*yin(  29)*zin(   2)+xin(  35)*yin(  61)*zin(  34)+xin(  67)*yin(  93)*zin(  66)+xin(  99)*yin( 125)*zin(  98))
          eri_value(   40)=eri_value(   40)+d13bra(  5)*d11ket(  4)*(xin(   2)*yin(  31)*zin(   1)+xin(  34)*yin(  63)*zin(  33)+xin(  66)*yin(  95)*zin(  65)+xin(  98)*yin( 127)*zin(  97))
          eri_value(   41)=eri_value(   41)+d13bra(  5)*d11ket(  5)*(xin(   1)*yin(  32)*zin(   1)+xin(  33)*yin(  64)*zin(  33)+xin(  65)*yin(  96)*zin(  65)+xin(  97)*yin( 128)*zin(  97))
          eri_value(   42)=eri_value(   42)+d13bra(  5)*d11ket(  6)*(xin(   1)*yin(  31)*zin(   2)+xin(  33)*yin(  63)*zin(  34)+xin(  65)*yin(  95)*zin(  66)+xin(  97)*yin( 127)*zin(  98))
          eri_value(   43)=eri_value(   43)+d13bra(  5)*d11ket(  7)*(xin(   2)*yin(  29)*zin(   3)+xin(  34)*yin(  61)*zin(  35)+xin(  66)*yin(  93)*zin(  67)+xin(  98)*yin( 125)*zin(  99))
          eri_value(   44)=eri_value(   44)+d13bra(  5)*d11ket(  8)*(xin(   1)*yin(  30)*zin(   3)+xin(  33)*yin(  62)*zin(  35)+xin(  65)*yin(  94)*zin(  67)+xin(  97)*yin( 126)*zin(  99))
          eri_value(   45)=eri_value(   45)+d13bra(  5)*d11ket(  9)*(xin(   1)*yin(  29)*zin(   4)+xin(  33)*yin(  61)*zin(  36)+xin(  65)*yin(  93)*zin(  68)+xin(  97)*yin( 125)*zin( 100))
          eri_value(   46)=eri_value(   46)+d13bra(  6)*d11ket(  1)*(xin(   4)*yin(  25)*zin(   5)+xin(  36)*yin(  57)*zin(  37)+xin(  68)*yin(  89)*zin(  69)+xin( 100)*yin( 121)*zin( 101))
          eri_value(   47)=eri_value(   47)+d13bra(  6)*d11ket(  2)*(xin(   3)*yin(  26)*zin(   5)+xin(  35)*yin(  58)*zin(  37)+xin(  67)*yin(  90)*zin(  69)+xin(  99)*yin( 122)*zin( 101))
          eri_value(   48)=eri_value(   48)+d13bra(  6)*d11ket(  3)*(xin(   3)*yin(  25)*zin(   6)+xin(  35)*yin(  57)*zin(  38)+xin(  67)*yin(  89)*zin(  70)+xin(  99)*yin( 121)*zin( 102))
          eri_value(   49)=eri_value(   49)+d13bra(  6)*d11ket(  4)*(xin(   2)*yin(  27)*zin(   5)+xin(  34)*yin(  59)*zin(  37)+xin(  66)*yin(  91)*zin(  69)+xin(  98)*yin( 123)*zin( 101))
          eri_value(   50)=eri_value(   50)+d13bra(  6)*d11ket(  5)*(xin(   1)*yin(  28)*zin(   5)+xin(  33)*yin(  60)*zin(  37)+xin(  65)*yin(  92)*zin(  69)+xin(  97)*yin( 124)*zin( 101))
          eri_value(   51)=eri_value(   51)+d13bra(  6)*d11ket(  6)*(xin(   1)*yin(  27)*zin(   6)+xin(  33)*yin(  59)*zin(  38)+xin(  65)*yin(  91)*zin(  70)+xin(  97)*yin( 123)*zin( 102))
          eri_value(   52)=eri_value(   52)+d13bra(  6)*d11ket(  7)*(xin(   2)*yin(  25)*zin(   7)+xin(  34)*yin(  57)*zin(  39)+xin(  66)*yin(  89)*zin(  71)+xin(  98)*yin( 121)*zin( 103))
          eri_value(   53)=eri_value(   53)+d13bra(  6)*d11ket(  8)*(xin(   1)*yin(  26)*zin(   7)+xin(  33)*yin(  58)*zin(  39)+xin(  65)*yin(  90)*zin(  71)+xin(  97)*yin( 122)*zin( 103))
          eri_value(   54)=eri_value(   54)+d13bra(  6)*d11ket(  9)*(xin(   1)*yin(  25)*zin(   8)+xin(  33)*yin(  57)*zin(  40)+xin(  65)*yin(  89)*zin(  72)+xin(  97)*yin( 121)*zin( 104))
          eri_value(   55)=eri_value(   55)+d13bra(  7)*d11ket(  1)*(xin(   8)*yin(   1)*zin(  25)+xin(  40)*yin(  33)*zin(  57)+xin(  72)*yin(  65)*zin(  89)+xin( 104)*yin(  97)*zin( 121))
          eri_value(   56)=eri_value(   56)+d13bra(  7)*d11ket(  2)*(xin(   7)*yin(   2)*zin(  25)+xin(  39)*yin(  34)*zin(  57)+xin(  71)*yin(  66)*zin(  89)+xin( 103)*yin(  98)*zin( 121))
          eri_value(   57)=eri_value(   57)+d13bra(  7)*d11ket(  3)*(xin(   7)*yin(   1)*zin(  26)+xin(  39)*yin(  33)*zin(  58)+xin(  71)*yin(  65)*zin(  90)+xin( 103)*yin(  97)*zin( 122))
          eri_value(   58)=eri_value(   58)+d13bra(  7)*d11ket(  4)*(xin(   6)*yin(   3)*zin(  25)+xin(  38)*yin(  35)*zin(  57)+xin(  70)*yin(  67)*zin(  89)+xin( 102)*yin(  99)*zin( 121))
          eri_value(   59)=eri_value(   59)+d13bra(  7)*d11ket(  5)*(xin(   5)*yin(   4)*zin(  25)+xin(  37)*yin(  36)*zin(  57)+xin(  69)*yin(  68)*zin(  89)+xin( 101)*yin( 100)*zin( 121))
          eri_value(   60)=eri_value(   60)+d13bra(  7)*d11ket(  6)*(xin(   5)*yin(   3)*zin(  26)+xin(  37)*yin(  35)*zin(  58)+xin(  69)*yin(  67)*zin(  90)+xin( 101)*yin(  99)*zin( 122))
          eri_value(   61)=eri_value(   61)+d13bra(  7)*d11ket(  7)*(xin(   6)*yin(   1)*zin(  27)+xin(  38)*yin(  33)*zin(  59)+xin(  70)*yin(  65)*zin(  91)+xin( 102)*yin(  97)*zin( 123))
          eri_value(   62)=eri_value(   62)+d13bra(  7)*d11ket(  8)*(xin(   5)*yin(   2)*zin(  27)+xin(  37)*yin(  34)*zin(  59)+xin(  69)*yin(  66)*zin(  91)+xin( 101)*yin(  98)*zin( 123))
          eri_value(   63)=eri_value(   63)+d13bra(  7)*d11ket(  9)*(xin(   5)*yin(   1)*zin(  28)+xin(  37)*yin(  33)*zin(  60)+xin(  69)*yin(  65)*zin(  92)+xin( 101)*yin(  97)*zin( 124))
          eri_value(   64)=eri_value(   64)+d13bra(  8)*d11ket(  1)*(xin(   4)*yin(   5)*zin(  25)+xin(  36)*yin(  37)*zin(  57)+xin(  68)*yin(  69)*zin(  89)+xin( 100)*yin( 101)*zin( 121))
          eri_value(   65)=eri_value(   65)+d13bra(  8)*d11ket(  2)*(xin(   3)*yin(   6)*zin(  25)+xin(  35)*yin(  38)*zin(  57)+xin(  67)*yin(  70)*zin(  89)+xin(  99)*yin( 102)*zin( 121))
          eri_value(   66)=eri_value(   66)+d13bra(  8)*d11ket(  3)*(xin(   3)*yin(   5)*zin(  26)+xin(  35)*yin(  37)*zin(  58)+xin(  67)*yin(  69)*zin(  90)+xin(  99)*yin( 101)*zin( 122))
          eri_value(   67)=eri_value(   67)+d13bra(  8)*d11ket(  4)*(xin(   2)*yin(   7)*zin(  25)+xin(  34)*yin(  39)*zin(  57)+xin(  66)*yin(  71)*zin(  89)+xin(  98)*yin( 103)*zin( 121))
          eri_value(   68)=eri_value(   68)+d13bra(  8)*d11ket(  5)*(xin(   1)*yin(   8)*zin(  25)+xin(  33)*yin(  40)*zin(  57)+xin(  65)*yin(  72)*zin(  89)+xin(  97)*yin( 104)*zin( 121))
          eri_value(   69)=eri_value(   69)+d13bra(  8)*d11ket(  6)*(xin(   1)*yin(   7)*zin(  26)+xin(  33)*yin(  39)*zin(  58)+xin(  65)*yin(  71)*zin(  90)+xin(  97)*yin( 103)*zin( 122))
          eri_value(   70)=eri_value(   70)+d13bra(  8)*d11ket(  7)*(xin(   2)*yin(   5)*zin(  27)+xin(  34)*yin(  37)*zin(  59)+xin(  66)*yin(  69)*zin(  91)+xin(  98)*yin( 101)*zin( 123))
          eri_value(   71)=eri_value(   71)+d13bra(  8)*d11ket(  8)*(xin(   1)*yin(   6)*zin(  27)+xin(  33)*yin(  38)*zin(  59)+xin(  65)*yin(  70)*zin(  91)+xin(  97)*yin( 102)*zin( 123))
          eri_value(   72)=eri_value(   72)+d13bra(  8)*d11ket(  9)*(xin(   1)*yin(   5)*zin(  28)+xin(  33)*yin(  37)*zin(  60)+xin(  65)*yin(  69)*zin(  92)+xin(  97)*yin( 101)*zin( 124))
          eri_value(   73)=eri_value(   73)+d13bra(  9)*d11ket(  1)*(xin(   4)*yin(   1)*zin(  29)+xin(  36)*yin(  33)*zin(  61)+xin(  68)*yin(  65)*zin(  93)+xin( 100)*yin(  97)*zin( 125))
          eri_value(   74)=eri_value(   74)+d13bra(  9)*d11ket(  2)*(xin(   3)*yin(   2)*zin(  29)+xin(  35)*yin(  34)*zin(  61)+xin(  67)*yin(  66)*zin(  93)+xin(  99)*yin(  98)*zin( 125))
          eri_value(   75)=eri_value(   75)+d13bra(  9)*d11ket(  3)*(xin(   3)*yin(   1)*zin(  30)+xin(  35)*yin(  33)*zin(  62)+xin(  67)*yin(  65)*zin(  94)+xin(  99)*yin(  97)*zin( 126))
          eri_value(   76)=eri_value(   76)+d13bra(  9)*d11ket(  4)*(xin(   2)*yin(   3)*zin(  29)+xin(  34)*yin(  35)*zin(  61)+xin(  66)*yin(  67)*zin(  93)+xin(  98)*yin(  99)*zin( 125))
          eri_value(   77)=eri_value(   77)+d13bra(  9)*d11ket(  5)*(xin(   1)*yin(   4)*zin(  29)+xin(  33)*yin(  36)*zin(  61)+xin(  65)*yin(  68)*zin(  93)+xin(  97)*yin( 100)*zin( 125))
          eri_value(   78)=eri_value(   78)+d13bra(  9)*d11ket(  6)*(xin(   1)*yin(   3)*zin(  30)+xin(  33)*yin(  35)*zin(  62)+xin(  65)*yin(  67)*zin(  94)+xin(  97)*yin(  99)*zin( 126))
          eri_value(   79)=eri_value(   79)+d13bra(  9)*d11ket(  7)*(xin(   2)*yin(   1)*zin(  31)+xin(  34)*yin(  33)*zin(  63)+xin(  66)*yin(  65)*zin(  95)+xin(  98)*yin(  97)*zin( 127))
          eri_value(   80)=eri_value(   80)+d13bra(  9)*d11ket(  8)*(xin(   1)*yin(   2)*zin(  31)+xin(  33)*yin(  34)*zin(  63)+xin(  65)*yin(  66)*zin(  95)+xin(  97)*yin(  98)*zin( 127))
          eri_value(   81)=eri_value(   81)+d13bra(  9)*d11ket(  9)*(xin(   1)*yin(   1)*zin(  32)+xin(  33)*yin(  33)*zin(  64)+xin(  65)*yin(  65)*zin(  96)+xin(  97)*yin(  97)*zin( 128))
          eri_value(   82)=eri_value(   82)+d13bra( 10)*d11ket(  1)*(xin(  24)*yin(   9)*zin(   1)+xin(  56)*yin(  41)*zin(  33)+xin(  88)*yin(  73)*zin(  65)+xin( 120)*yin( 105)*zin(  97))
          eri_value(   83)=eri_value(   83)+d13bra( 10)*d11ket(  2)*(xin(  23)*yin(  10)*zin(   1)+xin(  55)*yin(  42)*zin(  33)+xin(  87)*yin(  74)*zin(  65)+xin( 119)*yin( 106)*zin(  97))
          eri_value(   84)=eri_value(   84)+d13bra( 10)*d11ket(  3)*(xin(  23)*yin(   9)*zin(   2)+xin(  55)*yin(  41)*zin(  34)+xin(  87)*yin(  73)*zin(  66)+xin( 119)*yin( 105)*zin(  98))
          eri_value(   85)=eri_value(   85)+d13bra( 10)*d11ket(  4)*(xin(  22)*yin(  11)*zin(   1)+xin(  54)*yin(  43)*zin(  33)+xin(  86)*yin(  75)*zin(  65)+xin( 118)*yin( 107)*zin(  97))
          eri_value(   86)=eri_value(   86)+d13bra( 10)*d11ket(  5)*(xin(  21)*yin(  12)*zin(   1)+xin(  53)*yin(  44)*zin(  33)+xin(  85)*yin(  76)*zin(  65)+xin( 117)*yin( 108)*zin(  97))
          eri_value(   87)=eri_value(   87)+d13bra( 10)*d11ket(  6)*(xin(  21)*yin(  11)*zin(   2)+xin(  53)*yin(  43)*zin(  34)+xin(  85)*yin(  75)*zin(  66)+xin( 117)*yin( 107)*zin(  98))
          eri_value(   88)=eri_value(   88)+d13bra( 10)*d11ket(  7)*(xin(  22)*yin(   9)*zin(   3)+xin(  54)*yin(  41)*zin(  35)+xin(  86)*yin(  73)*zin(  67)+xin( 118)*yin( 105)*zin(  99))
          eri_value(   89)=eri_value(   89)+d13bra( 10)*d11ket(  8)*(xin(  21)*yin(  10)*zin(   3)+xin(  53)*yin(  42)*zin(  35)+xin(  85)*yin(  74)*zin(  67)+xin( 117)*yin( 106)*zin(  99))
          eri_value(   90)=eri_value(   90)+d13bra( 10)*d11ket(  9)*(xin(  21)*yin(   9)*zin(   4)+xin(  53)*yin(  41)*zin(  36)+xin(  85)*yin(  73)*zin(  68)+xin( 117)*yin( 105)*zin( 100))
          eri_value(   91)=eri_value(   91)+d13bra( 11)*d11ket(  1)*(xin(  20)*yin(  13)*zin(   1)+xin(  52)*yin(  45)*zin(  33)+xin(  84)*yin(  77)*zin(  65)+xin( 116)*yin( 109)*zin(  97))
          eri_value(   92)=eri_value(   92)+d13bra( 11)*d11ket(  2)*(xin(  19)*yin(  14)*zin(   1)+xin(  51)*yin(  46)*zin(  33)+xin(  83)*yin(  78)*zin(  65)+xin( 115)*yin( 110)*zin(  97))
          eri_value(   93)=eri_value(   93)+d13bra( 11)*d11ket(  3)*(xin(  19)*yin(  13)*zin(   2)+xin(  51)*yin(  45)*zin(  34)+xin(  83)*yin(  77)*zin(  66)+xin( 115)*yin( 109)*zin(  98))
          eri_value(   94)=eri_value(   94)+d13bra( 11)*d11ket(  4)*(xin(  18)*yin(  15)*zin(   1)+xin(  50)*yin(  47)*zin(  33)+xin(  82)*yin(  79)*zin(  65)+xin( 114)*yin( 111)*zin(  97))
          eri_value(   95)=eri_value(   95)+d13bra( 11)*d11ket(  5)*(xin(  17)*yin(  16)*zin(   1)+xin(  49)*yin(  48)*zin(  33)+xin(  81)*yin(  80)*zin(  65)+xin( 113)*yin( 112)*zin(  97))
          eri_value(   96)=eri_value(   96)+d13bra( 11)*d11ket(  6)*(xin(  17)*yin(  15)*zin(   2)+xin(  49)*yin(  47)*zin(  34)+xin(  81)*yin(  79)*zin(  66)+xin( 113)*yin( 111)*zin(  98))
          eri_value(   97)=eri_value(   97)+d13bra( 11)*d11ket(  7)*(xin(  18)*yin(  13)*zin(   3)+xin(  50)*yin(  45)*zin(  35)+xin(  82)*yin(  77)*zin(  67)+xin( 114)*yin( 109)*zin(  99))
          eri_value(   98)=eri_value(   98)+d13bra( 11)*d11ket(  8)*(xin(  17)*yin(  14)*zin(   3)+xin(  49)*yin(  46)*zin(  35)+xin(  81)*yin(  78)*zin(  67)+xin( 113)*yin( 110)*zin(  99))
          eri_value(   99)=eri_value(   99)+d13bra( 11)*d11ket(  9)*(xin(  17)*yin(  13)*zin(   4)+xin(  49)*yin(  45)*zin(  36)+xin(  81)*yin(  77)*zin(  68)+xin( 113)*yin( 109)*zin( 100))
          eri_value(  100)=eri_value(  100)+d13bra( 12)*d11ket(  1)*(xin(  20)*yin(   9)*zin(   5)+xin(  52)*yin(  41)*zin(  37)+xin(  84)*yin(  73)*zin(  69)+xin( 116)*yin( 105)*zin( 101))
          eri_value(  101)=eri_value(  101)+d13bra( 12)*d11ket(  2)*(xin(  19)*yin(  10)*zin(   5)+xin(  51)*yin(  42)*zin(  37)+xin(  83)*yin(  74)*zin(  69)+xin( 115)*yin( 106)*zin( 101))
          eri_value(  102)=eri_value(  102)+d13bra( 12)*d11ket(  3)*(xin(  19)*yin(   9)*zin(   6)+xin(  51)*yin(  41)*zin(  38)+xin(  83)*yin(  73)*zin(  70)+xin( 115)*yin( 105)*zin( 102))
          eri_value(  103)=eri_value(  103)+d13bra( 12)*d11ket(  4)*(xin(  18)*yin(  11)*zin(   5)+xin(  50)*yin(  43)*zin(  37)+xin(  82)*yin(  75)*zin(  69)+xin( 114)*yin( 107)*zin( 101))
          eri_value(  104)=eri_value(  104)+d13bra( 12)*d11ket(  5)*(xin(  17)*yin(  12)*zin(   5)+xin(  49)*yin(  44)*zin(  37)+xin(  81)*yin(  76)*zin(  69)+xin( 113)*yin( 108)*zin( 101))
          eri_value(  105)=eri_value(  105)+d13bra( 12)*d11ket(  6)*(xin(  17)*yin(  11)*zin(   6)+xin(  49)*yin(  43)*zin(  38)+xin(  81)*yin(  75)*zin(  70)+xin( 113)*yin( 107)*zin( 102))
          eri_value(  106)=eri_value(  106)+d13bra( 12)*d11ket(  7)*(xin(  18)*yin(   9)*zin(   7)+xin(  50)*yin(  41)*zin(  39)+xin(  82)*yin(  73)*zin(  71)+xin( 114)*yin( 105)*zin( 103))
          eri_value(  107)=eri_value(  107)+d13bra( 12)*d11ket(  8)*(xin(  17)*yin(  10)*zin(   7)+xin(  49)*yin(  42)*zin(  39)+xin(  81)*yin(  74)*zin(  71)+xin( 113)*yin( 106)*zin( 103))
          eri_value(  108)=eri_value(  108)+d13bra( 12)*d11ket(  9)*(xin(  17)*yin(   9)*zin(   8)+xin(  49)*yin(  41)*zin(  40)+xin(  81)*yin(  73)*zin(  72)+xin( 113)*yin( 105)*zin( 104))
          eri_value(  109)=eri_value(  109)+d13bra( 13)*d11ket(  1)*(xin(  24)*yin(   1)*zin(   9)+xin(  56)*yin(  33)*zin(  41)+xin(  88)*yin(  65)*zin(  73)+xin( 120)*yin(  97)*zin( 105))
          eri_value(  110)=eri_value(  110)+d13bra( 13)*d11ket(  2)*(xin(  23)*yin(   2)*zin(   9)+xin(  55)*yin(  34)*zin(  41)+xin(  87)*yin(  66)*zin(  73)+xin( 119)*yin(  98)*zin( 105))
          eri_value(  111)=eri_value(  111)+d13bra( 13)*d11ket(  3)*(xin(  23)*yin(   1)*zin(  10)+xin(  55)*yin(  33)*zin(  42)+xin(  87)*yin(  65)*zin(  74)+xin( 119)*yin(  97)*zin( 106))
          eri_value(  112)=eri_value(  112)+d13bra( 13)*d11ket(  4)*(xin(  22)*yin(   3)*zin(   9)+xin(  54)*yin(  35)*zin(  41)+xin(  86)*yin(  67)*zin(  73)+xin( 118)*yin(  99)*zin( 105))
          eri_value(  113)=eri_value(  113)+d13bra( 13)*d11ket(  5)*(xin(  21)*yin(   4)*zin(   9)+xin(  53)*yin(  36)*zin(  41)+xin(  85)*yin(  68)*zin(  73)+xin( 117)*yin( 100)*zin( 105))
          eri_value(  114)=eri_value(  114)+d13bra( 13)*d11ket(  6)*(xin(  21)*yin(   3)*zin(  10)+xin(  53)*yin(  35)*zin(  42)+xin(  85)*yin(  67)*zin(  74)+xin( 117)*yin(  99)*zin( 106))
          eri_value(  115)=eri_value(  115)+d13bra( 13)*d11ket(  7)*(xin(  22)*yin(   1)*zin(  11)+xin(  54)*yin(  33)*zin(  43)+xin(  86)*yin(  65)*zin(  75)+xin( 118)*yin(  97)*zin( 107))
          eri_value(  116)=eri_value(  116)+d13bra( 13)*d11ket(  8)*(xin(  21)*yin(   2)*zin(  11)+xin(  53)*yin(  34)*zin(  43)+xin(  85)*yin(  66)*zin(  75)+xin( 117)*yin(  98)*zin( 107))
          eri_value(  117)=eri_value(  117)+d13bra( 13)*d11ket(  9)*(xin(  21)*yin(   1)*zin(  12)+xin(  53)*yin(  33)*zin(  44)+xin(  85)*yin(  65)*zin(  76)+xin( 117)*yin(  97)*zin( 108))
          eri_value(  118)=eri_value(  118)+d13bra( 14)*d11ket(  1)*(xin(  20)*yin(   5)*zin(   9)+xin(  52)*yin(  37)*zin(  41)+xin(  84)*yin(  69)*zin(  73)+xin( 116)*yin( 101)*zin( 105))
          eri_value(  119)=eri_value(  119)+d13bra( 14)*d11ket(  2)*(xin(  19)*yin(   6)*zin(   9)+xin(  51)*yin(  38)*zin(  41)+xin(  83)*yin(  70)*zin(  73)+xin( 115)*yin( 102)*zin( 105))
          eri_value(  120)=eri_value(  120)+d13bra( 14)*d11ket(  3)*(xin(  19)*yin(   5)*zin(  10)+xin(  51)*yin(  37)*zin(  42)+xin(  83)*yin(  69)*zin(  74)+xin( 115)*yin( 101)*zin( 106))
          eri_value(  121)=eri_value(  121)+d13bra( 14)*d11ket(  4)*(xin(  18)*yin(   7)*zin(   9)+xin(  50)*yin(  39)*zin(  41)+xin(  82)*yin(  71)*zin(  73)+xin( 114)*yin( 103)*zin( 105))
          eri_value(  122)=eri_value(  122)+d13bra( 14)*d11ket(  5)*(xin(  17)*yin(   8)*zin(   9)+xin(  49)*yin(  40)*zin(  41)+xin(  81)*yin(  72)*zin(  73)+xin( 113)*yin( 104)*zin( 105))
          eri_value(  123)=eri_value(  123)+d13bra( 14)*d11ket(  6)*(xin(  17)*yin(   7)*zin(  10)+xin(  49)*yin(  39)*zin(  42)+xin(  81)*yin(  71)*zin(  74)+xin( 113)*yin( 103)*zin( 106))
          eri_value(  124)=eri_value(  124)+d13bra( 14)*d11ket(  7)*(xin(  18)*yin(   5)*zin(  11)+xin(  50)*yin(  37)*zin(  43)+xin(  82)*yin(  69)*zin(  75)+xin( 114)*yin( 101)*zin( 107))
          eri_value(  125)=eri_value(  125)+d13bra( 14)*d11ket(  8)*(xin(  17)*yin(   6)*zin(  11)+xin(  49)*yin(  38)*zin(  43)+xin(  81)*yin(  70)*zin(  75)+xin( 113)*yin( 102)*zin( 107))
          eri_value(  126)=eri_value(  126)+d13bra( 14)*d11ket(  9)*(xin(  17)*yin(   5)*zin(  12)+xin(  49)*yin(  37)*zin(  44)+xin(  81)*yin(  69)*zin(  76)+xin( 113)*yin( 101)*zin( 108))
          eri_value(  127)=eri_value(  127)+d13bra( 15)*d11ket(  1)*(xin(  20)*yin(   1)*zin(  13)+xin(  52)*yin(  33)*zin(  45)+xin(  84)*yin(  65)*zin(  77)+xin( 116)*yin(  97)*zin( 109))
          eri_value(  128)=eri_value(  128)+d13bra( 15)*d11ket(  2)*(xin(  19)*yin(   2)*zin(  13)+xin(  51)*yin(  34)*zin(  45)+xin(  83)*yin(  66)*zin(  77)+xin( 115)*yin(  98)*zin( 109))
          eri_value(  129)=eri_value(  129)+d13bra( 15)*d11ket(  3)*(xin(  19)*yin(   1)*zin(  14)+xin(  51)*yin(  33)*zin(  46)+xin(  83)*yin(  65)*zin(  78)+xin( 115)*yin(  97)*zin( 110))
          eri_value(  130)=eri_value(  130)+d13bra( 15)*d11ket(  4)*(xin(  18)*yin(   3)*zin(  13)+xin(  50)*yin(  35)*zin(  45)+xin(  82)*yin(  67)*zin(  77)+xin( 114)*yin(  99)*zin( 109))
          eri_value(  131)=eri_value(  131)+d13bra( 15)*d11ket(  5)*(xin(  17)*yin(   4)*zin(  13)+xin(  49)*yin(  36)*zin(  45)+xin(  81)*yin(  68)*zin(  77)+xin( 113)*yin( 100)*zin( 109))
          eri_value(  132)=eri_value(  132)+d13bra( 15)*d11ket(  6)*(xin(  17)*yin(   3)*zin(  14)+xin(  49)*yin(  35)*zin(  46)+xin(  81)*yin(  67)*zin(  78)+xin( 113)*yin(  99)*zin( 110))
          eri_value(  133)=eri_value(  133)+d13bra( 15)*d11ket(  7)*(xin(  18)*yin(   1)*zin(  15)+xin(  50)*yin(  33)*zin(  47)+xin(  82)*yin(  65)*zin(  79)+xin( 114)*yin(  97)*zin( 111))
          eri_value(  134)=eri_value(  134)+d13bra( 15)*d11ket(  8)*(xin(  17)*yin(   2)*zin(  15)+xin(  49)*yin(  34)*zin(  47)+xin(  81)*yin(  66)*zin(  79)+xin( 113)*yin(  98)*zin( 111))
          eri_value(  135)=eri_value(  135)+d13bra( 15)*d11ket(  9)*(xin(  17)*yin(   1)*zin(  16)+xin(  49)*yin(  33)*zin(  48)+xin(  81)*yin(  65)*zin(  80)+xin( 113)*yin(  97)*zin( 112))
          eri_value(  136)=eri_value(  136)+d13bra( 16)*d11ket(  1)*(xin(  16)*yin(  17)*zin(   1)+xin(  48)*yin(  49)*zin(  33)+xin(  80)*yin(  81)*zin(  65)+xin( 112)*yin( 113)*zin(  97))
          eri_value(  137)=eri_value(  137)+d13bra( 16)*d11ket(  2)*(xin(  15)*yin(  18)*zin(   1)+xin(  47)*yin(  50)*zin(  33)+xin(  79)*yin(  82)*zin(  65)+xin( 111)*yin( 114)*zin(  97))
          eri_value(  138)=eri_value(  138)+d13bra( 16)*d11ket(  3)*(xin(  15)*yin(  17)*zin(   2)+xin(  47)*yin(  49)*zin(  34)+xin(  79)*yin(  81)*zin(  66)+xin( 111)*yin( 113)*zin(  98))
          eri_value(  139)=eri_value(  139)+d13bra( 16)*d11ket(  4)*(xin(  14)*yin(  19)*zin(   1)+xin(  46)*yin(  51)*zin(  33)+xin(  78)*yin(  83)*zin(  65)+xin( 110)*yin( 115)*zin(  97))
          eri_value(  140)=eri_value(  140)+d13bra( 16)*d11ket(  5)*(xin(  13)*yin(  20)*zin(   1)+xin(  45)*yin(  52)*zin(  33)+xin(  77)*yin(  84)*zin(  65)+xin( 109)*yin( 116)*zin(  97))
          eri_value(  141)=eri_value(  141)+d13bra( 16)*d11ket(  6)*(xin(  13)*yin(  19)*zin(   2)+xin(  45)*yin(  51)*zin(  34)+xin(  77)*yin(  83)*zin(  66)+xin( 109)*yin( 115)*zin(  98))
          eri_value(  142)=eri_value(  142)+d13bra( 16)*d11ket(  7)*(xin(  14)*yin(  17)*zin(   3)+xin(  46)*yin(  49)*zin(  35)+xin(  78)*yin(  81)*zin(  67)+xin( 110)*yin( 113)*zin(  99))
          eri_value(  143)=eri_value(  143)+d13bra( 16)*d11ket(  8)*(xin(  13)*yin(  18)*zin(   3)+xin(  45)*yin(  50)*zin(  35)+xin(  77)*yin(  82)*zin(  67)+xin( 109)*yin( 114)*zin(  99))
          eri_value(  144)=eri_value(  144)+d13bra( 16)*d11ket(  9)*(xin(  13)*yin(  17)*zin(   4)+xin(  45)*yin(  49)*zin(  36)+xin(  77)*yin(  81)*zin(  68)+xin( 109)*yin( 113)*zin( 100))
          eri_value(  145)=eri_value(  145)+d13bra( 17)*d11ket(  1)*(xin(  12)*yin(  21)*zin(   1)+xin(  44)*yin(  53)*zin(  33)+xin(  76)*yin(  85)*zin(  65)+xin( 108)*yin( 117)*zin(  97))
          eri_value(  146)=eri_value(  146)+d13bra( 17)*d11ket(  2)*(xin(  11)*yin(  22)*zin(   1)+xin(  43)*yin(  54)*zin(  33)+xin(  75)*yin(  86)*zin(  65)+xin( 107)*yin( 118)*zin(  97))
          eri_value(  147)=eri_value(  147)+d13bra( 17)*d11ket(  3)*(xin(  11)*yin(  21)*zin(   2)+xin(  43)*yin(  53)*zin(  34)+xin(  75)*yin(  85)*zin(  66)+xin( 107)*yin( 117)*zin(  98))
          eri_value(  148)=eri_value(  148)+d13bra( 17)*d11ket(  4)*(xin(  10)*yin(  23)*zin(   1)+xin(  42)*yin(  55)*zin(  33)+xin(  74)*yin(  87)*zin(  65)+xin( 106)*yin( 119)*zin(  97))
          eri_value(  149)=eri_value(  149)+d13bra( 17)*d11ket(  5)*(xin(   9)*yin(  24)*zin(   1)+xin(  41)*yin(  56)*zin(  33)+xin(  73)*yin(  88)*zin(  65)+xin( 105)*yin( 120)*zin(  97))
          eri_value(  150)=eri_value(  150)+d13bra( 17)*d11ket(  6)*(xin(   9)*yin(  23)*zin(   2)+xin(  41)*yin(  55)*zin(  34)+xin(  73)*yin(  87)*zin(  66)+xin( 105)*yin( 119)*zin(  98))
          eri_value(  151)=eri_value(  151)+d13bra( 17)*d11ket(  7)*(xin(  10)*yin(  21)*zin(   3)+xin(  42)*yin(  53)*zin(  35)+xin(  74)*yin(  85)*zin(  67)+xin( 106)*yin( 117)*zin(  99))
          eri_value(  152)=eri_value(  152)+d13bra( 17)*d11ket(  8)*(xin(   9)*yin(  22)*zin(   3)+xin(  41)*yin(  54)*zin(  35)+xin(  73)*yin(  86)*zin(  67)+xin( 105)*yin( 118)*zin(  99))
          eri_value(  153)=eri_value(  153)+d13bra( 17)*d11ket(  9)*(xin(   9)*yin(  21)*zin(   4)+xin(  41)*yin(  53)*zin(  36)+xin(  73)*yin(  85)*zin(  68)+xin( 105)*yin( 117)*zin( 100))
          eri_value(  154)=eri_value(  154)+d13bra( 18)*d11ket(  1)*(xin(  12)*yin(  17)*zin(   5)+xin(  44)*yin(  49)*zin(  37)+xin(  76)*yin(  81)*zin(  69)+xin( 108)*yin( 113)*zin( 101))
          eri_value(  155)=eri_value(  155)+d13bra( 18)*d11ket(  2)*(xin(  11)*yin(  18)*zin(   5)+xin(  43)*yin(  50)*zin(  37)+xin(  75)*yin(  82)*zin(  69)+xin( 107)*yin( 114)*zin( 101))
          eri_value(  156)=eri_value(  156)+d13bra( 18)*d11ket(  3)*(xin(  11)*yin(  17)*zin(   6)+xin(  43)*yin(  49)*zin(  38)+xin(  75)*yin(  81)*zin(  70)+xin( 107)*yin( 113)*zin( 102))
          eri_value(  157)=eri_value(  157)+d13bra( 18)*d11ket(  4)*(xin(  10)*yin(  19)*zin(   5)+xin(  42)*yin(  51)*zin(  37)+xin(  74)*yin(  83)*zin(  69)+xin( 106)*yin( 115)*zin( 101))
          eri_value(  158)=eri_value(  158)+d13bra( 18)*d11ket(  5)*(xin(   9)*yin(  20)*zin(   5)+xin(  41)*yin(  52)*zin(  37)+xin(  73)*yin(  84)*zin(  69)+xin( 105)*yin( 116)*zin( 101))
          eri_value(  159)=eri_value(  159)+d13bra( 18)*d11ket(  6)*(xin(   9)*yin(  19)*zin(   6)+xin(  41)*yin(  51)*zin(  38)+xin(  73)*yin(  83)*zin(  70)+xin( 105)*yin( 115)*zin( 102))
          eri_value(  160)=eri_value(  160)+d13bra( 18)*d11ket(  7)*(xin(  10)*yin(  17)*zin(   7)+xin(  42)*yin(  49)*zin(  39)+xin(  74)*yin(  81)*zin(  71)+xin( 106)*yin( 113)*zin( 103))
          eri_value(  161)=eri_value(  161)+d13bra( 18)*d11ket(  8)*(xin(   9)*yin(  18)*zin(   7)+xin(  41)*yin(  50)*zin(  39)+xin(  73)*yin(  82)*zin(  71)+xin( 105)*yin( 114)*zin( 103))
          eri_value(  162)=eri_value(  162)+d13bra( 18)*d11ket(  9)*(xin(   9)*yin(  17)*zin(   8)+xin(  41)*yin(  49)*zin(  40)+xin(  73)*yin(  81)*zin(  72)+xin( 105)*yin( 113)*zin( 104))
          eri_value(  163)=eri_value(  163)+d13bra( 19)*d11ket(  1)*(xin(   8)*yin(  17)*zin(   9)+xin(  40)*yin(  49)*zin(  41)+xin(  72)*yin(  81)*zin(  73)+xin( 104)*yin( 113)*zin( 105))
          eri_value(  164)=eri_value(  164)+d13bra( 19)*d11ket(  2)*(xin(   7)*yin(  18)*zin(   9)+xin(  39)*yin(  50)*zin(  41)+xin(  71)*yin(  82)*zin(  73)+xin( 103)*yin( 114)*zin( 105))
          eri_value(  165)=eri_value(  165)+d13bra( 19)*d11ket(  3)*(xin(   7)*yin(  17)*zin(  10)+xin(  39)*yin(  49)*zin(  42)+xin(  71)*yin(  81)*zin(  74)+xin( 103)*yin( 113)*zin( 106))
          eri_value(  166)=eri_value(  166)+d13bra( 19)*d11ket(  4)*(xin(   6)*yin(  19)*zin(   9)+xin(  38)*yin(  51)*zin(  41)+xin(  70)*yin(  83)*zin(  73)+xin( 102)*yin( 115)*zin( 105))
          eri_value(  167)=eri_value(  167)+d13bra( 19)*d11ket(  5)*(xin(   5)*yin(  20)*zin(   9)+xin(  37)*yin(  52)*zin(  41)+xin(  69)*yin(  84)*zin(  73)+xin( 101)*yin( 116)*zin( 105))
          eri_value(  168)=eri_value(  168)+d13bra( 19)*d11ket(  6)*(xin(   5)*yin(  19)*zin(  10)+xin(  37)*yin(  51)*zin(  42)+xin(  69)*yin(  83)*zin(  74)+xin( 101)*yin( 115)*zin( 106))
          eri_value(  169)=eri_value(  169)+d13bra( 19)*d11ket(  7)*(xin(   6)*yin(  17)*zin(  11)+xin(  38)*yin(  49)*zin(  43)+xin(  70)*yin(  81)*zin(  75)+xin( 102)*yin( 113)*zin( 107))
          eri_value(  170)=eri_value(  170)+d13bra( 19)*d11ket(  8)*(xin(   5)*yin(  18)*zin(  11)+xin(  37)*yin(  50)*zin(  43)+xin(  69)*yin(  82)*zin(  75)+xin( 101)*yin( 114)*zin( 107))
          eri_value(  171)=eri_value(  171)+d13bra( 19)*d11ket(  9)*(xin(   5)*yin(  17)*zin(  12)+xin(  37)*yin(  49)*zin(  44)+xin(  69)*yin(  81)*zin(  76)+xin( 101)*yin( 113)*zin( 108))
          eri_value(  172)=eri_value(  172)+d13bra( 20)*d11ket(  1)*(xin(   4)*yin(  21)*zin(   9)+xin(  36)*yin(  53)*zin(  41)+xin(  68)*yin(  85)*zin(  73)+xin( 100)*yin( 117)*zin( 105))
          eri_value(  173)=eri_value(  173)+d13bra( 20)*d11ket(  2)*(xin(   3)*yin(  22)*zin(   9)+xin(  35)*yin(  54)*zin(  41)+xin(  67)*yin(  86)*zin(  73)+xin(  99)*yin( 118)*zin( 105))
          eri_value(  174)=eri_value(  174)+d13bra( 20)*d11ket(  3)*(xin(   3)*yin(  21)*zin(  10)+xin(  35)*yin(  53)*zin(  42)+xin(  67)*yin(  85)*zin(  74)+xin(  99)*yin( 117)*zin( 106))
          eri_value(  175)=eri_value(  175)+d13bra( 20)*d11ket(  4)*(xin(   2)*yin(  23)*zin(   9)+xin(  34)*yin(  55)*zin(  41)+xin(  66)*yin(  87)*zin(  73)+xin(  98)*yin( 119)*zin( 105))
          eri_value(  176)=eri_value(  176)+d13bra( 20)*d11ket(  5)*(xin(   1)*yin(  24)*zin(   9)+xin(  33)*yin(  56)*zin(  41)+xin(  65)*yin(  88)*zin(  73)+xin(  97)*yin( 120)*zin( 105))
          eri_value(  177)=eri_value(  177)+d13bra( 20)*d11ket(  6)*(xin(   1)*yin(  23)*zin(  10)+xin(  33)*yin(  55)*zin(  42)+xin(  65)*yin(  87)*zin(  74)+xin(  97)*yin( 119)*zin( 106))
          eri_value(  178)=eri_value(  178)+d13bra( 20)*d11ket(  7)*(xin(   2)*yin(  21)*zin(  11)+xin(  34)*yin(  53)*zin(  43)+xin(  66)*yin(  85)*zin(  75)+xin(  98)*yin( 117)*zin( 107))
          eri_value(  179)=eri_value(  179)+d13bra( 20)*d11ket(  8)*(xin(   1)*yin(  22)*zin(  11)+xin(  33)*yin(  54)*zin(  43)+xin(  65)*yin(  86)*zin(  75)+xin(  97)*yin( 118)*zin( 107))
          eri_value(  180)=eri_value(  180)+d13bra( 20)*d11ket(  9)*(xin(   1)*yin(  21)*zin(  12)+xin(  33)*yin(  53)*zin(  44)+xin(  65)*yin(  85)*zin(  76)+xin(  97)*yin( 117)*zin( 108))
          eri_value(  181)=eri_value(  181)+d13bra( 21)*d11ket(  1)*(xin(   4)*yin(  17)*zin(  13)+xin(  36)*yin(  49)*zin(  45)+xin(  68)*yin(  81)*zin(  77)+xin( 100)*yin( 113)*zin( 109))
          eri_value(  182)=eri_value(  182)+d13bra( 21)*d11ket(  2)*(xin(   3)*yin(  18)*zin(  13)+xin(  35)*yin(  50)*zin(  45)+xin(  67)*yin(  82)*zin(  77)+xin(  99)*yin( 114)*zin( 109))
          eri_value(  183)=eri_value(  183)+d13bra( 21)*d11ket(  3)*(xin(   3)*yin(  17)*zin(  14)+xin(  35)*yin(  49)*zin(  46)+xin(  67)*yin(  81)*zin(  78)+xin(  99)*yin( 113)*zin( 110))
          eri_value(  184)=eri_value(  184)+d13bra( 21)*d11ket(  4)*(xin(   2)*yin(  19)*zin(  13)+xin(  34)*yin(  51)*zin(  45)+xin(  66)*yin(  83)*zin(  77)+xin(  98)*yin( 115)*zin( 109))
          eri_value(  185)=eri_value(  185)+d13bra( 21)*d11ket(  5)*(xin(   1)*yin(  20)*zin(  13)+xin(  33)*yin(  52)*zin(  45)+xin(  65)*yin(  84)*zin(  77)+xin(  97)*yin( 116)*zin( 109))
          eri_value(  186)=eri_value(  186)+d13bra( 21)*d11ket(  6)*(xin(   1)*yin(  19)*zin(  14)+xin(  33)*yin(  51)*zin(  46)+xin(  65)*yin(  83)*zin(  78)+xin(  97)*yin( 115)*zin( 110))
          eri_value(  187)=eri_value(  187)+d13bra( 21)*d11ket(  7)*(xin(   2)*yin(  17)*zin(  15)+xin(  34)*yin(  49)*zin(  47)+xin(  66)*yin(  81)*zin(  79)+xin(  98)*yin( 113)*zin( 111))
          eri_value(  188)=eri_value(  188)+d13bra( 21)*d11ket(  8)*(xin(   1)*yin(  18)*zin(  15)+xin(  33)*yin(  50)*zin(  47)+xin(  65)*yin(  82)*zin(  79)+xin(  97)*yin( 114)*zin( 111))
          eri_value(  189)=eri_value(  189)+d13bra( 21)*d11ket(  9)*(xin(   1)*yin(  17)*zin(  16)+xin(  33)*yin(  49)*zin(  48)+xin(  65)*yin(  81)*zin(  80)+xin(  97)*yin( 113)*zin( 112))
          eri_value(  190)=eri_value(  190)+d13bra( 22)*d11ket(  1)*(xin(  16)*yin(   1)*zin(  17)+xin(  48)*yin(  33)*zin(  49)+xin(  80)*yin(  65)*zin(  81)+xin( 112)*yin(  97)*zin( 113))
          eri_value(  191)=eri_value(  191)+d13bra( 22)*d11ket(  2)*(xin(  15)*yin(   2)*zin(  17)+xin(  47)*yin(  34)*zin(  49)+xin(  79)*yin(  66)*zin(  81)+xin( 111)*yin(  98)*zin( 113))
          eri_value(  192)=eri_value(  192)+d13bra( 22)*d11ket(  3)*(xin(  15)*yin(   1)*zin(  18)+xin(  47)*yin(  33)*zin(  50)+xin(  79)*yin(  65)*zin(  82)+xin( 111)*yin(  97)*zin( 114))
          eri_value(  193)=eri_value(  193)+d13bra( 22)*d11ket(  4)*(xin(  14)*yin(   3)*zin(  17)+xin(  46)*yin(  35)*zin(  49)+xin(  78)*yin(  67)*zin(  81)+xin( 110)*yin(  99)*zin( 113))
          eri_value(  194)=eri_value(  194)+d13bra( 22)*d11ket(  5)*(xin(  13)*yin(   4)*zin(  17)+xin(  45)*yin(  36)*zin(  49)+xin(  77)*yin(  68)*zin(  81)+xin( 109)*yin( 100)*zin( 113))
          eri_value(  195)=eri_value(  195)+d13bra( 22)*d11ket(  6)*(xin(  13)*yin(   3)*zin(  18)+xin(  45)*yin(  35)*zin(  50)+xin(  77)*yin(  67)*zin(  82)+xin( 109)*yin(  99)*zin( 114))
          eri_value(  196)=eri_value(  196)+d13bra( 22)*d11ket(  7)*(xin(  14)*yin(   1)*zin(  19)+xin(  46)*yin(  33)*zin(  51)+xin(  78)*yin(  65)*zin(  83)+xin( 110)*yin(  97)*zin( 115))
          eri_value(  197)=eri_value(  197)+d13bra( 22)*d11ket(  8)*(xin(  13)*yin(   2)*zin(  19)+xin(  45)*yin(  34)*zin(  51)+xin(  77)*yin(  66)*zin(  83)+xin( 109)*yin(  98)*zin( 115))
          eri_value(  198)=eri_value(  198)+d13bra( 22)*d11ket(  9)*(xin(  13)*yin(   1)*zin(  20)+xin(  45)*yin(  33)*zin(  52)+xin(  77)*yin(  65)*zin(  84)+xin( 109)*yin(  97)*zin( 116))
          eri_value(  199)=eri_value(  199)+d13bra( 23)*d11ket(  1)*(xin(  12)*yin(   5)*zin(  17)+xin(  44)*yin(  37)*zin(  49)+xin(  76)*yin(  69)*zin(  81)+xin( 108)*yin( 101)*zin( 113))
          eri_value(  200)=eri_value(  200)+d13bra( 23)*d11ket(  2)*(xin(  11)*yin(   6)*zin(  17)+xin(  43)*yin(  38)*zin(  49)+xin(  75)*yin(  70)*zin(  81)+xin( 107)*yin( 102)*zin( 113))
          eri_value(  201)=eri_value(  201)+d13bra( 23)*d11ket(  3)*(xin(  11)*yin(   5)*zin(  18)+xin(  43)*yin(  37)*zin(  50)+xin(  75)*yin(  69)*zin(  82)+xin( 107)*yin( 101)*zin( 114))
          eri_value(  202)=eri_value(  202)+d13bra( 23)*d11ket(  4)*(xin(  10)*yin(   7)*zin(  17)+xin(  42)*yin(  39)*zin(  49)+xin(  74)*yin(  71)*zin(  81)+xin( 106)*yin( 103)*zin( 113))
          eri_value(  203)=eri_value(  203)+d13bra( 23)*d11ket(  5)*(xin(   9)*yin(   8)*zin(  17)+xin(  41)*yin(  40)*zin(  49)+xin(  73)*yin(  72)*zin(  81)+xin( 105)*yin( 104)*zin( 113))
          eri_value(  204)=eri_value(  204)+d13bra( 23)*d11ket(  6)*(xin(   9)*yin(   7)*zin(  18)+xin(  41)*yin(  39)*zin(  50)+xin(  73)*yin(  71)*zin(  82)+xin( 105)*yin( 103)*zin( 114))
          eri_value(  205)=eri_value(  205)+d13bra( 23)*d11ket(  7)*(xin(  10)*yin(   5)*zin(  19)+xin(  42)*yin(  37)*zin(  51)+xin(  74)*yin(  69)*zin(  83)+xin( 106)*yin( 101)*zin( 115))
          eri_value(  206)=eri_value(  206)+d13bra( 23)*d11ket(  8)*(xin(   9)*yin(   6)*zin(  19)+xin(  41)*yin(  38)*zin(  51)+xin(  73)*yin(  70)*zin(  83)+xin( 105)*yin( 102)*zin( 115))
          eri_value(  207)=eri_value(  207)+d13bra( 23)*d11ket(  9)*(xin(   9)*yin(   5)*zin(  20)+xin(  41)*yin(  37)*zin(  52)+xin(  73)*yin(  69)*zin(  84)+xin( 105)*yin( 101)*zin( 116))
          eri_value(  208)=eri_value(  208)+d13bra( 24)*d11ket(  1)*(xin(  12)*yin(   1)*zin(  21)+xin(  44)*yin(  33)*zin(  53)+xin(  76)*yin(  65)*zin(  85)+xin( 108)*yin(  97)*zin( 117))
          eri_value(  209)=eri_value(  209)+d13bra( 24)*d11ket(  2)*(xin(  11)*yin(   2)*zin(  21)+xin(  43)*yin(  34)*zin(  53)+xin(  75)*yin(  66)*zin(  85)+xin( 107)*yin(  98)*zin( 117))
          eri_value(  210)=eri_value(  210)+d13bra( 24)*d11ket(  3)*(xin(  11)*yin(   1)*zin(  22)+xin(  43)*yin(  33)*zin(  54)+xin(  75)*yin(  65)*zin(  86)+xin( 107)*yin(  97)*zin( 118))
          eri_value(  211)=eri_value(  211)+d13bra( 24)*d11ket(  4)*(xin(  10)*yin(   3)*zin(  21)+xin(  42)*yin(  35)*zin(  53)+xin(  74)*yin(  67)*zin(  85)+xin( 106)*yin(  99)*zin( 117))
          eri_value(  212)=eri_value(  212)+d13bra( 24)*d11ket(  5)*(xin(   9)*yin(   4)*zin(  21)+xin(  41)*yin(  36)*zin(  53)+xin(  73)*yin(  68)*zin(  85)+xin( 105)*yin( 100)*zin( 117))
          eri_value(  213)=eri_value(  213)+d13bra( 24)*d11ket(  6)*(xin(   9)*yin(   3)*zin(  22)+xin(  41)*yin(  35)*zin(  54)+xin(  73)*yin(  67)*zin(  86)+xin( 105)*yin(  99)*zin( 118))
          eri_value(  214)=eri_value(  214)+d13bra( 24)*d11ket(  7)*(xin(  10)*yin(   1)*zin(  23)+xin(  42)*yin(  33)*zin(  55)+xin(  74)*yin(  65)*zin(  87)+xin( 106)*yin(  97)*zin( 119))
          eri_value(  215)=eri_value(  215)+d13bra( 24)*d11ket(  8)*(xin(   9)*yin(   2)*zin(  23)+xin(  41)*yin(  34)*zin(  55)+xin(  73)*yin(  66)*zin(  87)+xin( 105)*yin(  98)*zin( 119))
          eri_value(  216)=eri_value(  216)+d13bra( 24)*d11ket(  9)*(xin(   9)*yin(   1)*zin(  24)+xin(  41)*yin(  33)*zin(  56)+xin(  73)*yin(  65)*zin(  88)+xin( 105)*yin(  97)*zin( 120))
          eri_value(  217)=eri_value(  217)+d13bra( 25)*d11ket(  1)*(xin(   8)*yin(   9)*zin(  17)+xin(  40)*yin(  41)*zin(  49)+xin(  72)*yin(  73)*zin(  81)+xin( 104)*yin( 105)*zin( 113))
          eri_value(  218)=eri_value(  218)+d13bra( 25)*d11ket(  2)*(xin(   7)*yin(  10)*zin(  17)+xin(  39)*yin(  42)*zin(  49)+xin(  71)*yin(  74)*zin(  81)+xin( 103)*yin( 106)*zin( 113))
          eri_value(  219)=eri_value(  219)+d13bra( 25)*d11ket(  3)*(xin(   7)*yin(   9)*zin(  18)+xin(  39)*yin(  41)*zin(  50)+xin(  71)*yin(  73)*zin(  82)+xin( 103)*yin( 105)*zin( 114))
          eri_value(  220)=eri_value(  220)+d13bra( 25)*d11ket(  4)*(xin(   6)*yin(  11)*zin(  17)+xin(  38)*yin(  43)*zin(  49)+xin(  70)*yin(  75)*zin(  81)+xin( 102)*yin( 107)*zin( 113))
          eri_value(  221)=eri_value(  221)+d13bra( 25)*d11ket(  5)*(xin(   5)*yin(  12)*zin(  17)+xin(  37)*yin(  44)*zin(  49)+xin(  69)*yin(  76)*zin(  81)+xin( 101)*yin( 108)*zin( 113))
          eri_value(  222)=eri_value(  222)+d13bra( 25)*d11ket(  6)*(xin(   5)*yin(  11)*zin(  18)+xin(  37)*yin(  43)*zin(  50)+xin(  69)*yin(  75)*zin(  82)+xin( 101)*yin( 107)*zin( 114))
          eri_value(  223)=eri_value(  223)+d13bra( 25)*d11ket(  7)*(xin(   6)*yin(   9)*zin(  19)+xin(  38)*yin(  41)*zin(  51)+xin(  70)*yin(  73)*zin(  83)+xin( 102)*yin( 105)*zin( 115))
          eri_value(  224)=eri_value(  224)+d13bra( 25)*d11ket(  8)*(xin(   5)*yin(  10)*zin(  19)+xin(  37)*yin(  42)*zin(  51)+xin(  69)*yin(  74)*zin(  83)+xin( 101)*yin( 106)*zin( 115))
          eri_value(  225)=eri_value(  225)+d13bra( 25)*d11ket(  9)*(xin(   5)*yin(   9)*zin(  20)+xin(  37)*yin(  41)*zin(  52)+xin(  69)*yin(  73)*zin(  84)+xin( 101)*yin( 105)*zin( 116))
          eri_value(  226)=eri_value(  226)+d13bra( 26)*d11ket(  1)*(xin(   4)*yin(  13)*zin(  17)+xin(  36)*yin(  45)*zin(  49)+xin(  68)*yin(  77)*zin(  81)+xin( 100)*yin( 109)*zin( 113))
          eri_value(  227)=eri_value(  227)+d13bra( 26)*d11ket(  2)*(xin(   3)*yin(  14)*zin(  17)+xin(  35)*yin(  46)*zin(  49)+xin(  67)*yin(  78)*zin(  81)+xin(  99)*yin( 110)*zin( 113))
          eri_value(  228)=eri_value(  228)+d13bra( 26)*d11ket(  3)*(xin(   3)*yin(  13)*zin(  18)+xin(  35)*yin(  45)*zin(  50)+xin(  67)*yin(  77)*zin(  82)+xin(  99)*yin( 109)*zin( 114))
          eri_value(  229)=eri_value(  229)+d13bra( 26)*d11ket(  4)*(xin(   2)*yin(  15)*zin(  17)+xin(  34)*yin(  47)*zin(  49)+xin(  66)*yin(  79)*zin(  81)+xin(  98)*yin( 111)*zin( 113))
          eri_value(  230)=eri_value(  230)+d13bra( 26)*d11ket(  5)*(xin(   1)*yin(  16)*zin(  17)+xin(  33)*yin(  48)*zin(  49)+xin(  65)*yin(  80)*zin(  81)+xin(  97)*yin( 112)*zin( 113))
          eri_value(  231)=eri_value(  231)+d13bra( 26)*d11ket(  6)*(xin(   1)*yin(  15)*zin(  18)+xin(  33)*yin(  47)*zin(  50)+xin(  65)*yin(  79)*zin(  82)+xin(  97)*yin( 111)*zin( 114))
          eri_value(  232)=eri_value(  232)+d13bra( 26)*d11ket(  7)*(xin(   2)*yin(  13)*zin(  19)+xin(  34)*yin(  45)*zin(  51)+xin(  66)*yin(  77)*zin(  83)+xin(  98)*yin( 109)*zin( 115))
          eri_value(  233)=eri_value(  233)+d13bra( 26)*d11ket(  8)*(xin(   1)*yin(  14)*zin(  19)+xin(  33)*yin(  46)*zin(  51)+xin(  65)*yin(  78)*zin(  83)+xin(  97)*yin( 110)*zin( 115))
          eri_value(  234)=eri_value(  234)+d13bra( 26)*d11ket(  9)*(xin(   1)*yin(  13)*zin(  20)+xin(  33)*yin(  45)*zin(  52)+xin(  65)*yin(  77)*zin(  84)+xin(  97)*yin( 109)*zin( 116))
          eri_value(  235)=eri_value(  235)+d13bra( 27)*d11ket(  1)*(xin(   4)*yin(   9)*zin(  21)+xin(  36)*yin(  41)*zin(  53)+xin(  68)*yin(  73)*zin(  85)+xin( 100)*yin( 105)*zin( 117))
          eri_value(  236)=eri_value(  236)+d13bra( 27)*d11ket(  2)*(xin(   3)*yin(  10)*zin(  21)+xin(  35)*yin(  42)*zin(  53)+xin(  67)*yin(  74)*zin(  85)+xin(  99)*yin( 106)*zin( 117))
          eri_value(  237)=eri_value(  237)+d13bra( 27)*d11ket(  3)*(xin(   3)*yin(   9)*zin(  22)+xin(  35)*yin(  41)*zin(  54)+xin(  67)*yin(  73)*zin(  86)+xin(  99)*yin( 105)*zin( 118))
          eri_value(  238)=eri_value(  238)+d13bra( 27)*d11ket(  4)*(xin(   2)*yin(  11)*zin(  21)+xin(  34)*yin(  43)*zin(  53)+xin(  66)*yin(  75)*zin(  85)+xin(  98)*yin( 107)*zin( 117))
          eri_value(  239)=eri_value(  239)+d13bra( 27)*d11ket(  5)*(xin(   1)*yin(  12)*zin(  21)+xin(  33)*yin(  44)*zin(  53)+xin(  65)*yin(  76)*zin(  85)+xin(  97)*yin( 108)*zin( 117))
          eri_value(  240)=eri_value(  240)+d13bra( 27)*d11ket(  6)*(xin(   1)*yin(  11)*zin(  22)+xin(  33)*yin(  43)*zin(  54)+xin(  65)*yin(  75)*zin(  86)+xin(  97)*yin( 107)*zin( 118))
          eri_value(  241)=eri_value(  241)+d13bra( 27)*d11ket(  7)*(xin(   2)*yin(   9)*zin(  23)+xin(  34)*yin(  41)*zin(  55)+xin(  66)*yin(  73)*zin(  87)+xin(  98)*yin( 105)*zin( 119))
          eri_value(  242)=eri_value(  242)+d13bra( 27)*d11ket(  8)*(xin(   1)*yin(  10)*zin(  23)+xin(  33)*yin(  42)*zin(  55)+xin(  65)*yin(  74)*zin(  87)+xin(  97)*yin( 106)*zin( 119))
          eri_value(  243)=eri_value(  243)+d13bra( 27)*d11ket(  9)*(xin(   1)*yin(   9)*zin(  24)+xin(  33)*yin(  41)*zin(  56)+xin(  65)*yin(  73)*zin(  88)+xin(  97)*yin( 105)*zin( 120))
          eri_value(  244)=eri_value(  244)+d13bra( 28)*d11ket(  1)*(xin(  16)*yin(   9)*zin(   9)+xin(  48)*yin(  41)*zin(  41)+xin(  80)*yin(  73)*zin(  73)+xin( 112)*yin( 105)*zin( 105))
          eri_value(  245)=eri_value(  245)+d13bra( 28)*d11ket(  2)*(xin(  15)*yin(  10)*zin(   9)+xin(  47)*yin(  42)*zin(  41)+xin(  79)*yin(  74)*zin(  73)+xin( 111)*yin( 106)*zin( 105))
          eri_value(  246)=eri_value(  246)+d13bra( 28)*d11ket(  3)*(xin(  15)*yin(   9)*zin(  10)+xin(  47)*yin(  41)*zin(  42)+xin(  79)*yin(  73)*zin(  74)+xin( 111)*yin( 105)*zin( 106))
          eri_value(  247)=eri_value(  247)+d13bra( 28)*d11ket(  4)*(xin(  14)*yin(  11)*zin(   9)+xin(  46)*yin(  43)*zin(  41)+xin(  78)*yin(  75)*zin(  73)+xin( 110)*yin( 107)*zin( 105))
          eri_value(  248)=eri_value(  248)+d13bra( 28)*d11ket(  5)*(xin(  13)*yin(  12)*zin(   9)+xin(  45)*yin(  44)*zin(  41)+xin(  77)*yin(  76)*zin(  73)+xin( 109)*yin( 108)*zin( 105))
          eri_value(  249)=eri_value(  249)+d13bra( 28)*d11ket(  6)*(xin(  13)*yin(  11)*zin(  10)+xin(  45)*yin(  43)*zin(  42)+xin(  77)*yin(  75)*zin(  74)+xin( 109)*yin( 107)*zin( 106))
          eri_value(  250)=eri_value(  250)+d13bra( 28)*d11ket(  7)*(xin(  14)*yin(   9)*zin(  11)+xin(  46)*yin(  41)*zin(  43)+xin(  78)*yin(  73)*zin(  75)+xin( 110)*yin( 105)*zin( 107))
          eri_value(  251)=eri_value(  251)+d13bra( 28)*d11ket(  8)*(xin(  13)*yin(  10)*zin(  11)+xin(  45)*yin(  42)*zin(  43)+xin(  77)*yin(  74)*zin(  75)+xin( 109)*yin( 106)*zin( 107))
          eri_value(  252)=eri_value(  252)+d13bra( 28)*d11ket(  9)*(xin(  13)*yin(   9)*zin(  12)+xin(  45)*yin(  41)*zin(  44)+xin(  77)*yin(  73)*zin(  76)+xin( 109)*yin( 105)*zin( 108))
          eri_value(  253)=eri_value(  253)+d13bra( 29)*d11ket(  1)*(xin(  12)*yin(  13)*zin(   9)+xin(  44)*yin(  45)*zin(  41)+xin(  76)*yin(  77)*zin(  73)+xin( 108)*yin( 109)*zin( 105))
          eri_value(  254)=eri_value(  254)+d13bra( 29)*d11ket(  2)*(xin(  11)*yin(  14)*zin(   9)+xin(  43)*yin(  46)*zin(  41)+xin(  75)*yin(  78)*zin(  73)+xin( 107)*yin( 110)*zin( 105))
          eri_value(  255)=eri_value(  255)+d13bra( 29)*d11ket(  3)*(xin(  11)*yin(  13)*zin(  10)+xin(  43)*yin(  45)*zin(  42)+xin(  75)*yin(  77)*zin(  74)+xin( 107)*yin( 109)*zin( 106))
          eri_value(  256)=eri_value(  256)+d13bra( 29)*d11ket(  4)*(xin(  10)*yin(  15)*zin(   9)+xin(  42)*yin(  47)*zin(  41)+xin(  74)*yin(  79)*zin(  73)+xin( 106)*yin( 111)*zin( 105))
          eri_value(  257)=eri_value(  257)+d13bra( 29)*d11ket(  5)*(xin(   9)*yin(  16)*zin(   9)+xin(  41)*yin(  48)*zin(  41)+xin(  73)*yin(  80)*zin(  73)+xin( 105)*yin( 112)*zin( 105))
          eri_value(  258)=eri_value(  258)+d13bra( 29)*d11ket(  6)*(xin(   9)*yin(  15)*zin(  10)+xin(  41)*yin(  47)*zin(  42)+xin(  73)*yin(  79)*zin(  74)+xin( 105)*yin( 111)*zin( 106))
          eri_value(  259)=eri_value(  259)+d13bra( 29)*d11ket(  7)*(xin(  10)*yin(  13)*zin(  11)+xin(  42)*yin(  45)*zin(  43)+xin(  74)*yin(  77)*zin(  75)+xin( 106)*yin( 109)*zin( 107))
          eri_value(  260)=eri_value(  260)+d13bra( 29)*d11ket(  8)*(xin(   9)*yin(  14)*zin(  11)+xin(  41)*yin(  46)*zin(  43)+xin(  73)*yin(  78)*zin(  75)+xin( 105)*yin( 110)*zin( 107))
          eri_value(  261)=eri_value(  261)+d13bra( 29)*d11ket(  9)*(xin(   9)*yin(  13)*zin(  12)+xin(  41)*yin(  45)*zin(  44)+xin(  73)*yin(  77)*zin(  76)+xin( 105)*yin( 109)*zin( 108))
          eri_value(  262)=eri_value(  262)+d13bra( 30)*d11ket(  1)*(xin(  12)*yin(   9)*zin(  13)+xin(  44)*yin(  41)*zin(  45)+xin(  76)*yin(  73)*zin(  77)+xin( 108)*yin( 105)*zin( 109))
          eri_value(  263)=eri_value(  263)+d13bra( 30)*d11ket(  2)*(xin(  11)*yin(  10)*zin(  13)+xin(  43)*yin(  42)*zin(  45)+xin(  75)*yin(  74)*zin(  77)+xin( 107)*yin( 106)*zin( 109))
          eri_value(  264)=eri_value(  264)+d13bra( 30)*d11ket(  3)*(xin(  11)*yin(   9)*zin(  14)+xin(  43)*yin(  41)*zin(  46)+xin(  75)*yin(  73)*zin(  78)+xin( 107)*yin( 105)*zin( 110))
          eri_value(  265)=eri_value(  265)+d13bra( 30)*d11ket(  4)*(xin(  10)*yin(  11)*zin(  13)+xin(  42)*yin(  43)*zin(  45)+xin(  74)*yin(  75)*zin(  77)+xin( 106)*yin( 107)*zin( 109))
          eri_value(  266)=eri_value(  266)+d13bra( 30)*d11ket(  5)*(xin(   9)*yin(  12)*zin(  13)+xin(  41)*yin(  44)*zin(  45)+xin(  73)*yin(  76)*zin(  77)+xin( 105)*yin( 108)*zin( 109))
          eri_value(  267)=eri_value(  267)+d13bra( 30)*d11ket(  6)*(xin(   9)*yin(  11)*zin(  14)+xin(  41)*yin(  43)*zin(  46)+xin(  73)*yin(  75)*zin(  78)+xin( 105)*yin( 107)*zin( 110))
          eri_value(  268)=eri_value(  268)+d13bra( 30)*d11ket(  7)*(xin(  10)*yin(   9)*zin(  15)+xin(  42)*yin(  41)*zin(  47)+xin(  74)*yin(  73)*zin(  79)+xin( 106)*yin( 105)*zin( 111))
          eri_value(  269)=eri_value(  269)+d13bra( 30)*d11ket(  8)*(xin(   9)*yin(  10)*zin(  15)+xin(  41)*yin(  42)*zin(  47)+xin(  73)*yin(  74)*zin(  79)+xin( 105)*yin( 106)*zin( 111))
          eri_value(  270)=eri_value(  270)+d13bra( 30)*d11ket(  9)*(xin(   9)*yin(   9)*zin(  16)+xin(  41)*yin(  41)*zin(  48)+xin(  73)*yin(  73)*zin(  80)+xin( 105)*yin( 105)*zin( 112))

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxl = 3
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*27 ! Stride between functions in i

                                    do j = 1, 3 ! # of cartesians in j

                                      maxl2 = maxl

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*9 + ip ! Add stride between functions in j

                                      do k = 1, 3 ! # of cartesians in k

                                        if (kandl) maxl2 = k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*3 + ijp ! Add stride between functions in k

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

                              deallocate (n13bra)
                              deallocate (xint13bra)
                              deallocate (n11ket)
                              deallocate (xint11ket)

                              end subroutine int3111
                              end submodule
