! The total angular momentum of this class is:           4
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3100_impl
contains
  module subroutine int3100(pf_pair, ss_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pf_pair, ss_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n13bra(:), n00ket(:)
    real(dp), allocatable :: xint13bra(:), xint00ket(:)
    integer(kind=int64) :: npfbra, nssket
    real(dp) :: scutpfbra, scutssket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
    integer(kind=int64) :: n, i1, i3, i4, i5, nm, nn, km, nj, ni
    real(dp) :: cp10, c10
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
    real(dp) :: xin(24), yin(24), zin(24)
    real(dp) :: eri_value(30)
    real(dp) :: d13bra(30), d00ket(1)
    integer(kind=int64) :: ix(10), jx(3), kx(1), lx(1)
    integer(kind=int64) :: iy(10), jy(3), ky(1), ly(1)
    integer(kind=int64) :: iz(10), jz(3), kz(1), lz(1)
    integer(kind=int64) :: in(5), in1(5), kn(1)
    integer(kind=int64) :: ijx(30), ijy(30), ijz(30)
    integer(kind=int64) :: klx(1), kly(1), klz(1)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 3
    in1(3) = 5
    in1(4) = 7
    in1(5) = 8

    kn(1) = 0

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 0

    jx(1) = 1
    jx(2) = 0
    jx(3) = 0

    ix(1) = 7
    ix(2) = 1
    ix(3) = 1
    ix(4) = 5
    ix(5) = 5
    ix(6) = 3
    ix(7) = 1
    ix(8) = 3
    ix(9) = 1
    ix(10) = 3

    ! y-arrays

    ly(1) = 0

    ky(1) = 0

    jy(1) = 0
    jy(2) = 1
    jy(3) = 0

    iy(1) = 1
    iy(2) = 7
    iy(3) = 1
    iy(4) = 3
    iy(5) = 1
    iy(6) = 5
    iy(7) = 5
    iy(8) = 1
    iy(9) = 3
    iy(10) = 3

    ! z-arrays

    lz(1) = 0

    kz(1) = 0

    jz(1) = 0
    jz(2) = 0
    jz(3) = 1

    iz(1) = 1
    iz(2) = 1
    iz(3) = 7
    iz(4) = 1
    iz(5) = 3
    iz(6) = 1
    iz(7) = 3
    iz(8) = 5
    iz(9) = 5
    iz(10) = 3

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 8
    ijx(2) = 7
    ijx(3) = 7
    ijx(4) = 2
    ijx(5) = 1
    ijx(6) = 1
    ijx(7) = 2
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 6
    ijx(11) = 5
    ijx(12) = 5
    ijx(13) = 6
    ijx(14) = 5
    ijx(15) = 5
    ijx(16) = 4
    ijx(17) = 3
    ijx(18) = 3
    ijx(19) = 2
    ijx(20) = 1
    ijx(21) = 1
    ijx(22) = 4
    ijx(23) = 3
    ijx(24) = 3
    ijx(25) = 2
    ijx(26) = 1
    ijx(27) = 1
    ijx(28) = 4
    ijx(29) = 3
    ijx(30) = 3

    ijy(1) = 1
    ijy(2) = 2
    ijy(3) = 1
    ijy(4) = 7
    ijy(5) = 8
    ijy(6) = 7
    ijy(7) = 1
    ijy(8) = 2
    ijy(9) = 1
    ijy(10) = 3
    ijy(11) = 4
    ijy(12) = 3
    ijy(13) = 1
    ijy(14) = 2
    ijy(15) = 1
    ijy(16) = 5
    ijy(17) = 6
    ijy(18) = 5
    ijy(19) = 5
    ijy(20) = 6
    ijy(21) = 5
    ijy(22) = 1
    ijy(23) = 2
    ijy(24) = 1
    ijy(25) = 3
    ijy(26) = 4
    ijy(27) = 3
    ijy(28) = 3
    ijy(29) = 4
    ijy(30) = 3

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 2
    ijz(4) = 1
    ijz(5) = 1
    ijz(6) = 2
    ijz(7) = 7
    ijz(8) = 7
    ijz(9) = 8
    ijz(10) = 1
    ijz(11) = 1
    ijz(12) = 2
    ijz(13) = 3
    ijz(14) = 3
    ijz(15) = 4
    ijz(16) = 1
    ijz(17) = 1
    ijz(18) = 2
    ijz(19) = 3
    ijz(20) = 3
    ijz(21) = 4
    ijz(22) = 5
    ijz(23) = 5
    ijz(24) = 6
    ijz(25) = 5
    ijz(26) = 5
    ijz(27) = 6
    ijz(28) = 3
    ijz(29) = 3
    ijz(30) = 4

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 0

    kly(1) = 0

    klz(1) = 0

    allocate (n13bra(res%n_p_shl*res%n_f_shl))
    allocate (xint13bra(res%n_p_shl*res%n_f_shl))
    allocate (n00ket(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00ket(res%n_s_shl*(res%n_s_shl + 1)/2))

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

    if ((npfbra*nssket) .le. nchunksize_int64) nchunksize_int64 = npfbra*nssket
    ntile = int(npfbra*nssket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = npfbra*nssket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, npfbra, xint13bra, n13bra, xint00ket, n00ket, pf_pair, ss_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d00ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d13bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,nm,nn,km,nj,ni) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, npfbra) + 1
              kl_tmp = (iquart - 1)/npfbra + 1

              test = xint13bra(ij_tmp)*xint00ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n13bra(ij_tmp)
                kl = n00ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_p_shl(jsh_tmp)
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

                                      ! i2 = in(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(3) = xc00
                                      yin(3) = yc00
                                      zin(3) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =    3

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    5
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(5) = c10*xin(1) + xc00*xin(3)
                                      yin(5) = c10*yin(1) + yc00*yin(3)
                                      zin(5) = c10*zin(1) + zc00*zin(3)

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    5

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    7
                                      ! i3 =    3
                                      ! i4 =    5

                                      xin(7) = c10*xin(3) + xc00*xin(5)
                                      yin(7) = c10*yin(3) + yc00*yin(5)
                                      zin(7) = c10*zin(3) + zc00*zin(5)

                                      ! i3 = i4 =    5
                                      ! i4 = i5 =    7

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =    8
                                      ! i3 =    5
                                      ! i4 =    7

                                      xin(8) = c10*xin(5) + xc00*xin(7)
                                      yin(8) = c10*yin(5) + yc00*yin(7)
                                      zin(8) = c10*zin(5) + zc00*zin(7)

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =    8

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =    8

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =    7

                                      xin(8) = xin(8) + dxij*xin(7)
                                      yin(8) = yin(8) + dyij*yin(7)
                                      zin(8) = zin(8) + dzij*zin(7)

                                      ! i3 = i4 =    7
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    2

                                      ! do nj = 1,    1

                                      ! i4 = i3 =    2

                                      ! do ni = 1,    3

                                      xin(2) = xin(3) + dxij*xin(1)
                                      yin(2) = yin(3) + dyij*yin(1)
                                      zin(2) = zin(3) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    4

                                      ! ni =    2

                                      xin(4) = xin(5) + dxij*xin(3)
                                      yin(4) = yin(5) + dyij*yin(3)
                                      zin(4) = zin(5) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    6

                                      ! ni =    3

                                      xin(6) = xin(7) + dxij*xin(5)
                                      yin(6) = yin(7) + dyij*yin(5)
                                      zin(6) = zin(7) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =    8

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    3

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =    8

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

                                      ! i1 = in(1) =    9

                                      xin(9) = 1.0_dp
                                      yin(9) = 1.0_dp
                                      zin(9) = f00

                                      ! i2 = in(2) =   11
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(11) = xc00
                                      yin(11) = yc00
                                      zin(11) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    9
                                      ! i4 = i2 =   11

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   13
                                      ! i3 =    9
                                      ! i4 =   11

                                      xin(13) = c10*xin(9) + xc00*xin(11)
                                      yin(13) = c10*yin(9) + yc00*yin(11)
                                      zin(13) = c10*zin(9) + zc00*zin(11)

                                      ! i3 = i4 =   11
                                      ! i4 = i5 =   13

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   15
                                      ! i3 =   11
                                      ! i4 =   13

                                      xin(15) = c10*xin(11) + xc00*xin(13)
                                      yin(15) = c10*yin(11) + yc00*yin(13)
                                      zin(15) = c10*zin(11) + zc00*zin(13)

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   15

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   16
                                      ! i3 =   13
                                      ! i4 =   15

                                      xin(16) = c10*xin(13) + xc00*xin(15)
                                      yin(16) = c10*yin(13) + yc00*yin(15)
                                      zin(16) = c10*zin(13) + zc00*zin(15)

                                      ! i3 = i4 =   15
                                      ! i4 = i5 =   16

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   16

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   16

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   15

                                      xin(16) = xin(16) + dxij*xin(15)
                                      yin(16) = yin(16) + dyij*yin(15)
                                      zin(16) = zin(16) + dzij*zin(15)

                                      ! i3 = i4 =   15
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   10

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(11) + dxij*xin(9)
                                      yin(10) = yin(11) + dyij*yin(9)
                                      zin(10) = zin(11) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   12

                                      ! ni =    2

                                      xin(12) = xin(13) + dxij*xin(11)
                                      yin(12) = yin(13) + dyij*yin(11)
                                      zin(12) = zin(13) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   14

                                      ! ni =    3

                                      xin(14) = xin(15) + dxij*xin(13)
                                      yin(14) = yin(15) + dyij*yin(13)
                                      zin(14) = zin(15) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   16

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   16

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

                                      ! i1 = in(1) =   17

                                      xin(17) = 1.0_dp
                                      yin(17) = 1.0_dp
                                      zin(17) = f00

                                      ! i2 = in(2) =   19
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(19) = xc00
                                      yin(19) = yc00
                                      zin(19) = zc00*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   17
                                      ! i4 = i2 =   19

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   21
                                      ! i3 =   17
                                      ! i4 =   19

                                      xin(21) = c10*xin(17) + xc00*xin(19)
                                      yin(21) = c10*yin(17) + yc00*yin(19)
                                      zin(21) = c10*zin(17) + zc00*zin(19)

                                      ! i3 = i4 =   19
                                      ! i4 = i5 =   21

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   23
                                      ! i3 =   19
                                      ! i4 =   21

                                      xin(23) = c10*xin(19) + xc00*xin(21)
                                      yin(23) = c10*yin(19) + yc00*yin(21)
                                      zin(23) = c10*zin(19) + zc00*zin(21)

                                      ! i3 = i4 =   21
                                      ! i4 = i5 =   23

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   24
                                      ! i3 =   21
                                      ! i4 =   23

                                      xin(24) = c10*xin(21) + xc00*xin(23)
                                      yin(24) = c10*yin(21) + yc00*yin(23)
                                      zin(24) = c10*zin(21) + zc00*zin(23)

                                      ! i3 = i4 =   23
                                      ! i4 = i5 =   24

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   24

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   24

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   23

                                      xin(24) = xin(24) + dxij*xin(23)
                                      yin(24) = yin(24) + dyij*yin(23)
                                      zin(24) = zin(24) + dzij*zin(23)

                                      ! i3 = i4 =   23
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   18

                                      ! do nj = 1,    1

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(19) + dxij*xin(17)
                                      yin(18) = yin(19) + dyij*yin(17)
                                      zin(18) = zin(19) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   20

                                      ! ni =    2

                                      xin(20) = xin(21) + dxij*xin(19)
                                      yin(20) = yin(21) + dyij*yin(19)
                                      zin(20) = zin(21) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    3

                                      xin(22) = xin(23) + dxij*xin(21)
                                      yin(22) = yin(23) + dyij*yin(21)
                                      zin(22) = zin(23) + dzij*zin(21)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   24

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

          eri_value(1) = eri_value(1) + d13bra(1)*d00ket(1)*(xin(8)*yin(1)*zin(1) + xin(16)*yin(9)*zin(9) + xin(24)*yin(17)*zin(17))
         eri_value(2) = eri_value(2) + d13bra(2)*d00ket(1)*(xin(7)*yin(2)*zin(1) + xin(15)*yin(10)*zin(9) + xin(23)*yin(18)*zin(17))
         eri_value(3) = eri_value(3) + d13bra(3)*d00ket(1)*(xin(7)*yin(1)*zin(2) + xin(15)*yin(9)*zin(10) + xin(23)*yin(17)*zin(18))
         eri_value(4) = eri_value(4) + d13bra(4)*d00ket(1)*(xin(2)*yin(7)*zin(1) + xin(10)*yin(15)*zin(9) + xin(18)*yin(23)*zin(17))
          eri_value(5) = eri_value(5) + d13bra(5)*d00ket(1)*(xin(1)*yin(8)*zin(1) + xin(9)*yin(16)*zin(9) + xin(17)*yin(24)*zin(17))
         eri_value(6) = eri_value(6) + d13bra(6)*d00ket(1)*(xin(1)*yin(7)*zin(2) + xin(9)*yin(15)*zin(10) + xin(17)*yin(23)*zin(18))
         eri_value(7) = eri_value(7) + d13bra(7)*d00ket(1)*(xin(2)*yin(1)*zin(7) + xin(10)*yin(9)*zin(15) + xin(18)*yin(17)*zin(23))
         eri_value(8) = eri_value(8) + d13bra(8)*d00ket(1)*(xin(1)*yin(2)*zin(7) + xin(9)*yin(10)*zin(15) + xin(17)*yin(18)*zin(23))
          eri_value(9) = eri_value(9) + d13bra(9)*d00ket(1)*(xin(1)*yin(1)*zin(8) + xin(9)*yin(9)*zin(16) + xin(17)*yin(17)*zin(24))
      eri_value(10) = eri_value(10) + d13bra(10)*d00ket(1)*(xin(6)*yin(3)*zin(1) + xin(14)*yin(11)*zin(9) + xin(22)*yin(19)*zin(17))
      eri_value(11) = eri_value(11) + d13bra(11)*d00ket(1)*(xin(5)*yin(4)*zin(1) + xin(13)*yin(12)*zin(9) + xin(21)*yin(20)*zin(17))
     eri_value(12) = eri_value(12) + d13bra(12)*d00ket(1)*(xin(5)*yin(3)*zin(2) + xin(13)*yin(11)*zin(10) + xin(21)*yin(19)*zin(18))
      eri_value(13) = eri_value(13) + d13bra(13)*d00ket(1)*(xin(6)*yin(1)*zin(3) + xin(14)*yin(9)*zin(11) + xin(22)*yin(17)*zin(19))
     eri_value(14) = eri_value(14) + d13bra(14)*d00ket(1)*(xin(5)*yin(2)*zin(3) + xin(13)*yin(10)*zin(11) + xin(21)*yin(18)*zin(19))
      eri_value(15) = eri_value(15) + d13bra(15)*d00ket(1)*(xin(5)*yin(1)*zin(4) + xin(13)*yin(9)*zin(12) + xin(21)*yin(17)*zin(20))
      eri_value(16) = eri_value(16) + d13bra(16)*d00ket(1)*(xin(4)*yin(5)*zin(1) + xin(12)*yin(13)*zin(9) + xin(20)*yin(21)*zin(17))
      eri_value(17) = eri_value(17) + d13bra(17)*d00ket(1)*(xin(3)*yin(6)*zin(1) + xin(11)*yin(14)*zin(9) + xin(19)*yin(22)*zin(17))
     eri_value(18) = eri_value(18) + d13bra(18)*d00ket(1)*(xin(3)*yin(5)*zin(2) + xin(11)*yin(13)*zin(10) + xin(19)*yin(21)*zin(18))
     eri_value(19) = eri_value(19) + d13bra(19)*d00ket(1)*(xin(2)*yin(5)*zin(3) + xin(10)*yin(13)*zin(11) + xin(18)*yin(21)*zin(19))
      eri_value(20) = eri_value(20) + d13bra(20)*d00ket(1)*(xin(1)*yin(6)*zin(3) + xin(9)*yin(14)*zin(11) + xin(17)*yin(22)*zin(19))
      eri_value(21) = eri_value(21) + d13bra(21)*d00ket(1)*(xin(1)*yin(5)*zin(4) + xin(9)*yin(13)*zin(12) + xin(17)*yin(21)*zin(20))
      eri_value(22) = eri_value(22) + d13bra(22)*d00ket(1)*(xin(4)*yin(1)*zin(5) + xin(12)*yin(9)*zin(13) + xin(20)*yin(17)*zin(21))
     eri_value(23) = eri_value(23) + d13bra(23)*d00ket(1)*(xin(3)*yin(2)*zin(5) + xin(11)*yin(10)*zin(13) + xin(19)*yin(18)*zin(21))
      eri_value(24) = eri_value(24) + d13bra(24)*d00ket(1)*(xin(3)*yin(1)*zin(6) + xin(11)*yin(9)*zin(14) + xin(19)*yin(17)*zin(22))
     eri_value(25) = eri_value(25) + d13bra(25)*d00ket(1)*(xin(2)*yin(3)*zin(5) + xin(10)*yin(11)*zin(13) + xin(18)*yin(19)*zin(21))
      eri_value(26) = eri_value(26) + d13bra(26)*d00ket(1)*(xin(1)*yin(4)*zin(5) + xin(9)*yin(12)*zin(13) + xin(17)*yin(20)*zin(21))
      eri_value(27) = eri_value(27) + d13bra(27)*d00ket(1)*(xin(1)*yin(3)*zin(6) + xin(9)*yin(11)*zin(14) + xin(17)*yin(19)*zin(22))
     eri_value(28) = eri_value(28) + d13bra(28)*d00ket(1)*(xin(4)*yin(3)*zin(3) + xin(12)*yin(11)*zin(11) + xin(20)*yin(19)*zin(19))
     eri_value(29) = eri_value(29) + d13bra(29)*d00ket(1)*(xin(3)*yin(4)*zin(3) + xin(11)*yin(12)*zin(11) + xin(19)*yin(20)*zin(19))
     eri_value(30) = eri_value(30) + d13bra(30)*d00ket(1)*(xin(3)*yin(3)*zin(4) + xin(11)*yin(11)*zin(12) + xin(19)*yin(19)*zin(20))

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
                                    ip = (i - 1)*3 ! Stride between functions in i

                                    do j = 1, 3 ! # of cartesians in j

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

                              deallocate (n13bra)
                              deallocate (xint13bra)
                              deallocate (n00ket)
                              deallocate (xint00ket)

                              end subroutine int3100
                              end submodule
