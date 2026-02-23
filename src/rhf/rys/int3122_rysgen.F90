! The total angular momentum of this class is:           8
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int2231_impl
contains
  module subroutine int2231(dd_pair, pf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: dd_pair, pf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n22bra(:), n13ket(:)
    real(dp), allocatable :: xint22bra(:), xint13ket(:)
    integer(kind=int64) :: nddbra, npfket
    real(dp) :: scutddbra, scutpfket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxj2
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
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(360), yin(360), zin(360)
    real(dp) :: eri_value(1080)
    real(dp) :: d22bra(36), d13ket(30)
    integer(kind=int64) :: ix(6), jx(6), kx(10), lx(3)
    integer(kind=int64) :: iy(6), jy(6), ky(10), ly(3)
    integer(kind=int64) :: iz(6), jz(6), kz(10), lz(3)
    integer(kind=int64) :: in(5), in1(5), kn(5)
    integer(kind=int64) :: ijx(36), ijy(36), ijz(36)
    integer(kind=int64) :: klx(30), kly(30), klz(30)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: iandj

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 25
    in1(3) = 49
    in1(4) = 57
    in1(5) = 65

    kn(1) = 0
    kn(2) = 2
    kn(3) = 4
    kn(4) = 6
    kn(5) = 7

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 1
    lx(2) = 0
    lx(3) = 0

    kx(1) = 6
    kx(2) = 0
    kx(3) = 0
    kx(4) = 4
    kx(5) = 4
    kx(6) = 2
    kx(7) = 0
    kx(8) = 2
    kx(9) = 0
    kx(10) = 2

    jx(1) = 16
    jx(2) = 0
    jx(3) = 0
    jx(4) = 8
    jx(5) = 8
    jx(6) = 0

    ix(1) = 49
    ix(2) = 1
    ix(3) = 1
    ix(4) = 25
    ix(5) = 25
    ix(6) = 1

    ! y-arrays

    ly(1) = 0
    ly(2) = 1
    ly(3) = 0

    ky(1) = 0
    ky(2) = 6
    ky(3) = 0
    ky(4) = 2
    ky(5) = 0
    ky(6) = 4
    ky(7) = 4
    ky(8) = 0
    ky(9) = 2
    ky(10) = 2

    jy(1) = 0
    jy(2) = 16
    jy(3) = 0
    jy(4) = 8
    jy(5) = 0
    jy(6) = 8

    iy(1) = 1
    iy(2) = 49
    iy(3) = 1
    iy(4) = 25
    iy(5) = 1
    iy(6) = 25

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 6
    kz(4) = 0
    kz(5) = 2
    kz(6) = 0
    kz(7) = 2
    kz(8) = 4
    kz(9) = 4
    kz(10) = 2

    jz(1) = 0
    jz(2) = 0
    jz(3) = 16
    jz(4) = 0
    jz(5) = 8
    jz(6) = 8

    iz(1) = 1
    iz(2) = 1
    iz(3) = 49
    iz(4) = 1
    iz(5) = 25
    iz(6) = 25

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 65
    ijx(2) = 49
    ijx(3) = 49
    ijx(4) = 57
    ijx(5) = 57
    ijx(6) = 49
    ijx(7) = 17
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 9
    ijx(11) = 9
    ijx(12) = 1
    ijx(13) = 17
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 9
    ijx(17) = 9
    ijx(18) = 1
    ijx(19) = 41
    ijx(20) = 25
    ijx(21) = 25
    ijx(22) = 33
    ijx(23) = 33
    ijx(24) = 25
    ijx(25) = 41
    ijx(26) = 25
    ijx(27) = 25
    ijx(28) = 33
    ijx(29) = 33
    ijx(30) = 25
    ijx(31) = 17
    ijx(32) = 1
    ijx(33) = 1
    ijx(34) = 9
    ijx(35) = 9
    ijx(36) = 1

    ijy(1) = 1
    ijy(2) = 17
    ijy(3) = 1
    ijy(4) = 9
    ijy(5) = 1
    ijy(6) = 9
    ijy(7) = 49
    ijy(8) = 65
    ijy(9) = 49
    ijy(10) = 57
    ijy(11) = 49
    ijy(12) = 57
    ijy(13) = 1
    ijy(14) = 17
    ijy(15) = 1
    ijy(16) = 9
    ijy(17) = 1
    ijy(18) = 9
    ijy(19) = 25
    ijy(20) = 41
    ijy(21) = 25
    ijy(22) = 33
    ijy(23) = 25
    ijy(24) = 33
    ijy(25) = 1
    ijy(26) = 17
    ijy(27) = 1
    ijy(28) = 9
    ijy(29) = 1
    ijy(30) = 9
    ijy(31) = 25
    ijy(32) = 41
    ijy(33) = 25
    ijy(34) = 33
    ijy(35) = 25
    ijy(36) = 33

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 17
    ijz(4) = 1
    ijz(5) = 9
    ijz(6) = 9
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 17
    ijz(10) = 1
    ijz(11) = 9
    ijz(12) = 9
    ijz(13) = 49
    ijz(14) = 49
    ijz(15) = 65
    ijz(16) = 49
    ijz(17) = 57
    ijz(18) = 57
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 17
    ijz(22) = 1
    ijz(23) = 9
    ijz(24) = 9
    ijz(25) = 25
    ijz(26) = 25
    ijz(27) = 41
    ijz(28) = 25
    ijz(29) = 33
    ijz(30) = 33
    ijz(31) = 25
    ijz(32) = 25
    ijz(33) = 41
    ijz(34) = 25
    ijz(35) = 33
    ijz(36) = 33

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 7
    klx(2) = 6
    klx(3) = 6
    klx(4) = 1
    klx(5) = 0
    klx(6) = 0
    klx(7) = 1
    klx(8) = 0
    klx(9) = 0
    klx(10) = 5
    klx(11) = 4
    klx(12) = 4
    klx(13) = 5
    klx(14) = 4
    klx(15) = 4
    klx(16) = 3
    klx(17) = 2
    klx(18) = 2
    klx(19) = 1
    klx(20) = 0
    klx(21) = 0
    klx(22) = 3
    klx(23) = 2
    klx(24) = 2
    klx(25) = 1
    klx(26) = 0
    klx(27) = 0
    klx(28) = 3
    klx(29) = 2
    klx(30) = 2

    kly(1) = 0
    kly(2) = 1
    kly(3) = 0
    kly(4) = 6
    kly(5) = 7
    kly(6) = 6
    kly(7) = 0
    kly(8) = 1
    kly(9) = 0
    kly(10) = 2
    kly(11) = 3
    kly(12) = 2
    kly(13) = 0
    kly(14) = 1
    kly(15) = 0
    kly(16) = 4
    kly(17) = 5
    kly(18) = 4
    kly(19) = 4
    kly(20) = 5
    kly(21) = 4
    kly(22) = 0
    kly(23) = 1
    kly(24) = 0
    kly(25) = 2
    kly(26) = 3
    kly(27) = 2
    kly(28) = 2
    kly(29) = 3
    kly(30) = 2

    klz(1) = 0
    klz(2) = 0
    klz(3) = 1
    klz(4) = 0
    klz(5) = 0
    klz(6) = 1
    klz(7) = 6
    klz(8) = 6
    klz(9) = 7
    klz(10) = 0
    klz(11) = 0
    klz(12) = 1
    klz(13) = 2
    klz(14) = 2
    klz(15) = 3
    klz(16) = 0
    klz(17) = 0
    klz(18) = 1
    klz(19) = 2
    klz(20) = 2
    klz(21) = 3
    klz(22) = 4
    klz(23) = 4
    klz(24) = 5
    klz(25) = 4
    klz(26) = 4
    klz(27) = 5
    klz(28) = 2
    klz(29) = 2
    klz(30) = 3

    allocate (n22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22bra(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (n13ket(res%n_p_shl*res%n_f_shl))
    allocate (xint13ket(res%n_p_shl*res%n_f_shl))

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

    scutpfket = cutoff_schwarz/maxval(pf_pair%xints)
    npfket = 0
    do ij = 1, res%n_p_shl*res%n_f_shl
      if (pf_pair%xints(ij) .ge. scutpfket) then
        npfket = npfket + 1
        xint13ket(npfket) = pf_pair%xints(ij)
        n13ket(npfket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((nddbra*npfket) .le. nchunksize_int64) nchunksize_int64 = nddbra*npfket
    ntile = int(nddbra*npfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nddbra*npfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, nddbra, xint22bra, n22bra, xint13ket, n13ket, dd_pair, pf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d13ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d22bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxj2,iandj)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, nddbra) + 1
              kl_tmp = (iquart - 1)/nddbra + 1

              test = xint22bra(ij_tmp)*xint13ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n22bra(ij_tmp)
                kl = n13ket(kl_tmp)

                ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
                jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_d_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
                ksh = res%i_f_shl(ksh_tmp)
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

                  t_expon_cd = pf_pair%t_expon_ab(pf_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = pf_pair%expon_b(pf_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = pf_pair%expon_a(pf_pair%pair_loc(kl) + ket_loop)
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

                  d13ket(1) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(2) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(3) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(4) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(5) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(6) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(7) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(8) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(9) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d13ket(10) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(11) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(12) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(13) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(14) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(15) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(16) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(17) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(18) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(19) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(20) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(21) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(22) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(23) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(24) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(25) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(26) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(27) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d13ket(28) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                  d13ket(29) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                  d13ket(30) = pf_pair%d_coeff_alt(pf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3

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

                    if (xx .ge. 55.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.1175813202117781D+00
                      rts(2) = factr*0.1074562012436902D+01
                      rts(3) = factr*0.3085937443717543D+01
                      rts(4) = factr*0.6414729733662035D+01
                      rts(5) = factr*0.1180718948997174D+02

                      wts(1) = factw*0.6108626337353259D+00
                      wts(2) = factw*0.2401386110823149D+00
                      wts(3) = factw*0.3387439445548109D-01
                      wts(4) = factw*0.1343645746781227D-02
                      wts(5) = factw*0.7640432855232614D-05

                    else ! "regular" evaluation

                      rgrid(1) = 0.1314956323727580D-05
                      rgrid(2) = 0.3638644488856184D-04
                      rgrid(3) = 0.2184836363610052D-03
                      rgrid(4) = 0.7467882061060832D-03
                      rgrid(5) = 0.1898594940807536D-02
                      rgrid(6) = 0.4018347564842567D-02
                      rgrid(7) = 0.7503899366138534D-02
                      rgrid(8) = 0.1279045049262753D-01
                      rgrid(9) = 0.2033269213059855D-01
                      rgrid(10) = 0.3058574876834499D-01
                      rgrid(11) = 0.4398555195277359D-01
                      rgrid(12) = 0.6092930105499368D-01
                      rgrid(13) = 0.8175666785532240D-01
                      rgrid(14) = 0.1067323821680571D+00
                      rgrid(15) = 0.1360307958356441D+00
                      rgrid(16) = 0.1697229634514230D+00
                      rgrid(17) = 0.2077667019402563D+00
                      rgrid(18) = 0.2499999999999998D+00
                      rgrid(19) = 0.2961380452159158D+00
                      rgrid(20) = 0.3457740246174122D+00
                      rgrid(21) = 0.3983837370449402D+00
                      rgrid(22) = 0.4533339365988709D+00
                      rgrid(23) = 0.5098942093731368D+00
                      rgrid(24) = 0.5672520742964827D+00
                      rgrid(25) = 0.6245308967025374D+00
                      rgrid(26) = 0.6808101134342349D+00
                      rgrid(27) = 0.7351471936872270D+00
                      rgrid(28) = 0.7866007027795398D+00
                      rgrid(29) = 0.8342537984583634D+00
                      rgrid(30) = 0.8772374725900651D+00
                      rgrid(31) = 0.9147528563001249D+00
                      rgrid(32) = 0.9460919364139329D+00
                      rgrid(33) = 0.9706560996755904D+00
                      rgrid(34) = 0.9879721508887392D+00
                      rgrid(35) = 0.9977078840559234D+00

                      wgrid(1) = 0.2941716710221880D-02*exp(-xx*0.1314956323727580D-05)
                      wgrid(2) = 0.6825414174180503D-02*exp(-xx*0.3638644488856184D-04)
                      wgrid(3) = 0.1066148995574242D-01*exp(-xx*0.2184836363610052D-03)
                      wgrid(4) = 0.1441463005444643D-01*exp(-xx*0.7467882061060832D-03)
                      wgrid(5) = 0.1805505793173160D-01*exp(-xx*0.1898594940807536D-02)
                      wgrid(6) = 0.2155421116308512D-01*exp(-xx*0.4018347564842567D-02)
                      wgrid(7) = 0.2488468520067677D-01*exp(-xx*0.7503899366138534D-02)
                      wgrid(8) = 0.2802040810618503D-01*exp(-xx*0.1279045049262753D-01)
                      wgrid(9) = 0.3093683598304008D-01*exp(-xx*0.2033269213059855D-01)
                      wgrid(10) = 0.3361114263454359D-01*exp(-xx*0.3058574876834499D-01)
                      wgrid(11) = 0.3602239738628011D-01*exp(-xx*0.4398555195277359D-01)
                      wgrid(12) = 0.3815172857772094D-01*exp(-xx*0.6092930105499368D-01)
                      wgrid(13) = 0.3998247112116213D-01*exp(-xx*0.8175666785532240D-01)
                      wgrid(14) = 0.4150029686442836D-01*exp(-xx*0.1067323821680571D+00)
                      wgrid(15) = 0.4269332669604951D-01*exp(-xx*0.1360307958356441D+00)
                      wgrid(16) = 0.4355222349859141D-01*exp(-xx*0.1697229634514230D+00)
                      wgrid(17) = 0.4407026521513794D-01*exp(-xx*0.2077667019402563D+00)
                      wgrid(18) = 0.4424339745355207D-01*exp(-xx*0.2499999999999998D+00)
                      wgrid(19) = 0.4407026521513754D-01*exp(-xx*0.2961380452159158D+00)
                      wgrid(20) = 0.4355222349859152D-01*exp(-xx*0.3457740246174122D+00)
                      wgrid(21) = 0.4269332669604975D-01*exp(-xx*0.3983837370449402D+00)
                      wgrid(22) = 0.4150029686442849D-01*exp(-xx*0.4533339365988709D+00)
                      wgrid(23) = 0.3998247112116278D-01*exp(-xx*0.5098942093731368D+00)
                      wgrid(24) = 0.3815172857772100D-01*exp(-xx*0.5672520742964827D+00)
                      wgrid(25) = 0.3602239738628051D-01*exp(-xx*0.6245308967025374D+00)
                      wgrid(26) = 0.3361114263454346D-01*exp(-xx*0.6808101134342349D+00)
                      wgrid(27) = 0.3093683598304035D-01*exp(-xx*0.7351471936872270D+00)
                      wgrid(28) = 0.2802040810618515D-01*exp(-xx*0.7866007027795398D+00)
                      wgrid(29) = 0.2488468520067681D-01*exp(-xx*0.8342537984583634D+00)
                      wgrid(30) = 0.2155421116308484D-01*exp(-xx*0.8772374725900651D+00)
                      wgrid(31) = 0.1805505793173228D-01*exp(-xx*0.9147528563001249D+00)
                      wgrid(32) = 0.1441463005444692D-01*exp(-xx*0.9460919364139329D+00)
                      wgrid(33) = 0.1066148995574185D-01*exp(-xx*0.9706560996755904D+00)
                      wgrid(34) = 0.6825414174181070D-02*exp(-xx*0.9879721508887392D+00)
                      wgrid(35) = 0.2941716710221409D-02*exp(-xx*0.9977078840559234D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 35
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 35
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 4

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 35

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
                        wrk(5) = 0.0D+00
                        do 100 kk = 2, 5

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 5

                            jj = 0

105                         do 110 m = l, 5
                              if (m .eq. 5) go to 120
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

                                do 300 ii = 2, 5

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 5
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

                                    do 310 kk = 1, 5
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 5
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

                                      ! i2 = in(2) =   25
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(25) = xc00
                                      yin(25) = yc00
                                      zin(25) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    3

                                      xin(3) = xcp00
                                      yin(3) = ycp00
                                      zin(3) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   27
                                      ! i2 =   25

                                      xin(27) = xcp00*xin(25) + cp10
                                      yin(27) = ycp00*yin(25) + cp10
                                      zin(27) = zcp00*zin(25) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   49
                                      ! i3 =    1
                                      ! i4 =   25

                                      xin(49) = c10*xin(1) + xc00*xin(25)
                                      yin(49) = c10*yin(1) + yc00*yin(25)
                                      zin(49) = c10*zin(1) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   51
                                      ! i5 =   49
                                      ! i4 =   25

                                      xin(51) = xcp00*xin(49) + cp10*xin(25)
                                      yin(51) = ycp00*yin(49) + cp10*yin(25)
                                      zin(51) = zcp00*zin(49) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   57
                                      ! i3 =   25
                                      ! i4 =   49

                                      xin(57) = c10*xin(25) + xc00*xin(49)
                                      yin(57) = c10*yin(25) + yc00*yin(49)
                                      zin(57) = c10*zin(25) + zc00*zin(49)

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

                                      ! i5 = in(n+1) =   65
                                      ! i3 =   49
                                      ! i4 =   57

                                      xin(65) = c10*xin(49) + xc00*xin(57)
                                      yin(65) = c10*yin(49) + yc00*yin(57)
                                      zin(65) = c10*zin(49) + zc00*zin(57)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   67
                                      ! i5 =   65
                                      ! i4 =   57

                                      xin(67) = xcp00*xin(65) + cp10*xin(57)
                                      yin(67) = ycp00*yin(65) + cp10*yin(57)
                                      zin(67) = zcp00*zin(65) + cp10*zin(57)

                                      ! ------------------

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    3

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    5
                                      ! i3 =    1
                                      ! i4 =    3

                                      xin(5) = cp01*xin(1) + xcp00*xin(3)
                                      yin(5) = cp01*yin(1) + ycp00*yin(3)
                                      zin(5) = cp01*zin(1) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   29

                                      xin(29) = xc00*xin(5) + c01*xin(3)
                                      yin(29) = yc00*yin(5) + c01*yin(3)
                                      zin(29) = zc00*zin(5) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    5

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    7
                                      ! i3 =    3
                                      ! i4 =    5

                                      xin(7) = cp01*xin(3) + xcp00*xin(5)
                                      yin(7) = cp01*yin(3) + ycp00*yin(5)
                                      zin(7) = cp01*zin(3) + zcp00*zin(5)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   31

                                      xin(31) = xc00*xin(7) + c01*xin(5)
                                      yin(31) = yc00*yin(7) + c01*yin(5)
                                      zin(31) = zc00*zin(7) + c01*zin(5)

                                      ! ------------------

                                      ! i3 = i4 =    5
                                      ! i4 = i5 =    7

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    8
                                      ! i3 =    5
                                      ! i4 =    7

                                      xin(8) = cp01*xin(5) + xcp00*xin(7)
                                      yin(8) = cp01*yin(5) + ycp00*yin(7)
                                      zin(8) = cp01*zin(5) + zcp00*zin(7)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   32

                                      xin(32) = xc00*xin(8) + c01*xin(7)
                                      yin(32) = yc00*yin(8) + c01*yin(7)
                                      zin(32) = zc00*zin(8) + c01*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   49

                                      xin(53) = c10*xin(5) + xc00*xin(29) + c01*xin(27)
                                      yin(53) = c10*yin(5) + yc00*yin(29) + c01*yin(27)
                                      zin(53) = c10*zin(5) + zc00*zin(29) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   57

                                      xin(61) = c10*xin(29) + xc00*xin(53) + c01*xin(51)
                                      yin(61) = c10*yin(29) + yc00*yin(53) + c01*yin(51)
                                      zin(61) = c10*zin(29) + zc00*zin(53) + c01*zin(51)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   57

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   65

                                      xin(69) = c10*xin(53) + xc00*xin(61) + c01*xin(59)
                                      yin(69) = c10*yin(53) + yc00*yin(61) + c01*yin(59)
                                      zin(69) = c10*zin(53) + zc00*zin(61) + c01*zin(59)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   49

                                      xin(55) = c10*xin(7) + xc00*xin(31) + c01*xin(29)
                                      yin(55) = c10*yin(7) + yc00*yin(31) + c01*yin(29)
                                      zin(55) = c10*zin(7) + zc00*zin(31) + c01*zin(29)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   57

                                      xin(63) = c10*xin(31) + xc00*xin(55) + c01*xin(53)
                                      yin(63) = c10*yin(31) + yc00*yin(55) + c01*yin(53)
                                      zin(63) = c10*zin(31) + zc00*zin(55) + c01*zin(53)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   57

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   65

                                      xin(71) = c10*xin(55) + xc00*xin(63) + c01*xin(61)
                                      yin(71) = c10*yin(55) + yc00*yin(63) + c01*yin(61)
                                      zin(71) = c10*zin(55) + zc00*zin(63) + c01*zin(61)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   25

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =   49

                                      xin(56) = c10*xin(8) + xc00*xin(32) + c01*xin(31)
                                      yin(56) = c10*yin(8) + yc00*yin(32) + c01*yin(31)
                                      zin(56) = c10*zin(8) + zc00*zin(32) + c01*zin(31)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   49

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   57

                                      xin(64) = c10*xin(32) + xc00*xin(56) + c01*xin(55)
                                      yin(64) = c10*yin(32) + yc00*yin(56) + c01*yin(55)
                                      zin(64) = c10*zin(32) + zc00*zin(56) + c01*zin(55)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   49
                                      ! i4 = i5 =   57

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   65

                                      xin(72) = c10*xin(56) + xc00*xin(64) + c01*xin(63)
                                      yin(72) = c10*yin(56) + yc00*yin(64) + c01*yin(63)
                                      zin(72) = c10*zin(56) + zc00*zin(64) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   57
                                      ! i4 = i5 =   65

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   65

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   65

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(65) = xin(65) + dxij*xin(57)
                                      yin(65) = yin(65) + dyij*yin(57)
                                      zin(65) = zin(65) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   49

                                      xin(57) = xin(57) + dxij*xin(49)
                                      yin(57) = yin(57) + dyij*yin(49)
                                      zin(57) = zin(57) + dzij*zin(49)

                                      ! i3 = i4 =   49
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   65

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   57

                                      xin(65) = xin(65) + dxij*xin(57)
                                      yin(65) = yin(65) + dyij*yin(57)
                                      zin(65) = zin(65) + dzij*zin(57)

                                      ! i3 = i4 =   57
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    9

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    2

                                      xin(9) = xin(25) + dxij*xin(1)
                                      yin(9) = yin(25) + dyij*yin(1)
                                      zin(9) = zin(25) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    2

                                      xin(33) = xin(49) + dxij*xin(25)
                                      yin(33) = yin(49) + dyij*yin(25)
                                      zin(33) = zin(49) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   57

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   17

                                      ! nj =    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    2

                                      xin(17) = xin(33) + dxij*xin(9)
                                      yin(17) = yin(33) + dyij*yin(9)
                                      zin(17) = zin(33) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    2

                                      xin(41) = xin(57) + dxij*xin(33)
                                      yin(41) = yin(57) + dyij*yin(33)
                                      zin(41) = zin(57) + dzij*zin(33)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   67

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(67) = xin(67) + dxij*xin(59)
                                      yin(67) = yin(67) + dyij*yin(59)
                                      zin(67) = zin(67) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   51

                                      xin(59) = xin(59) + dxij*xin(51)
                                      yin(59) = yin(59) + dyij*yin(51)
                                      zin(59) = zin(59) + dzij*zin(51)

                                      ! i3 = i4 =   51
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   67

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   59

                                      xin(67) = xin(67) + dxij*xin(59)
                                      yin(67) = yin(67) + dyij*yin(59)
                                      zin(67) = zin(67) + dzij*zin(59)

                                      ! i3 = i4 =   59
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   11

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    2

                                      xin(11) = xin(27) + dxij*xin(3)
                                      yin(11) = yin(27) + dyij*yin(3)
                                      zin(11) = zin(27) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    2

                                      xin(35) = xin(51) + dxij*xin(27)
                                      yin(35) = yin(51) + dyij*yin(27)
                                      zin(35) = zin(51) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    2

                                      xin(19) = xin(35) + dxij*xin(11)
                                      yin(19) = yin(35) + dyij*yin(11)
                                      zin(19) = zin(35) + dzij*zin(11)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    2

                                      xin(43) = xin(59) + dxij*xin(35)
                                      yin(43) = yin(59) + dyij*yin(35)
                                      zin(43) = zin(59) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   61

                                      xin(69) = xin(69) + dxij*xin(61)
                                      yin(69) = yin(69) + dyij*yin(61)
                                      zin(69) = zin(69) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   53

                                      xin(61) = xin(61) + dxij*xin(53)
                                      yin(61) = yin(61) + dyij*yin(53)
                                      zin(61) = zin(61) + dzij*zin(53)

                                      ! i3 = i4 =   53
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   69

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   61

                                      xin(69) = xin(69) + dxij*xin(61)
                                      yin(69) = yin(69) + dyij*yin(61)
                                      zin(69) = zin(69) + dzij*zin(61)

                                      ! i3 = i4 =   61
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   13

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    2

                                      xin(13) = xin(29) + dxij*xin(5)
                                      yin(13) = yin(29) + dyij*yin(5)
                                      zin(13) = zin(29) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! ni =    2

                                      xin(37) = xin(53) + dxij*xin(29)
                                      yin(37) = yin(53) + dyij*yin(29)
                                      zin(37) = zin(53) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   21

                                      ! nj =    2

                                      ! i4 = i3 =   21

                                      ! do ni = 1,    2

                                      xin(21) = xin(37) + dxij*xin(13)
                                      yin(21) = yin(37) + dyij*yin(13)
                                      zin(21) = zin(37) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    2

                                      xin(45) = xin(61) + dxij*xin(37)
                                      yin(45) = yin(61) + dyij*yin(37)
                                      zin(45) = zin(61) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   29

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(71) = xin(71) + dxij*xin(63)
                                      yin(71) = yin(71) + dyij*yin(63)
                                      zin(71) = zin(71) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   55

                                      xin(63) = xin(63) + dxij*xin(55)
                                      yin(63) = yin(63) + dyij*yin(55)
                                      zin(63) = zin(63) + dzij*zin(55)

                                      ! i3 = i4 =   55
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   71

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   63

                                      xin(71) = xin(71) + dxij*xin(63)
                                      yin(71) = yin(71) + dyij*yin(63)
                                      zin(71) = zin(71) + dzij*zin(63)

                                      ! i3 = i4 =   63
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   15

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   15

                                      ! do ni = 1,    2

                                      xin(15) = xin(31) + dxij*xin(7)
                                      yin(15) = yin(31) + dyij*yin(7)
                                      zin(15) = zin(31) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   39

                                      ! ni =    2

                                      xin(39) = xin(55) + dxij*xin(31)
                                      yin(39) = yin(55) + dyij*yin(31)
                                      zin(39) = zin(55) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   63

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   23

                                      ! nj =    2

                                      ! i4 = i3 =   23

                                      ! do ni = 1,    2

                                      xin(23) = xin(39) + dxij*xin(15)
                                      yin(23) = yin(39) + dyij*yin(15)
                                      zin(23) = zin(39) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    2

                                      xin(47) = xin(63) + dxij*xin(39)
                                      yin(47) = yin(63) + dyij*yin(39)
                                      zin(47) = zin(63) + dzij*zin(39)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   31

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   64

                                      xin(72) = xin(72) + dxij*xin(64)
                                      yin(72) = yin(72) + dyij*yin(64)
                                      zin(72) = zin(72) + dzij*zin(64)

                                      ! i3 = i4 =   64
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =   56

                                      xin(64) = xin(64) + dxij*xin(56)
                                      yin(64) = yin(64) + dyij*yin(56)
                                      zin(64) = zin(64) + dzij*zin(56)

                                      ! i3 = i4 =   56
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =   72

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   64

                                      xin(72) = xin(72) + dxij*xin(64)
                                      yin(72) = yin(72) + dyij*yin(64)
                                      zin(72) = zin(72) + dzij*zin(64)

                                      ! i3 = i4 =   64
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   16

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   16

                                      ! do ni = 1,    2

                                      xin(16) = xin(32) + dxij*xin(8)
                                      yin(16) = yin(32) + dyij*yin(8)
                                      zin(16) = zin(32) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    2

                                      xin(40) = xin(56) + dxij*xin(32)
                                      yin(40) = yin(56) + dyij*yin(32)
                                      zin(40) = zin(56) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   24

                                      ! nj =    2

                                      ! i4 = i3 =   24

                                      ! do ni = 1,    2

                                      xin(24) = xin(40) + dxij*xin(16)
                                      yin(24) = yin(40) + dyij*yin(16)
                                      zin(24) = zin(40) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    2

                                      xin(48) = xin(64) + dxij*xin(40)
                                      yin(48) = yin(64) + dyij*yin(40)
                                      zin(48) = zin(64) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   32

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    8

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    1

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    3

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

                                      xin(6) = xin(7) + dxkl*xin(5)
                                      yin(6) = yin(7) + dykl*yin(5)
                                      zin(6) = zin(7) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =    9

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   16

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   15

                                      xin(16) = xin(16) + dxkl*xin(15)
                                      yin(16) = yin(16) + dykl*yin(15)
                                      zin(16) = zin(16) + dzkl*zin(15)

                                      ! i3 = i4 =   15
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   10

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   10

                                      ! do nk = 1,    3

                                      xin(10) = xin(11) + dxkl*xin(9)
                                      yin(10) = yin(11) + dykl*yin(9)
                                      zin(10) = zin(11) + dzkl*zin(9)
                                      ! i4 = i4 + lang+1 =   12

                                      ! nk =    2

                                      xin(12) = xin(13) + dxkl*xin(11)
                                      yin(12) = yin(13) + dykl*yin(11)
                                      zin(12) = zin(13) + dzkl*zin(11)
                                      ! i4 = i4 + lang+1 =   14

                                      ! nk =    3

                                      xin(14) = xin(15) + dxkl*xin(13)
                                      yin(14) = yin(15) + dykl*yin(13)
                                      zin(14) = zin(15) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   16

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   11

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   17

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   24

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   23

                                      xin(24) = xin(24) + dxkl*xin(23)
                                      yin(24) = yin(24) + dykl*yin(23)
                                      zin(24) = zin(24) + dzkl*zin(23)

                                      ! i3 = i4 =   23
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   18

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   18

                                      ! do nk = 1,    3

                                      xin(18) = xin(19) + dxkl*xin(17)
                                      yin(18) = yin(19) + dykl*yin(17)
                                      zin(18) = zin(19) + dzkl*zin(17)
                                      ! i4 = i4 + lang+1 =   20

                                      ! nk =    2

                                      xin(20) = xin(21) + dxkl*xin(19)
                                      yin(20) = yin(21) + dykl*yin(19)
                                      zin(20) = zin(21) + dzkl*zin(19)
                                      ! i4 = i4 + lang+1 =   22

                                      ! nk =    3

                                      xin(22) = xin(23) + dxkl*xin(21)
                                      yin(22) = yin(23) + dykl*yin(21)
                                      zin(22) = zin(23) + dzkl*zin(21)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   19

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   25

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   25

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   32

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   31

                                      xin(32) = xin(32) + dxkl*xin(31)
                                      yin(32) = yin(32) + dykl*yin(31)
                                      zin(32) = zin(32) + dzkl*zin(31)

                                      ! i3 = i4 =   31
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   26

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   26

                                      ! do nk = 1,    3

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

                                      xin(30) = xin(31) + dxkl*xin(29)
                                      yin(30) = yin(31) + dykl*yin(29)
                                      zin(30) = zin(31) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   27

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   33

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   40

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   39

                                      xin(40) = xin(40) + dxkl*xin(39)
                                      yin(40) = yin(40) + dykl*yin(39)
                                      zin(40) = zin(40) + dzkl*zin(39)

                                      ! i3 = i4 =   39
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   34

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   34

                                      ! do nk = 1,    3

                                      xin(34) = xin(35) + dxkl*xin(33)
                                      yin(34) = yin(35) + dykl*yin(33)
                                      zin(34) = zin(35) + dzkl*zin(33)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    2

                                      xin(36) = xin(37) + dxkl*xin(35)
                                      yin(36) = yin(37) + dykl*yin(35)
                                      zin(36) = zin(37) + dzkl*zin(35)
                                      ! i4 = i4 + lang+1 =   38

                                      ! nk =    3

                                      xin(38) = xin(39) + dxkl*xin(37)
                                      yin(38) = yin(39) + dykl*yin(37)
                                      zin(38) = zin(39) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   40

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   35

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   41

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   48

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   47

                                      xin(48) = xin(48) + dxkl*xin(47)
                                      yin(48) = yin(48) + dykl*yin(47)
                                      zin(48) = zin(48) + dzkl*zin(47)

                                      ! i3 = i4 =   47
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   42

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   42

                                      ! do nk = 1,    3

                                      xin(42) = xin(43) + dxkl*xin(41)
                                      yin(42) = yin(43) + dykl*yin(41)
                                      zin(42) = zin(43) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    2

                                      xin(44) = xin(45) + dxkl*xin(43)
                                      yin(44) = yin(45) + dykl*yin(43)
                                      zin(44) = zin(45) + dzkl*zin(43)
                                      ! i4 = i4 + lang+1 =   46

                                      ! nk =    3

                                      xin(46) = xin(47) + dxkl*xin(45)
                                      yin(46) = yin(47) + dykl*yin(45)
                                      zin(46) = zin(47) + dzkl*zin(45)
                                      ! i4 = i4 + lang+1 =   48

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   43

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   49

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   56

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   55

                                      xin(56) = xin(56) + dxkl*xin(55)
                                      yin(56) = yin(56) + dykl*yin(55)
                                      zin(56) = zin(56) + dzkl*zin(55)

                                      ! i3 = i4 =   55
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   50

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   50

                                      ! do nk = 1,    3

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

                                      xin(54) = xin(55) + dxkl*xin(53)
                                      yin(54) = yin(55) + dykl*yin(53)
                                      zin(54) = zin(55) + dzkl*zin(53)
                                      ! i4 = i4 + lang+1 =   56

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   51

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   57

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   64

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   63

                                      xin(64) = xin(64) + dxkl*xin(63)
                                      yin(64) = yin(64) + dykl*yin(63)
                                      zin(64) = zin(64) + dzkl*zin(63)

                                      ! i3 = i4 =   63
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   58

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   58

                                      ! do nk = 1,    3

                                      xin(58) = xin(59) + dxkl*xin(57)
                                      yin(58) = yin(59) + dykl*yin(57)
                                      zin(58) = zin(59) + dzkl*zin(57)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      xin(60) = xin(61) + dxkl*xin(59)
                                      yin(60) = yin(61) + dykl*yin(59)
                                      zin(60) = zin(61) + dzkl*zin(59)
                                      ! i4 = i4 + lang+1 =   62

                                      ! nk =    3

                                      xin(62) = xin(63) + dxkl*xin(61)
                                      yin(62) = yin(63) + dykl*yin(61)
                                      zin(62) = zin(63) + dzkl*zin(61)
                                      ! i4 = i4 + lang+1 =   64

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   59

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   65

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   66

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    3

                                      xin(66) = xin(67) + dxkl*xin(65)
                                      yin(66) = yin(67) + dykl*yin(65)
                                      zin(66) = zin(67) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      xin(68) = xin(69) + dxkl*xin(67)
                                      yin(68) = yin(69) + dykl*yin(67)
                                      zin(68) = zin(69) + dzkl*zin(67)
                                      ! i4 = i4 + lang+1 =   70

                                      ! nk =    3

                                      xin(70) = xin(71) + dxkl*xin(69)
                                      yin(70) = yin(71) + dykl*yin(69)
                                      zin(70) = zin(71) + dzkl*zin(69)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   72

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

                                      ! i1 = in(1) =   73

                                      xin(73) = 1.0_dp
                                      yin(73) = 1.0_dp
                                      zin(73) = f00

                                      ! i2 = in(2) =   97
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(97) = xc00
                                      yin(97) = yc00
                                      zin(97) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   75

                                      xin(75) = xcp00
                                      yin(75) = ycp00
                                      zin(75) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   99
                                      ! i2 =   97

                                      xin(99) = xcp00*xin(97) + cp10
                                      yin(99) = ycp00*yin(97) + cp10
                                      zin(99) = zcp00*zin(97) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   97

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =   73
                                      ! i4 =   97

                                      xin(121) = c10*xin(73) + xc00*xin(97)
                                      yin(121) = c10*yin(73) + yc00*yin(97)
                                      zin(121) = c10*zin(73) + zc00*zin(97)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  123
                                      ! i5 =  121
                                      ! i4 =   97

                                      xin(123) = xcp00*xin(121) + cp10*xin(97)
                                      yin(123) = ycp00*yin(121) + cp10*yin(97)
                                      zin(123) = zcp00*zin(121) + cp10*zin(97)

                                      ! ------------------

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  121

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  129
                                      ! i3 =   97
                                      ! i4 =  121

                                      xin(129) = c10*xin(97) + xc00*xin(121)
                                      yin(129) = c10*yin(97) + yc00*yin(121)
                                      zin(129) = c10*zin(97) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  131
                                      ! i5 =  129
                                      ! i4 =  121

                                      xin(131) = xcp00*xin(129) + cp10*xin(121)
                                      yin(131) = ycp00*yin(129) + cp10*yin(121)
                                      zin(131) = zcp00*zin(129) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  129

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  137
                                      ! i3 =  121
                                      ! i4 =  129

                                      xin(137) = c10*xin(121) + xc00*xin(129)
                                      yin(137) = c10*yin(121) + yc00*yin(129)
                                      zin(137) = c10*zin(121) + zc00*zin(129)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  139
                                      ! i5 =  137
                                      ! i4 =  129

                                      xin(139) = xcp00*xin(137) + cp10*xin(129)
                                      yin(139) = ycp00*yin(137) + cp10*yin(129)
                                      zin(139) = zcp00*zin(137) + cp10*zin(129)

                                      ! ------------------

                                      ! i3 = i4 =  129
                                      ! i4 = i5 =  137

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   73
                                      ! i4 = i1+k2 =   75

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   77
                                      ! i3 =   73
                                      ! i4 =   75

                                      xin(77) = cp01*xin(73) + xcp00*xin(75)
                                      yin(77) = cp01*yin(73) + ycp00*yin(75)
                                      zin(77) = cp01*zin(73) + zcp00*zin(75)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  101

                                      xin(101) = xc00*xin(77) + c01*xin(75)
                                      yin(101) = yc00*yin(77) + c01*yin(75)
                                      zin(101) = zc00*zin(77) + c01*zin(75)

                                      ! ------------------

                                      ! i3 = i4 =   75
                                      ! i4 = i5 =   77

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   79
                                      ! i3 =   75
                                      ! i4 =   77

                                      xin(79) = cp01*xin(75) + xcp00*xin(77)
                                      yin(79) = cp01*yin(75) + ycp00*yin(77)
                                      zin(79) = cp01*zin(75) + zcp00*zin(77)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  103

                                      xin(103) = xc00*xin(79) + c01*xin(77)
                                      yin(103) = yc00*yin(79) + c01*yin(77)
                                      zin(103) = zc00*zin(79) + c01*zin(77)

                                      ! ------------------

                                      ! i3 = i4 =   77
                                      ! i4 = i5 =   79

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   80
                                      ! i3 =   77
                                      ! i4 =   79

                                      xin(80) = cp01*xin(77) + xcp00*xin(79)
                                      yin(80) = cp01*yin(77) + ycp00*yin(79)
                                      zin(80) = cp01*zin(77) + zcp00*zin(79)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  104

                                      xin(104) = xc00*xin(80) + c01*xin(79)
                                      yin(104) = yc00*yin(80) + c01*yin(79)
                                      zin(104) = zc00*zin(80) + c01*zin(79)

                                      ! ------------------

                                      ! i3 = i4 =   79
                                      ! i4 = i5 =   80

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   97

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  121

                                      xin(125) = c10*xin(77) + xc00*xin(101) + c01*xin(99)
                                      yin(125) = c10*yin(77) + yc00*yin(101) + c01*yin(99)
                                      zin(125) = c10*zin(77) + zc00*zin(101) + c01*zin(99)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  129

                                      xin(133) = c10*xin(101) + xc00*xin(125) + c01*xin(123)
                                      yin(133) = c10*yin(101) + yc00*yin(125) + c01*yin(123)
                                      zin(133) = c10*zin(101) + zc00*zin(125) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  129

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(141) = c10*xin(125) + xc00*xin(133) + c01*xin(131)
                                      yin(141) = c10*yin(125) + yc00*yin(133) + c01*yin(131)
                                      zin(141) = c10*zin(125) + zc00*zin(133) + c01*zin(131)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  129
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   97

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  121

                                      xin(127) = c10*xin(79) + xc00*xin(103) + c01*xin(101)
                                      yin(127) = c10*yin(79) + yc00*yin(103) + c01*yin(101)
                                      zin(127) = c10*zin(79) + zc00*zin(103) + c01*zin(101)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  129

                                      xin(135) = c10*xin(103) + xc00*xin(127) + c01*xin(125)
                                      yin(135) = c10*yin(103) + yc00*yin(127) + c01*yin(125)
                                      zin(135) = c10*zin(103) + zc00*zin(127) + c01*zin(125)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  129

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(143) = c10*xin(127) + xc00*xin(135) + c01*xin(133)
                                      yin(143) = c10*yin(127) + yc00*yin(135) + c01*yin(133)
                                      zin(143) = c10*zin(127) + zc00*zin(135) + c01*zin(133)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  129
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =   73
                                      ! i4 = i2 =   97

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  121

                                      xin(128) = c10*xin(80) + xc00*xin(104) + c01*xin(103)
                                      yin(128) = c10*yin(80) + yc00*yin(104) + c01*yin(103)
                                      zin(128) = c10*zin(80) + zc00*zin(104) + c01*zin(103)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   97
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  129

                                      xin(136) = c10*xin(104) + xc00*xin(128) + c01*xin(127)
                                      yin(136) = c10*yin(104) + yc00*yin(128) + c01*yin(127)
                                      zin(136) = c10*zin(104) + zc00*zin(128) + c01*zin(127)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  129

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(144) = c10*xin(128) + xc00*xin(136) + c01*xin(135)
                                      yin(144) = c10*yin(128) + yc00*yin(136) + c01*yin(135)
                                      zin(144) = c10*zin(128) + zc00*zin(136) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  129
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  137

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  137

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  129

                                      xin(137) = xin(137) + dxij*xin(129)
                                      yin(137) = yin(137) + dyij*yin(129)
                                      zin(137) = zin(137) + dzij*zin(129)

                                      ! i3 = i4 =  129
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  121

                                      xin(129) = xin(129) + dxij*xin(121)
                                      yin(129) = yin(129) + dyij*yin(121)
                                      zin(129) = zin(129) + dzij*zin(121)

                                      ! i3 = i4 =  121
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  137

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  129

                                      xin(137) = xin(137) + dxij*xin(129)
                                      yin(137) = yin(137) + dyij*yin(129)
                                      zin(137) = zin(137) + dzij*zin(129)

                                      ! i3 = i4 =  129
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   81

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   81

                                      ! do ni = 1,    2

                                      xin(81) = xin(97) + dxij*xin(73)
                                      yin(81) = yin(97) + dyij*yin(73)
                                      zin(81) = zin(97) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  105

                                      ! ni =    2

                                      xin(105) = xin(121) + dxij*xin(97)
                                      yin(105) = yin(121) + dyij*yin(97)
                                      zin(105) = zin(121) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   89

                                      ! nj =    2

                                      ! i4 = i3 =   89

                                      ! do ni = 1,    2

                                      xin(89) = xin(105) + dxij*xin(81)
                                      yin(89) = yin(105) + dyij*yin(81)
                                      zin(89) = zin(105) + dzij*zin(81)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! ni =    2

                                      xin(113) = xin(129) + dxij*xin(105)
                                      yin(113) = yin(129) + dyij*yin(105)
                                      zin(113) = zin(129) + dzij*zin(105)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   97

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  131

                                      xin(139) = xin(139) + dxij*xin(131)
                                      yin(139) = yin(139) + dyij*yin(131)
                                      zin(139) = zin(139) + dzij*zin(131)

                                      ! i3 = i4 =  131
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  123

                                      xin(131) = xin(131) + dxij*xin(123)
                                      yin(131) = yin(131) + dyij*yin(123)
                                      zin(131) = zin(131) + dzij*zin(123)

                                      ! i3 = i4 =  123
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  139

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  131

                                      xin(139) = xin(139) + dxij*xin(131)
                                      yin(139) = yin(139) + dyij*yin(131)
                                      zin(139) = zin(139) + dzij*zin(131)

                                      ! i3 = i4 =  131
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   83

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   83

                                      ! do ni = 1,    2

                                      xin(83) = xin(99) + dxij*xin(75)
                                      yin(83) = yin(99) + dyij*yin(75)
                                      zin(83) = zin(99) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    2

                                      xin(107) = xin(123) + dxij*xin(99)
                                      yin(107) = yin(123) + dyij*yin(99)
                                      zin(107) = zin(123) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   91

                                      ! nj =    2

                                      ! i4 = i3 =   91

                                      ! do ni = 1,    2

                                      xin(91) = xin(107) + dxij*xin(83)
                                      yin(91) = yin(107) + dyij*yin(83)
                                      zin(91) = zin(107) + dzij*zin(83)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    2

                                      xin(115) = xin(131) + dxij*xin(107)
                                      yin(115) = yin(131) + dyij*yin(107)
                                      zin(115) = zin(131) + dzij*zin(107)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   99

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(141) = xin(141) + dxij*xin(133)
                                      yin(141) = yin(141) + dyij*yin(133)
                                      zin(141) = zin(141) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  125

                                      xin(133) = xin(133) + dxij*xin(125)
                                      yin(133) = yin(133) + dyij*yin(125)
                                      zin(133) = zin(133) + dzij*zin(125)

                                      ! i3 = i4 =  125
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  133

                                      xin(141) = xin(141) + dxij*xin(133)
                                      yin(141) = yin(141) + dyij*yin(133)
                                      zin(141) = zin(141) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   85

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   85

                                      ! do ni = 1,    2

                                      xin(85) = xin(101) + dxij*xin(77)
                                      yin(85) = yin(101) + dyij*yin(77)
                                      zin(85) = zin(101) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! ni =    2

                                      xin(109) = xin(125) + dxij*xin(101)
                                      yin(109) = yin(125) + dyij*yin(101)
                                      zin(109) = zin(125) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   93

                                      ! nj =    2

                                      ! i4 = i3 =   93

                                      ! do ni = 1,    2

                                      xin(93) = xin(109) + dxij*xin(85)
                                      yin(93) = yin(109) + dyij*yin(85)
                                      zin(93) = zin(109) + dzij*zin(85)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    2

                                      xin(117) = xin(133) + dxij*xin(109)
                                      yin(117) = yin(133) + dyij*yin(109)
                                      zin(117) = zin(133) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  101

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(143) = xin(143) + dxij*xin(135)
                                      yin(143) = yin(143) + dyij*yin(135)
                                      zin(143) = zin(143) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  127

                                      xin(135) = xin(135) + dxij*xin(127)
                                      yin(135) = yin(135) + dyij*yin(127)
                                      zin(135) = zin(135) + dzij*zin(127)

                                      ! i3 = i4 =  127
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  135

                                      xin(143) = xin(143) + dxij*xin(135)
                                      yin(143) = yin(143) + dyij*yin(135)
                                      zin(143) = zin(143) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   87

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   87

                                      ! do ni = 1,    2

                                      xin(87) = xin(103) + dxij*xin(79)
                                      yin(87) = yin(103) + dyij*yin(79)
                                      zin(87) = zin(103) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  111

                                      ! ni =    2

                                      xin(111) = xin(127) + dxij*xin(103)
                                      yin(111) = yin(127) + dyij*yin(103)
                                      zin(111) = zin(127) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  135

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   95

                                      ! nj =    2

                                      ! i4 = i3 =   95

                                      ! do ni = 1,    2

                                      xin(95) = xin(111) + dxij*xin(87)
                                      yin(95) = yin(111) + dyij*yin(87)
                                      zin(95) = zin(111) + dzij*zin(87)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    2

                                      xin(119) = xin(135) + dxij*xin(111)
                                      yin(119) = yin(135) + dyij*yin(111)
                                      zin(119) = zin(135) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  103

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  136

                                      xin(144) = xin(144) + dxij*xin(136)
                                      yin(144) = yin(144) + dyij*yin(136)
                                      zin(144) = zin(144) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  128

                                      xin(136) = xin(136) + dxij*xin(128)
                                      yin(136) = yin(136) + dyij*yin(128)
                                      zin(136) = zin(136) + dzij*zin(128)

                                      ! i3 = i4 =  128
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  136

                                      xin(144) = xin(144) + dxij*xin(136)
                                      yin(144) = yin(144) + dyij*yin(136)
                                      zin(144) = zin(144) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   88

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   88

                                      ! do ni = 1,    2

                                      xin(88) = xin(104) + dxij*xin(80)
                                      yin(88) = yin(104) + dyij*yin(80)
                                      zin(88) = zin(104) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  112

                                      ! ni =    2

                                      xin(112) = xin(128) + dxij*xin(104)
                                      yin(112) = yin(128) + dyij*yin(104)
                                      zin(112) = zin(128) + dzij*zin(104)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  136

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   96

                                      ! nj =    2

                                      ! i4 = i3 =   96

                                      ! do ni = 1,    2

                                      xin(96) = xin(112) + dxij*xin(88)
                                      yin(96) = yin(112) + dyij*yin(88)
                                      zin(96) = zin(112) + dzij*zin(88)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    2

                                      xin(120) = xin(136) + dxij*xin(112)
                                      yin(120) = yin(136) + dyij*yin(112)
                                      zin(120) = zin(136) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  104

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =   73

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   80

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    3

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

                                      xin(78) = xin(79) + dxkl*xin(77)
                                      yin(78) = yin(79) + dykl*yin(77)
                                      zin(78) = zin(79) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   81

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   88

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   87

                                      xin(88) = xin(88) + dxkl*xin(87)
                                      yin(88) = yin(88) + dykl*yin(87)
                                      zin(88) = zin(88) + dzkl*zin(87)

                                      ! i3 = i4 =   87
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   82

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   82

                                      ! do nk = 1,    3

                                      xin(82) = xin(83) + dxkl*xin(81)
                                      yin(82) = yin(83) + dykl*yin(81)
                                      zin(82) = zin(83) + dzkl*zin(81)
                                      ! i4 = i4 + lang+1 =   84

                                      ! nk =    2

                                      xin(84) = xin(85) + dxkl*xin(83)
                                      yin(84) = yin(85) + dykl*yin(83)
                                      zin(84) = zin(85) + dzkl*zin(83)
                                      ! i4 = i4 + lang+1 =   86

                                      ! nk =    3

                                      xin(86) = xin(87) + dxkl*xin(85)
                                      yin(86) = yin(87) + dykl*yin(85)
                                      zin(86) = zin(87) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   88

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   83

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   89

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   96

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   95

                                      xin(96) = xin(96) + dxkl*xin(95)
                                      yin(96) = yin(96) + dykl*yin(95)
                                      zin(96) = zin(96) + dzkl*zin(95)

                                      ! i3 = i4 =   95
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   90

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   90

                                      ! do nk = 1,    3

                                      xin(90) = xin(91) + dxkl*xin(89)
                                      yin(90) = yin(91) + dykl*yin(89)
                                      zin(90) = zin(91) + dzkl*zin(89)
                                      ! i4 = i4 + lang+1 =   92

                                      ! nk =    2

                                      xin(92) = xin(93) + dxkl*xin(91)
                                      yin(92) = yin(93) + dykl*yin(91)
                                      zin(92) = zin(93) + dzkl*zin(91)
                                      ! i4 = i4 + lang+1 =   94

                                      ! nk =    3

                                      xin(94) = xin(95) + dxkl*xin(93)
                                      yin(94) = yin(95) + dykl*yin(93)
                                      zin(94) = zin(95) + dzkl*zin(93)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   91

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   97

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  104

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  103

                                      xin(104) = xin(104) + dxkl*xin(103)
                                      yin(104) = yin(104) + dykl*yin(103)
                                      zin(104) = zin(104) + dzkl*zin(103)

                                      ! i3 = i4 =  103
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   98

                                      ! do nl = 1,    1

                                      ! i4 = i3 =   98

                                      ! do nk = 1,    3

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

                                      xin(102) = xin(103) + dxkl*xin(101)
                                      yin(102) = yin(103) + dykl*yin(101)
                                      zin(102) = zin(103) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =   99

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  105

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  112

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  111

                                      xin(112) = xin(112) + dxkl*xin(111)
                                      yin(112) = yin(112) + dykl*yin(111)
                                      zin(112) = zin(112) + dzkl*zin(111)

                                      ! i3 = i4 =  111
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  106

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  106

                                      ! do nk = 1,    3

                                      xin(106) = xin(107) + dxkl*xin(105)
                                      yin(106) = yin(107) + dykl*yin(105)
                                      zin(106) = zin(107) + dzkl*zin(105)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    2

                                      xin(108) = xin(109) + dxkl*xin(107)
                                      yin(108) = yin(109) + dykl*yin(107)
                                      zin(108) = zin(109) + dzkl*zin(107)
                                      ! i4 = i4 + lang+1 =  110

                                      ! nk =    3

                                      xin(110) = xin(111) + dxkl*xin(109)
                                      yin(110) = yin(111) + dykl*yin(109)
                                      zin(110) = zin(111) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  112

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  107

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  113

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  120

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  119

                                      xin(120) = xin(120) + dxkl*xin(119)
                                      yin(120) = yin(120) + dykl*yin(119)
                                      zin(120) = zin(120) + dzkl*zin(119)

                                      ! i3 = i4 =  119
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  114

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  114

                                      ! do nk = 1,    3

                                      xin(114) = xin(115) + dxkl*xin(113)
                                      yin(114) = yin(115) + dykl*yin(113)
                                      zin(114) = zin(115) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    2

                                      xin(116) = xin(117) + dxkl*xin(115)
                                      yin(116) = yin(117) + dykl*yin(115)
                                      zin(116) = zin(117) + dzkl*zin(115)
                                      ! i4 = i4 + lang+1 =  118

                                      ! nk =    3

                                      xin(118) = xin(119) + dxkl*xin(117)
                                      yin(118) = yin(119) + dykl*yin(117)
                                      zin(118) = zin(119) + dzkl*zin(117)
                                      ! i4 = i4 + lang+1 =  120

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  115

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  121

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  121

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  128

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  127

                                      xin(128) = xin(128) + dxkl*xin(127)
                                      yin(128) = yin(128) + dykl*yin(127)
                                      zin(128) = zin(128) + dzkl*zin(127)

                                      ! i3 = i4 =  127
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  122

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  122

                                      ! do nk = 1,    3

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

                                      xin(126) = xin(127) + dxkl*xin(125)
                                      yin(126) = yin(127) + dykl*yin(125)
                                      zin(126) = zin(127) + dzkl*zin(125)
                                      ! i4 = i4 + lang+1 =  128

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  123

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  129

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  136

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  135

                                      xin(136) = xin(136) + dxkl*xin(135)
                                      yin(136) = yin(136) + dykl*yin(135)
                                      zin(136) = zin(136) + dzkl*zin(135)

                                      ! i3 = i4 =  135
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  130

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  130

                                      ! do nk = 1,    3

                                      xin(130) = xin(131) + dxkl*xin(129)
                                      yin(130) = yin(131) + dykl*yin(129)
                                      zin(130) = zin(131) + dzkl*zin(129)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    2

                                      xin(132) = xin(133) + dxkl*xin(131)
                                      yin(132) = yin(133) + dykl*yin(131)
                                      zin(132) = zin(133) + dzkl*zin(131)
                                      ! i4 = i4 + lang+1 =  134

                                      ! nk =    3

                                      xin(134) = xin(135) + dxkl*xin(133)
                                      yin(134) = yin(135) + dykl*yin(133)
                                      zin(134) = zin(135) + dzkl*zin(133)
                                      ! i4 = i4 + lang+1 =  136

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  131

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  137

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  138

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  138

                                      ! do nk = 1,    3

                                      xin(138) = xin(139) + dxkl*xin(137)
                                      yin(138) = yin(139) + dykl*yin(137)
                                      zin(138) = zin(139) + dzkl*zin(137)
                                      ! i4 = i4 + lang+1 =  140

                                      ! nk =    2

                                      xin(140) = xin(141) + dxkl*xin(139)
                                      yin(140) = yin(141) + dykl*yin(139)
                                      zin(140) = zin(141) + dzkl*zin(139)
                                      ! i4 = i4 + lang+1 =  142

                                      ! nk =    3

                                      xin(142) = xin(143) + dxkl*xin(141)
                                      yin(142) = yin(143) + dykl*yin(141)
                                      zin(142) = zin(143) + dzkl*zin(141)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  139

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

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

                                      ! i1 = in(1) =  145

                                      xin(145) = 1.0_dp
                                      yin(145) = 1.0_dp
                                      zin(145) = f00

                                      ! i2 = in(2) =  169
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(169) = xc00
                                      yin(169) = yc00
                                      zin(169) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  147

                                      xin(147) = xcp00
                                      yin(147) = ycp00
                                      zin(147) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  171
                                      ! i2 =  169

                                      xin(171) = xcp00*xin(169) + cp10
                                      yin(171) = ycp00*yin(169) + cp10
                                      zin(171) = zcp00*zin(169) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  169

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  193
                                      ! i3 =  145
                                      ! i4 =  169

                                      xin(193) = c10*xin(145) + xc00*xin(169)
                                      yin(193) = c10*yin(145) + yc00*yin(169)
                                      zin(193) = c10*zin(145) + zc00*zin(169)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  195
                                      ! i5 =  193
                                      ! i4 =  169

                                      xin(195) = xcp00*xin(193) + cp10*xin(169)
                                      yin(195) = ycp00*yin(193) + cp10*yin(169)
                                      zin(195) = zcp00*zin(193) + cp10*zin(169)

                                      ! ------------------

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  193

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  201
                                      ! i3 =  169
                                      ! i4 =  193

                                      xin(201) = c10*xin(169) + xc00*xin(193)
                                      yin(201) = c10*yin(169) + yc00*yin(193)
                                      zin(201) = c10*zin(169) + zc00*zin(193)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  203
                                      ! i5 =  201
                                      ! i4 =  193

                                      xin(203) = xcp00*xin(201) + cp10*xin(193)
                                      yin(203) = ycp00*yin(201) + cp10*yin(193)
                                      zin(203) = zcp00*zin(201) + cp10*zin(193)

                                      ! ------------------

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  201

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  209
                                      ! i3 =  193
                                      ! i4 =  201

                                      xin(209) = c10*xin(193) + xc00*xin(201)
                                      yin(209) = c10*yin(193) + yc00*yin(201)
                                      zin(209) = c10*zin(193) + zc00*zin(201)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  211
                                      ! i5 =  209
                                      ! i4 =  201

                                      xin(211) = xcp00*xin(209) + cp10*xin(201)
                                      yin(211) = ycp00*yin(209) + cp10*yin(201)
                                      zin(211) = zcp00*zin(209) + cp10*zin(201)

                                      ! ------------------

                                      ! i3 = i4 =  201
                                      ! i4 = i5 =  209

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  147

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  149
                                      ! i3 =  145
                                      ! i4 =  147

                                      xin(149) = cp01*xin(145) + xcp00*xin(147)
                                      yin(149) = cp01*yin(145) + ycp00*yin(147)
                                      zin(149) = cp01*zin(145) + zcp00*zin(147)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  173

                                      xin(173) = xc00*xin(149) + c01*xin(147)
                                      yin(173) = yc00*yin(149) + c01*yin(147)
                                      zin(173) = zc00*zin(149) + c01*zin(147)

                                      ! ------------------

                                      ! i3 = i4 =  147
                                      ! i4 = i5 =  149

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  151
                                      ! i3 =  147
                                      ! i4 =  149

                                      xin(151) = cp01*xin(147) + xcp00*xin(149)
                                      yin(151) = cp01*yin(147) + ycp00*yin(149)
                                      zin(151) = cp01*zin(147) + zcp00*zin(149)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  175

                                      xin(175) = xc00*xin(151) + c01*xin(149)
                                      yin(175) = yc00*yin(151) + c01*yin(149)
                                      zin(175) = zc00*zin(151) + c01*zin(149)

                                      ! ------------------

                                      ! i3 = i4 =  149
                                      ! i4 = i5 =  151

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  152
                                      ! i3 =  149
                                      ! i4 =  151

                                      xin(152) = cp01*xin(149) + xcp00*xin(151)
                                      yin(152) = cp01*yin(149) + ycp00*yin(151)
                                      zin(152) = cp01*zin(149) + zcp00*zin(151)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  176

                                      xin(176) = xc00*xin(152) + c01*xin(151)
                                      yin(176) = yc00*yin(152) + c01*yin(151)
                                      zin(176) = zc00*zin(152) + c01*zin(151)

                                      ! ------------------

                                      ! i3 = i4 =  151
                                      ! i4 = i5 =  152

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  169

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  193

                                      xin(197) = c10*xin(149) + xc00*xin(173) + c01*xin(171)
                                      yin(197) = c10*yin(149) + yc00*yin(173) + c01*yin(171)
                                      zin(197) = c10*zin(149) + zc00*zin(173) + c01*zin(171)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  201

                                      xin(205) = c10*xin(173) + xc00*xin(197) + c01*xin(195)
                                      yin(205) = c10*yin(173) + yc00*yin(197) + c01*yin(195)
                                      zin(205) = c10*zin(173) + zc00*zin(197) + c01*zin(195)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  201

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  209

                                      xin(213) = c10*xin(197) + xc00*xin(205) + c01*xin(203)
                                      yin(213) = c10*yin(197) + yc00*yin(205) + c01*yin(203)
                                      zin(213) = c10*zin(197) + zc00*zin(205) + c01*zin(203)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  201
                                      ! i4 = i5 =  209

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  169

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  193

                                      xin(199) = c10*xin(151) + xc00*xin(175) + c01*xin(173)
                                      yin(199) = c10*yin(151) + yc00*yin(175) + c01*yin(173)
                                      zin(199) = c10*zin(151) + zc00*zin(175) + c01*zin(173)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  201

                                      xin(207) = c10*xin(175) + xc00*xin(199) + c01*xin(197)
                                      yin(207) = c10*yin(175) + yc00*yin(199) + c01*yin(197)
                                      zin(207) = c10*zin(175) + zc00*zin(199) + c01*zin(197)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  201

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  209

                                      xin(215) = c10*xin(199) + xc00*xin(207) + c01*xin(205)
                                      yin(215) = c10*yin(199) + yc00*yin(207) + c01*yin(205)
                                      zin(215) = c10*zin(199) + zc00*zin(207) + c01*zin(205)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  201
                                      ! i4 = i5 =  209

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  169

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  193

                                      xin(200) = c10*xin(152) + xc00*xin(176) + c01*xin(175)
                                      yin(200) = c10*yin(152) + yc00*yin(176) + c01*yin(175)
                                      zin(200) = c10*zin(152) + zc00*zin(176) + c01*zin(175)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  193

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  201

                                      xin(208) = c10*xin(176) + xc00*xin(200) + c01*xin(199)
                                      yin(208) = c10*yin(176) + yc00*yin(200) + c01*yin(199)
                                      zin(208) = c10*zin(176) + zc00*zin(200) + c01*zin(199)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  193
                                      ! i4 = i5 =  201

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  209

                                      xin(216) = c10*xin(200) + xc00*xin(208) + c01*xin(207)
                                      yin(216) = c10*yin(200) + yc00*yin(208) + c01*yin(207)
                                      zin(216) = c10*zin(200) + zc00*zin(208) + c01*zin(207)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  201
                                      ! i4 = i5 =  209

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  209

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  209

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  201

                                      xin(209) = xin(209) + dxij*xin(201)
                                      yin(209) = yin(209) + dyij*yin(201)
                                      zin(209) = zin(209) + dzij*zin(201)

                                      ! i3 = i4 =  201
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  193

                                      xin(201) = xin(201) + dxij*xin(193)
                                      yin(201) = yin(201) + dyij*yin(193)
                                      zin(201) = zin(201) + dzij*zin(193)

                                      ! i3 = i4 =  193
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  209

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  201

                                      xin(209) = xin(209) + dxij*xin(201)
                                      yin(209) = yin(209) + dyij*yin(201)
                                      zin(209) = zin(209) + dzij*zin(201)

                                      ! i3 = i4 =  201
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  153

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    2

                                      xin(153) = xin(169) + dxij*xin(145)
                                      yin(153) = yin(169) + dyij*yin(145)
                                      zin(153) = zin(169) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    2

                                      xin(177) = xin(193) + dxij*xin(169)
                                      yin(177) = yin(193) + dyij*yin(169)
                                      zin(177) = zin(193) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  201

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  161

                                      ! nj =    2

                                      ! i4 = i3 =  161

                                      ! do ni = 1,    2

                                      xin(161) = xin(177) + dxij*xin(153)
                                      yin(161) = yin(177) + dyij*yin(153)
                                      zin(161) = zin(177) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    2

                                      xin(185) = xin(201) + dxij*xin(177)
                                      yin(185) = yin(201) + dyij*yin(177)
                                      zin(185) = zin(201) + dzij*zin(177)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  169

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  203

                                      xin(211) = xin(211) + dxij*xin(203)
                                      yin(211) = yin(211) + dyij*yin(203)
                                      zin(211) = zin(211) + dzij*zin(203)

                                      ! i3 = i4 =  203
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  195

                                      xin(203) = xin(203) + dxij*xin(195)
                                      yin(203) = yin(203) + dyij*yin(195)
                                      zin(203) = zin(203) + dzij*zin(195)

                                      ! i3 = i4 =  195
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  203

                                      xin(211) = xin(211) + dxij*xin(203)
                                      yin(211) = yin(211) + dyij*yin(203)
                                      zin(211) = zin(211) + dzij*zin(203)

                                      ! i3 = i4 =  203
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  155

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    2

                                      xin(155) = xin(171) + dxij*xin(147)
                                      yin(155) = yin(171) + dyij*yin(147)
                                      zin(155) = zin(171) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    2

                                      xin(179) = xin(195) + dxij*xin(171)
                                      yin(179) = yin(195) + dyij*yin(171)
                                      zin(179) = zin(195) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  203

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  163

                                      ! nj =    2

                                      ! i4 = i3 =  163

                                      ! do ni = 1,    2

                                      xin(163) = xin(179) + dxij*xin(155)
                                      yin(163) = yin(179) + dyij*yin(155)
                                      zin(163) = zin(179) + dzij*zin(155)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    2

                                      xin(187) = xin(203) + dxij*xin(179)
                                      yin(187) = yin(203) + dyij*yin(179)
                                      zin(187) = zin(203) + dzij*zin(179)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  171

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(213) = xin(213) + dxij*xin(205)
                                      yin(213) = yin(213) + dyij*yin(205)
                                      zin(213) = zin(213) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  197

                                      xin(205) = xin(205) + dxij*xin(197)
                                      yin(205) = yin(205) + dyij*yin(197)
                                      zin(205) = zin(205) + dzij*zin(197)

                                      ! i3 = i4 =  197
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  213

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(213) = xin(213) + dxij*xin(205)
                                      yin(213) = yin(213) + dyij*yin(205)
                                      zin(213) = zin(213) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  157

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  157

                                      ! do ni = 1,    2

                                      xin(157) = xin(173) + dxij*xin(149)
                                      yin(157) = yin(173) + dyij*yin(149)
                                      zin(157) = zin(173) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! ni =    2

                                      xin(181) = xin(197) + dxij*xin(173)
                                      yin(181) = yin(197) + dyij*yin(173)
                                      zin(181) = zin(197) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  165

                                      ! nj =    2

                                      ! i4 = i3 =  165

                                      ! do ni = 1,    2

                                      xin(165) = xin(181) + dxij*xin(157)
                                      yin(165) = yin(181) + dyij*yin(157)
                                      zin(165) = zin(181) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    2

                                      xin(189) = xin(205) + dxij*xin(181)
                                      yin(189) = yin(205) + dyij*yin(181)
                                      zin(189) = zin(205) + dzij*zin(181)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  173

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(215) = xin(215) + dxij*xin(207)
                                      yin(215) = yin(215) + dyij*yin(207)
                                      zin(215) = zin(215) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  199

                                      xin(207) = xin(207) + dxij*xin(199)
                                      yin(207) = yin(207) + dyij*yin(199)
                                      zin(207) = zin(207) + dzij*zin(199)

                                      ! i3 = i4 =  199
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(215) = xin(215) + dxij*xin(207)
                                      yin(215) = yin(215) + dyij*yin(207)
                                      zin(215) = zin(215) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  159

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  159

                                      ! do ni = 1,    2

                                      xin(159) = xin(175) + dxij*xin(151)
                                      yin(159) = yin(175) + dyij*yin(151)
                                      zin(159) = zin(175) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  183

                                      ! ni =    2

                                      xin(183) = xin(199) + dxij*xin(175)
                                      yin(183) = yin(199) + dyij*yin(175)
                                      zin(183) = zin(199) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  167

                                      ! nj =    2

                                      ! i4 = i3 =  167

                                      ! do ni = 1,    2

                                      xin(167) = xin(183) + dxij*xin(159)
                                      yin(167) = yin(183) + dyij*yin(159)
                                      zin(167) = zin(183) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    2

                                      xin(191) = xin(207) + dxij*xin(183)
                                      yin(191) = yin(207) + dyij*yin(183)
                                      zin(191) = zin(207) + dzij*zin(183)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  175

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  208

                                      xin(216) = xin(216) + dxij*xin(208)
                                      yin(216) = yin(216) + dyij*yin(208)
                                      zin(216) = zin(216) + dzij*zin(208)

                                      ! i3 = i4 =  208
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  200

                                      xin(208) = xin(208) + dxij*xin(200)
                                      yin(208) = yin(208) + dyij*yin(200)
                                      zin(208) = zin(208) + dzij*zin(200)

                                      ! i3 = i4 =  200
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  208

                                      xin(216) = xin(216) + dxij*xin(208)
                                      yin(216) = yin(216) + dyij*yin(208)
                                      zin(216) = zin(216) + dzij*zin(208)

                                      ! i3 = i4 =  208
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  160

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  160

                                      ! do ni = 1,    2

                                      xin(160) = xin(176) + dxij*xin(152)
                                      yin(160) = yin(176) + dyij*yin(152)
                                      zin(160) = zin(176) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  184

                                      ! ni =    2

                                      xin(184) = xin(200) + dxij*xin(176)
                                      yin(184) = yin(200) + dyij*yin(176)
                                      zin(184) = zin(200) + dzij*zin(176)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  168

                                      ! nj =    2

                                      ! i4 = i3 =  168

                                      ! do ni = 1,    2

                                      xin(168) = xin(184) + dxij*xin(160)
                                      yin(168) = yin(184) + dyij*yin(160)
                                      zin(168) = zin(184) + dzij*zin(160)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    2

                                      xin(192) = xin(208) + dxij*xin(184)
                                      yin(192) = yin(208) + dyij*yin(184)
                                      zin(192) = zin(208) + dzij*zin(184)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  176

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  145

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  152

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  151

                                      xin(152) = xin(152) + dxkl*xin(151)
                                      yin(152) = yin(152) + dykl*yin(151)
                                      zin(152) = zin(152) + dzkl*zin(151)

                                      ! i3 = i4 =  151
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    3

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

                                      xin(150) = xin(151) + dxkl*xin(149)
                                      yin(150) = yin(151) + dykl*yin(149)
                                      zin(150) = zin(151) + dzkl*zin(149)
                                      ! i4 = i4 + lang+1 =  152

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  153

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  160

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  159

                                      xin(160) = xin(160) + dxkl*xin(159)
                                      yin(160) = yin(160) + dykl*yin(159)
                                      zin(160) = zin(160) + dzkl*zin(159)

                                      ! i3 = i4 =  159
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  154

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  154

                                      ! do nk = 1,    3

                                      xin(154) = xin(155) + dxkl*xin(153)
                                      yin(154) = yin(155) + dykl*yin(153)
                                      zin(154) = zin(155) + dzkl*zin(153)
                                      ! i4 = i4 + lang+1 =  156

                                      ! nk =    2

                                      xin(156) = xin(157) + dxkl*xin(155)
                                      yin(156) = yin(157) + dykl*yin(155)
                                      zin(156) = zin(157) + dzkl*zin(155)
                                      ! i4 = i4 + lang+1 =  158

                                      ! nk =    3

                                      xin(158) = xin(159) + dxkl*xin(157)
                                      yin(158) = yin(159) + dykl*yin(157)
                                      zin(158) = zin(159) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  160

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  155

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  161

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  168

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  167

                                      xin(168) = xin(168) + dxkl*xin(167)
                                      yin(168) = yin(168) + dykl*yin(167)
                                      zin(168) = zin(168) + dzkl*zin(167)

                                      ! i3 = i4 =  167
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  162

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  162

                                      ! do nk = 1,    3

                                      xin(162) = xin(163) + dxkl*xin(161)
                                      yin(162) = yin(163) + dykl*yin(161)
                                      zin(162) = zin(163) + dzkl*zin(161)
                                      ! i4 = i4 + lang+1 =  164

                                      ! nk =    2

                                      xin(164) = xin(165) + dxkl*xin(163)
                                      yin(164) = yin(165) + dykl*yin(163)
                                      zin(164) = zin(165) + dzkl*zin(163)
                                      ! i4 = i4 + lang+1 =  166

                                      ! nk =    3

                                      xin(166) = xin(167) + dxkl*xin(165)
                                      yin(166) = yin(167) + dykl*yin(165)
                                      zin(166) = zin(167) + dzkl*zin(165)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  163

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  169

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  169

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  176

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  175

                                      xin(176) = xin(176) + dxkl*xin(175)
                                      yin(176) = yin(176) + dykl*yin(175)
                                      zin(176) = zin(176) + dzkl*zin(175)

                                      ! i3 = i4 =  175
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  170

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  170

                                      ! do nk = 1,    3

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

                                      xin(174) = xin(175) + dxkl*xin(173)
                                      yin(174) = yin(175) + dykl*yin(173)
                                      zin(174) = zin(175) + dzkl*zin(173)
                                      ! i4 = i4 + lang+1 =  176

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  171

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  177

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  184

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  183

                                      xin(184) = xin(184) + dxkl*xin(183)
                                      yin(184) = yin(184) + dykl*yin(183)
                                      zin(184) = zin(184) + dzkl*zin(183)

                                      ! i3 = i4 =  183
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  178

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  178

                                      ! do nk = 1,    3

                                      xin(178) = xin(179) + dxkl*xin(177)
                                      yin(178) = yin(179) + dykl*yin(177)
                                      zin(178) = zin(179) + dzkl*zin(177)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    2

                                      xin(180) = xin(181) + dxkl*xin(179)
                                      yin(180) = yin(181) + dykl*yin(179)
                                      zin(180) = zin(181) + dzkl*zin(179)
                                      ! i4 = i4 + lang+1 =  182

                                      ! nk =    3

                                      xin(182) = xin(183) + dxkl*xin(181)
                                      yin(182) = yin(183) + dykl*yin(181)
                                      zin(182) = zin(183) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  184

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  179

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  185

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  192

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  191

                                      xin(192) = xin(192) + dxkl*xin(191)
                                      yin(192) = yin(192) + dykl*yin(191)
                                      zin(192) = zin(192) + dzkl*zin(191)

                                      ! i3 = i4 =  191
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  186

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  186

                                      ! do nk = 1,    3

                                      xin(186) = xin(187) + dxkl*xin(185)
                                      yin(186) = yin(187) + dykl*yin(185)
                                      zin(186) = zin(187) + dzkl*zin(185)
                                      ! i4 = i4 + lang+1 =  188

                                      ! nk =    2

                                      xin(188) = xin(189) + dxkl*xin(187)
                                      yin(188) = yin(189) + dykl*yin(187)
                                      zin(188) = zin(189) + dzkl*zin(187)
                                      ! i4 = i4 + lang+1 =  190

                                      ! nk =    3

                                      xin(190) = xin(191) + dxkl*xin(189)
                                      yin(190) = yin(191) + dykl*yin(189)
                                      zin(190) = zin(191) + dzkl*zin(189)
                                      ! i4 = i4 + lang+1 =  192

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  187

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  193

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  193

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  200

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  199

                                      xin(200) = xin(200) + dxkl*xin(199)
                                      yin(200) = yin(200) + dykl*yin(199)
                                      zin(200) = zin(200) + dzkl*zin(199)

                                      ! i3 = i4 =  199
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  194

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  194

                                      ! do nk = 1,    3

                                      xin(194) = xin(195) + dxkl*xin(193)
                                      yin(194) = yin(195) + dykl*yin(193)
                                      zin(194) = zin(195) + dzkl*zin(193)
                                      ! i4 = i4 + lang+1 =  196

                                      ! nk =    2

                                      xin(196) = xin(197) + dxkl*xin(195)
                                      yin(196) = yin(197) + dykl*yin(195)
                                      zin(196) = zin(197) + dzkl*zin(195)
                                      ! i4 = i4 + lang+1 =  198

                                      ! nk =    3

                                      xin(198) = xin(199) + dxkl*xin(197)
                                      yin(198) = yin(199) + dykl*yin(197)
                                      zin(198) = zin(199) + dzkl*zin(197)
                                      ! i4 = i4 + lang+1 =  200

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  195

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  201

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  208

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  207

                                      xin(208) = xin(208) + dxkl*xin(207)
                                      yin(208) = yin(208) + dykl*yin(207)
                                      zin(208) = zin(208) + dzkl*zin(207)

                                      ! i3 = i4 =  207
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  202

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  202

                                      ! do nk = 1,    3

                                      xin(202) = xin(203) + dxkl*xin(201)
                                      yin(202) = yin(203) + dykl*yin(201)
                                      zin(202) = zin(203) + dzkl*zin(201)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    2

                                      xin(204) = xin(205) + dxkl*xin(203)
                                      yin(204) = yin(205) + dykl*yin(203)
                                      zin(204) = zin(205) + dzkl*zin(203)
                                      ! i4 = i4 + lang+1 =  206

                                      ! nk =    3

                                      xin(206) = xin(207) + dxkl*xin(205)
                                      yin(206) = yin(207) + dykl*yin(205)
                                      zin(206) = zin(207) + dzkl*zin(205)
                                      ! i4 = i4 + lang+1 =  208

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  203

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  209

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  210

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  210

                                      ! do nk = 1,    3

                                      xin(210) = xin(211) + dxkl*xin(209)
                                      yin(210) = yin(211) + dykl*yin(209)
                                      zin(210) = zin(211) + dzkl*zin(209)
                                      ! i4 = i4 + lang+1 =  212

                                      ! nk =    2

                                      xin(212) = xin(213) + dxkl*xin(211)
                                      yin(212) = yin(213) + dykl*yin(211)
                                      zin(212) = zin(213) + dzkl*zin(211)
                                      ! i4 = i4 + lang+1 =  214

                                      ! nk =    3

                                      xin(214) = xin(215) + dxkl*xin(213)
                                      yin(214) = yin(215) + dykl*yin(213)
                                      zin(214) = zin(215) + dzkl*zin(213)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  211

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  216

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

                                      ! i1 = in(1) =  217

                                      xin(217) = 1.0_dp
                                      yin(217) = 1.0_dp
                                      zin(217) = f00

                                      ! i2 = in(2) =  241
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(241) = xc00
                                      yin(241) = yc00
                                      zin(241) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  219

                                      xin(219) = xcp00
                                      yin(219) = ycp00
                                      zin(219) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  243
                                      ! i2 =  241

                                      xin(243) = xcp00*xin(241) + cp10
                                      yin(243) = ycp00*yin(241) + cp10
                                      zin(243) = zcp00*zin(241) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  241

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  265
                                      ! i3 =  217
                                      ! i4 =  241

                                      xin(265) = c10*xin(217) + xc00*xin(241)
                                      yin(265) = c10*yin(217) + yc00*yin(241)
                                      zin(265) = c10*zin(217) + zc00*zin(241)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  267
                                      ! i5 =  265
                                      ! i4 =  241

                                      xin(267) = xcp00*xin(265) + cp10*xin(241)
                                      yin(267) = ycp00*yin(265) + cp10*yin(241)
                                      zin(267) = zcp00*zin(265) + cp10*zin(241)

                                      ! ------------------

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  273
                                      ! i3 =  241
                                      ! i4 =  265

                                      xin(273) = c10*xin(241) + xc00*xin(265)
                                      yin(273) = c10*yin(241) + yc00*yin(265)
                                      zin(273) = c10*zin(241) + zc00*zin(265)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  275
                                      ! i5 =  273
                                      ! i4 =  265

                                      xin(275) = xcp00*xin(273) + cp10*xin(265)
                                      yin(275) = ycp00*yin(273) + cp10*yin(265)
                                      zin(275) = zcp00*zin(273) + cp10*zin(265)

                                      ! ------------------

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  273

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  281
                                      ! i3 =  265
                                      ! i4 =  273

                                      xin(281) = c10*xin(265) + xc00*xin(273)
                                      yin(281) = c10*yin(265) + yc00*yin(273)
                                      zin(281) = c10*zin(265) + zc00*zin(273)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  283
                                      ! i5 =  281
                                      ! i4 =  273

                                      xin(283) = xcp00*xin(281) + cp10*xin(273)
                                      yin(283) = ycp00*yin(281) + cp10*yin(273)
                                      zin(283) = zcp00*zin(281) + cp10*zin(273)

                                      ! ------------------

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  281

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  217
                                      ! i4 = i1+k2 =  219

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  221
                                      ! i3 =  217
                                      ! i4 =  219

                                      xin(221) = cp01*xin(217) + xcp00*xin(219)
                                      yin(221) = cp01*yin(217) + ycp00*yin(219)
                                      zin(221) = cp01*zin(217) + zcp00*zin(219)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  245

                                      xin(245) = xc00*xin(221) + c01*xin(219)
                                      yin(245) = yc00*yin(221) + c01*yin(219)
                                      zin(245) = zc00*zin(221) + c01*zin(219)

                                      ! ------------------

                                      ! i3 = i4 =  219
                                      ! i4 = i5 =  221

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  223
                                      ! i3 =  219
                                      ! i4 =  221

                                      xin(223) = cp01*xin(219) + xcp00*xin(221)
                                      yin(223) = cp01*yin(219) + ycp00*yin(221)
                                      zin(223) = cp01*zin(219) + zcp00*zin(221)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  247

                                      xin(247) = xc00*xin(223) + c01*xin(221)
                                      yin(247) = yc00*yin(223) + c01*yin(221)
                                      zin(247) = zc00*zin(223) + c01*zin(221)

                                      ! ------------------

                                      ! i3 = i4 =  221
                                      ! i4 = i5 =  223

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  224
                                      ! i3 =  221
                                      ! i4 =  223

                                      xin(224) = cp01*xin(221) + xcp00*xin(223)
                                      yin(224) = cp01*yin(221) + ycp00*yin(223)
                                      zin(224) = cp01*zin(221) + zcp00*zin(223)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  248

                                      xin(248) = xc00*xin(224) + c01*xin(223)
                                      yin(248) = yc00*yin(224) + c01*yin(223)
                                      zin(248) = zc00*zin(224) + c01*zin(223)

                                      ! ------------------

                                      ! i3 = i4 =  223
                                      ! i4 = i5 =  224

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  241

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  265

                                      xin(269) = c10*xin(221) + xc00*xin(245) + c01*xin(243)
                                      yin(269) = c10*yin(221) + yc00*yin(245) + c01*yin(243)
                                      zin(269) = c10*zin(221) + zc00*zin(245) + c01*zin(243)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  273

                                      xin(277) = c10*xin(245) + xc00*xin(269) + c01*xin(267)
                                      yin(277) = c10*yin(245) + yc00*yin(269) + c01*yin(267)
                                      zin(277) = c10*zin(245) + zc00*zin(269) + c01*zin(267)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  273

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  281

                                      xin(285) = c10*xin(269) + xc00*xin(277) + c01*xin(275)
                                      yin(285) = c10*yin(269) + yc00*yin(277) + c01*yin(275)
                                      zin(285) = c10*zin(269) + zc00*zin(277) + c01*zin(275)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  281

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  241

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  265

                                      xin(271) = c10*xin(223) + xc00*xin(247) + c01*xin(245)
                                      yin(271) = c10*yin(223) + yc00*yin(247) + c01*yin(245)
                                      zin(271) = c10*zin(223) + zc00*zin(247) + c01*zin(245)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  273

                                      xin(279) = c10*xin(247) + xc00*xin(271) + c01*xin(269)
                                      yin(279) = c10*yin(247) + yc00*yin(271) + c01*yin(269)
                                      zin(279) = c10*zin(247) + zc00*zin(271) + c01*zin(269)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  273

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  281

                                      xin(287) = c10*xin(271) + xc00*xin(279) + c01*xin(277)
                                      yin(287) = c10*yin(271) + yc00*yin(279) + c01*yin(277)
                                      zin(287) = c10*zin(271) + zc00*zin(279) + c01*zin(277)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  281

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  241

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  265

                                      xin(272) = c10*xin(224) + xc00*xin(248) + c01*xin(247)
                                      yin(272) = c10*yin(224) + yc00*yin(248) + c01*yin(247)
                                      zin(272) = c10*zin(224) + zc00*zin(248) + c01*zin(247)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  241
                                      ! i4 = i5 =  265

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  273

                                      xin(280) = c10*xin(248) + xc00*xin(272) + c01*xin(271)
                                      yin(280) = c10*yin(248) + yc00*yin(272) + c01*yin(271)
                                      zin(280) = c10*zin(248) + zc00*zin(272) + c01*zin(271)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  265
                                      ! i4 = i5 =  273

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  281

                                      xin(288) = c10*xin(272) + xc00*xin(280) + c01*xin(279)
                                      yin(288) = c10*yin(272) + yc00*yin(280) + c01*yin(279)
                                      zin(288) = c10*zin(272) + zc00*zin(280) + c01*zin(279)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  273
                                      ! i4 = i5 =  281

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  281

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  281

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  273

                                      xin(281) = xin(281) + dxij*xin(273)
                                      yin(281) = yin(281) + dyij*yin(273)
                                      zin(281) = zin(281) + dzij*zin(273)

                                      ! i3 = i4 =  273
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  265

                                      xin(273) = xin(273) + dxij*xin(265)
                                      yin(273) = yin(273) + dyij*yin(265)
                                      zin(273) = zin(273) + dzij*zin(265)

                                      ! i3 = i4 =  265
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  281

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  273

                                      xin(281) = xin(281) + dxij*xin(273)
                                      yin(281) = yin(281) + dyij*yin(273)
                                      zin(281) = zin(281) + dzij*zin(273)

                                      ! i3 = i4 =  273
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  225

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  225

                                      ! do ni = 1,    2

                                      xin(225) = xin(241) + dxij*xin(217)
                                      yin(225) = yin(241) + dyij*yin(217)
                                      zin(225) = zin(241) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  249

                                      ! ni =    2

                                      xin(249) = xin(265) + dxij*xin(241)
                                      yin(249) = yin(265) + dyij*yin(241)
                                      zin(249) = zin(265) + dzij*zin(241)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  273

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  233

                                      ! nj =    2

                                      ! i4 = i3 =  233

                                      ! do ni = 1,    2

                                      xin(233) = xin(249) + dxij*xin(225)
                                      yin(233) = yin(249) + dyij*yin(225)
                                      zin(233) = zin(249) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  257

                                      ! ni =    2

                                      xin(257) = xin(273) + dxij*xin(249)
                                      yin(257) = yin(273) + dyij*yin(249)
                                      zin(257) = zin(273) + dzij*zin(249)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  281

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  241

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  275

                                      xin(283) = xin(283) + dxij*xin(275)
                                      yin(283) = yin(283) + dyij*yin(275)
                                      zin(283) = zin(283) + dzij*zin(275)

                                      ! i3 = i4 =  275
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  267

                                      xin(275) = xin(275) + dxij*xin(267)
                                      yin(275) = yin(275) + dyij*yin(267)
                                      zin(275) = zin(275) + dzij*zin(267)

                                      ! i3 = i4 =  267
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  283

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  275

                                      xin(283) = xin(283) + dxij*xin(275)
                                      yin(283) = yin(283) + dyij*yin(275)
                                      zin(283) = zin(283) + dzij*zin(275)

                                      ! i3 = i4 =  275
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  227

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  227

                                      ! do ni = 1,    2

                                      xin(227) = xin(243) + dxij*xin(219)
                                      yin(227) = yin(243) + dyij*yin(219)
                                      zin(227) = zin(243) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                      ! ni =    2

                                      xin(251) = xin(267) + dxij*xin(243)
                                      yin(251) = yin(267) + dyij*yin(243)
                                      zin(251) = zin(267) + dzij*zin(243)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  275

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  235

                                      ! nj =    2

                                      ! i4 = i3 =  235

                                      ! do ni = 1,    2

                                      xin(235) = xin(251) + dxij*xin(227)
                                      yin(235) = yin(251) + dyij*yin(227)
                                      zin(235) = zin(251) + dzij*zin(227)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  259

                                      ! ni =    2

                                      xin(259) = xin(275) + dxij*xin(251)
                                      yin(259) = yin(275) + dyij*yin(251)
                                      zin(259) = zin(275) + dzij*zin(251)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  243

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(285) = xin(285) + dxij*xin(277)
                                      yin(285) = yin(285) + dyij*yin(277)
                                      zin(285) = zin(285) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  269

                                      xin(277) = xin(277) + dxij*xin(269)
                                      yin(277) = yin(277) + dyij*yin(269)
                                      zin(277) = zin(277) + dzij*zin(269)

                                      ! i3 = i4 =  269
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  285

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  277

                                      xin(285) = xin(285) + dxij*xin(277)
                                      yin(285) = yin(285) + dyij*yin(277)
                                      zin(285) = zin(285) + dzij*zin(277)

                                      ! i3 = i4 =  277
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  229

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  229

                                      ! do ni = 1,    2

                                      xin(229) = xin(245) + dxij*xin(221)
                                      yin(229) = yin(245) + dyij*yin(221)
                                      zin(229) = zin(245) + dzij*zin(221)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! ni =    2

                                      xin(253) = xin(269) + dxij*xin(245)
                                      yin(253) = yin(269) + dyij*yin(245)
                                      zin(253) = zin(269) + dzij*zin(245)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  237

                                      ! nj =    2

                                      ! i4 = i3 =  237

                                      ! do ni = 1,    2

                                      xin(237) = xin(253) + dxij*xin(229)
                                      yin(237) = yin(253) + dyij*yin(229)
                                      zin(237) = zin(253) + dzij*zin(229)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  261

                                      ! ni =    2

                                      xin(261) = xin(277) + dxij*xin(253)
                                      yin(261) = yin(277) + dyij*yin(253)
                                      zin(261) = zin(277) + dzij*zin(253)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  285

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  245

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(287) = xin(287) + dxij*xin(279)
                                      yin(287) = yin(287) + dyij*yin(279)
                                      zin(287) = zin(287) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  271

                                      xin(279) = xin(279) + dxij*xin(271)
                                      yin(279) = yin(279) + dyij*yin(271)
                                      zin(279) = zin(279) + dzij*zin(271)

                                      ! i3 = i4 =  271
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  287

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  279

                                      xin(287) = xin(287) + dxij*xin(279)
                                      yin(287) = yin(287) + dyij*yin(279)
                                      zin(287) = zin(287) + dzij*zin(279)

                                      ! i3 = i4 =  279
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  231

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  231

                                      ! do ni = 1,    2

                                      xin(231) = xin(247) + dxij*xin(223)
                                      yin(231) = yin(247) + dyij*yin(223)
                                      zin(231) = zin(247) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  255

                                      ! ni =    2

                                      xin(255) = xin(271) + dxij*xin(247)
                                      yin(255) = yin(271) + dyij*yin(247)
                                      zin(255) = zin(271) + dzij*zin(247)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  279

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  239

                                      ! nj =    2

                                      ! i4 = i3 =  239

                                      ! do ni = 1,    2

                                      xin(239) = xin(255) + dxij*xin(231)
                                      yin(239) = yin(255) + dyij*yin(231)
                                      zin(239) = zin(255) + dzij*zin(231)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  263

                                      ! ni =    2

                                      xin(263) = xin(279) + dxij*xin(255)
                                      yin(263) = yin(279) + dyij*yin(255)
                                      zin(263) = zin(279) + dzij*zin(255)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  247

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  280

                                      xin(288) = xin(288) + dxij*xin(280)
                                      yin(288) = yin(288) + dyij*yin(280)
                                      zin(288) = zin(288) + dzij*zin(280)

                                      ! i3 = i4 =  280
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  272

                                      xin(280) = xin(280) + dxij*xin(272)
                                      yin(280) = yin(280) + dyij*yin(272)
                                      zin(280) = zin(280) + dzij*zin(272)

                                      ! i3 = i4 =  272
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  288

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  280

                                      xin(288) = xin(288) + dxij*xin(280)
                                      yin(288) = yin(288) + dyij*yin(280)
                                      zin(288) = zin(288) + dzij*zin(280)

                                      ! i3 = i4 =  280
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  232

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  232

                                      ! do ni = 1,    2

                                      xin(232) = xin(248) + dxij*xin(224)
                                      yin(232) = yin(248) + dyij*yin(224)
                                      zin(232) = zin(248) + dzij*zin(224)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  256

                                      ! ni =    2

                                      xin(256) = xin(272) + dxij*xin(248)
                                      yin(256) = yin(272) + dyij*yin(248)
                                      zin(256) = zin(272) + dzij*zin(248)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  280

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  240

                                      ! nj =    2

                                      ! i4 = i3 =  240

                                      ! do ni = 1,    2

                                      xin(240) = xin(256) + dxij*xin(232)
                                      yin(240) = yin(256) + dyij*yin(232)
                                      zin(240) = zin(256) + dzij*zin(232)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  264

                                      ! ni =    2

                                      xin(264) = xin(280) + dxij*xin(256)
                                      yin(264) = yin(280) + dyij*yin(256)
                                      zin(264) = zin(280) + dzij*zin(256)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  248

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  217

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  224

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  223

                                      xin(224) = xin(224) + dxkl*xin(223)
                                      yin(224) = yin(224) + dykl*yin(223)
                                      zin(224) = zin(224) + dzkl*zin(223)

                                      ! i3 = i4 =  223
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    3

                                      xin(218) = xin(219) + dxkl*xin(217)
                                      yin(218) = yin(219) + dykl*yin(217)
                                      zin(218) = zin(219) + dzkl*zin(217)
                                      ! i4 = i4 + lang+1 =  220

                                      ! nk =    2

                                      xin(220) = xin(221) + dxkl*xin(219)
                                      yin(220) = yin(221) + dykl*yin(219)
                                      zin(220) = zin(221) + dzkl*zin(219)
                                      ! i4 = i4 + lang+1 =  222

                                      ! nk =    3

                                      xin(222) = xin(223) + dxkl*xin(221)
                                      yin(222) = yin(223) + dykl*yin(221)
                                      zin(222) = zin(223) + dzkl*zin(221)
                                      ! i4 = i4 + lang+1 =  224

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  225

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  232

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  231

                                      xin(232) = xin(232) + dxkl*xin(231)
                                      yin(232) = yin(232) + dykl*yin(231)
                                      zin(232) = zin(232) + dzkl*zin(231)

                                      ! i3 = i4 =  231
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  226

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  226

                                      ! do nk = 1,    3

                                      xin(226) = xin(227) + dxkl*xin(225)
                                      yin(226) = yin(227) + dykl*yin(225)
                                      zin(226) = zin(227) + dzkl*zin(225)
                                      ! i4 = i4 + lang+1 =  228

                                      ! nk =    2

                                      xin(228) = xin(229) + dxkl*xin(227)
                                      yin(228) = yin(229) + dykl*yin(227)
                                      zin(228) = zin(229) + dzkl*zin(227)
                                      ! i4 = i4 + lang+1 =  230

                                      ! nk =    3

                                      xin(230) = xin(231) + dxkl*xin(229)
                                      yin(230) = yin(231) + dykl*yin(229)
                                      zin(230) = zin(231) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  232

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  227

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  233

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  240

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  239

                                      xin(240) = xin(240) + dxkl*xin(239)
                                      yin(240) = yin(240) + dykl*yin(239)
                                      zin(240) = zin(240) + dzkl*zin(239)

                                      ! i3 = i4 =  239
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  234

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  234

                                      ! do nk = 1,    3

                                      xin(234) = xin(235) + dxkl*xin(233)
                                      yin(234) = yin(235) + dykl*yin(233)
                                      zin(234) = zin(235) + dzkl*zin(233)
                                      ! i4 = i4 + lang+1 =  236

                                      ! nk =    2

                                      xin(236) = xin(237) + dxkl*xin(235)
                                      yin(236) = yin(237) + dykl*yin(235)
                                      zin(236) = zin(237) + dzkl*zin(235)
                                      ! i4 = i4 + lang+1 =  238

                                      ! nk =    3

                                      xin(238) = xin(239) + dxkl*xin(237)
                                      yin(238) = yin(239) + dykl*yin(237)
                                      zin(238) = zin(239) + dzkl*zin(237)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  235

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  241

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  241

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  248

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  247

                                      xin(248) = xin(248) + dxkl*xin(247)
                                      yin(248) = yin(248) + dykl*yin(247)
                                      zin(248) = zin(248) + dzkl*zin(247)

                                      ! i3 = i4 =  247
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  242

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  242

                                      ! do nk = 1,    3

                                      xin(242) = xin(243) + dxkl*xin(241)
                                      yin(242) = yin(243) + dykl*yin(241)
                                      zin(242) = zin(243) + dzkl*zin(241)
                                      ! i4 = i4 + lang+1 =  244

                                      ! nk =    2

                                      xin(244) = xin(245) + dxkl*xin(243)
                                      yin(244) = yin(245) + dykl*yin(243)
                                      zin(244) = zin(245) + dzkl*zin(243)
                                      ! i4 = i4 + lang+1 =  246

                                      ! nk =    3

                                      xin(246) = xin(247) + dxkl*xin(245)
                                      yin(246) = yin(247) + dykl*yin(245)
                                      zin(246) = zin(247) + dzkl*zin(245)
                                      ! i4 = i4 + lang+1 =  248

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  243

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  249

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  256

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  255

                                      xin(256) = xin(256) + dxkl*xin(255)
                                      yin(256) = yin(256) + dykl*yin(255)
                                      zin(256) = zin(256) + dzkl*zin(255)

                                      ! i3 = i4 =  255
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  250

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  250

                                      ! do nk = 1,    3

                                      xin(250) = xin(251) + dxkl*xin(249)
                                      yin(250) = yin(251) + dykl*yin(249)
                                      zin(250) = zin(251) + dzkl*zin(249)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    2

                                      xin(252) = xin(253) + dxkl*xin(251)
                                      yin(252) = yin(253) + dykl*yin(251)
                                      zin(252) = zin(253) + dzkl*zin(251)
                                      ! i4 = i4 + lang+1 =  254

                                      ! nk =    3

                                      xin(254) = xin(255) + dxkl*xin(253)
                                      yin(254) = yin(255) + dykl*yin(253)
                                      zin(254) = zin(255) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  256

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  251

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  257

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  264

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  263

                                      xin(264) = xin(264) + dxkl*xin(263)
                                      yin(264) = yin(264) + dykl*yin(263)
                                      zin(264) = zin(264) + dzkl*zin(263)

                                      ! i3 = i4 =  263
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  258

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  258

                                      ! do nk = 1,    3

                                      xin(258) = xin(259) + dxkl*xin(257)
                                      yin(258) = yin(259) + dykl*yin(257)
                                      zin(258) = zin(259) + dzkl*zin(257)
                                      ! i4 = i4 + lang+1 =  260

                                      ! nk =    2

                                      xin(260) = xin(261) + dxkl*xin(259)
                                      yin(260) = yin(261) + dykl*yin(259)
                                      zin(260) = zin(261) + dzkl*zin(259)
                                      ! i4 = i4 + lang+1 =  262

                                      ! nk =    3

                                      xin(262) = xin(263) + dxkl*xin(261)
                                      yin(262) = yin(263) + dykl*yin(261)
                                      zin(262) = zin(263) + dzkl*zin(261)
                                      ! i4 = i4 + lang+1 =  264

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  259

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  265

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  265

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  272

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  271

                                      xin(272) = xin(272) + dxkl*xin(271)
                                      yin(272) = yin(272) + dykl*yin(271)
                                      zin(272) = zin(272) + dzkl*zin(271)

                                      ! i3 = i4 =  271
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  266

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  266

                                      ! do nk = 1,    3

                                      xin(266) = xin(267) + dxkl*xin(265)
                                      yin(266) = yin(267) + dykl*yin(265)
                                      zin(266) = zin(267) + dzkl*zin(265)
                                      ! i4 = i4 + lang+1 =  268

                                      ! nk =    2

                                      xin(268) = xin(269) + dxkl*xin(267)
                                      yin(268) = yin(269) + dykl*yin(267)
                                      zin(268) = zin(269) + dzkl*zin(267)
                                      ! i4 = i4 + lang+1 =  270

                                      ! nk =    3

                                      xin(270) = xin(271) + dxkl*xin(269)
                                      yin(270) = yin(271) + dykl*yin(269)
                                      zin(270) = zin(271) + dzkl*zin(269)
                                      ! i4 = i4 + lang+1 =  272

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  267

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  273

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  280

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  279

                                      xin(280) = xin(280) + dxkl*xin(279)
                                      yin(280) = yin(280) + dykl*yin(279)
                                      zin(280) = zin(280) + dzkl*zin(279)

                                      ! i3 = i4 =  279
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  274

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  274

                                      ! do nk = 1,    3

                                      xin(274) = xin(275) + dxkl*xin(273)
                                      yin(274) = yin(275) + dykl*yin(273)
                                      zin(274) = zin(275) + dzkl*zin(273)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    2

                                      xin(276) = xin(277) + dxkl*xin(275)
                                      yin(276) = yin(277) + dykl*yin(275)
                                      zin(276) = zin(277) + dzkl*zin(275)
                                      ! i4 = i4 + lang+1 =  278

                                      ! nk =    3

                                      xin(278) = xin(279) + dxkl*xin(277)
                                      yin(278) = yin(279) + dykl*yin(277)
                                      zin(278) = zin(279) + dzkl*zin(277)
                                      ! i4 = i4 + lang+1 =  280

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  275

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  281

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  282

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  282

                                      ! do nk = 1,    3

                                      xin(282) = xin(283) + dxkl*xin(281)
                                      yin(282) = yin(283) + dykl*yin(281)
                                      zin(282) = zin(283) + dzkl*zin(281)
                                      ! i4 = i4 + lang+1 =  284

                                      ! nk =    2

                                      xin(284) = xin(285) + dxkl*xin(283)
                                      yin(284) = yin(285) + dykl*yin(283)
                                      zin(284) = zin(285) + dzkl*zin(283)
                                      ! i4 = i4 + lang+1 =  286

                                      ! nk =    3

                                      xin(286) = xin(287) + dxkl*xin(285)
                                      yin(286) = yin(287) + dykl*yin(285)
                                      zin(286) = zin(287) + dzkl*zin(285)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  283

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  288

                                      u2 = roots(5)*rho
                                      f00 = expe*wghts(5)

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

                                      ! i1 = in(1) =  289

                                      xin(289) = 1.0_dp
                                      yin(289) = 1.0_dp
                                      zin(289) = f00

                                      ! i2 = in(2) =  313
                                      ! k2 = kn(2) =    2
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(313) = xc00
                                      yin(313) = yc00
                                      zin(313) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  291

                                      xin(291) = xcp00
                                      yin(291) = ycp00
                                      zin(291) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  315
                                      ! i2 =  313

                                      xin(315) = xcp00*xin(313) + cp10
                                      yin(315) = ycp00*yin(313) + cp10
                                      zin(315) = zcp00*zin(313) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      ! do n = 2,   4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  337
                                      ! i3 =  289
                                      ! i4 =  313

                                      xin(337) = c10*xin(289) + xc00*xin(313)
                                      yin(337) = c10*yin(289) + yc00*yin(313)
                                      zin(337) = c10*zin(289) + zc00*zin(313)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  339
                                      ! i5 =  337
                                      ! i4 =  313

                                      xin(339) = xcp00*xin(337) + cp10*xin(313)
                                      yin(339) = ycp00*yin(337) + cp10*yin(313)
                                      zin(339) = zcp00*zin(337) + cp10*zin(313)

                                      ! ------------------

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  345
                                      ! i3 =  313
                                      ! i4 =  337

                                      xin(345) = c10*xin(313) + xc00*xin(337)
                                      yin(345) = c10*yin(313) + yc00*yin(337)
                                      zin(345) = c10*zin(313) + zc00*zin(337)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  347
                                      ! i5 =  345
                                      ! i4 =  337

                                      xin(347) = xcp00*xin(345) + cp10*xin(337)
                                      yin(347) = ycp00*yin(345) + cp10*yin(337)
                                      zin(347) = zcp00*zin(345) + cp10*zin(337)

                                      ! ------------------

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  345

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  353
                                      ! i3 =  337
                                      ! i4 =  345

                                      xin(353) = c10*xin(337) + xc00*xin(345)
                                      yin(353) = c10*yin(337) + yc00*yin(345)
                                      zin(353) = c10*zin(337) + zc00*zin(345)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  355
                                      ! i5 =  353
                                      ! i4 =  345

                                      xin(355) = xcp00*xin(353) + cp10*xin(345)
                                      yin(355) = ycp00*yin(353) + cp10*yin(345)
                                      zin(355) = zcp00*zin(353) + cp10*zin(345)

                                      ! ------------------

                                      ! i3 = i4 =  345
                                      ! i4 = i5 =  353

                                      ! n =    5

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  289
                                      ! i4 = i1+k2 =  291

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  293
                                      ! i3 =  289
                                      ! i4 =  291

                                      xin(293) = cp01*xin(289) + xcp00*xin(291)
                                      yin(293) = cp01*yin(289) + ycp00*yin(291)
                                      zin(293) = cp01*zin(289) + zcp00*zin(291)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  317

                                      xin(317) = xc00*xin(293) + c01*xin(291)
                                      yin(317) = yc00*yin(293) + c01*yin(291)
                                      zin(317) = zc00*zin(293) + c01*zin(291)

                                      ! ------------------

                                      ! i3 = i4 =  291
                                      ! i4 = i5 =  293

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  295
                                      ! i3 =  291
                                      ! i4 =  293

                                      xin(295) = cp01*xin(291) + xcp00*xin(293)
                                      yin(295) = cp01*yin(291) + ycp00*yin(293)
                                      zin(295) = cp01*zin(291) + zcp00*zin(293)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  319

                                      xin(319) = xc00*xin(295) + c01*xin(293)
                                      yin(319) = yc00*yin(295) + c01*yin(293)
                                      zin(319) = zc00*zin(295) + c01*zin(293)

                                      ! ------------------

                                      ! i3 = i4 =  293
                                      ! i4 = i5 =  295

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  296
                                      ! i3 =  293
                                      ! i4 =  295

                                      xin(296) = cp01*xin(293) + xcp00*xin(295)
                                      yin(296) = cp01*yin(293) + ycp00*yin(295)
                                      zin(296) = cp01*zin(293) + zcp00*zin(295)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  320

                                      xin(320) = xc00*xin(296) + c01*xin(295)
                                      yin(320) = yc00*yin(296) + c01*yin(295)
                                      zin(320) = zc00*zin(296) + c01*zin(295)

                                      ! ------------------

                                      ! i3 = i4 =  295
                                      ! i4 = i5 =  296

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    2

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    4
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  337

                                      xin(341) = c10*xin(293) + xc00*xin(317) + c01*xin(315)
                                      yin(341) = c10*yin(293) + yc00*yin(317) + c01*yin(315)
                                      zin(341) = c10*zin(293) + zc00*zin(317) + c01*zin(315)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  345

                                      xin(349) = c10*xin(317) + xc00*xin(341) + c01*xin(339)
                                      yin(349) = c10*yin(317) + yc00*yin(341) + c01*yin(339)
                                      zin(349) = c10*zin(317) + zc00*zin(341) + c01*zin(339)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  345

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  353

                                      xin(357) = c10*xin(341) + xc00*xin(349) + c01*xin(347)
                                      yin(357) = c10*yin(341) + yc00*yin(349) + c01*yin(347)
                                      zin(357) = c10*zin(341) + zc00*zin(349) + c01*zin(347)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  345
                                      ! i4 = i5 =  353

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   4

                                      ! n =    3

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  337

                                      xin(343) = c10*xin(295) + xc00*xin(319) + c01*xin(317)
                                      yin(343) = c10*yin(295) + yc00*yin(319) + c01*yin(317)
                                      zin(343) = c10*zin(295) + zc00*zin(319) + c01*zin(317)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  345

                                      xin(351) = c10*xin(319) + xc00*xin(343) + c01*xin(341)
                                      yin(351) = c10*yin(319) + yc00*yin(343) + c01*yin(341)
                                      zin(351) = c10*zin(319) + zc00*zin(343) + c01*zin(341)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  345

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  353

                                      xin(359) = c10*xin(343) + xc00*xin(351) + c01*xin(349)
                                      yin(359) = c10*yin(343) + yc00*yin(351) + c01*yin(349)
                                      zin(359) = c10*zin(343) + zc00*zin(351) + c01*zin(349)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  345
                                      ! i4 = i5 =  353

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    4

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  289
                                      ! i4 = i2 =  313

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    4

                                      ! i5 = in(nn+1) =  337

                                      xin(344) = c10*xin(296) + xc00*xin(320) + c01*xin(319)
                                      yin(344) = c10*yin(296) + yc00*yin(320) + c01*yin(319)
                                      zin(344) = c10*zin(296) + zc00*zin(320) + c01*zin(319)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  313
                                      ! i4 = i5 =  337

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  345

                                      xin(352) = c10*xin(320) + xc00*xin(344) + c01*xin(343)
                                      yin(352) = c10*yin(320) + yc00*yin(344) + c01*yin(343)
                                      zin(352) = c10*zin(320) + zc00*zin(344) + c01*zin(343)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  337
                                      ! i4 = i5 =  345

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  353

                                      xin(360) = c10*xin(344) + xc00*xin(352) + c01*xin(351)
                                      yin(360) = c10*yin(344) + yc00*yin(352) + c01*yin(351)
                                      zin(360) = c10*zin(344) + zc00*zin(352) + c01*zin(351)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  345
                                      ! i4 = i5 =  353

                                      ! nn =    5

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  353

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  353

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  345

                                      xin(353) = xin(353) + dxij*xin(345)
                                      yin(353) = yin(353) + dyij*yin(345)
                                      zin(353) = zin(353) + dzij*zin(345)

                                      ! i3 = i4 =  345
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  337

                                      xin(345) = xin(345) + dxij*xin(337)
                                      yin(345) = yin(345) + dyij*yin(337)
                                      zin(345) = zin(345) + dzij*zin(337)

                                      ! i3 = i4 =  337
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  353

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  345

                                      xin(353) = xin(353) + dxij*xin(345)
                                      yin(353) = yin(353) + dyij*yin(345)
                                      zin(353) = zin(353) + dzij*zin(345)

                                      ! i3 = i4 =  345
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  297

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  297

                                      ! do ni = 1,    2

                                      xin(297) = xin(313) + dxij*xin(289)
                                      yin(297) = yin(313) + dyij*yin(289)
                                      zin(297) = zin(313) + dzij*zin(289)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  321

                                      ! ni =    2

                                      xin(321) = xin(337) + dxij*xin(313)
                                      yin(321) = yin(337) + dyij*yin(313)
                                      zin(321) = zin(337) + dzij*zin(313)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  345

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  305

                                      ! nj =    2

                                      ! i4 = i3 =  305

                                      ! do ni = 1,    2

                                      xin(305) = xin(321) + dxij*xin(297)
                                      yin(305) = yin(321) + dyij*yin(297)
                                      zin(305) = zin(321) + dzij*zin(297)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  329

                                      ! ni =    2

                                      xin(329) = xin(345) + dxij*xin(321)
                                      yin(329) = yin(345) + dyij*yin(321)
                                      zin(329) = zin(345) + dzij*zin(321)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  353

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  313

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  355

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  347

                                      xin(355) = xin(355) + dxij*xin(347)
                                      yin(355) = yin(355) + dyij*yin(347)
                                      zin(355) = zin(355) + dzij*zin(347)

                                      ! i3 = i4 =  347
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  339

                                      xin(347) = xin(347) + dxij*xin(339)
                                      yin(347) = yin(347) + dyij*yin(339)
                                      zin(347) = zin(347) + dzij*zin(339)

                                      ! i3 = i4 =  339
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  355

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  347

                                      xin(355) = xin(355) + dxij*xin(347)
                                      yin(355) = yin(355) + dyij*yin(347)
                                      zin(355) = zin(355) + dzij*zin(347)

                                      ! i3 = i4 =  347
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  299

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  299

                                      ! do ni = 1,    2

                                      xin(299) = xin(315) + dxij*xin(291)
                                      yin(299) = yin(315) + dyij*yin(291)
                                      zin(299) = zin(315) + dzij*zin(291)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  323

                                      ! ni =    2

                                      xin(323) = xin(339) + dxij*xin(315)
                                      yin(323) = yin(339) + dyij*yin(315)
                                      zin(323) = zin(339) + dzij*zin(315)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  347

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  307

                                      ! nj =    2

                                      ! i4 = i3 =  307

                                      ! do ni = 1,    2

                                      xin(307) = xin(323) + dxij*xin(299)
                                      yin(307) = yin(323) + dyij*yin(299)
                                      zin(307) = zin(323) + dzij*zin(299)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  331

                                      ! ni =    2

                                      xin(331) = xin(347) + dxij*xin(323)
                                      yin(331) = yin(347) + dyij*yin(323)
                                      zin(331) = zin(347) + dzij*zin(323)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  355

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  315

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    4

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  357

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  349

                                      xin(357) = xin(357) + dxij*xin(349)
                                      yin(357) = yin(357) + dyij*yin(349)
                                      zin(357) = zin(357) + dzij*zin(349)

                                      ! i3 = i4 =  349
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  341

                                      xin(349) = xin(349) + dxij*xin(341)
                                      yin(349) = yin(349) + dyij*yin(341)
                                      zin(349) = zin(349) + dzij*zin(341)

                                      ! i3 = i4 =  341
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  357

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  349

                                      xin(357) = xin(357) + dxij*xin(349)
                                      yin(357) = yin(357) + dyij*yin(349)
                                      zin(357) = zin(357) + dzij*zin(349)

                                      ! i3 = i4 =  349
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  301

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  301

                                      ! do ni = 1,    2

                                      xin(301) = xin(317) + dxij*xin(293)
                                      yin(301) = yin(317) + dyij*yin(293)
                                      zin(301) = zin(317) + dzij*zin(293)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  325

                                      ! ni =    2

                                      xin(325) = xin(341) + dxij*xin(317)
                                      yin(325) = yin(341) + dyij*yin(317)
                                      zin(325) = zin(341) + dzij*zin(317)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  309

                                      ! nj =    2

                                      ! i4 = i3 =  309

                                      ! do ni = 1,    2

                                      xin(309) = xin(325) + dxij*xin(301)
                                      yin(309) = yin(325) + dyij*yin(301)
                                      zin(309) = zin(325) + dzij*zin(301)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  333

                                      ! ni =    2

                                      xin(333) = xin(349) + dxij*xin(325)
                                      yin(333) = yin(349) + dyij*yin(325)
                                      zin(333) = zin(349) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  357

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  317

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  359

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  351

                                      xin(359) = xin(359) + dxij*xin(351)
                                      yin(359) = yin(359) + dyij*yin(351)
                                      zin(359) = zin(359) + dzij*zin(351)

                                      ! i3 = i4 =  351
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  343

                                      xin(351) = xin(351) + dxij*xin(343)
                                      yin(351) = yin(351) + dyij*yin(343)
                                      zin(351) = zin(351) + dzij*zin(343)

                                      ! i3 = i4 =  343
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  359

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  351

                                      xin(359) = xin(359) + dxij*xin(351)
                                      yin(359) = yin(359) + dyij*yin(351)
                                      zin(359) = zin(359) + dzij*zin(351)

                                      ! i3 = i4 =  351
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  303

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  303

                                      ! do ni = 1,    2

                                      xin(303) = xin(319) + dxij*xin(295)
                                      yin(303) = yin(319) + dyij*yin(295)
                                      zin(303) = zin(319) + dzij*zin(295)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  327

                                      ! ni =    2

                                      xin(327) = xin(343) + dxij*xin(319)
                                      yin(327) = yin(343) + dyij*yin(319)
                                      zin(327) = zin(343) + dzij*zin(319)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  351

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  311

                                      ! nj =    2

                                      ! i4 = i3 =  311

                                      ! do ni = 1,    2

                                      xin(311) = xin(327) + dxij*xin(303)
                                      yin(311) = yin(327) + dyij*yin(303)
                                      zin(311) = zin(327) + dzij*zin(303)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  335

                                      ! ni =    2

                                      xin(335) = xin(351) + dxij*xin(327)
                                      yin(335) = yin(351) + dyij*yin(327)
                                      zin(335) = zin(351) + dzij*zin(327)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  359

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  319

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  360

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  352

                                      xin(360) = xin(360) + dxij*xin(352)
                                      yin(360) = yin(360) + dyij*yin(352)
                                      zin(360) = zin(360) + dzij*zin(352)

                                      ! i3 = i4 =  352
                                      ! nn = nn-1 =    3

                                      ! i4 = in(nn)+km =  344

                                      xin(352) = xin(352) + dxij*xin(344)
                                      yin(352) = yin(352) + dyij*yin(344)
                                      zin(352) = zin(352) + dzij*zin(344)

                                      ! i3 = i4 =  344
                                      ! nn = nn-1 =    2

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    4

                                      ! i3 = i5 + km =  360

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  352

                                      xin(360) = xin(360) + dxij*xin(352)
                                      yin(360) = yin(360) + dyij*yin(352)
                                      zin(360) = zin(360) + dzij*zin(352)

                                      ! i3 = i4 =  352
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  304

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  304

                                      ! do ni = 1,    2

                                      xin(304) = xin(320) + dxij*xin(296)
                                      yin(304) = yin(320) + dyij*yin(296)
                                      zin(304) = zin(320) + dzij*zin(296)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  328

                                      ! ni =    2

                                      xin(328) = xin(344) + dxij*xin(320)
                                      yin(328) = yin(344) + dyij*yin(320)
                                      zin(328) = zin(344) + dzij*zin(320)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  352

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  312

                                      ! nj =    2

                                      ! i4 = i3 =  312

                                      ! do ni = 1,    2

                                      xin(312) = xin(328) + dxij*xin(304)
                                      yin(312) = yin(328) + dyij*yin(304)
                                      zin(312) = zin(328) + dzij*zin(304)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  336

                                      ! ni =    2

                                      xin(336) = xin(352) + dxij*xin(328)
                                      yin(336) = yin(352) + dyij*yin(328)
                                      zin(336) = zin(352) + dzij*zin(328)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  360

                                      ! ni =    3

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  320

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    7

                                      ! iaa = i1 =  289

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  296

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  295

                                      xin(296) = xin(296) + dxkl*xin(295)
                                      yin(296) = yin(296) + dykl*yin(295)
                                      zin(296) = zin(296) + dzkl*zin(295)

                                      ! i3 = i4 =  295
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    3

                                      xin(290) = xin(291) + dxkl*xin(289)
                                      yin(290) = yin(291) + dykl*yin(289)
                                      zin(290) = zin(291) + dzkl*zin(289)
                                      ! i4 = i4 + lang+1 =  292

                                      ! nk =    2

                                      xin(292) = xin(293) + dxkl*xin(291)
                                      yin(292) = yin(293) + dykl*yin(291)
                                      zin(292) = zin(293) + dzkl*zin(291)
                                      ! i4 = i4 + lang+1 =  294

                                      ! nk =    3

                                      xin(294) = xin(295) + dxkl*xin(293)
                                      yin(294) = yin(295) + dykl*yin(293)
                                      zin(294) = zin(295) + dzkl*zin(293)
                                      ! i4 = i4 + lang+1 =  296

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  297

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  304

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  303

                                      xin(304) = xin(304) + dxkl*xin(303)
                                      yin(304) = yin(304) + dykl*yin(303)
                                      zin(304) = zin(304) + dzkl*zin(303)

                                      ! i3 = i4 =  303
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  298

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  298

                                      ! do nk = 1,    3

                                      xin(298) = xin(299) + dxkl*xin(297)
                                      yin(298) = yin(299) + dykl*yin(297)
                                      zin(298) = zin(299) + dzkl*zin(297)
                                      ! i4 = i4 + lang+1 =  300

                                      ! nk =    2

                                      xin(300) = xin(301) + dxkl*xin(299)
                                      yin(300) = yin(301) + dykl*yin(299)
                                      zin(300) = zin(301) + dzkl*zin(299)
                                      ! i4 = i4 + lang+1 =  302

                                      ! nk =    3

                                      xin(302) = xin(303) + dxkl*xin(301)
                                      yin(302) = yin(303) + dykl*yin(301)
                                      zin(302) = zin(303) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  304

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  299

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  305

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  312

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  311

                                      xin(312) = xin(312) + dxkl*xin(311)
                                      yin(312) = yin(312) + dykl*yin(311)
                                      zin(312) = zin(312) + dzkl*zin(311)

                                      ! i3 = i4 =  311
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  306

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  306

                                      ! do nk = 1,    3

                                      xin(306) = xin(307) + dxkl*xin(305)
                                      yin(306) = yin(307) + dykl*yin(305)
                                      zin(306) = zin(307) + dzkl*zin(305)
                                      ! i4 = i4 + lang+1 =  308

                                      ! nk =    2

                                      xin(308) = xin(309) + dxkl*xin(307)
                                      yin(308) = yin(309) + dykl*yin(307)
                                      zin(308) = zin(309) + dzkl*zin(307)
                                      ! i4 = i4 + lang+1 =  310

                                      ! nk =    3

                                      xin(310) = xin(311) + dxkl*xin(309)
                                      yin(310) = yin(311) + dykl*yin(309)
                                      zin(310) = zin(311) + dzkl*zin(309)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  307

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  313

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  320

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  319

                                      xin(320) = xin(320) + dxkl*xin(319)
                                      yin(320) = yin(320) + dykl*yin(319)
                                      zin(320) = zin(320) + dzkl*zin(319)

                                      ! i3 = i4 =  319
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  314

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  314

                                      ! do nk = 1,    3

                                      xin(314) = xin(315) + dxkl*xin(313)
                                      yin(314) = yin(315) + dykl*yin(313)
                                      zin(314) = zin(315) + dzkl*zin(313)
                                      ! i4 = i4 + lang+1 =  316

                                      ! nk =    2

                                      xin(316) = xin(317) + dxkl*xin(315)
                                      yin(316) = yin(317) + dykl*yin(315)
                                      zin(316) = zin(317) + dzkl*zin(315)
                                      ! i4 = i4 + lang+1 =  318

                                      ! nk =    3

                                      xin(318) = xin(319) + dxkl*xin(317)
                                      yin(318) = yin(319) + dykl*yin(317)
                                      zin(318) = zin(319) + dzkl*zin(317)
                                      ! i4 = i4 + lang+1 =  320

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  315

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  321

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  328

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  327

                                      xin(328) = xin(328) + dxkl*xin(327)
                                      yin(328) = yin(328) + dykl*yin(327)
                                      zin(328) = zin(328) + dzkl*zin(327)

                                      ! i3 = i4 =  327
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  322

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  322

                                      ! do nk = 1,    3

                                      xin(322) = xin(323) + dxkl*xin(321)
                                      yin(322) = yin(323) + dykl*yin(321)
                                      zin(322) = zin(323) + dzkl*zin(321)
                                      ! i4 = i4 + lang+1 =  324

                                      ! nk =    2

                                      xin(324) = xin(325) + dxkl*xin(323)
                                      yin(324) = yin(325) + dykl*yin(323)
                                      zin(324) = zin(325) + dzkl*zin(323)
                                      ! i4 = i4 + lang+1 =  326

                                      ! nk =    3

                                      xin(326) = xin(327) + dxkl*xin(325)
                                      yin(326) = yin(327) + dykl*yin(325)
                                      zin(326) = zin(327) + dzkl*zin(325)
                                      ! i4 = i4 + lang+1 =  328

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  323

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  329

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  336

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  335

                                      xin(336) = xin(336) + dxkl*xin(335)
                                      yin(336) = yin(336) + dykl*yin(335)
                                      zin(336) = zin(336) + dzkl*zin(335)

                                      ! i3 = i4 =  335
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  330

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  330

                                      ! do nk = 1,    3

                                      xin(330) = xin(331) + dxkl*xin(329)
                                      yin(330) = yin(331) + dykl*yin(329)
                                      zin(330) = zin(331) + dzkl*zin(329)
                                      ! i4 = i4 + lang+1 =  332

                                      ! nk =    2

                                      xin(332) = xin(333) + dxkl*xin(331)
                                      yin(332) = yin(333) + dykl*yin(331)
                                      zin(332) = zin(333) + dzkl*zin(331)
                                      ! i4 = i4 + lang+1 =  334

                                      ! nk =    3

                                      xin(334) = xin(335) + dxkl*xin(333)
                                      yin(334) = yin(335) + dykl*yin(333)
                                      zin(334) = zin(335) + dzkl*zin(333)
                                      ! i4 = i4 + lang+1 =  336

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  331

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  337

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  337

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  344

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  343

                                      xin(344) = xin(344) + dxkl*xin(343)
                                      yin(344) = yin(344) + dykl*yin(343)
                                      zin(344) = zin(344) + dzkl*zin(343)

                                      ! i3 = i4 =  343
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  338

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  338

                                      ! do nk = 1,    3

                                      xin(338) = xin(339) + dxkl*xin(337)
                                      yin(338) = yin(339) + dykl*yin(337)
                                      zin(338) = zin(339) + dzkl*zin(337)
                                      ! i4 = i4 + lang+1 =  340

                                      ! nk =    2

                                      xin(340) = xin(341) + dxkl*xin(339)
                                      yin(340) = yin(341) + dykl*yin(339)
                                      zin(340) = zin(341) + dzkl*zin(339)
                                      ! i4 = i4 + lang+1 =  342

                                      ! nk =    3

                                      xin(342) = xin(343) + dxkl*xin(341)
                                      yin(342) = yin(343) + dykl*yin(341)
                                      zin(342) = zin(343) + dzkl*zin(341)
                                      ! i4 = i4 + lang+1 =  344

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  339

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  345

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  352

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  351

                                      xin(352) = xin(352) + dxkl*xin(351)
                                      yin(352) = yin(352) + dykl*yin(351)
                                      zin(352) = zin(352) + dzkl*zin(351)

                                      ! i3 = i4 =  351
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  346

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  346

                                      ! do nk = 1,    3

                                      xin(346) = xin(347) + dxkl*xin(345)
                                      yin(346) = yin(347) + dykl*yin(345)
                                      zin(346) = zin(347) + dzkl*zin(345)
                                      ! i4 = i4 + lang+1 =  348

                                      ! nk =    2

                                      xin(348) = xin(349) + dxkl*xin(347)
                                      yin(348) = yin(349) + dykl*yin(347)
                                      zin(348) = zin(349) + dzkl*zin(347)
                                      ! i4 = i4 + lang+1 =  350

                                      ! nk =    3

                                      xin(350) = xin(351) + dxkl*xin(349)
                                      yin(350) = yin(351) + dykl*yin(349)
                                      zin(350) = zin(351) + dzkl*zin(349)
                                      ! i4 = i4 + lang+1 =  352

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  347

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  353

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  354

                                      ! do nl = 1,    1

                                      ! i4 = i3 =  354

                                      ! do nk = 1,    3

                                      xin(354) = xin(355) + dxkl*xin(353)
                                      yin(354) = yin(355) + dykl*yin(353)
                                      zin(354) = zin(355) + dzkl*zin(353)
                                      ! i4 = i4 + lang+1 =  356

                                      ! nk =    2

                                      xin(356) = xin(357) + dxkl*xin(355)
                                      yin(356) = yin(357) + dykl*yin(355)
                                      zin(356) = zin(357) + dzkl*zin(355)
                                      ! i4 = i4 + lang+1 =  358

                                      ! nk =    3

                                      xin(358) = xin(359) + dxkl*xin(357)
                                      yin(358) = yin(359) + dykl*yin(357)
                                      zin(358) = zin(359) + dzkl*zin(357)
                                      ! i4 = i4 + lang+1 =  360

                                      ! nk =    4

                                      ! end do

                                      ! i3 = i3 + 1 =  355

                                      ! nl =    2

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  361

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  360

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 1080! loop over all integrals

                                        l = n - 30*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d22bra(j)*d13ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 72)*yin(my + 72)*zin(mz + 72) & ! root  2
                                                        + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  3
                                                        + xin(mx + 216)*yin(my + 216)*zin(mz + 216) & ! root  4
                                                        + xin(mx + 288)*yin(my + 288)*zin(mz + 288)) ! root  5

                                        j = int(n/30) + 1 ! index for the next bra cartesian pair

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
                                    ip = (i - 1)*180 ! Stride between functions in i

                                    do j = 1, maxj2

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*30 + ip ! Add stride between functions in j

                                      do k = 1, 10 ! # of cartesians in k

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

                              deallocate (n22bra)
                              deallocate (xint22bra)
                              deallocate (n13ket)
                              deallocate (xint13ket)

                              end subroutine int2231
                              end submodule
