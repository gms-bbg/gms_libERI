submodule(rot_axis_kernels) int0122_impl
contains
  module subroutine int0122(sp_pair, dd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sp_pair, dd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n01bra(:), n22ket(:)
    real(dp), allocatable :: xint01bra(:), xint22ket(:)
    integer(kind=int64) :: nspbra, nddket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, nquarts, nquart_thrd(42), nthreads, thread_id
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutspbra, scutddket
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d01p, d22p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
  real(dp) ::  zero_m_0(1), zero_m_1(4), zero_m_2(6), zero_m_3(8), zero_m_4(10), zero_m_5(6), boys0, boys1, boys2, boys3, boys4, boys5, boys6, boys7
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21, y03
    !for r-integrals
    real(dp) :: r00(5), r01(27), r02(48), r03(60), r04(45), r05(21)
    !eri_value
    real(dp) :: eri_value(108), trans(6), qx, qz
    !digestion
    logical :: iandj, kandl, same
  integer(kind=int64) :: ii1,kk1,nij,maxl2,jj1,i2,j2,ijp,nkl,itmp,ijklp,ll2,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ll1,ijkp,ij_index,ijk_index,ijkl_index
    integer(kind=int64) :: maxj2, loci, locj, lock, locl, ip
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: nchunksize_int64
    integer(kind=int64) :: istart, iend, itile, ntile, ijkl_collapsed
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp
    real(dp) :: kernel_full1, kernel_full2, kernel_only1, kernel_only2, first_screen1, first_screen2
    integer :: shp_thresh
    mini = 1  !0
    maxi = 1  !1
    minj = 1  !2
    maxj = 3  !2
    mink = 1
    maxk = 6
    minl = 1
    maxl = 6
    allocate (n01bra(res%n_s_shl*res%n_p_shl))
    allocate (xint01bra(res%n_s_shl*res%n_p_shl))
    allocate (n22ket(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22ket(res%n_d_shl*(res%n_d_shl + 1)/2))
    !start screening
    first_screen1 = omp_get_wtime()
    scutspbra = cutoff_schwarz/maxval(dd_pair%xints)
    nspbra = 0
    do ij = 1, res%n_s_shl*res%n_p_shl
    if (sp_pair%xints(ij) .ge. scutspbra) then
      nspbra = nspbra + 1
      xint01bra(nspbra) = sp_pair%xints(ij)
      n01bra(nspbra) = ij
    end if
    end do
    scutddket = cutoff_schwarz/maxval(sp_pair%xints)
    nddket = 0
    do kl = 1, res%n_d_shl*(res%n_d_shl + 1)/2
    if (dd_pair%xints(kl) .ge. scutddket) then
      nddket = nddket + 1
      xint22ket(nddket) = dd_pair%xints(kl)
      n22ket(nddket) = kl
    end if
    end do
    first_screen2 = omp_get_wtime()

    nchunksize_int64 = 375000000
    kernel_full1 = omp_get_wtime()
    if ((nspbra*nddket) .le. nchunksize_int64) nchunksize_int64 = nspbra*nddket
    ntile = int(nspbra*nddket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nspbra*nddket

      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()
 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, dd_pair, sp_pair) &
 !$omp shared(nquart_start, nquart_end, xint01bra, xint22ket) &
 !$omp shared(nddket, n01bra, n22ket) &
 !$omp private(shp_thresh, test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(tx21, r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d01p, t_alpha, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, zero_m_4, zero_m_5, r00, r01, r02, r03, r04, r05, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, boys4, boys5, boys6, boys7, t_inverse, t_new ) &
 !$omp private(qx, qz, eri_value, kandl, same, maxl2, ijk_index, nij, nkl, ijkp, ijkl_index, ijklp, mini, maxi, minj, maxj, mink, maxk, minl, maxl) &
 !$omp private(ip, ijp, itmp, trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l, iandj, maxj2, loci, locj, lock, locl)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nddket + 1
        kl_tmp = mod(iquart - 1, nddket) + 1

        test = xint01bra(ij_tmp)*xint22ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then

          ij = n01bra(ij_tmp)
          kl = n22ket(kl_tmp)
          ish_tmp = (ij - 1)/res%n_p_shl + 1
          jsh_tmp = mod(ij - 1, res%n_p_shl) + 1
          ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
          lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

          ii = res%i_s_shl(ish_tmp)
          jj = res%i_p_shl(jsh_tmp)
          kk = res%i_d_shl(ksh_tmp)
          ll = res%i_d_shl(lsh_tmp)

          ish = ii
          jsh = jj
          ksh = kk
          lsh = ll

          !get distances r12 and r34
          r12 = 0.0_dp
          r12 = r12 + (res%coord_sh(jsh, 1) - res%coord_sh(ish, 1))*(res%coord_sh(jsh, 1) - res%coord_sh(ish, 1))
          r12 = r12 + (res%coord_sh(jsh, 2) - res%coord_sh(ish, 2))*(res%coord_sh(jsh, 2) - res%coord_sh(ish, 2))
          r12 = r12 + (res%coord_sh(jsh, 3) - res%coord_sh(ish, 3))*(res%coord_sh(jsh, 3) - res%coord_sh(ish, 3))
          r34 = 0.0_dp
          r34 = r34 + (res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1))*(res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1))
          r34 = r34 + (res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2))*(res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2))
          r34 = r34 + (res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3))*(res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3))

          !rotate axis
          buff(1) = 0.0_dp
          buff(2) = 0.0_dp
          buff(3) = 1.0_dp
          rab = 0.0_dp
          if (r12 .ne. 0.0_dp) then
            rab = sqrt(r12)
            tmp = 1.0_dp/rab
            buff(1) = (res%coord_sh(jsh, 1) - res%coord_sh(ish, 1))*tmp
            buff(2) = (res%coord_sh(jsh, 2) - res%coord_sh(ish, 2))*tmp
            buff(3) = (res%coord_sh(jsh, 3) - res%coord_sh(ish, 3))*tmp
          end if

          buff(4) = 0.0_dp
          buff(5) = 0.0_dp
          buff(6) = 1.0_dp
          rcd = 0.0_dp
          if (r34 .ne. 0.0_dp) then
            rcd = sqrt(r34)
            tmp = 1.0_dp/rcd
            buff(4) = (res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1))*tmp
            buff(5) = (res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2))*tmp
            buff(6) = (res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3))*tmp
          end if

          cosg = buff(4)*buff(1) + buff(5)*buff(2) + buff(6)*buff(3)
          cosg = min(1.0_dp, cosg)
          cosg = max(-1.0_dp, cosg)

          buff(7) = buff(6)*buff(2) - buff(5)*buff(3)
          buff(8) = buff(4)*buff(3) - buff(6)*buff(1)
          buff(9) = buff(5)*buff(1) - buff(4)*buff(2)
          if (abs(cosg) .gt. 0.9_dp) then
            sing = sqrt(buff(7)*buff(7) + buff(8)*buff(8) + buff(9)*buff(9))
          else
            sing = sqrt(1.0_dp - cosg*cosg)
          end if
          if (abs(cosg) .le. 0.9_dp .or. sing .ge. tenm12) then
            tmp = 1.0_dp/sing
            buff(7) = buff(7)*tmp
            buff(8) = buff(8)*tmp
            buff(9) = buff(9)*tmp
          else
            tmp = buff(3)*buff(3)
            if (abs(buff(1)) .le. 0.7_dp) tmp = buff(1)*buff(1)
            tmp = min(1.0_dp, tmp)
            tmp = sqrt(1.0_dp - tmp)
            if (tmp .ne. 0.0_dp) tmp = 1.0_dp/tmp
            if (abs(buff(1)) .le. 0.7_dp) then
              buff(7) = 0.0_dp
              buff(8) = buff(3)*tmp
              buff(9) = -buff(2)*tmp
            else
              buff(7) = buff(2)*tmp
              buff(8) = -buff(1)*tmp
              buff(9) = 0.0_dp
            end if
          end if
          !find the coordinates of c relative to local axes at a
  acx = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* (buff(8)*buff(3) - buff(9)*buff(2))+(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* (buff(9)*buff(1) - buff(7)*buff(3)) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* (buff(7)*buff(2) - buff(8)*buff(1))
  acy = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* buff(7) +(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* buff(8) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* buff(9)
  acz = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* buff(1) +(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* buff(2) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* buff(3)
          if (abs(acy) .le. acycut) then
            acy = 0.0_dp
            acy2 = 0.0_dp
          else
            acy2 = acy*acy
          end if
          r00 = 0.0_dp

          ! r01 = 0.0_dp
          ! r02 = 0.0_dp
          ! r03 = 0.0_dp
          ! r04 = 0.0_dp
          ! r05 = 0.0_dp
          r01(1) = 0.0_dp
          r01(2) = 0.0_dp
          r01(3) = 0.0_dp
          r01(4) = 0.0_dp
          r01(5) = 0.0_dp
          r01(6) = 0.0_dp
          r01(7) = 0.0_dp
          r01(8) = 0.0_dp
          r01(9) = 0.0_dp
          r01(10) = 0.0_dp
          r01(11) = 0.0_dp
          r01(12) = 0.0_dp
          r01(13) = 0.0_dp
          r01(14) = 0.0_dp
          r01(15) = 0.0_dp
          r01(16) = 0.0_dp
          r01(17) = 0.0_dp
          r01(18) = 0.0_dp
          r01(19) = 0.0_dp
          r01(20) = 0.0_dp
          r01(21) = 0.0_dp
          r01(22) = 0.0_dp
          r01(23) = 0.0_dp
          r01(24) = 0.0_dp
          r01(25) = 0.0_dp
          r01(26) = 0.0_dp
          r01(27) = 0.0_dp

          r02(1) = 0.0_dp
          r02(2) = 0.0_dp
          r02(3) = 0.0_dp
          r02(4) = 0.0_dp
          r02(5) = 0.0_dp
          r02(6) = 0.0_dp
          r02(7) = 0.0_dp
          r02(8) = 0.0_dp
          r02(9) = 0.0_dp
          r02(10) = 0.0_dp
          r02(11) = 0.0_dp
          r02(12) = 0.0_dp
          r02(13) = 0.0_dp
          r02(14) = 0.0_dp
          r02(15) = 0.0_dp
          r02(16) = 0.0_dp
          r02(17) = 0.0_dp
          r02(18) = 0.0_dp
          r02(19) = 0.0_dp
          r02(20) = 0.0_dp
          r02(21) = 0.0_dp
          r02(22) = 0.0_dp
          r02(23) = 0.0_dp
          r02(24) = 0.0_dp
          r02(25) = 0.0_dp
          r02(26) = 0.0_dp
          r02(27) = 0.0_dp
          r02(28) = 0.0_dp
          r02(29) = 0.0_dp
          r02(30) = 0.0_dp
          r02(31) = 0.0_dp
          r02(32) = 0.0_dp
          r02(33) = 0.0_dp
          r02(34) = 0.0_dp
          r02(35) = 0.0_dp
          r02(36) = 0.0_dp
          r02(37) = 0.0_dp
          r02(38) = 0.0_dp
          r02(39) = 0.0_dp
          r02(40) = 0.0_dp
          r02(41) = 0.0_dp
          r02(42) = 0.0_dp
          r02(43) = 0.0_dp
          r02(44) = 0.0_dp
          r02(45) = 0.0_dp
          r02(46) = 0.0_dp
          r02(47) = 0.0_dp
          r02(48) = 0.0_dp

          r03(1) = 0.0_dp
          r03(2) = 0.0_dp
          r03(3) = 0.0_dp
          r03(4) = 0.0_dp
          r03(5) = 0.0_dp
          r03(6) = 0.0_dp
          r03(7) = 0.0_dp
          r03(8) = 0.0_dp
          r03(9) = 0.0_dp
          r03(10) = 0.0_dp
          r03(11) = 0.0_dp
          r03(12) = 0.0_dp
          r03(13) = 0.0_dp
          r03(14) = 0.0_dp
          r03(15) = 0.0_dp
          r03(16) = 0.0_dp
          r03(17) = 0.0_dp
          r03(18) = 0.0_dp
          r03(19) = 0.0_dp
          r03(20) = 0.0_dp
          r03(21) = 0.0_dp
          r03(22) = 0.0_dp
          r03(23) = 0.0_dp
          r03(24) = 0.0_dp
          r03(25) = 0.0_dp
          r03(26) = 0.0_dp
          r03(27) = 0.0_dp
          r03(28) = 0.0_dp
          r03(29) = 0.0_dp
          r03(30) = 0.0_dp
          r03(31) = 0.0_dp
          r03(32) = 0.0_dp
          r03(33) = 0.0_dp
          r03(34) = 0.0_dp
          r03(35) = 0.0_dp
          r03(36) = 0.0_dp
          r03(37) = 0.0_dp
          r03(38) = 0.0_dp
          r03(39) = 0.0_dp
          r03(40) = 0.0_dp
          r03(41) = 0.0_dp
          r03(42) = 0.0_dp
          r03(43) = 0.0_dp
          r03(44) = 0.0_dp
          r03(45) = 0.0_dp
          r03(46) = 0.0_dp
          r03(47) = 0.0_dp
          r03(48) = 0.0_dp
          r03(49) = 0.0_dp
          r03(50) = 0.0_dp
          r03(51) = 0.0_dp
          r03(52) = 0.0_dp
          r03(53) = 0.0_dp
          r03(54) = 0.0_dp
          r03(55) = 0.0_dp
          r03(56) = 0.0_dp
          r03(57) = 0.0_dp
          r03(58) = 0.0_dp
          r03(59) = 0.0_dp
          r03(60) = 0.0_dp

          r04(1) = 0.0_dp
          r04(2) = 0.0_dp
          r04(3) = 0.0_dp
          r04(4) = 0.0_dp
          r04(5) = 0.0_dp
          r04(6) = 0.0_dp
          r04(7) = 0.0_dp
          r04(8) = 0.0_dp
          r04(9) = 0.0_dp
          r04(10) = 0.0_dp
          r04(11) = 0.0_dp
          r04(12) = 0.0_dp
          r04(13) = 0.0_dp
          r04(14) = 0.0_dp
          r04(15) = 0.0_dp
          r04(16) = 0.0_dp
          r04(17) = 0.0_dp
          r04(18) = 0.0_dp
          r04(19) = 0.0_dp
          r04(20) = 0.0_dp
          r04(21) = 0.0_dp
          r04(22) = 0.0_dp
          r04(23) = 0.0_dp
          r04(24) = 0.0_dp
          r04(25) = 0.0_dp
          r04(26) = 0.0_dp
          r04(27) = 0.0_dp
          r04(28) = 0.0_dp
          r04(29) = 0.0_dp
          r04(30) = 0.0_dp
          r04(31) = 0.0_dp
          r04(32) = 0.0_dp
          r04(33) = 0.0_dp
          r04(34) = 0.0_dp
          r04(35) = 0.0_dp
          r04(36) = 0.0_dp
          r04(37) = 0.0_dp
          r04(38) = 0.0_dp
          r04(39) = 0.0_dp
          r04(40) = 0.0_dp
          r04(41) = 0.0_dp
          r04(42) = 0.0_dp
          r04(43) = 0.0_dp
          r04(44) = 0.0_dp
          r04(45) = 0.0_dp

          r05(1) = 0.0_dp
          r05(2) = 0.0_dp
          r05(3) = 0.0_dp
          r05(4) = 0.0_dp
          r05(5) = 0.0_dp
          r05(6) = 0.0_dp
          r05(7) = 0.0_dp
          r05(8) = 0.0_dp
          r05(9) = 0.0_dp
          r05(10) = 0.0_dp
          r05(11) = 0.0_dp
          r05(12) = 0.0_dp
          r05(13) = 0.0_dp
          r05(14) = 0.0_dp
          r05(15) = 0.0_dp
          r05(16) = 0.0_dp
          r05(17) = 0.0_dp
          r05(18) = 0.0_dp
          r05(19) = 0.0_dp
          r05(20) = 0.0_dp
          r05(21) = 0.0_dp

          ket_loop = 0
          do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1
            t_expon_cd = dd_pair%t_expon_ab(dd_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_sp
            t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
            y03 = t_inverse_expon_cd*dd_pair%expon_a(dd_pair%pair_loc(kl) + ket_loop) !y03
            sq = dd_pair%sq(dd_pair%pair_loc(kl) + ket_loop)
            cq = dd_pair%t_alpha(dd_pair%pair_loc(kl) + ket_loop)
            cqx = cq*sing
            cqz = cq*cosg
            aqx = acx + cqx
            aqx2 = aqx*aqx
            aqxy = aqx*acy
            aqz = acz + cqz
            qps = aqx2 + acy2

            zero_m_0 = 0.0_dp
            zero_m_1 = 0.0_dp
            zero_m_2 = 0.0_dp
            zero_m_3 = 0.0_dp
            zero_m_4 = 0.0_dp
            zero_m_5 = 0.0_dp
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              bra_loop = bra_loop + 1
              !get bra shell pair info
              shp_thresh = sp_pair%ismlp(sp_pair%pair_loc(ij) + bra_loop) + dd_pair%ismlp(dd_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              t_expon_ab = sp_pair%t_expon_ab(sp_pair%pair_loc(ij) + bra_loop) !tx12
              t_alpha = sp_pair%t_alpha(sp_pair%pair_loc(ij) + bra_loop) !ty02
              t_beta = sp_pair%t_beta(sp_pair%pair_loc(ij) + bra_loop) !ty01
              tx21 = sp_pair%t_inverse_expon_ab(sp_pair%pair_loc(ij) + bra_loop) !tx21
              d01p = sp_pair%d_coeff(sp_pair%pair_loc(ij) + bra_loop) !d01p
              expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd) !x41
              pqr = t_alpha - aqz
              pqs = pqr*pqr
              rho = t_expon_ab*t_expon_cd*expon_abcd_inverse
              t = (pqs + qps)*rho
              rho = rho + rho
              if (t .le. t_max) then
                !evaluate boys function with recursion
                t_new = t*1.4323099396562323d+01 !f_increment(6)
                t_int = nint(t_new)
                fmt = boys_grid_zero((5*451*5) + (t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((5*451*5) + (t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((5*451*5) + (t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((5*451*5) + (t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((5*451*5) + (t_int*5) + 1)
                !use extrapolation for exp(-t)
                t_new = t*m_increment !rxinc
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)
                boys7 = ((t + t)*fmt + expt)*6.66666666666666666d-002!*rmr(8)
                ! boys6 = ((t + t)*boys7 + expt)* 7.6923076923076927d-002!*rmr(7)
                boys6 = ((t + t)*boys7 + expt)*7.692307692307692307692d-2
                boys5 = ((t + t)*boys6 + expt)*9.090909090909090909091d-2!*rmr(6)
                boys4 = ((t + t)*boys5 + expt)*0.111111111111111111d00!*rmr(5)
                boys3 = ((t + t)*boys4 + expt)*0.142857142857142857d00 !*rmr(4)
                boys2 = ((t + t)*boys3 + expt)*0.2000000000000000d00 !*rmr(3)
                boys1 = ((t + t)*boys2 + expt)*0.333333333333333333d00 !*rmr(2)
                boys0 = ((t + t)*boys1 + expt)!*rmr(1)
                boys0 = boys0*sqrt(expon_abcd_inverse)*d01p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d01p*rho
                boys2 = boys2*sqrt(expon_abcd_inverse)*d01p*rho*rho
                boys3 = boys3*sqrt(expon_abcd_inverse)*d01p*rho*rho*rho
                boys4 = boys4*sqrt(expon_abcd_inverse)*d01p*rho*rho*rho*rho
                boys5 = boys5*sqrt(expon_abcd_inverse)*d01p*rho*rho*rho*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d01p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
                boys3 = boys2*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho)
                boys4 = boys3*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho + t_inverse*rho)
                boys5 = boys4*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho + t_inverse*rho + t_inverse*rho)
              end if
              ! write(*,*) t_beta = , t_beta
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0*t_beta

              zero_m_1(1) = zero_m_1(1) + boys1*t_beta
              zero_m_1(2) = zero_m_1(2) + boys1*tx21
              zero_m_1(3) = zero_m_1(3) + boys1*t_beta*pqr
              zero_m_1(4) = zero_m_1(4) + boys1*tx21*pqr

              zero_m_2(1) = zero_m_2(1) + boys2*t_beta
              zero_m_2(2) = zero_m_2(2) + boys2*tx21
              zero_m_2(3) = zero_m_2(3) + boys2*t_beta*pqr
              zero_m_2(4) = zero_m_2(4) + boys2*tx21*pqr
              zero_m_2(5) = zero_m_2(5) + boys2*t_beta*pqs
              zero_m_2(6) = zero_m_2(6) + boys2*tx21*pqs

              zero_m_3(1) = zero_m_3(1) + boys3*t_beta !check ?
              zero_m_3(2) = zero_m_3(2) + boys3*t_beta*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*t_beta*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*t_beta*pqr*pqs
              zero_m_3(5) = zero_m_3(5) + boys3*tx21
              zero_m_3(6) = zero_m_3(6) + boys3*tx21*pqr
              zero_m_3(7) = zero_m_3(7) + boys3*tx21*pqs
              zero_m_3(8) = zero_m_3(8) + boys3*tx21*pqr*pqs

              zero_m_4(1) = zero_m_4(1) + boys4*t_beta
              zero_m_4(2) = zero_m_4(2) + boys4*t_beta*pqr
              zero_m_4(3) = zero_m_4(3) + boys4*t_beta*pqs
              zero_m_4(4) = zero_m_4(4) + boys4*t_beta*pqr*pqs
              zero_m_4(5) = zero_m_4(5) + boys4*t_beta*pqs*pqs
              zero_m_4(6) = zero_m_4(6) + boys4*tx21
              zero_m_4(7) = zero_m_4(7) + boys4*tx21*pqr
              zero_m_4(8) = zero_m_4(8) + boys4*tx21*pqs
              zero_m_4(9) = zero_m_4(9) + boys4*tx21*pqr*pqs
              zero_m_4(10) = zero_m_4(10) + boys4*tx21*pqs*pqs

              zero_m_5(1) = zero_m_5(1) + boys5*tx21
              zero_m_5(2) = zero_m_5(2) + boys5*tx21*pqr
              zero_m_5(3) = zero_m_5(3) + boys5*tx21*pqs
              zero_m_5(4) = zero_m_5(4) + boys5*tx21*pqs*pqr
              zero_m_5(5) = zero_m_5(5) + boys5*tx21*pqs*pqs
              zero_m_5(6) = zero_m_5(6) + boys5*tx21*pqs*pqs*pqr
            end do
            ! write(*,*) "zero_m_0(1) =" , zero_m_0(1) !good
            ! write(*,*) "zero_m_1(1) =" , zero_m_1(1)
            ! write(*,*) "zero_m_1(2) =" , zero_m_1(2)
            ! write(*,*) "zero_m_1(3) =" , zero_m_1(3)
            ! write(*,*) "zero_m_1(4) =" , zero_m_1(4)
            ! write(*,*) "zero_m_2(1) =" , zero_m_2(1)
            ! write(*,*) "zero_m_2(2) =" , zero_m_2(2)
            ! write(*,*) "zero_m_2(3) =" , zero_m_2(3)
            ! write(*,*) "zero_m_2(4) =" , zero_m_2(4)
            ! write(*,*) "zero_m_2(5) =" , zero_m_2(5)
            ! write(*,*) "zero_m_2(6) =" , zero_m_2(6)
            ! write(*,*) "zero_m_3(1) =" , zero_m_3(1)
            ! write(*,*) "zero_m_3(2) =" , zero_m_3(2)
            ! write(*,*) "zero_m_3(3) =" , zero_m_3(3)
            ! write(*,*) "zero_m_3(4) =" , zero_m_3(4)
            ! write(*,*) "zero_m_3(5) =" , zero_m_3(5)
            ! write(*,*) "zero_m_3(6) =" , zero_m_3(6)
            ! write(*,*) "zero_m_3(7) =" , zero_m_3(7)
            ! write(*,*) "zero_m_3(8) =" , zero_m_3(8)
            ! write(*,*) "zero_m_4(1) =" , zero_m_4(1)
            ! write(*,*) "zero_m_4(2) =" , zero_m_4(2)
            ! write(*,*) "zero_m_4(3) =" , zero_m_4(3)
            ! write(*,*) "zero_m_4(4) =" , zero_m_4(4)
            ! write(*,*) "zero_m_4(5) =" , zero_m_4(5)
            ! write(*,*) "zero_m_4(6) =" , zero_m_4(6)
            ! write(*,*) "zero_m_4(7) =" , zero_m_4(7)
            ! write(*,*) "zero_m_4(8) =" , zero_m_4(8)
            ! write(*,*) "zero_m_4(9) =" , zero_m_4(9)
            ! write(*,*) "zero_m_4(10) =" , zero_m_4(10)
            ! write(*,*) "zero_m_5(1) =" , zero_m_5(1)
            ! write(*,*) "zero_m_5(2) =" , zero_m_5(2)
            ! write(*,*) "zero_m_5(3) =" , zero_m_5(3)
            ! write(*,*) "zero_m_5(4) =" , zero_m_5(4)
            ! write(*,*) "zero_m_5(5) =" , zero_m_5(5)
            ! write(*,*) "zero_m_5(6) =" , zero_m_5(6)
            !r integrals here
            ! write(*,*) x43 = , t_inverse_expon_cd
            ! write(*,*) xmd4 = , t_inverse_expon_cd*0.5_dp*t_inverse_expon_cd*0.5_dp*sq

            r00(1) = r00(1) + zero_m_0(1)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r00(2) = r00(2) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r00(3) = r00(3) - zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            r00(4) = r00(4) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq*(1 - y03)*(1 - y03)
            r00(5) = r00(5) + zero_m_0(1)*y03*y03*(1 - y03)*(1 - y03)*sq !ok
            ! do i = 1, 5
            !     write(*,*) "rd0", r00(i)
            ! enddo

            r01(1) = r01(1) + zero_m_1(1)*aqx*t_inverse_expon_cd*sq*t_inverse_expon_cd*y03*0.25_dp !5 !ok
            r01(2) = r01(2) + zero_m_1(1)*acy*t_inverse_expon_cd*sq*t_inverse_expon_cd*y03*0.25_dp
            r01(3) = r01(3) - zero_m_1(3)*t_inverse_expon_cd*sq*t_inverse_expon_cd*y03*0.25_dp
            r01(4) = r01(4) - zero_m_1(1)*aqx*t_inverse_expon_cd*sq*t_inverse_expon_cd*(1 - y03)*0.25_dp !6
            r01(5) = r01(5) - zero_m_1(1)*acy*t_inverse_expon_cd*sq*t_inverse_expon_cd*(1 - y03)*0.25_dp
            r01(6) = r01(6) + zero_m_1(3)*t_inverse_expon_cd*sq*t_inverse_expon_cd*(1 - y03)*0.25_dp
            r01(7) = r01(7) - zero_m_1(1)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03*y03*(1 - y03) !7
            r01(8) = r01(8) - zero_m_1(1)*acy*t_inverse_expon_cd*0.5_dp*sq*y03*y03*(1 - y03)
            r01(9) = r01(9) + zero_m_1(3)*t_inverse_expon_cd*0.5_dp*sq*y03*y03*(1 - y03)
            r01(10) = r01(10) - zero_m_1(1)*aqx*-t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)*(1 - y03) !8
            r01(11) = r01(11) - zero_m_1(1)*acy*-t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)*(1 - y03)
            r01(12) = r01(12) - zero_m_1(3)*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)*(1 - y03)
            r01(13) = r01(13) - zero_m_1(2)*aqx*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp !9
            r01(14) = r01(14) - zero_m_1(2)*acy*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r01(15) = r01(15) + zero_m_1(4)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r01(16) = r01(16) - zero_m_1(2)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03*y03 !10
            r01(17) = r01(17) - zero_m_1(2)*acy*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r01(18) = r01(18) + zero_m_1(4)*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r01(19) = r01(19) + zero_m_1(2)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03) !11
            r01(20) = r01(20) + zero_m_1(2)*acy*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            r01(21) = r01(21) - zero_m_1(4)*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            r01(22) = r01(22) - zero_m_1(2)*aqx*t_inverse_expon_cd*0.5_dp*sq*(1 - y03)*(1 - y03) !12
            r01(23) = r01(23) - zero_m_1(2)*acy*t_inverse_expon_cd*0.5_dp*sq*(1 - y03)*(1 - y03)
            r01(24) = r01(24) + zero_m_1(4)*t_inverse_expon_cd*0.5_dp*sq*(1 - y03)*(1 - y03)
            r01(25) = r01(25) - zero_m_1(2)*aqx*y03*y03*(1 - y03)*(1 - y03)*sq !13
            r01(26) = r01(26) - zero_m_1(2)*acy*y03*y03*(1 - y03)*(1 - y03)*sq
            r01(27) = r01(27) + zero_m_1(4)*y03*y03*(1 - y03)*(1 - y03)*sq
            ! do i = 1, 27
            !     write(*,*) "rd1", r01(i)
            ! enddo

            buff(10) = t_inverse_expon_cd*t_inverse_expon_cd*t_inverse_expon_cd*sq*0.125_dp
            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*buff(10) !5, worint64
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*buff(10)
            r02(3) = r02(3) + (zero_m_2(5) - zero_m_1(1))*buff(10)
            r02(4) = r02(4) + zero_m_2(1)*aqxy*buff(10)
            r02(5) = r02(5) - zero_m_2(3)*aqx*buff(10)
            r02(6) = r02(6) - zero_m_2(3)*acy*buff(10)

            buff(10) = t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq
            r02(7) = r02(7) + (zero_m_2(1)*aqx2 - zero_m_1(1))*buff(10)*y03*y03 !6, work11
            r02(8) = r02(8) + (zero_m_2(1)*acy2 - zero_m_1(1))*buff(10)*y03*y03
            r02(9) = r02(9) + (zero_m_2(5) - zero_m_1(1))*buff(10)*y03*y03
            r02(10) = r02(10) + zero_m_2(1)*aqxy*buff(10)*y03*y03
            r02(11) = r02(11) - zero_m_2(3)*aqx*buff(10)*y03*y03
            r02(12) = r02(12) - zero_m_2(3)*acy*buff(10)*y03*y03

            r02(13) = r02(13) - (zero_m_2(1)*aqx2 - zero_m_1(1))*buff(10)*y03*(1 - y03) !7, work12
            r02(14) = r02(14) - (zero_m_2(1)*acy2 - zero_m_1(1))*buff(10)*y03*(1 - y03)
            r02(15) = r02(15) - (zero_m_2(5) - zero_m_1(1))*buff(10)*y03*(1 - y03)
            r02(16) = r02(16) - zero_m_2(1)*aqxy*buff(10)*y03*(1 - y03)
            r02(17) = r02(17) + zero_m_2(3)*aqx*buff(10)*y03*(1 - y03)
            r02(18) = r02(18) + zero_m_2(3)*acy*buff(10)*y03*(1 - y03)

            r02(19) = r02(19) + (zero_m_2(1)*aqx2 - zero_m_1(1))*buff(10)*(1 - y03)*(1 - y03) !8, work13
            r02(20) = r02(20) + (zero_m_2(1)*acy2 - zero_m_1(1))*buff(10)*(1 - y03)*(1 - y03)
            r02(21) = r02(21) + (zero_m_2(5) - zero_m_1(1))*buff(10)*(1 - y03)*(1 - y03)
            r02(22) = r02(22) + zero_m_2(1)*aqxy*buff(10)*(1 - y03)*(1 - y03)
            r02(23) = r02(23) - zero_m_2(3)*aqx*buff(10)*(1 - y03)*(1 - y03)
            r02(24) = r02(24) - zero_m_2(3)*acy*buff(10)*(1 - y03)*(1 - y03)

            r02(25) = r02(25) - (zero_m_2(2)*aqx2 - zero_m_1(2))*buff(10)*y03 !9, work6
            r02(26) = r02(26) - (zero_m_2(2)*acy2 - zero_m_1(2))*buff(10)*y03
            r02(27) = r02(27) - (zero_m_2(6) - zero_m_1(2))*buff(10)*y03
            r02(28) = r02(28) - zero_m_2(2)*aqxy*buff(10)*y03
            r02(29) = r02(29) + zero_m_2(4)*aqx*buff(10)*y03
            r02(30) = r02(30) + zero_m_2(4)*acy*buff(10)*y03

            r02(31) = r02(31) + (zero_m_2(2)*aqx2 - zero_m_1(2))*buff(10)*(1 - y03) !10, work7
            r02(32) = r02(32) + (zero_m_2(2)*acy2 - zero_m_1(2))*buff(10)*(1 - y03)
            r02(33) = r02(33) + (zero_m_2(6) - zero_m_1(2))*buff(10)*(1 - y03)
            r02(34) = r02(34) + zero_m_2(2)*aqxy*buff(10)*(1 - y03)
            r02(35) = r02(35) - zero_m_2(4)*aqx*buff(10)*(1 - y03)
            r02(36) = r02(36) - zero_m_2(4)*acy*buff(10)*(1 - y03)

            buff(10) = t_inverse_expon_cd*0.5_dp*sq
            r02(37) = r02(37) + (zero_m_2(2)*aqx2 - zero_m_1(2))*buff(10)*y03*y03*(1 - y03) !11, work8
            r02(38) = r02(38) + (zero_m_2(2)*acy2 - zero_m_1(2))*buff(10)*y03*y03*(1 - y03)
            r02(39) = r02(39) + (zero_m_2(6) - zero_m_1(2))*buff(10)*y03*y03*(1 - y03)
            r02(40) = r02(40) + zero_m_2(2)*aqxy*buff(10)*y03*y03*(1 - y03)
            r02(41) = r02(41) - zero_m_2(4)*aqx*buff(10)*y03*y03*(1 - y03)
            r02(42) = r02(42) - zero_m_2(4)*acy*buff(10)*y03*y03*(1 - y03)

            r02(43) = r02(43) - (zero_m_2(2)*aqx2 - zero_m_1(2))*buff(10)*y03*(1 - y03)*(1 - y03) !12, work9
            r02(44) = r02(44) - (zero_m_2(2)*acy2 - zero_m_1(2))*buff(10)*y03*(1 - y03)*(1 - y03)
            r02(45) = r02(45) - (zero_m_2(6) - zero_m_1(2))*buff(10)*y03*(1 - y03)*(1 - y03)
            r02(46) = r02(46) - zero_m_2(2)*aqxy*buff(10)*y03*(1 - y03)*(1 - y03)
            r02(47) = r02(47) + zero_m_2(4)*aqx*buff(10)*y03*(1 - y03)*(1 - y03)
            r02(48) = r02(48) + zero_m_2(4)*acy*buff(10)*y03*(1 - y03)*(1 - y03)
            ! do i = 1, 48
            !     write(*,*) "rd2", r02(i)
            ! enddo !ok

            buff(10) = t_inverse_expon_cd*t_inverse_expon_cd*t_inverse_expon_cd*sq*0.125_dp
            r03(1) = r03(1) + (zero_m_3(1)*aqx2 - zero_m_2(1)*3.0_dp)*aqx*buff(10)*y03 !3, work14
            r03(2) = r03(2) + (zero_m_3(1)*aqx2 - zero_m_2(1))*acy*buff(10)*y03
            r03(3) = r03(3) - (zero_m_3(2)*aqx2 - zero_m_2(3))*buff(10)*y03
            r03(4) = r03(4) + (zero_m_3(1)*acy2 - zero_m_2(1))*aqx*buff(10)*y03
            r03(5) = r03(5) - zero_m_3(2)*aqxy*buff(10)*y03
            r03(6) = r03(6) + (zero_m_3(3) - zero_m_2(1))*aqx*buff(10)*y03
            r03(7) = r03(7) + (zero_m_3(1)*acy2 - zero_m_2(1)*3.0_dp)*acy*buff(10)*y03
            r03(8) = r03(8) - (zero_m_3(2)*acy2 - zero_m_2(3))*buff(10)*y03
            r03(9) = r03(9) + (zero_m_3(3) - zero_m_2(1))*acy*buff(10)*y03
            r03(10) = r03(10) - (zero_m_3(4) - zero_m_2(3)*3.0_dp)*buff(10)*y03

            r03(11) = r03(11) - (zero_m_3(1)*aqx2 - zero_m_2(1)*3.0_dp)*aqx*buff(10)*(1 - y03) !4, work15
            r03(12) = r03(12) - (zero_m_3(1)*aqx2 - zero_m_2(1))*acy*buff(10)*(1 - y03)
            r03(13) = r03(13) + (zero_m_3(2)*aqx2 - zero_m_2(3))*buff(10)*(1 - y03)
            r03(14) = r03(14) - (zero_m_3(1)*acy2 - zero_m_2(1))*aqx*buff(10)*(1 - y03)
            r03(15) = r03(15) + zero_m_3(2)*aqxy*buff(10)*(1 - y03)
            r03(16) = r03(16) - (zero_m_3(3) - zero_m_2(1))*aqx*buff(10)*(1 - y03)
            r03(17) = r03(17) - (zero_m_3(1)*acy2 - zero_m_2(1)*3.0_dp)*acy*buff(10)*(1 - y03)
            r03(18) = r03(18) + (zero_m_3(2)*acy2 - zero_m_2(3))*buff(10)*(1 - y03)
            r03(19) = r03(19) - (zero_m_3(3) - zero_m_2(1))*acy*buff(10)*(1 - y03)
            r03(20) = r03(20) + (zero_m_3(4) - zero_m_2(3)*3.0_dp)*buff(10)*(1 - y03)

            r03(21) = r03(21) - (zero_m_3(5)*aqx2 - zero_m_2(2)*3.0_dp)*aqx*buff(10) !5,worint64
            r03(22) = r03(22) - (zero_m_3(5)*aqx2 - zero_m_2(2))*acy*buff(10)
            r03(23) = r03(23) + (zero_m_3(6)*aqx2 - zero_m_2(4))*buff(10)
            r03(24) = r03(24) - (zero_m_3(5)*acy2 - zero_m_2(2))*aqx*buff(10)
            r03(25) = r03(25) + zero_m_3(6)*aqxy*buff(10)
            r03(26) = r03(26) - (zero_m_3(7) - zero_m_2(2))*aqx*buff(10)
            r03(27) = r03(27) - (zero_m_3(5)*acy2 - zero_m_2(2)*3.0_dp)*acy*buff(10)
            r03(28) = r03(28) + (zero_m_3(6)*acy2 - zero_m_2(4))*buff(10)
            r03(29) = r03(29) - (zero_m_3(7) - zero_m_2(2))*acy*buff(10)
            r03(30) = r03(30) + (zero_m_3(8) - zero_m_2(4)*3.0_dp)*buff(10)

            buff(10) = t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq
            r03(31) = r03(31) - (zero_m_3(5)*aqx2 - zero_m_2(2)*3.0_dp)*aqx*buff(10)*y03*y03 !6,work11
            r03(32) = r03(32) - (zero_m_3(5)*aqx2 - zero_m_2(2))*acy*buff(10)*y03*y03
            r03(33) = r03(33) + (zero_m_3(6)*aqx2 - zero_m_2(4))*buff(10)*y03*y03 !here
            r03(34) = r03(34) - (zero_m_3(5)*acy2 - zero_m_2(2))*aqx*buff(10)*y03*y03
            r03(35) = r03(35) + zero_m_3(6)*aqxy*buff(10)*y03*y03
            r03(36) = r03(36) - (zero_m_3(7) - zero_m_2(2))*aqx*buff(10)*y03*y03
            r03(37) = r03(37) - (zero_m_3(5)*acy2 - zero_m_2(2)*3.0_dp)*acy*buff(10)*y03*y03
            r03(38) = r03(38) + (zero_m_3(6)*acy2 - zero_m_2(4))*buff(10)*y03*y03
            r03(39) = r03(39) - (zero_m_3(7) - zero_m_2(2))*acy*buff(10)*y03*y03
            r03(40) = r03(40) + (zero_m_3(8) - zero_m_2(4)*3.0_dp)*buff(10)*y03*y03

            r03(41) = r03(41) + (zero_m_3(5)*aqx2 - zero_m_2(2)*3.0_dp)*aqx*buff(10)*y03*(1 - y03) !7,work12
            r03(42) = r03(42) + (zero_m_3(5)*aqx2 - zero_m_2(2))*acy*buff(10)*y03*(1 - y03)
            r03(43) = r03(43) - (zero_m_3(6)*aqx2 - zero_m_2(4))*buff(10)*y03*(1 - y03)
            r03(44) = r03(44) + (zero_m_3(5)*acy2 - zero_m_2(2))*aqx*buff(10)*y03*(1 - y03)
            r03(45) = r03(45) - zero_m_3(6)*aqxy*buff(10)*y03*(1 - y03)
            r03(46) = r03(46) + (zero_m_3(7) - zero_m_2(2))*aqx*buff(10)*y03*(1 - y03)
            r03(47) = r03(47) + (zero_m_3(5)*acy2 - zero_m_2(2)*3.0_dp)*acy*buff(10)*y03*(1 - y03)
            r03(48) = r03(48) - (zero_m_3(6)*acy2 - zero_m_2(4))*buff(10)*y03*(1 - y03)
            r03(49) = r03(49) + (zero_m_3(7) - zero_m_2(2))*acy*buff(10)*y03*(1 - y03)
            r03(50) = r03(50) - (zero_m_3(8) - zero_m_2(4)*3.0_dp)*buff(10)*y03*(1 - y03)

            r03(51) = r03(51) - (zero_m_3(5)*aqx2 - zero_m_2(2)*3.0_dp)*aqx*buff(10)*(1 - y03)*(1 - y03) !8 !,work13
            r03(52) = r03(52) - (zero_m_3(5)*aqx2 - zero_m_2(2))*acy*buff(10)*(1 - y03)*(1 - y03)
            r03(53) = r03(53) + (zero_m_3(6)*aqx2 - zero_m_2(4))*buff(10)*(1 - y03)*(1 - y03)
            r03(54) = r03(54) - (zero_m_3(5)*acy2 - zero_m_2(2))*aqx*buff(10)*(1 - y03)*(1 - y03)
            r03(55) = r03(55) + zero_m_3(6)*aqxy*buff(10)*(1 - y03)*(1 - y03)
            r03(56) = r03(56) - (zero_m_3(7) - zero_m_2(2))*aqx*buff(10)*(1 - y03)*(1 - y03)
            r03(57) = r03(57) - (zero_m_3(5)*acy2 - zero_m_2(2)*3.0_dp)*acy*buff(10)*(1 - y03)*(1 - y03)
            r03(58) = r03(58) + (zero_m_3(6)*acy2 - zero_m_2(4))*buff(10)*(1 - y03)*(1 - y03)
            r03(59) = r03(59) - (zero_m_3(7) - zero_m_2(2))*acy*buff(10)*(1 - y03)*(1 - y03)
            r03(60) = r03(60) + (zero_m_3(8) - zero_m_2(4)*3.0_dp)*buff(10)*(1 - y03)*(1 - y03)
            ! do i = 1, 60
            !     write(*,*) "rd3", r03(i)
            ! enddo !ok

            buff(4) = t_inverse_expon_cd*sq*t_inverse_expon_cd*t_inverse_expon_cd*t_inverse_expon_cd*0.0625_dp !ok
            r04(1) = r04(1) + (zero_m_4(1)*aqx2*aqx2 - zero_m_3(1)*6.0_dp*aqx2 + zero_m_2(1)*3.0_dp)*buff(4) !2
            r04(2) = r04(2) + (zero_m_4(1)*aqx2 - zero_m_3(1)*3.0_dp)*buff(4)*aqxy
            r04(3) = r04(3) - (zero_m_4(2)*aqx2 - zero_m_3(2)*3.0_dp)*buff(4)*aqx
            r04(4) = r04(4) + (zero_m_4(1)*aqx2*acy2 - zero_m_3(1)*(aqx2 + acy2) + zero_m_2(1))*buff(4)
            r04(5) = r04(5) - (zero_m_4(2)*aqx2 - zero_m_3(2))*buff(4)*acy
            r04(6) = r04(6) + (zero_m_4(3)*aqx2 - zero_m_3(1)*aqx2 - zero_m_3(3) + zero_m_2(1))*buff(4)
            r04(7) = r04(7) + (zero_m_4(1)*acy2 - zero_m_3(1)*3.0_dp)*buff(4)*aqxy
            r04(8) = r04(8) - (zero_m_4(2)*acy2 - zero_m_3(2))*buff(4)*aqx
            r04(9) = r04(9) + (zero_m_4(3) - zero_m_3(1))*buff(4)*aqxy
            r04(10) = r04(10) - (zero_m_4(4) - zero_m_3(2)*3.0_dp)*buff(4)*aqx
            r04(11) = r04(11) + (zero_m_4(1)*acy2*acy2 - zero_m_3(1)*6.0_dp*acy2 + zero_m_2(1)*3.0_dp)*buff(4)
            r04(12) = r04(12) - (zero_m_4(2)*acy2 - zero_m_3(2)*3.0_dp)*buff(4)*acy
            r04(13) = r04(13) + (zero_m_4(3)*acy2 - zero_m_3(1)*acy2 - zero_m_3(3) + zero_m_2(1))*buff(4)
            r04(14) = r04(14) - (zero_m_4(4) - zero_m_3(2)*3.0_dp)*buff(4)*acy
            r04(15) = r04(15) + (zero_m_4(5) - zero_m_3(3)*6.0_dp + zero_m_2(1)*3.0_dp)*buff(4)

            buff(5) = t_inverse_expon_cd*sq*t_inverse_expon_cd*t_inverse_expon_cd*0.125_dp
            r04(16) = r04(16) - (zero_m_4(6)*aqx2*aqx2 - zero_m_3(5)*6.0_dp*aqx2 + zero_m_2(2)*3.0_dp)*buff(5)*y03 !3, work(14)
            r04(17) = r04(17) - (zero_m_4(6)*aqx2 - zero_m_3(5)*3.0_dp)*aqxy*buff(5)*y03
            r04(18) = r04(18) + (zero_m_4(7)*aqx2 - zero_m_3(6)*3.0_dp)*aqx*buff(5)*y03
            r04(19) = r04(19) - (zero_m_4(6)*aqx2*acy2 - zero_m_3(5)*(aqx2 + acy2) + zero_m_2(2))*buff(5)*y03
            r04(20) = r04(20) + (zero_m_4(7)*aqx2 - zero_m_3(6))*acy*buff(5)*y03
            r04(21) = r04(21) - (zero_m_4(8)*aqx2 - zero_m_3(5)*aqx2 - zero_m_3(7) + zero_m_2(2))*buff(5)*y03
            r04(22) = r04(22) - (zero_m_4(6)*acy2 - zero_m_3(5)*3.0_dp)*aqxy*buff(5)*y03
            r04(23) = r04(23) + (zero_m_4(7)*acy2 - zero_m_3(6))*aqx*buff(5)*y03
            r04(24) = r04(24) - (zero_m_4(8) - zero_m_3(5))*aqxy*buff(5)*y03
            r04(25) = r04(25) + (zero_m_4(9) - zero_m_3(6)*3.0_dp)*aqx*buff(5)*y03
            r04(26) = r04(26) - (zero_m_4(6)*acy2*acy2 - zero_m_3(5)*6.0_dp*acy2 + zero_m_2(2)*3.0_dp)*buff(5)*y03
            r04(27) = r04(27) + (zero_m_4(7)*acy2 - zero_m_3(6)*3.0_dp)*acy*buff(5)*y03
            r04(28) = r04(28) - (zero_m_4(8)*acy2 - zero_m_3(5)*acy2 - zero_m_3(7) + zero_m_2(2))*buff(5)*y03
            r04(29) = r04(29) + (zero_m_4(9) - zero_m_3(6)*3.0_dp)*acy*buff(5)*y03
            r04(30) = r04(30) - (zero_m_4(10) - zero_m_3(7)*6.0_dp + zero_m_2(2)*3.0_dp)*buff(5)*y03

            r04(31) = r04(31) + (zero_m_4(6)*aqx2*aqx2 - zero_m_3(5)*6.0_dp*aqx2 + zero_m_2(2)*3.0_dp)*buff(5)*(1 - y03) !4
            r04(32) = r04(32) + (zero_m_4(6)*aqx2 - zero_m_3(5)*3.0_dp)*aqxy*buff(5)*(1 - y03)
            r04(33) = r04(33) - (zero_m_4(7)*aqx2 - zero_m_3(6)*3.0_dp)*aqx*buff(5)*(1 - y03)
            r04(34) = r04(34) + (zero_m_4(6)*aqx2*acy2 - zero_m_3(5)*(aqx2 + acy2) + zero_m_2(2))*buff(5)*(1 - y03)
            r04(35) = r04(35) - (zero_m_4(7)*aqx2 - zero_m_3(6))*acy*buff(5)*(1 - y03)
            r04(36) = r04(36) + (zero_m_4(8)*aqx2 - zero_m_3(5)*aqx2 - zero_m_3(7) + zero_m_2(2))*buff(5)*(1 - y03)
            r04(37) = r04(37) + (zero_m_4(6)*acy2 - zero_m_3(5)*3.0_dp)*aqxy*buff(5)*(1 - y03)
            r04(38) = r04(38) - (zero_m_4(7)*acy2 - zero_m_3(6))*aqx*buff(5)*(1 - y03)
            r04(39) = r04(39) + (zero_m_4(8) - zero_m_3(5))*aqxy*buff(5)*(1 - y03)
            r04(40) = r04(40) - (zero_m_4(9) - zero_m_3(6)*3.0_dp)*aqx*buff(5)*(1 - y03)
            r04(41) = r04(41) + (zero_m_4(6)*acy2*acy2 - zero_m_3(5)*6.0_dp*acy2 + zero_m_2(2)*3.0_dp)*buff(5)*(1 - y03)
            r04(42) = r04(42) - (zero_m_4(7)*acy2 - zero_m_3(6)*3.0_dp)*acy*buff(5)*(1 - y03)
            r04(43) = r04(43) + (zero_m_4(8)*acy2 - zero_m_3(5)*acy2 - zero_m_3(7) + zero_m_2(2))*buff(5)*(1 - y03)
            r04(44) = r04(44) - (zero_m_4(9) - zero_m_3(6)*3.0_dp)*acy*buff(5)*(1 - y03)
            r04(45) = r04(45) + (zero_m_4(10) - zero_m_3(7)*6.0_dp + zero_m_2(2)*3.0_dp)*buff(5)*(1 - y03)
            ! write(*,*) fwk(15,3) = , (zero_m_4(10)-zero_m_3(7)*6.0_dp+zero_m_2(2)*3.0_dp)
            ! write(*,*) work = , t_inverse_expon_cd *sq*t_inverse_expon_cd *(1-y03)*0.25_dp
            ! write(*,*) r04(45), r04(45)
            ! do i = 1, 45
            !     write(*,*) "rd4", r04(i)
            ! enddo
            r05(1) = r05(1) - (zero_m_5(1)*aqx2*aqx2 - zero_m_4(6)*10.0_dp*aqx2 + zero_m_3(5)*15.0_dp)*buff(4)*aqx
            r05(2) = r05(2) - (zero_m_5(1)*aqx2*aqx2 - zero_m_4(6)*6.0_dp*aqx2 + zero_m_3(5)*3.0_dp)*buff(4)*acy
            r05(3) = r05(3) + (zero_m_5(2)*aqx2*aqx2 - zero_m_4(7)*6.0_dp*aqx2 + zero_m_3(6)*3.0_dp)*buff(4)
            r05(4) = r05(4) - (zero_m_5(1)*aqx2*acy2 - zero_m_4(6)*aqx2 - zero_m_4(6)*3.0_dp*acy2 + zero_m_3(5)*3.0_dp)*buff(4)*aqx
            r05(5) = r05(5) + (zero_m_5(2)*aqx2 - zero_m_4(7)*3.0_dp)*buff(4)*aqxy
            r05(6) = r05(6) - (zero_m_5(3)*aqx2 - zero_m_4(6)*aqx2 - zero_m_4(8)*3.0_dp + zero_m_3(5)*3.0_dp)*buff(4)*aqx
            r05(7) = r05(7) - (zero_m_5(1)*aqx2*acy2 - zero_m_4(6)*3.0_dp*aqx2 - zero_m_4(6)*acy2 + zero_m_3(5)*3.0_dp)*buff(4)*acy
            r05(8) = r05(8) + (zero_m_5(2)*aqx2*acy2 - zero_m_4(7)*(aqx2 + acy2) + zero_m_3(6))*buff(4)
            r05(9) = r05(9) - (zero_m_5(3)*aqx2 - zero_m_4(6)*aqx2 - zero_m_4(8) + zero_m_3(5))*buff(4)*acy
            r05(10) = r05(10) + (zero_m_5(4)*aqx2 - zero_m_4(7)*3.0_dp*aqx2 - zero_m_4(9) + zero_m_3(6)*3.0_dp)*buff(4)
            r05(11) = r05(11) - (zero_m_5(1)*acy2*acy2 - zero_m_4(6)*6.0_dp*acy2 + zero_m_3(5)*3.0_dp)*buff(4)*aqx
            r05(12) = r05(12) + (zero_m_5(2)*acy2 - zero_m_4(7)*3.0_dp)*buff(4)*aqxy
            r05(13) = r05(13) - (zero_m_5(3)*acy2 - zero_m_4(6)*acy2 - zero_m_4(8) + zero_m_3(5))*buff(4)*aqx
            r05(14) = r05(14) + (zero_m_5(4) - zero_m_4(7)*3.0_dp)*buff(4)*aqxy
            r05(15) = r05(15) - (zero_m_5(5) - zero_m_4(8)*6.0_dp + zero_m_3(5)*3.0_dp)*buff(4)*aqx
            r05(16) = r05(16) - (zero_m_5(1)*acy2*acy2 - zero_m_4(6)*10.0_dp*acy2 + zero_m_3(5)*15.0_dp)*buff(4)*acy
            r05(17) = r05(17) + (zero_m_5(2)*acy2*acy2 - zero_m_4(7)*6.0_dp*acy2 + zero_m_3(6)*3.0_dp)*buff(4)
            r05(18) = r05(18) - (zero_m_5(3)*acy2 - zero_m_4(6)*acy2 - zero_m_4(8)*3.0_dp + zero_m_3(5)*3.0_dp)*buff(4)*acy
            r05(19) = r05(19) + (zero_m_5(4)*acy2 - zero_m_4(7)*3.0_dp*acy2 - zero_m_4(9) + zero_m_3(6)*3.0_dp)*buff(4)
            r05(20) = r05(20) - (zero_m_5(5) - zero_m_4(8)*6.0_dp + zero_m_3(5)*3.0_dp)*buff(4)*acy
            r05(21) = r05(21) + (zero_m_5(6) - zero_m_4(9)*10.0_dp + zero_m_3(6)*15.0_dp)*buff(4)
            ! do i = 1, 21
            !    write(*,*) "rd5", r05(i)
            ! enddo

          end do !end ket loop
          qx = rcd*sing
          qz = rcd*cosg
  eri_value(1) = -(3.0_dp*r01(13)) - qx*qx*qx*qx*r01(25) + 2.0_dp*qx*qx*qx*(-r02(37) - r02(43)) - 6.0_dp*r03(21) + qx*qx*(-r01(16) - 4.0_dp*r01(19) - r01(22) - r03(31) - 4.0_dp*r03(41) - r03(51)) + 2.0_dp*qx*(3.0_dp*(-r02(25) - r02(31)) - r04(16) - r04(31)) - r05(1)
          eri_value(2) = -r01(13) - r03(21) - r03(24) + qx*qx*(-r01(22) - r03(54)) + 2.0_dp*qx*(-r02(31) - r04(34)) - r05(4)
  eri_value(3) = -r01(13) - qx*qx*qz*qz*r01(25) - 2.0_dp*qx*qz*qz*r02(37) - 2.0_dp*qx*qx*qz*r02(47) - r03(21) - r03(26) + qz*qz*(-r01(16) - r03(31)) - 4*qx*qz*r03(43) + qx*qx*(-r01(22) - r03(56)) + 2.0_dp*qz*(-r02(29) - r04(18)) + 2.0_dp*qx*(-r02(31) - r04(36)) - r05(6)
  eri_value(4) = -(qx*qx*qx*r02(46)) - 3.0_dp*r03(22) + qx*qx*(-(2.0_dp*r03(42)) - r03(52)) + qx*(-r02(28) - 2.0_dp*r02(34) - r04(17) - 2.0_dp*r04(32)) - r05(2)
  eri_value(5) = -(qx*qx*qx*qz*r01(25)) + qx*qx*qz*(-(2.0_dp*r02(37)) - r02(43)) - qx*qx*qx*r02(47) - 3.0_dp*r03(23) + qx*qz*(-r01(16) - 2.0_dp*r01(19) - r03(31) - 2.0_dp*r03(41)) + qx*qx*(-(2.0_dp*r03(43)) - r03(53)) + qz*(-(3.0_dp*r02(25)) - r04(16)) + qx*(-r02(29) - 2.0_dp*r02(35) - r04(18) - 2.0_dp*r04(33)) - r05(3)
eri_value(6) = -(qx*qx*qz*r02(46)) - r03(25) - 2.0_dp*qx*qz*r03(42) - qx*qx*r03(55) + qz*(-r02(28) - r04(17)) - 2.0_dp*qx*r04(35) - r05(5)
          eri_value(7) = -r01(13) - r03(21) - r03(24) + qx*qx*(-r01(16) - r03(34)) + 2.0_dp*qx*(-r02(25) - r04(19)) - r05(4)
          eri_value(8) = -(3.0_dp*r01(13)) - 6.0_dp*r03(24) - r05(11)
          eri_value(9) = -r01(13) - r03(24) - r03(26) + qz*qz*(-r01(16) - r03(34)) + 2.0_dp*qz*(-r02(29) - r04(23)) - r05(13)
          eri_value(10) = -(3.0_dp*r03(22)) + qx*(-(3.0_dp*r02(28)) - r04(22)) - r05(7)
          eri_value(11) = -r03(23) + qx*qz*(-r01(16) - r03(34)) + qz*(-r02(25) - r04(19)) + qx*(-r02(29) - r04(23)) - r05(8)
          eri_value(12) = -(3.0_dp*r03(25)) + qz*(-(3.0_dp*r02(28)) - r04(22)) - r05(12)
  eri_value(13) = -r01(13) - qx*qx*qz*qz*r01(25) - 2.0_dp*qx*qx*qz*r02(41) - 2.0_dp*qx*qz*qz*r02(43) - r03(21) - r03(26) + qx*qx*(-r01(16) - r03(36)) - 4*qx*qz*r03(43) + qz*qz*(-r01(22) - r03(51)) + 2.0_dp*qx*(-r02(25) - r04(21)) + 2.0_dp*qz*(-r02(35) - r04(33)) - r05(6)
          eri_value(14) = -r01(13) - r03(24) - r03(26) + qz*qz*(-r01(22) - r03(54)) + 2.0_dp*qz*(-r02(35) - r04(38)) - r05(13)
  eri_value(15) = -(3.0_dp*r01(13)) - qz*qz*qz*qz*r01(25) + 2.0_dp*qz*qz*qz*(-r02(41) - r02(47)) - 6.0_dp*r03(26) + qz*qz*(-r01(16) - 4.0_dp*r01(19) - r01(22) - r03(36) - 4.0_dp*r03(46) - r03(56)) + 2.0_dp*qz*(3.0_dp*(-r02(29) - r02(35)) - r04(25) - r04(40)) - r05(15)
  eri_value(16) = -(qx*qz*qz*r02(46)) - r03(22) - 2.0_dp*qx*qz*r03(45) - qz*qz*r03(52) + qx*(-r02(28) - r04(24)) - 2.0_dp*qz*r04(35) - r05(9)
  eri_value(17) = -(qx*qz*qz*qz*r01(25)) - qz*qz*qz*r02(43) + qx*qz*qz*(-(2.0_dp*r02(41)) - r02(47)) - 3.0_dp*r03(23) + qx*qz*(-r01(16) - 2.0_dp*r01(19) - r03(36) - 2.0_dp*r03(46)) + qz*qz*(-(2.0_dp*r03(43)) - r03(53)) + qx*(-(3.0_dp*r02(29)) - r04(25)) + qz*(-r02(25) - 2.0_dp*r02(31) - r04(21) - 2.0_dp*r04(36)) - r05(10)
  eri_value(18) = -(qz*qz*qz*r02(46)) - 3.0_dp*r03(25) + qz*qz*(-(2.0_dp*r03(45)) - r03(55)) + qz*(-r02(28) - 2.0_dp*r02(34) - r04(24) - 2.0_dp*r04(39)) - r05(14)
  eri_value(19) = -(qx*qx*qx*r02(40)) - 3.0_dp*r03(22) + qx*qx*(-r03(32) - 2.0_dp*r03(42)) + qx*(-(2.0_dp*r02(28)) - r02(34) - 2.0_dp*r04(17) - r04(32)) - r05(2)
          eri_value(20) = -(3.0_dp*r03(22)) + qx*(-(3.0_dp*r02(34)) - r04(37)) - r05(7)
  eri_value(21) = -(qx*qz*qz*r02(40)) - r03(22) - qz*qz*r03(32) - 2.0_dp*qx*qz*r03(45) - 2.0_dp*qz*r04(20) + qx*(-r02(34) - r04(39)) - r05(9)
    eri_value(22) = -r01(13) - r03(21) - r03(24) + qx*qx*(-r01(19) - r03(44)) + qx*(-r02(25) - r02(31) - r04(19) - r04(34)) - r05(4)
  eri_value(23) = -(qx*qx*qz*r02(40)) - r03(25) + qx*qz*(-r03(32) - r03(42)) - qx*qx*r03(45) + qz*(-r02(28) - r04(17)) + qx*(-r04(20) - r04(35)) - r05(5)
          eri_value(24) = -r03(23) + qx*qz*(-r01(19) - r03(44)) + qz*(-r02(25) - r04(19)) + qx*(-r02(35) - r04(38)) - r05(8)
  eri_value(25) = -(qx*qx*qx*qz*r01(25)) - qx*qx*qx*r02(41) + qx*qx*qz*(-r02(37) - 2.0_dp*r02(43)) - 3.0_dp*r03(23) + qx*qx*(-r03(33) - 2.0_dp*r03(43)) + qx*qz*(-(2.0_dp*r01(19)) - r01(22) - 2.0_dp*r03(41) - r03(51)) + qz*(-(3.0_dp*r02(31)) - r04(31)) + qx*(-(2.0_dp*r02(29)) - r02(35) - 2.0_dp*r04(18) - r04(33)) - r05(3)
          eri_value(26) = -r03(23) + qx*qz*(-r01(22) - r03(54)) + qz*(-r02(31) - r04(34)) + qx*(-r02(35) - r04(38)) - r05(8)
  eri_value(27) = -(qx*qz*qz*qz*r01(25)) - qz*qz*qz*r02(37) + qx*qz*qz*(-r02(41) - 2.0_dp*r02(47)) - 3.0_dp*r03(23) + qz*qz*(-r03(33) - 2.0_dp*r03(43)) + qx*qz*(-(2.0_dp*r01(19)) - r01(22) - 2.0_dp*r03(46) - r03(56)) + qz*(-(2.0_dp*r02(25)) - r02(31) - 2.0_dp*r04(21) - r04(36)) + qx*(-(3.0_dp*r02(35)) - r04(40)) - r05(10)
  eri_value(28) = -(qx*qx*qz*r02(46)) - r03(25) - qx*qx*r03(45) + qx*qz*(-r03(42) - r03(52)) + qz*(-r02(34) - r04(32)) + qx*(-r04(20) - r04(35)) - r05(5)
  eri_value(29) = -r01(13) - qx*qx*qz*qz*r01(25) + qx*qz*qz*(-r02(37) - r02(43)) + qx*qx*qz*(-r02(41) - r02(47)) - r03(21) - r03(26) + qz*qz*(-r01(19) - r03(41)) + qx*qx*(-r01(19) - r03(46)) + qx*qz*(-r03(33) - 2.0_dp*r03(43) - r03(53)) + qz*(-r02(29) - r02(35) - r04(18) - r04(33)) + qx*(-r02(25) - r02(31) - r04(21) - r04(36)) - r05(6)
  eri_value(30) = -(qx*qz*qz*r02(46)) - r03(22) - qz*qz*r03(42) + qx*qz*(-r03(45) - r03(55)) + qz*(-r04(20) - r04(35)) + qx*(-r02(34) - r04(39)) - r05(9)
  eri_value(31) = -(qx*qx*qz*r02(40)) - r03(25) - qx*qx*r03(35) - 2.0_dp*qx*qz*r03(42) - 2.0_dp*qx*r04(20) + qz*(-r02(34) - r04(32)) - r05(5)
          eri_value(32) = -(3.0_dp*r03(25)) + qz*(-(3.0_dp*r02(34)) - r04(37)) - r05(12)
  eri_value(33) = -(qz*qz*qz*r02(40)) - 3.0_dp*r03(25) + qz*qz*(-r03(35) - 2.0_dp*r03(45)) + qz*(-(2.0_dp*r02(28)) - r02(34) - 2.0_dp*r04(24) - r04(39)) - r05(14)
          eri_value(34) = -r03(23) + qx*qz*(-r01(19) - r03(44)) + qx*(-r02(29) - r04(23)) + qz*(-r02(31) - r04(34)) - r05(8)
  eri_value(35) = -(qx*qz*qz*r02(40)) - r03(22) - qz*qz*r03(42) + qx*qz*(-r03(35) - r03(45)) + qx*(-r02(28) - r04(24)) + qz*(-r04(20) - r04(35)) - r05(9)
   eri_value(36) = -r01(13) - r03(24) - r03(26) + qz*qz*(-r01(19) - r03(44)) + qz*(-r02(29) - r02(35) - r04(23) - r04(38)) - r05(13)
  eri_value(37) = -(3.0_dp*r01(14)) - qx*qx*qx*qx*r01(26) + 2.0_dp*qx*qx*qx*(-r02(40) - r02(46)) - 6.0_dp*r03(22) + qx*qx*(-r01(17) - 4.0_dp*r01(20) - r01(23) - r03(32) - 4.0_dp*r03(42) - r03(52)) + 2.0_dp*qx*(3.0_dp*(-r02(28) - r02(34)) - r04(17) - r04(32)) - r05(2)
          eri_value(38) = -r01(14) - r03(22) - r03(27) + qx*qx*(-r01(23) - r03(57)) + 2.0_dp*qx*(-r02(34) - r04(37)) - r05(7)
  eri_value(39) = -r01(14) - qx*qx*qz*qz*r01(26) - 2.0_dp*qx*qz*qz*r02(40) - 2.0_dp*qx*qx*qz*r02(48) - r03(22) - r03(29) + qz*qz*(-r01(17) - r03(32)) - 4*qx*qz*r03(45) + qx*qx*(-r01(23) - r03(59)) + 2.0_dp*qz*(-r02(30) - r04(20)) + 2.0_dp*qx*(-r02(34) - r04(39)) - r05(9)
  eri_value(40) = -(qx*qx*qx*r02(44)) - 3.0_dp*r03(24) + qx*qx*(-(2.0_dp*r03(44)) - r03(54)) + qx*(-r02(26) - 2.0_dp*r02(32) - r04(19) - 2.0_dp*r04(34)) - r05(4)
  eri_value(41) = -(qx*qx*qx*qz*r01(26)) + qx*qx*qz*(-(2.0_dp*r02(40)) - r02(46)) - qx*qx*qx*r02(48) - 3.0_dp*r03(25) + qx*qz*(-r01(17) - 2.0_dp*r01(20) - r03(32) - 2.0_dp*r03(42)) + qx*qx*(-(2.0_dp*r03(45)) - r03(55)) + qz*(-(3.0_dp*r02(28)) - r04(17)) + qx*(-r02(30) - 2.0_dp*r02(36) - r04(20) - 2.0_dp*r04(35)) - r05(5)
  eri_value(42) = -(qx*qx*qz*r02(44)) - r03(28) - 2.0_dp*qx*qz*r03(44) - qx*qx*r03(58) + qz*(-r02(26) - r04(19)) - 2.0_dp*qx*r04(38) - r05(8)
          eri_value(43) = -r01(14) - r03(22) - r03(27) + qx*qx*(-r01(17) - r03(37)) + 2.0_dp*qx*(-r02(28) - r04(22)) - r05(7)
          eri_value(44) = -(3.0_dp*r01(14)) - 6.0_dp*r03(27) - r05(16)
          eri_value(45) = -r01(14) - r03(27) - r03(29) + qz*qz*(-r01(17) - r03(37)) + 2.0_dp*qz*(-r02(30) - r04(27)) - r05(18)
          eri_value(46) = -(3.0_dp*r03(24)) + qx*(-(3.0_dp*r02(26)) - r04(26)) - r05(11)
          eri_value(47) = -r03(25) + qx*qz*(-r01(17) - r03(37)) + qz*(-r02(28) - r04(22)) + qx*(-r02(30) - r04(27)) - r05(12)
          eri_value(48) = -(3.0_dp*r03(28)) + qz*(-(3.0_dp*r02(26)) - r04(26)) - r05(17)
  eri_value(49) = -r01(14) - qx*qx*qz*qz*r01(26) - 2.0_dp*qx*qx*qz*r02(42) - 2.0_dp*qx*qz*qz*r02(46) - r03(22) - r03(29) + qx*qx*(-r01(17) - r03(39)) - 4*qx*qz*r03(45) + qz*qz*(-r01(23) - r03(52)) + 2.0_dp*qx*(-r02(28) - r04(24)) + 2.0_dp*qz*(-r02(36) - r04(35)) - r05(9)
          eri_value(50) = -r01(14) - r03(27) - r03(29) + qz*qz*(-r01(23) - r03(57)) + 2.0_dp*qz*(-r02(36) - r04(42)) - r05(18)
  eri_value(51) = -(3.0_dp*r01(14)) - qz*qz*qz*qz*r01(26) + 2.0_dp*qz*qz*qz*(-r02(42) - r02(48)) - 6.0_dp*r03(29) + qz*qz*(-r01(17) - 4.0_dp*r01(20) - r01(23) - r03(39) - 4.0_dp*r03(49) - r03(59)) + 2.0_dp*qz*(3.0_dp*(-r02(30) - r02(36)) - r04(29) - r04(44)) - r05(20)
  eri_value(52) = -(qx*qz*qz*r02(44)) - r03(24) - 2.0_dp*qx*qz*r03(48) - qz*qz*r03(54) + qx*(-r02(26) - r04(28)) - 2.0_dp*qz*r04(38) - r05(13)
  eri_value(53) = -(qx*qz*qz*qz*r01(26)) - qz*qz*qz*r02(46) + qx*qz*qz*(-(2.0_dp*r02(42)) - r02(48)) - 3.0_dp*r03(25) + qx*qz*(-r01(17) - 2.0_dp*r01(20) - r03(39) - 2.0_dp*r03(49)) + qz*qz*(-(2.0_dp*r03(45)) - r03(55)) + qx*(-(3.0_dp*r02(30)) - r04(29)) + qz*(-r02(28) - 2.0_dp*r02(34) - r04(24) - 2.0_dp*r04(39)) - r05(14)
  eri_value(54) = -(qz*qz*qz*r02(44)) - 3.0_dp*r03(28) + qz*qz*(-(2.0_dp*r03(48)) - r03(58)) + qz*(-r02(26) - 2.0_dp*r02(32) - r04(28) - 2.0_dp*r04(43)) - r05(19)
  eri_value(55) = -(qx*qx*qx*r02(38)) - 3.0_dp*r03(24) + qx*qx*(-r03(34) - 2.0_dp*r03(44)) + qx*(-(2.0_dp*r02(26)) - r02(32) - 2.0_dp*r04(19) - r04(34)) - r05(4)
          eri_value(56) = -(3.0_dp*r03(24)) + qx*(-(3.0_dp*r02(32)) - r04(41)) - r05(11)
  eri_value(57) = -(qx*qz*qz*r02(38)) - r03(24) - qz*qz*r03(34) - 2.0_dp*qx*qz*r03(48) - 2.0_dp*qz*r04(23) + qx*(-r02(32) - r04(43)) - r05(13)
    eri_value(58) = -r01(14) - r03(22) - r03(27) + qx*qx*(-r01(20) - r03(47)) + qx*(-r02(28) - r02(34) - r04(22) - r04(37)) - r05(7)
  eri_value(59) = -(qx*qx*qz*r02(38)) - r03(28) + qx*qz*(-r03(34) - r03(44)) - qx*qx*r03(48) + qz*(-r02(26) - r04(19)) + qx*(-r04(23) - r04(38)) - r05(8)
          eri_value(60) = -r03(25) + qx*qz*(-r01(20) - r03(47)) + qz*(-r02(28) - r04(22)) + qx*(-r02(36) - r04(42)) - r05(12)
  eri_value(61) = -(qx*qx*qx*qz*r01(26)) - qx*qx*qx*r02(42) + qx*qx*qz*(-r02(40) - 2.0_dp*r02(46)) - 3.0_dp*r03(25) + qx*qx*(-r03(35) - 2.0_dp*r03(45)) + qx*qz*(-(2.0_dp*r01(20)) - r01(23) - 2.0_dp*r03(42) - r03(52)) + qz*(-(3.0_dp*r02(34)) - r04(32)) + qx*(-(2.0_dp*r02(30)) - r02(36) - 2.0_dp*r04(20) - r04(35)) - r05(5)
          eri_value(62) = -r03(25) + qx*qz*(-r01(23) - r03(57)) + qz*(-r02(34) - r04(37)) + qx*(-r02(36) - r04(42)) - r05(12)
  eri_value(63) = -(qx*qz*qz*qz*r01(26)) - qz*qz*qz*r02(40) + qx*qz*qz*(-r02(42) - 2.0_dp*r02(48)) - 3.0_dp*r03(25) + qz*qz*(-r03(35) - 2.0_dp*r03(45)) + qx*qz*(-(2.0_dp*r01(20)) - r01(23) - 2.0_dp*r03(49) - r03(59)) + qz*(-(2.0_dp*r02(28)) - r02(34) - 2.0_dp*r04(24) - r04(39)) + qx*(-(3.0_dp*r02(36)) - r04(44)) - r05(14)
  eri_value(64) = -(qx*qx*qz*r02(44)) - r03(28) - qx*qx*r03(48) + qx*qz*(-r03(44) - r03(54)) + qz*(-r02(32) - r04(34)) + qx*(-r04(23) - r04(38)) - r05(8)
  eri_value(65) = -r01(14) - qx*qx*qz*qz*r01(26) + qx*qz*qz*(-r02(40) - r02(46)) + qx*qx*qz*(-r02(42) - r02(48)) - r03(22) - r03(29) + qz*qz*(-r01(20) - r03(42)) + qx*qx*(-r01(20) - r03(49)) + qx*qz*(-r03(35) - 2.0_dp*r03(45) - r03(55)) + qz*(-r02(30) - r02(36) - r04(20) - r04(35)) + qx*(-r02(28) - r02(34) - r04(24) - r04(39)) - r05(9)
  eri_value(66) = -(qx*qz*qz*r02(44)) - r03(24) - qz*qz*r03(44) + qx*qz*(-r03(48) - r03(58)) + qz*(-r04(23) - r04(38)) + qx*(-r02(32) - r04(43)) - r05(13)
  eri_value(67) = -(qx*qx*qz*r02(38)) - r03(28) - qx*qx*r03(38) - 2.0_dp*qx*qz*r03(44) - 2.0_dp*qx*r04(23) + qz*(-r02(32) - r04(34)) - r05(8)
          eri_value(68) = -(3.0_dp*r03(28)) + qz*(-(3.0_dp*r02(32)) - r04(41)) - r05(17)
  eri_value(69) = -(qz*qz*qz*r02(38)) - 3.0_dp*r03(28) + qz*qz*(-r03(38) - 2.0_dp*r03(48)) + qz*(-(2.0_dp*r02(26)) - r02(32) - 2.0_dp*r04(28) - r04(43)) - r05(19)
          eri_value(70) = -r03(25) + qx*qz*(-r01(20) - r03(47)) + qx*(-r02(30) - r04(27)) + qz*(-r02(34) - r04(37)) - r05(12)
  eri_value(71) = -(qx*qz*qz*r02(38)) - r03(24) - qz*qz*r03(44) + qx*qz*(-r03(38) - r03(48)) + qx*(-r02(26) - r04(28)) + qz*(-r04(23) - r04(38)) - r05(13)
   eri_value(72) = -r01(14) - r03(27) - r03(29) + qz*qz*(-r01(20) - r03(47)) + qz*(-r02(30) - r02(36) - r04(27) - r04(42)) - r05(18)
  eri_value(73) = 3.0_dp*(r00(1) - r01(15)) + qx*qx*qx*qx*(r00(5) - r01(27)) + 2.0_dp*qx*qx*qx*(r01(7) + r01(10) - r02(41) - r02(47)) + 6.0_dp*(r02(1) - r03(23)) + qx*qx*(r00(2) + r00(4) - r01(18) + 4.0_dp*(r00(3) - r01(21)) - r01(24) + r02(7) + r02(19) - r03(33) + 4.0_dp*(r02(13) - r03(43)) - r03(53)) + r04(1) + 2.0_dp*qx*(3.0_dp*(r01(1) + r01(4) - r02(29) - r02(35)) + r03(1) + r03(11) - r04(18) - r04(33)) - r05(3)
  eri_value(74) = r00(1) - r01(15) + r02(1) + r02(2) - r03(23) - r03(28) + qx*qx*(r00(4) - r01(24) + r02(20) - r03(58)) + r04(4) + 2.0_dp*qx*(r01(4) - r02(35) + r03(14) - r04(38)) - r05(8)
  eri_value(75) = r00(1) - r01(15) + qx*qx*qz*qz*(r00(5) - r01(27)) + r02(1) + r02(3) + 2.0_dp*qx*qz*qz*(r01(7) - r02(41)) + 2.0_dp*qx*qx*qz*(r01(12) - r02(45)) - r03(23) - r03(30) + qz*qz*(r00(2) - r01(18) + r02(7) - r03(33)) + 4*qx*qz*(r02(17) - r03(46)) + qx*qx*(r00(4) - r01(24) + r02(21) - r03(60)) + r04(6) + 2.0_dp*qz*(r01(3) - r02(27) + r03(3) - r04(21)) + 2.0_dp*qx*(r01(4) - r02(35) + r03(16) - r04(40)) - r05(10)
  eri_value(76) = qx*qx*qx*(r01(11) - r02(48)) + 3.0_dp*(r02(4) - r03(25)) + qx*qx*(r02(22) + 2.0_dp*(r02(16) - r03(45)) - r03(55)) + r04(2) + qx*(r01(2) - r02(30) + 2.0_dp*(r01(5) - r02(36)) + r03(2) - r04(20) + 2.0_dp*(r03(12) - r04(35))) - r05(5)
  eri_value(77) = qx*qx*qx*qz*(r00(5) - r01(27)) + qx*qx*qx*(r01(12) - r02(45)) + qx*qx*qz*(r01(10) + 2.0_dp*(r01(7) - r02(41)) - r02(47)) + 3.0_dp*(r02(5) - r03(26)) + qx*qz*(r00(2) - r01(18) + 2.0_dp*(r00(3) - r01(21)) + r02(7) - r03(33) + 2.0_dp*(r02(13) - r03(43))) + qx*qx*(r02(23) + 2.0_dp*(r02(17) - r03(46)) - r03(56)) + r04(3) + qz*(3.0_dp*(r01(1) - r02(29)) + r03(1) - r04(18)) + qx*(r01(3) - r02(27) + 2.0_dp*(r01(6) - r02(33)) + r03(3) - r04(21) + 2.0_dp*(r03(13) - r04(36))) - r05(6)
  eri_value(78) = r02(6) + qx*qx*qz*(r01(11) - r02(48)) - r03(29) + 2.0_dp*qx*qz*(r02(16) - r03(45)) + qx*qx*(r02(24) - r03(59)) + r04(5) + qz*(r01(2) - r02(30) + r03(2) - r04(20)) + 2.0_dp*qx*(r03(15) - r04(39)) - r05(9)
  eri_value(79) = r00(1) - r01(15) + r02(1) + r02(2) - r03(23) - r03(28) + qx*qx*(r00(2) - r01(18) + r02(8) - r03(38)) + r04(4) + 2.0_dp*qx*(r01(1) - r02(29) + r03(4) - r04(23)) - r05(8)
          eri_value(80) = 3.0_dp*(r00(1) - r01(15)) + 6.0_dp*(r02(2) - r03(28)) + r04(11) - r05(17)
  eri_value(81) = r00(1) - r01(15) + r02(2) + r02(3) - r03(28) - r03(30) + qz*qz*(r00(2) - r01(18) + r02(8) - r03(38)) + r04(13) + 2.0_dp*qz*(r01(3) - r02(27) + r03(8) - r04(28)) - r05(19)
          eri_value(82) = 3.0_dp*(r02(4) - r03(25)) + r04(7) + qx*(3.0_dp*(r01(2) - r02(30)) + r03(7) - r04(27)) - r05(12)
  eri_value(83) = r02(5) - r03(26) + qx*qz*(r00(2) - r01(18) + r02(8) - r03(38)) + r04(8) + qz*(r01(1) - r02(29) + r03(4) - r04(23)) + qx*(r01(3) - r02(27) + r03(8) - r04(28)) - r05(13)
          eri_value(84) = 3.0_dp*(r02(6) - r03(29)) + r04(12) + qz*(3.0_dp*(r01(2) - r02(30)) + r03(7) - r04(27)) - r05(18)
  eri_value(85) = r00(1) - r01(15) + qx*qx*qz*qz*(r00(5) - r01(27)) + r02(1) + r02(3) + 2.0_dp*qx*qx*qz*(r01(9) - r02(39)) + 2.0_dp*qx*qz*qz*(r01(10) - r02(47)) - r03(23) - r03(30) + qx*qx*(r00(2) - r01(18) + r02(9) - r03(40)) + 4*qx*qz*(r02(17) - r03(46)) + qz*qz*(r00(4) - r01(24) + r02(19) - r03(53)) + r04(6) + 2.0_dp*qx*(r01(1) - r02(29) + r03(6) - r04(25)) + 2.0_dp*qz*(r01(6) - r02(33) + r03(13) - r04(36)) - r05(10)
  eri_value(86) = r00(1) - r01(15) + r02(2) + r02(3) - r03(28) - r03(30) + qz*qz*(r00(4) - r01(24) + r02(20) - r03(58)) + r04(13) + 2.0_dp*qz*(r01(6) - r02(33) + r03(18) - r04(43)) - r05(19)
  eri_value(87) = 3.0_dp*(r00(1) - r01(15)) + qz*qz*qz*qz*(r00(5) - r01(27)) + 2.0_dp*qz*qz*qz*(r01(9) + r01(12) - r02(39) - r02(45)) + 6.0_dp*(r02(3) - r03(30)) + qz*qz*(r00(2) + r00(4) - r01(18) + 4.0_dp*(r00(3) - r01(21)) - r01(24) + r02(9) + r02(21) - r03(40) + 4.0_dp*(r02(15) - r03(50)) - r03(60)) + r04(15) + 2.0_dp*qz*(3.0_dp*(r01(3) + r01(6) - r02(27) - r02(33)) + r03(10) + r03(20) - r04(30) - r04(45)) - r05(21)
  eri_value(88) = r02(4) + qx*qz*qz*(r01(11) - r02(48)) - r03(25) + 2.0_dp*qx*qz*(r02(18) - r03(49)) + qz*qz*(r02(22) - r03(55)) + r04(9) + qx*(r01(2) - r02(30) + r03(9) - r04(29)) + 2.0_dp*qz*(r03(15) - r04(39)) - r05(14)
  eri_value(89) = qx*qz*qz*qz*(r00(5) - r01(27)) + qx*qz*qz*(r01(12) + 2.0_dp*(r01(9) - r02(39)) - r02(45)) + qz*qz*qz*(r01(10) - r02(47)) + 3.0_dp*(r02(5) - r03(26)) + qx*qz*(r00(2) - r01(18) + 2.0_dp*(r00(3) - r01(21)) + r02(9) - r03(40) + 2.0_dp*(r02(15) - r03(50))) + qz*qz*(r02(23) + 2.0_dp*(r02(17) - r03(46)) - r03(56)) + r04(10) + qx*(3.0_dp*(r01(3) - r02(27)) + r03(10) - r04(30)) + qz*(r01(1) - r02(29) + 2.0_dp*(r01(4) - r02(35)) + r03(6) - r04(25) + 2.0_dp*(r03(16) - r04(40))) - r05(15)
  eri_value(90) = qz*qz*qz*(r01(11) - r02(48)) + 3.0_dp*(r02(6) - r03(29)) + qz*qz*(r02(24) + 2.0_dp*(r02(18) - r03(49)) - r03(59)) + r04(14) + qz*(r01(2) - r02(30) + 2.0_dp*(r01(5) - r02(36)) + r03(9) - r04(29) + 2.0_dp*(r03(19) - r04(44))) - r05(20)
  eri_value(91) = qx*qx*qx*(r01(8) - r02(42)) + 3.0_dp*(r02(4) - r03(25)) + qx*qx*(r02(10) - r03(35) + 2.0_dp*(r02(16) - r03(45))) + r04(2) + qx*(r01(5) + 2.0_dp*(r01(2) - r02(30)) - r02(36) + r03(12) + 2.0_dp*(r03(2) - r04(20)) - r04(35)) - r05(5)
          eri_value(92) = 3.0_dp*(r02(4) - r03(25)) + r04(7) + qx*(3.0_dp*(r01(5) - r02(36)) + r03(17) - r04(42)) - r05(12)
  eri_value(93) = r02(4) + qx*qz*qz*(r01(8) - r02(42)) - r03(25) + qz*qz*(r02(10) - r03(35)) + 2.0_dp*qx*qz*(r02(18) - r03(49)) + r04(9) + 2.0_dp*qz*(r03(5) - r04(24)) + qx*(r01(5) - r02(36) + r03(19) - r04(44)) - r05(14)
  eri_value(94) = r00(1) - r01(15) + r02(1) + r02(2) - r03(23) - r03(28) + qx*qx*(r00(3) - r01(21) + r02(14) - r03(48)) + r04(4) + qx*(r01(1) + r01(4) - r02(29) - r02(35) + r03(4) + r03(14) - r04(23) - r04(38)) - r05(8)
  eri_value(95) = r02(6) + qx*qx*qz*(r01(8) - r02(42)) - r03(29) + qx*qz*(r02(10) + r02(16) - r03(35) - r03(45)) + qx*qx*(r02(18) - r03(49)) + r04(5) + qz*(r01(2) - r02(30) + r03(2) - r04(20)) + qx*(r03(5) + r03(15) - r04(24) - r04(39)) - r05(9)
  eri_value(96) = r02(5) - r03(26) + qx*qz*(r00(3) - r01(21) + r02(14) - r03(48)) + r04(8) + qz*(r01(1) - r02(29) + r03(4) - r04(23)) + qx*(r01(6) - r02(33) + r03(18) - r04(43)) - r05(13)
  eri_value(97) = qx*qx*qx*qz*(r00(5) - r01(27)) + qx*qx*qx*(r01(9) - r02(39)) + qx*qx*qz*(r01(7) - r02(41) + 2.0_dp*(r01(10) - r02(47))) + 3.0_dp*(r02(5) - r03(26)) + qx*qx*(r02(11) - r03(36) + 2.0_dp*(r02(17) - r03(46))) + qx*qz*(r00(4) + 2.0_dp*(r00(3) - r01(21)) - r01(24) + r02(19) + 2.0_dp*(r02(13) - r03(43)) - r03(53)) + r04(3) + qz*(3.0_dp*(r01(4) - r02(35)) + r03(11) - r04(33)) + qx*(r01(6) + 2.0_dp*(r01(3) - r02(27)) - r02(33) + r03(13) + 2.0_dp*(r03(3) - r04(21)) - r04(36)) - r05(6)
  eri_value(98) = r02(5) - r03(26) + qx*qz*(r00(4) - r01(24) + r02(20) - r03(58)) + r04(8) + qz*(r01(4) - r02(35) + r03(14) - r04(38)) + qx*(r01(6) - r02(33) + r03(18) - r04(43)) - r05(13)
  eri_value(99) = qx*qz*qz*qz*(r00(5) - r01(27)) + qz*qz*qz*(r01(7) - r02(41)) + qx*qz*qz*(r01(9) - r02(39) + 2.0_dp*(r01(12) - r02(45))) + 3.0_dp*(r02(5) - r03(26)) + qz*qz*(r02(11) - r03(36) + 2.0_dp*(r02(17) - r03(46))) + qx*qz*(r00(4) + 2.0_dp*(r00(3) - r01(21)) - r01(24) + r02(21) + 2.0_dp*(r02(15) - r03(50)) - r03(60)) + r04(10) + qz*(r01(4) + 2.0_dp*(r01(1) - r02(29)) - r02(35) + r03(16) + 2.0_dp*(r03(6) - r04(25)) - r04(40)) + qx*(3.0_dp*(r01(6) - r02(33)) + r03(20) - r04(45)) - r05(15)
  eri_value(100) = r02(6) + qx*qx*qz*(r01(11) - r02(48)) - r03(29) + qx*qx*(r02(18) - r03(49)) + qx*qz*(r02(16) + r02(22) - r03(45) - r03(55)) + r04(5) + qz*(r01(5) - r02(36) + r03(12) - r04(35)) + qx*(r03(5) + r03(15) - r04(24) - r04(39)) - r05(9)
  eri_value(101) = r00(1) - r01(15) + qx*qx*qz*qz*(r00(5) - r01(27)) + r02(1) + r02(3) + qx*qx*qz*(r01(9) + r01(12) - r02(39) - r02(45)) + qx*qz*qz*(r01(7) + r01(10) - r02(41) - r02(47)) - r03(23) - r03(30) + qz*qz*(r00(3) - r01(21) + r02(13) - r03(43)) + qx*qx*(r00(3) - r01(21) + r02(15) - r03(50)) + qx*qz*(r02(11) + r02(23) - r03(36) + 2.0_dp*(r02(17) - r03(46)) - r03(56)) + r04(6) + qz*(r01(3) + r01(6) - r02(27) - r02(33) + r03(3) + r03(13) - r04(21) - r04(36)) + qx*(r01(1) + r01(4) - r02(29) - r02(35) + r03(6) + r03(16) - r04(25) - r04(40)) - r05(10)
  eri_value(102) = r02(4) + qx*qz*qz*(r01(11) - r02(48)) - r03(25) + qz*qz*(r02(16) - r03(45)) + qx*qz*(r02(18) + r02(24) - r03(49) - r03(59)) + r04(9) + qz*(r03(5) + r03(15) - r04(24) - r04(39)) + qx*(r01(5) - r02(36) + r03(19) - r04(44)) - r05(14)
  eri_value(103) = r02(6) + qx*qx*qz*(r01(8) - r02(42)) - r03(29) + qx*qx*(r02(12) - r03(39)) + 2.0_dp*qx*qz*(r02(16) - r03(45)) + r04(5) + 2.0_dp*qx*(r03(5) - r04(24)) + qz*(r01(5) - r02(36) + r03(12) - r04(35)) - r05(9)
          eri_value(104) = 3.0_dp*(r02(6) - r03(29)) + r04(12) + qz*(3.0_dp*(r01(5) - r02(36)) + r03(17) - r04(42)) - r05(18)
  eri_value(105) = qz*qz*qz*(r01(8) - r02(42)) + 3.0_dp*(r02(6) - r03(29)) + qz*qz*(r02(12) - r03(39) + 2.0_dp*(r02(18) - r03(49))) + r04(14) + qz*(r01(5) + 2.0_dp*(r01(2) - r02(30)) - r02(36) + r03(19) + 2.0_dp*(r03(9) - r04(29)) - r04(44)) - r05(20)
  eri_value(106) = r02(5) - r03(26) + qx*qz*(r00(3) - r01(21) + r02(14) - r03(48)) + r04(8) + qx*(r01(3) - r02(27) + r03(8) - r04(28)) + qz*(r01(4) - r02(35) + r03(14) - r04(38)) - r05(13)
  eri_value(107) = r02(4) + qx*qz*qz*(r01(8) - r02(42)) - r03(25) + qz*qz*(r02(16) - r03(45)) + qx*qz*(r02(12) + r02(18) - r03(39) - r03(49)) + r04(9) + qx*(r01(2) - r02(30) + r03(9) - r04(29)) + qz*(r03(5) + r03(15) - r04(24) - r04(39)) - r05(14)
  eri_value(108) = r00(1) - r01(15) + r02(2) + r02(3) - r03(28) - r03(30) + qz*qz*(r00(3) - r01(21) + r02(14) - r03(48)) + r04(13) + qz*(r01(3) + r01(6) - r02(27) - r02(33) + r03(8) + r03(18) - r04(28) - r04(43)) - r05(19)
          ! write(*,*) eri value 8, eri_value(8)
          ! ! write(*,*) r05(11), r05(11)
          ! write(*,*) r01(5), r01(13)
          ! write(*,*) r03(21), r03(24)
          ! eri_value(8) = -(3.0_dp*r01(5)) - 6.0_dp*r03(21) - r05(11)
          ! do i = 1,108
          !     write(*,*) eri_value(i)
          ! enddo
          buff(10) = (buff(8)*buff(3) - buff(9)*buff(2))
          buff(11) = (buff(9)*buff(1) - buff(7)*buff(3))
          buff(12) = (buff(7)*buff(2) - buff(8)*buff(1))

          trans(1) = eri_value(1)
          trans(2) = eri_value(37)
          trans(3) = eri_value(73)
          eri_value(1) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(37) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(73) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(2)
          trans(2) = eri_value(38)
          trans(3) = eri_value(74)
          eri_value(2) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(38) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(74) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(3)
          trans(2) = eri_value(39)
          trans(3) = eri_value(75)
          eri_value(3) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(39) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(75) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(4)
          trans(2) = eri_value(40)
          trans(3) = eri_value(76)
          eri_value(4) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(40) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(76) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(5)
          trans(2) = eri_value(41)
          trans(3) = eri_value(77)
          eri_value(5) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(41) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(77) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(6)
          trans(2) = eri_value(42)
          trans(3) = eri_value(78)
          eri_value(6) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(42) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(78) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(7)
          trans(2) = eri_value(43)
          trans(3) = eri_value(79)
          eri_value(7) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(43) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(79) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(8)
          trans(2) = eri_value(44)
          trans(3) = eri_value(80)
          eri_value(8) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(44) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(80) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(9)
          trans(2) = eri_value(45)
          trans(3) = eri_value(81)
          eri_value(9) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(45) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(81) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(10)
          trans(2) = eri_value(46)
          trans(3) = eri_value(82)
          eri_value(10) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(46) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(82) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(11)
          trans(2) = eri_value(47)
          trans(3) = eri_value(83)
          eri_value(11) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(47) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(83) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(12)
          trans(2) = eri_value(48)
          trans(3) = eri_value(84)
          eri_value(12) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(48) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(84) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(13)
          trans(2) = eri_value(49)
          trans(3) = eri_value(85)
          eri_value(13) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(49) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(85) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(14)
          trans(2) = eri_value(50)
          trans(3) = eri_value(86)
          eri_value(14) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(50) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(86) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(15)
          trans(2) = eri_value(51)
          trans(3) = eri_value(87)
          eri_value(15) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(51) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(87) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(16)
          trans(2) = eri_value(52)
          trans(3) = eri_value(88)
          eri_value(16) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(52) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(88) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(17)
          trans(2) = eri_value(53)
          trans(3) = eri_value(89)
          eri_value(17) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(53) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(89) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(18)
          trans(2) = eri_value(54)
          trans(3) = eri_value(90)
          eri_value(18) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(54) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(90) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(19)
          trans(2) = eri_value(55)
          trans(3) = eri_value(91)
          eri_value(19) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(55) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(91) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(20)
          trans(2) = eri_value(56)
          trans(3) = eri_value(92)
          eri_value(20) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(56) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(92) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(21)
          trans(2) = eri_value(57)
          trans(3) = eri_value(93)
          eri_value(21) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(57) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(93) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(22)
          trans(2) = eri_value(58)
          trans(3) = eri_value(94)
          eri_value(22) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(58) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(94) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(23)
          trans(2) = eri_value(59)
          trans(3) = eri_value(95)
          eri_value(23) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(59) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(95) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(24)
          trans(2) = eri_value(60)
          trans(3) = eri_value(96)
          eri_value(24) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(60) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(96) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(25)
          trans(2) = eri_value(61)
          trans(3) = eri_value(97)
          eri_value(25) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(61) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(97) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(26)
          trans(2) = eri_value(62)
          trans(3) = eri_value(98)
          eri_value(26) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(62) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(98) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(27)
          trans(2) = eri_value(63)
          trans(3) = eri_value(99)
          eri_value(27) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(63) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(99) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(28)
          trans(2) = eri_value(64)
          trans(3) = eri_value(100)
          eri_value(28) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(64) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(100) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(29)
          trans(2) = eri_value(65)
          trans(3) = eri_value(101)
          eri_value(29) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(65) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(101) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(30)
          trans(2) = eri_value(66)
          trans(3) = eri_value(102)
          eri_value(30) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(66) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(102) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(31)
          trans(2) = eri_value(67)
          trans(3) = eri_value(103)
          eri_value(31) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(67) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(103) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(32)
          trans(2) = eri_value(68)
          trans(3) = eri_value(104)
          eri_value(32) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(68) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(104) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(33)
          trans(2) = eri_value(69)
          trans(3) = eri_value(105)
          eri_value(33) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(69) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(105) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(34)
          trans(2) = eri_value(70)
          trans(3) = eri_value(106)
          eri_value(34) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(70) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(106) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(35)
          trans(2) = eri_value(71)
          trans(3) = eri_value(107)
          eri_value(35) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(71) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(107) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(36)
          trans(2) = eri_value(72)
          trans(3) = eri_value(108)
          eri_value(36) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(72) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(108) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)

          trans(1) = eri_value(1)
          trans(2) = eri_value(7)
          trans(3) = eri_value(13)
          trans(4) = eri_value(19)
          trans(5) = eri_value(25)
          trans(6) = eri_value(31)
          ! write(*,*) "trans", trans(1)
          ! write(*,*) "trans", trans(2)
          ! write(*,*) "trans", trans(3)
          ! write(*,*) "trans", trans(4)
          ! write(*,*) "trans", trans(5)
          ! write(*,*) "trans", trans(6)
  eri_value(1)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(7)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(13)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(19)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(25)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(31)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          ! write(*,*) "eri_value", eri_value(1)
          ! write(*,*) "eri_value", eri_value(7)
          ! write(*,*) "eri_value", eri_value(13)
          ! write(*,*) "eri_value", eri_value(19)
          ! write(*,*) "eri_value", eri_value(25)
          ! write(*,*) "eri_value", eri_value(31)
          trans(1) = eri_value(2)
          trans(2) = eri_value(8)
          trans(3) = eri_value(14)
          trans(4) = eri_value(20)
          trans(5) = eri_value(26)
          trans(6) = eri_value(32)
          ! write(*,*) "trans", trans(1)
          ! write(*,*) "trans", trans(2)
          ! write(*,*) "trans", trans(3)
          ! write(*,*) "trans", trans(4)
          ! write(*,*) "trans", trans(5)
          ! write(*,*) "trans", trans(6)
  eri_value(2)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(8)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(14)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(20)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(26)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(32)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          ! write(*,*) "eri_value", eri_value(2)
          ! write(*,*) "eri_value", eri_value(8)
          ! write(*,*) "eri_value", eri_value(14)
          ! write(*,*) "eri_value", eri_value(20)
          ! write(*,*) "eri_value", eri_value(26)
          ! write(*,*) "eri_value", eri_value(32)
          trans(1) = eri_value(3)
          trans(2) = eri_value(9)
          trans(3) = eri_value(15)
          trans(4) = eri_value(21)
          trans(5) = eri_value(27)
          trans(6) = eri_value(33)
          ! write(*,*) "trans", trans(1)
          ! write(*,*) "trans", trans(2)
          ! write(*,*) "trans", trans(3)
          ! write(*,*) "trans", trans(4)
          ! write(*,*) "trans", trans(5)
          ! write(*,*) "trans", trans(6)
  eri_value(3)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(9)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(15)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(21)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(27)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(33)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          ! write(*,*) "eri_value", eri_value(3)
          ! write(*,*) "eri_value", eri_value(9)
          ! write(*,*) "eri_value", eri_value(15)
          ! write(*,*) "eri_value", eri_value(22)
          ! write(*,*) "eri_value", eri_value(27)
          ! write(*,*) "eri_value", eri_value(33)
          trans(1) = eri_value(4)
          trans(2) = eri_value(10)
          trans(3) = eri_value(16)
          trans(4) = eri_value(22)
          trans(5) = eri_value(28)
          trans(6) = eri_value(34)
  eri_value(4)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(10)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(16)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(22)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(28)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(34)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(5)
          trans(2) = eri_value(11)
          trans(3) = eri_value(17)
          trans(4) = eri_value(23)
          trans(5) = eri_value(29)
          trans(6) = eri_value(35)
  eri_value(5)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(11)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(17)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(23)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(29)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(35)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(6)
          trans(2) = eri_value(12)
          trans(3) = eri_value(18)
          trans(4) = eri_value(24)
          trans(5) = eri_value(30)
          trans(6) = eri_value(36)
  eri_value(6)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(12)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(18)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(24)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(30)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(36)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(37)
          trans(2) = eri_value(43)
          trans(3) = eri_value(49)
          trans(4) = eri_value(55)
          trans(5) = eri_value(61)
          trans(6) = eri_value(67)
  eri_value(37)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(43)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(49)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(55)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(61)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(67)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(38)
          trans(2) = eri_value(44)
          trans(3) = eri_value(50)
          trans(4) = eri_value(56)
          trans(5) = eri_value(62)
          trans(6) = eri_value(68)
  eri_value(38)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(44)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(50)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(56)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(62)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(68)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(39)
          trans(2) = eri_value(45)
          trans(3) = eri_value(51)
          trans(4) = eri_value(57)
          trans(5) = eri_value(63)
          trans(6) = eri_value(69)
  eri_value(39)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(45)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(51)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(57)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(63)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(69)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(40)
          trans(2) = eri_value(46)
          trans(3) = eri_value(52)
          trans(4) = eri_value(58)
          trans(5) = eri_value(64)
          trans(6) = eri_value(70)
  eri_value(40)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(46)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(52)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(58)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(64)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(70)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(41)
          trans(2) = eri_value(47)
          trans(3) = eri_value(53)
          trans(4) = eri_value(59)
          trans(5) = eri_value(65)
          trans(6) = eri_value(71)
  eri_value(41)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(47)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(53)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(59)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(65)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(71)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(42)
          trans(2) = eri_value(48)
          trans(3) = eri_value(54)
          trans(4) = eri_value(60)
          trans(5) = eri_value(66)
          trans(6) = eri_value(72)
  eri_value(42)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(48)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(54)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(60)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(66)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(72)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(73)
          trans(2) = eri_value(79)
          trans(3) = eri_value(85)
          trans(4) = eri_value(91)
          trans(5) = eri_value(97)
          trans(6) = eri_value(103)
  eri_value(73)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(79)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(85)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(91)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(97)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(103)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(74)
          trans(2) = eri_value(80)
          trans(3) = eri_value(86)
          trans(4) = eri_value(92)
          trans(5) = eri_value(98)
          trans(6) = eri_value(104)
  eri_value(74)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(80)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(86)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(92)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(98)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(104)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(75)
          trans(2) = eri_value(81)
          trans(3) = eri_value(87)
          trans(4) = eri_value(93)
          trans(5) = eri_value(99)
          trans(6) = eri_value(105)
  eri_value(75)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(81)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(87)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(93)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(99)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(105)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(76)
          trans(2) = eri_value(82)
          trans(3) = eri_value(88)
          trans(4) = eri_value(94)
          trans(5) = eri_value(100)
          trans(6) = eri_value(106)
  eri_value(76)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(82)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(88)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(94)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(100)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(106)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(77)
          trans(2) = eri_value(83)
          trans(3) = eri_value(89)
          trans(4) = eri_value(95)
          trans(5) = eri_value(101)
          trans(6) = eri_value(107)
  eri_value(77)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(83)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(89)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(95)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(101)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(107)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(78)
          trans(2) = eri_value(84)
          trans(3) = eri_value(90)
          trans(4) = eri_value(96)
          trans(5) = eri_value(102)
          trans(6) = eri_value(108)
  eri_value(78)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(84)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(90)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(96)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(102)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(108)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          ! do i = 1, 108
          !    write(*,*) "eri_value2", eri_value(i)
          ! enddo
          trans(1) = eri_value(1)
          trans(2) = eri_value(2)
          trans(3) = eri_value(3)
          trans(4) = eri_value(4)
          trans(5) = eri_value(5)
          trans(6) = eri_value(6)
  eri_value(1)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(2)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(3)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(4)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(5)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(6)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(7)
          trans(2) = eri_value(8)
          trans(3) = eri_value(9)
          trans(4) = eri_value(10)
          trans(5) = eri_value(11)
          trans(6) = eri_value(12)
  eri_value(7)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(8)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(9)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(10)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(11)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(12)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(13)
          trans(2) = eri_value(14)
          trans(3) = eri_value(15)
          trans(4) = eri_value(16)
          trans(5) = eri_value(17)
          trans(6) = eri_value(18)
  eri_value(13)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(14)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(15)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(16)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(17)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(18)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(19)
          trans(2) = eri_value(20)
          trans(3) = eri_value(21)
          trans(4) = eri_value(22)
          trans(5) = eri_value(23)
          trans(6) = eri_value(24)
  eri_value(19)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(20)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(21)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(22)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(23)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(24)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(25)
          trans(2) = eri_value(26)
          trans(3) = eri_value(27)
          trans(4) = eri_value(28)
          trans(5) = eri_value(29)
          trans(6) = eri_value(30)
  eri_value(25)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(26)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(27)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(28)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(29)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(30)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(31)
          trans(2) = eri_value(32)
          trans(3) = eri_value(33)
          trans(4) = eri_value(34)
          trans(5) = eri_value(35)
          trans(6) = eri_value(36)
  eri_value(31)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(32)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(33)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(34)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(35)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(36)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(37)
          trans(2) = eri_value(38)
          trans(3) = eri_value(39)
          trans(4) = eri_value(40)
          trans(5) = eri_value(41)
          trans(6) = eri_value(42)
  eri_value(37)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(38)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(39)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(40)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(41)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(42)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(43)
          trans(2) = eri_value(44)
          trans(3) = eri_value(45)
          trans(4) = eri_value(46)
          trans(5) = eri_value(47)
          trans(6) = eri_value(48)
  eri_value(43)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(44)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(45)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(46)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(47)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(48)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(49)
          trans(2) = eri_value(50)
          trans(3) = eri_value(51)
          trans(4) = eri_value(52)
          trans(5) = eri_value(53)
          trans(6) = eri_value(54)
  eri_value(49)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(50)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(51)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(52)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(53)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(54)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(55)
          trans(2) = eri_value(56)
          trans(3) = eri_value(57)
          trans(4) = eri_value(58)
          trans(5) = eri_value(59)
          trans(6) = eri_value(60)
  eri_value(55)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(56)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(57)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(58)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(59)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(60)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(61)
          trans(2) = eri_value(62)
          trans(3) = eri_value(63)
          trans(4) = eri_value(64)
          trans(5) = eri_value(65)
          trans(6) = eri_value(66)
  eri_value(61)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(62)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(63)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(64)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(65)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(66)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(67)
          trans(2) = eri_value(68)
          trans(3) = eri_value(69)
          trans(4) = eri_value(70)
          trans(5) = eri_value(71)
          trans(6) = eri_value(72)
  eri_value(67)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(68)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(69)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(70)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(71)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(72)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(73)
          trans(2) = eri_value(74)
          trans(3) = eri_value(75)
          trans(4) = eri_value(76)
          trans(5) = eri_value(77)
          trans(6) = eri_value(78)
  eri_value(73)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(74)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(75)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(76)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(77)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(78)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(79)
          trans(2) = eri_value(80)
          trans(3) = eri_value(81)
          trans(4) = eri_value(82)
          trans(5) = eri_value(83)
          trans(6) = eri_value(84)
  eri_value(79)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(80)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(81)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(82)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(83)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(84)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(85)
          trans(2) = eri_value(86)
          trans(3) = eri_value(87)
          trans(4) = eri_value(88)
          trans(5) = eri_value(89)
          trans(6) = eri_value(90)
  eri_value(85)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(86)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(87)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(88)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(89)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(90)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(91)
          trans(2) = eri_value(92)
          trans(3) = eri_value(93)
          trans(4) = eri_value(94)
          trans(5) = eri_value(95)
          trans(6) = eri_value(96)
  eri_value(91)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(92)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(93)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(94)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(95)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(96)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(97)
          trans(2) = eri_value(98)
          trans(3) = eri_value(99)
          trans(4) = eri_value(100)
          trans(5) = eri_value(101)
          trans(6) = eri_value(102)
  eri_value(97)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(98)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(99)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(100)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(101)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(102)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(103)
          trans(2) = eri_value(104)
          trans(3) = eri_value(105)
          trans(4) = eri_value(106)
          trans(5) = eri_value(107)
          trans(6) = eri_value(108)
  eri_value(103)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(104)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(105)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(106)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(107)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(108)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          ! do i = 1, 108
          !   write(*,*) "eri_value", eri_value(i)
          ! enddo

          mini = 1  !0
          maxi = 1  !1
          minj = 1  !2
          maxj = 3  !2
          mink = 1
          maxk = 6
          minl = 1
          maxl = 6
          same = ish .eq. ksh .and. jsh .eq. lsh
          iandj = ish .eq. jsh
          kandl = ksh .eq. lsh

          maxl2 = maxl
          maxj2 = maxj

          ii1 = res%atom_loc(ish)
          jj1 = res%atom_loc(jsh)
          kk1 = res%atom_loc(ksh)
          loci = res%atom_loc(ish) - 1
          locj = res%atom_loc(jsh) - 1
          locl = res%atom_loc(lsh) - 1
          lock = res%atom_loc(ksh) - 1
          nij = 0
          do i = mini, maxi
            if (iandj) maxj2 = i
            ii1 = i + loci
            ip = (i - 1)*0 + 1 !gpople index 64

            do j = minj, maxj2
              nij = nij + 1
              maxl2 = maxl
              jj1 = j + locj
              i2 = ii1
              j2 = jj1
              if (ii1 .lt. jj1) then ! sort <ij|
                i2 = jj1
                j2 = ii1
              end if
              ijp = (j - 1)*36 + ip !gpople index 16
              nkl = nij

              do k = mink, maxk
                if (kandl) maxl2 = k
                kk1 = k + lock
                if (same) then ! account for non-unique permutations
                  itmp = min(maxl2 - minl + 1, nkl)
                  if (itmp .eq. 0) cycle
                  maxl2 = minl + itmp - 1
                  nkl = nkl - itmp
                end if
                ijkp = (k - 1)*6 + ijp !gpople index 4

                do l = minl, maxl2
                  ijklp = (l - 1)*1 + ijkp !gpople index 1

                  buff(1) = eri_value(ijklp)

                  if (abs(buff(1)) .lt. 5.0d-11) cycle ! goto 300

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
                  !
                  !       account for identical permutations.
                  !
                  if (ii .eq. jj) buff(1) = buff(1)*0.5_dp
                  if (kk .eq. ll) buff(1) = buff(1)*0.5_dp
                  if (ii .eq. kk .and. jj .eq. ll) buff(1) = buff(1)*0.5_dp
                  buff(2) = buff(1)*1.0_dp
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
                end do
              end do
            end do
          end do
        end if
        ! endif
      end do !main loop
      !$omp end target teams distribute parallel do
      kernel_only2 = omp_get_wtime()
      ! write(*,*) "time in 0122 kernel", kernel_only2-kernel_only1
    end do !tiles
    kernel_full2 = omp_get_wtime()
    ! write(*,*) "time in 0122 full", kernel_full2-kernel_full1
    ! do i = 1, size_of_matrix
    !         if(fock(i).gt.1d-5) write(*,*) "fa final", i, fock(i)
    ! end do
    deallocate (n01bra)
    deallocate (xint01bra)
    deallocate (n22ket)
    deallocate (xint22ket)
  end subroutine int0122
end submodule
