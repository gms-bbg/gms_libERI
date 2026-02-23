submodule(rot_axis_kernels) int1102_impl
contains
  module subroutine int1102(pp_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pp_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n11bra(:), n02ket(:)
    real(dp), allocatable :: xint11bra(:), xint02ket(:)
    integer(kind=int64) :: nppbra, nsdket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutppbra, scutsdket
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d11p, d02p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) :: zero_m_0(2), zero_m_1(7), zero_m_2(12), zero_m_3(12), zero_m_4(5), boys0, boys1, boys2, boys3, boys4
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21
    !r-integrals
    real(dp) :: r00(4), r01(18), r02(36), r03(30), r04(15)
    real(dp) :: y03
    !eri_value
    real(dp) :: eri_value(54), trans(6), qx, qz
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
    mini = 1  !1
    maxi = 3  !1
    minj = 1  !0
    maxj = 3  !2
    mink = 1
    maxk = 1
    minl = 1
    maxl = 6
    allocate (n11bra(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11bra(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))
    !start screening
    first_screen1 = omp_get_wtime()
    scutppbra = cutoff_schwarz/maxval(sd_pair%xints)
    nppbra = 0
    do ij = 1, res%n_p_shl*(res%n_p_shl + 1)/2
    if (pp_pair%xints(ij) .ge. scutppbra) then
      nppbra = nppbra + 1
      xint11bra(nppbra) = pp_pair%xints(ij)
      n11bra(nppbra) = ij
    end if
    end do
    scutsdket = cutoff_schwarz/maxval(pp_pair%xints)
    nsdket = 0
    do kl = 1, res%n_s_shl*res%n_d_shl
    if (sd_pair%xints(kl) .ge. scutsdket) then
      nsdket = nsdket + 1
      xint02ket(nsdket) = sd_pair%xints(kl)
      n02ket(nsdket) = kl
    end if
    end do
    first_screen2 = omp_get_wtime()

    nchunksize_int64 = 375000000
    kernel_full1 = omp_get_wtime()
    if ((nppbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = nppbra*nsdket
    ntile = int(nppbra*nsdket/nchunksize_int64) ! nint rounds to the nearest integer, int rounds down
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nppbra*nsdket

      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()
 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, mini, maxi, minj, maxj, mink, maxk, minl, maxl, pp_pair, sd_pair) &
 !$omp shared(boys_grid_zero, exponent_grid) &
 !$omp shared(nquart_start, nquart_end, xint11bra, xint02ket) &
 !$omp shared(nsdket, n11bra, n02ket) &
 !$omp private(shp_thresh, test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d11p, t_alpha, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, zero_m_4, r00, r01, r02, r03, r04, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, boys4, t_inverse, t_new ) &
 !$omp private(qx, qz, eri_value, kandl, same, maxl2, ijk_index, nij, nkl, ijkp, ijkl_index, ijklp) &
 !$omp private(tx21, loci, locj, locl, lock, ip, ijp, itmp, iandj, maxj2)&
 !$omp private(trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsdket + 1
        kl_tmp = mod(iquart - 1, nsdket) + 1

        test = xint11bra(ij_tmp)*xint02ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then

          ij = n11bra(ij_tmp)
          kl = n02ket(kl_tmp)
          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_d_shl + 1
          lsh_tmp = mod(kl - 1, res%n_d_shl) + 1

          ii = res%i_p_shl(ish_tmp)
          jj = res%i_p_shl(jsh_tmp)
          kk = res%i_s_shl(ksh_tmp)
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
          r04 = 0.0_dp
          ket_loop = 0
          do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1
            t_expon_cd = sd_pair%t_expon_ab(sd_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_sp
            t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
            y03 = t_inverse_expon_cd*sd_pair%expon_a(sd_pair%pair_loc(kl) + ket_loop) !y03
            sq = sd_pair%sq(sd_pair%pair_loc(kl) + ket_loop)
            cq = sd_pair%t_alpha(sd_pair%pair_loc(kl) + ket_loop)
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
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              !get bra shell pair info
              bra_loop = bra_loop + 1
              shp_thresh = pp_pair%ismlp(pp_pair%pair_loc(ij) + bra_loop) + sd_pair%ismlp(sd_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              t_expon_ab = pp_pair%t_expon_ab(pp_pair%pair_loc(ij) + bra_loop) !tx12
              t_alpha = pp_pair%t_alpha(pp_pair%pair_loc(ij) + bra_loop) !ty02
              t_beta = pp_pair%t_beta(pp_pair%pair_loc(ij) + bra_loop) !ty01
              tx21 = pp_pair%t_inverse_expon_ab(pp_pair%pair_loc(ij) + bra_loop) !tx21
              !    write(*,*) "t_beta",t_beta
              !    write(*,*) "t_alpha",t_alpha
              !    write(*,*) "tx21",tx21
              d11p = pp_pair%d_coeff(pp_pair%pair_loc(ij) + bra_loop) !d11p
              expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd) !x41
              pqr = t_alpha - aqz
              pqs = pqr*pqr
              !    write(*,*) "pqr",pqr
              !    write(*,*) "pqs",pqs
              rho = t_expon_ab*t_expon_cd*expon_abcd_inverse
              t = (pqs + qps)*rho
              rho = rho + rho
              !    write(*,*) "t",t
              if (t .le. t_max) then
                !evaluate boys function with recursion
                t_new = t*1.5365936651378012d+01 !f_increment(5)
                t_int = nint(t_new)
                fmt = boys_grid_zero((4*451*5) + (t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((4*451*5) + (t_int*5) + 1)
                !use extrapolation for exp(-t)
                t_new = t*m_increment !rxinc
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)
                boys3 = ((t + t)*fmt + expt)*0.142857142857142857d00 !*rmr(4)
                boys2 = ((t + t)*boys3 + expt)*0.2000000000000000d00 !*rmr(3)
                boys1 = ((t + t)*boys2 + expt)*0.3333333333333333d00 !*rmr(2)
                boys0 = ((t + t)*boys1 + expt)!*rmr(1)
                boys0 = boys0*sqrt(expon_abcd_inverse)*d11p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d11p*rho
                boys2 = boys2*sqrt(expon_abcd_inverse)*d11p*rho*rho
                boys3 = boys3*sqrt(expon_abcd_inverse)*d11p*rho*rho*rho
                boys4 = fmt*sqrt(expon_abcd_inverse)*d11p*rho*rho*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d11p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
                boys3 = boys2*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho)
                boys4 = boys3*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho + t_inverse*rho)
              end if
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0*t_beta*t_alpha
              zero_m_0(2) = zero_m_0(2) + boys0*tx21
              zero_m_1(1) = zero_m_1(1) + boys1*t_beta*t_alpha
              zero_m_1(2) = zero_m_1(2) + boys1*tx21
              zero_m_1(3) = zero_m_1(3) + boys1*t_beta*tx21
              zero_m_1(4) = zero_m_1(4) + boys1*tx21*tx21
              zero_m_1(5) = zero_m_1(5) + boys1*t_beta*t_alpha*pqr
              zero_m_1(6) = zero_m_1(6) + boys1*tx21*pqr
              zero_m_1(7) = zero_m_1(7) + boys1*t_beta*tx21*pqr
              zero_m_2(1) = zero_m_2(1) + boys2*t_beta*t_alpha
              zero_m_2(2) = zero_m_2(2) + boys2*tx21
              zero_m_2(3) = zero_m_2(3) + boys2*t_beta*tx21
              zero_m_2(4) = zero_m_2(4) + boys2*tx21*tx21
              zero_m_2(5) = zero_m_2(5) + boys2*t_beta*t_alpha*pqr
              zero_m_2(6) = zero_m_2(6) + boys2*tx21*pqr
              zero_m_2(7) = zero_m_2(7) + boys2*tx21*t_beta*pqr
              zero_m_2(8) = zero_m_2(8) + boys2*tx21*tx21*pqr
              zero_m_2(9) = zero_m_2(9) + boys2*t_beta*t_alpha*pqs
              zero_m_2(10) = zero_m_2(10) + boys2*tx21*pqs
              zero_m_2(11) = zero_m_2(11) + boys2*tx21*t_beta*pqs
              zero_m_2(12) = zero_m_2(12) + boys2*tx21*tx21*pqs

              zero_m_3(1) = zero_m_3(1) + boys3*tx21
              zero_m_3(2) = zero_m_3(2) + boys3*tx21*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*tx21*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*tx21*pqr*pqs
              zero_m_3(5) = zero_m_3(5) + boys3*tx21*t_beta
              zero_m_3(6) = zero_m_3(6) + boys3*tx21*t_beta*pqr
              zero_m_3(7) = zero_m_3(7) + boys3*tx21*t_beta*pqs
              zero_m_3(8) = zero_m_3(8) + boys3*tx21*t_beta*pqs*pqr
              zero_m_3(9) = zero_m_3(9) + boys3*tx21*tx21
              zero_m_3(10) = zero_m_3(10) + boys3*tx21*tx21*pqr
              zero_m_3(11) = zero_m_3(11) + boys3*tx21*tx21*pqs
              zero_m_3(12) = zero_m_3(12) + boys3*tx21*tx21*pqs*pqr

              zero_m_4(1) = zero_m_4(1) + boys4*tx21*tx21
              zero_m_4(2) = zero_m_4(2) + boys4*tx21*tx21*pqr
              zero_m_4(3) = zero_m_4(3) + boys4*tx21*tx21*pqs
              zero_m_4(4) = zero_m_4(4) + boys4*tx21*tx21*pqr*pqs
              zero_m_4(5) = zero_m_4(5) + boys4*tx21*tx21*pqr*pqr*pqs
            end do
            ! write(*,*) "zero_m_0(1) = ", zero_m_0(1)
            ! write(*,*) "zero_m_0(2) = ", zero_m_0(2)
            ! write(*,*) "zero_m_1(1) = ", zero_m_1(1)
            ! write(*,*) "zero_m_1(2) = ", zero_m_1(2)
            ! write(*,*) "zero_m_1(3) = ", zero_m_1(3)
            ! write(*,*) "zero_m_1(4) = ", zero_m_1(4)
            ! write(*,*) "zero_m_1(5) = ", zero_m_1(5)
            ! write(*,*) "zero_m_1(6) = ", zero_m_1(6)
            ! write(*,*) "zero_m_1(7) = ", zero_m_1(7)
            ! write(*,*) "zero_m_2(1 ) = ", zero_m_2(1 )
            ! write(*,*) "zero_m_2(2 ) = ", zero_m_2(2 )
            ! write(*,*) "zero_m_2(3 ) = ", zero_m_2(3 )
            ! write(*,*) "zero_m_2(4 ) = ", zero_m_2(4 )
            ! write(*,*) "zero_m_2(5 ) = ", zero_m_2(5 )
            ! write(*,*) "zero_m_2(6 ) = ", zero_m_2(6 )
            ! write(*,*) "zero_m_2(7 ) = ", zero_m_2(7 )
            ! write(*,*) "zero_m_2(8 ) = ", zero_m_2(8 )
            ! write(*,*) "zero_m_2(9 ) = ", zero_m_2(9 )
            ! write(*,*) "zero_m_2(10) = ", zero_m_2(10)
            ! write(*,*) "zero_m_2(11) = ", zero_m_2(11)
            ! write(*,*) "zero_m_2(12) = ", zero_m_2(12)
            ! write(*,*) "zero_m_3(1 ) = ", zero_m_3(1 )
            ! write(*,*) "zero_m_3(2 ) = ", zero_m_3(2 )
            ! write(*,*) "zero_m_3(3 ) = ", zero_m_3(3 )
            ! write(*,*) "zero_m_3(4 ) = ", zero_m_3(4 )
            ! write(*,*) "zero_m_3(5 ) = ", zero_m_3(5 )
            ! write(*,*) "zero_m_3(6 ) = ", zero_m_3(6 )
            ! write(*,*) "zero_m_3(7 ) = ", zero_m_3(7 )
            ! write(*,*) "zero_m_3(8 ) = ", zero_m_3(8 )
            ! write(*,*) "zero_m_3(9 ) = ", zero_m_3(9 )
            ! write(*,*) "zero_m_3(10) = ", zero_m_3(10)
            ! write(*,*) "zero_m_3(11) = ", zero_m_3(11)
            ! write(*,*) "zero_m_3(12) = ", zero_m_3(12)
            ! write(*,*) "zero_m_4(1 ) = ", zero_m_4(1 )
            ! write(*,*) "zero_m_4(2 ) = ", zero_m_4(2 )
            ! write(*,*) "zero_m_4(3 ) = ", zero_m_4(3 )
            ! write(*,*) "zero_m_4(4 ) = ", zero_m_4(4 )
            ! write(*,*) "zero_m_4(5 ) = ", zero_m_4(5 )
            !r integrals here
            buff(10) = t_inverse_expon_cd*0.5_dp*sq
            r00(1) = r00(1) + zero_m_0(1)*buff(10)
            r00(2) = r00(2) + zero_m_0(2)*buff(10)
            r00(3) = r00(3) + zero_m_0(1)*y03*y03*sq
            r00(4) = r00(4) + zero_m_0(2)*y03*y03*sq
            ! write(*,*) "r00(1) = ", r00(1)
            ! write(*,*) "r00(2) = ", r00(2)
            ! write(*,*) "r00(3) = ", r00(3)
            ! write(*,*) "r00(4) = ", r00(4)
            buff(4) = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            buff(5) = t_inverse_expon_cd*0.5_dp*sq*y03

            r01(1) = r01(1) + zero_m_1(1)*buff(5)*aqx !4
            ! write(*,*) "r01", r01(1)
            r01(2) = r01(2) + zero_m_1(1)*buff(5)*acy
            r01(3) = r01(3) - zero_m_1(5)*buff(5)
            r01(4) = r01(4) + zero_m_1(2)*buff(5)*aqx !5
            r01(5) = r01(5) + zero_m_1(2)*buff(5)*acy
            r01(6) = r01(6) - zero_m_1(6)*buff(5)
            ! r01(7 ) = r01(7 ) - zero_m_1(3) * aqx * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp !10
            ! r01(8 ) = r01(8 ) - zero_m_1(3) * acy * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp
            ! r01(9 ) = r01(9 ) + zero_m_1(7) * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp
            r01(7) = r01(7) - zero_m_1(3)*aqx*buff(10) !10
            r01(8) = r01(8) - zero_m_1(3)*acy*buff(10)
            r01(9) = r01(9) + zero_m_1(7)*buff(10)
            ! write(*,*) "r01 9", r01(9)
            r01(10) = r01(10) - zero_m_1(3)*aqx*y03*y03*sq !11
            r01(11) = r01(11) - zero_m_1(3)*acy*y03*y03*sq
            r01(12) = r01(12) + zero_m_1(7)*y03*y03*sq
            ! r01(13) = r01(13) - (zero_m_1(2)*rab +zero_m_1(3))*aqx * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp !12
            ! r01(14) = r01(14) - (zero_m_1(2)*rab +zero_m_1(3))*acy * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp
            ! r01(15) = r01(15) + (zero_m_1(6)*rab +zero_m_1(7)) * t_inverse_expon_cd *sq*t_inverse_expon_cd *0.25_dp
            r01(13) = r01(13) - (zero_m_1(2)*rab + zero_m_1(3))*aqx*buff(10) !12
            r01(14) = r01(14) - (zero_m_1(2)*rab + zero_m_1(3))*acy*buff(10)
            r01(15) = r01(15) + (zero_m_1(6)*rab + zero_m_1(7))*buff(10)
            ! write(*,*) "r01 15", r01(15)
            r01(16) = r01(16) - (zero_m_1(2)*rab + zero_m_1(3))*aqx*y03*y03*sq !13
            r01(17) = r01(17) - (zero_m_1(2)*rab + zero_m_1(3))*acy*y03*y03*sq
            r01(18) = r01(18) + (zero_m_1(6)*rab + zero_m_1(7))*y03*y03*sq

            ! do i = 1, 18
            !     write(*,*) "rd1", r01(i)
            ! enddo
            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*buff(4) !4
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*buff(4)
            r02(3) = r02(3) + (zero_m_2(9) - zero_m_1(1))*buff(4)
            r02(4) = r02(4) + zero_m_2(1)*buff(4)*aqxy
            r02(5) = r02(5) - zero_m_2(5)*buff(4)*aqx
            r02(6) = r02(6) - zero_m_2(5)*buff(4)*acy

            r02(7) = r02(7) + (zero_m_2(2)*aqx2 - zero_m_1(2))*buff(4) !5
            r02(8) = r02(8) + (zero_m_2(2)*acy2 - zero_m_1(2))*buff(4)
            r02(9) = r02(9) + (zero_m_2(10) - zero_m_1(2))*buff(4)
            r02(10) = r02(10) + zero_m_2(2)*buff(4)*aqxy
            r02(11) = r02(11) - zero_m_2(6)*buff(4)*aqx
            r02(12) = r02(12) - zero_m_2(6)*buff(4)*acy

            r02(13) = r02(13) - (zero_m_2(3)*aqx2 - zero_m_1(3))*buff(5) !8
            r02(14) = r02(14) - (zero_m_2(3)*acy2 - zero_m_1(3))*buff(5)
            r02(15) = r02(15) - (zero_m_2(11) - zero_m_1(3))*buff(5)
            r02(16) = r02(16) - zero_m_2(3)*buff(5)*aqxy
            r02(17) = r02(17) + zero_m_2(7)*buff(5)*aqx
            r02(18) = r02(18) + zero_m_2(7)*buff(5)*acy
            ! write(*,*) "r02(13), r02(13)

            zero_m_2(2) = zero_m_2(2)*rab + zero_m_2(3)
            zero_m_2(6) = zero_m_2(6)*rab + zero_m_2(7)
            zero_m_2(10) = zero_m_2(10)*rab + zero_m_2(11)

            r02(19) = r02(19) - (zero_m_2(2)*aqx2 - (zero_m_1(2)*rab + zero_m_1(3)))*buff(5) !9
            r02(20) = r02(20) - (zero_m_2(2)*acy2 - (zero_m_1(2)*rab + zero_m_1(3)))*buff(5)
            r02(21) = r02(21) - (zero_m_2(10) - (zero_m_1(2)*rab + zero_m_1(3)))*buff(5)
            r02(22) = r02(22) - zero_m_2(2)*buff(5)*aqxy
            r02(23) = r02(23) + zero_m_2(6)*buff(5)*aqx
            r02(24) = r02(24) + zero_m_2(6)*buff(5)*acy

            r02(25) = r02(25) + (zero_m_2(4)*aqx2 - zero_m_1(4))*buff(10) !10
            r02(26) = r02(26) + (zero_m_2(4)*acy2 - zero_m_1(4))*buff(10)
            r02(27) = r02(27) + (zero_m_2(12) - zero_m_1(4))*buff(10)
            r02(28) = r02(28) + zero_m_2(4)*aqxy*buff(10)
            r02(29) = r02(29) - zero_m_2(8)*aqx*buff(10)
            r02(30) = r02(30) - zero_m_2(8)*acy*buff(10)

            r02(31) = r02(31) + (zero_m_2(4)*aqx2 - zero_m_1(4))*y03*y03*sq !11
            r02(32) = r02(32) + (zero_m_2(4)*acy2 - zero_m_1(4))*y03*y03*sq
            r02(33) = r02(33) + (zero_m_2(12) - zero_m_1(4))*y03*y03*sq
            r02(34) = r02(34) + zero_m_2(4)*aqxy*y03*y03*sq
            r02(35) = r02(35) - zero_m_2(8)*aqx*y03*y03*sq
            r02(36) = r02(36) - zero_m_2(8)*acy*y03*y03*sq
            ! do i = 1, 36
            !    write(*,*) "rd2", r02(i)
            ! enddo
            !add r03
            zero_m_3(1) = zero_m_3(1)*rab + zero_m_3(5)
            zero_m_3(2) = zero_m_3(2)*rab + zero_m_3(6)
            zero_m_3(3) = zero_m_3(3)*rab + zero_m_3(7)
            zero_m_3(4) = zero_m_3(4)*rab + zero_m_3(8)
            ! write(*,*) "zero_m 2", zero_m_2(2)
            ! write(*,*) "zero_m 6", zero_m_2(6)

            r03(1) = r03(1) - (zero_m_3(1)*aqx2 - zero_m_2(2)*3.0_dp)*buff(4)*aqx !1
            r03(2) = r03(2) - (zero_m_3(1)*aqx2 - zero_m_2(2))*buff(4)*acy
            r03(3) = r03(3) + (zero_m_3(2)*aqx2 - zero_m_2(6))*buff(4)
            r03(4) = r03(4) - (zero_m_3(1)*acy2 - zero_m_2(2))*buff(4)*aqx
            r03(5) = r03(5) + zero_m_3(2)*buff(4)*aqxy
            r03(6) = r03(6) - (zero_m_3(3) - zero_m_2(2))*buff(4)*aqx
            r03(7) = r03(7) - (zero_m_3(1)*acy2 - zero_m_2(2)*3.0_dp)*buff(4)*acy
            r03(8) = r03(8) + (zero_m_3(2)*acy2 - zero_m_2(6))*buff(4)
            r03(9) = r03(9) - (zero_m_3(3) - zero_m_2(2))*buff(4)*acy
            r03(10) = r03(10) + (zero_m_3(4) - zero_m_2(6)*3.0_dp)*buff(4)

            r03(11) = r03(11) - (zero_m_3(5)*aqx2 - zero_m_2(3)*3.0_dp)*buff(4)*aqx !4
            r03(12) = r03(12) - (zero_m_3(5)*aqx2 - zero_m_2(3))*buff(4)*acy
            r03(13) = r03(13) + (zero_m_3(6)*aqx2 - zero_m_2(7))*buff(4)
            r03(14) = r03(14) - (zero_m_3(5)*acy2 - zero_m_2(3))*buff(4)*aqx
            r03(15) = r03(15) + zero_m_3(6)*buff(4)*aqxy
            r03(16) = r03(16) - (zero_m_3(7) - zero_m_2(3))*buff(4)*aqx
            r03(17) = r03(17) - (zero_m_3(5)*acy2 - zero_m_2(3)*3.0_dp)*buff(4)*acy
            r03(18) = r03(18) + (zero_m_3(6)*acy2 - zero_m_2(7))*buff(4)
            r03(19) = r03(19) - (zero_m_3(7) - zero_m_2(3))*buff(4)*acy !this guy
            r03(20) = r03(20) + (zero_m_3(8) - zero_m_2(7)*3.0_dp)*buff(4)

            r03(21) = r03(21) + (zero_m_3(9)*aqx2 - zero_m_2(4)*3.0_dp)*buff(5)*aqx !5
            r03(22) = r03(22) + (zero_m_3(9)*aqx2 - zero_m_2(4))*buff(5)*acy
            r03(23) = r03(23) - (zero_m_3(10)*aqx2 - zero_m_2(8))*buff(5)
            r03(24) = r03(24) + (zero_m_3(9)*acy2 - zero_m_2(4))*buff(5)*aqx
            r03(25) = r03(25) - zero_m_3(10)*buff(5)*aqxy
            r03(26) = r03(26) + (zero_m_3(11) - zero_m_2(4))*buff(5)*aqx
            r03(27) = r03(27) + (zero_m_3(9)*acy2 - zero_m_2(4)*3.0_dp)*buff(5)*acy
            r03(28) = r03(28) - (zero_m_3(10)*acy2 - zero_m_2(8))*buff(5)
            r03(29) = r03(29) + (zero_m_3(11) - zero_m_2(4))*buff(5)*acy
            r03(30) = r03(30) - (zero_m_3(12) - zero_m_2(8)*3.0_dp)*buff(5)
            ! do i = 1, 30
            !    write(*,*) "rd3", r03(i)
            ! enddo
            ! write(*,*) "rd3(1,1), r03(1)
            ! write(*,*) "rd3(10,1), r03(10)
            ! write(*,*) "rd3(1,4), r03(11)
            ! write(*,*) "rd3(10,4), r03(20)
            ! write(*,*) "rd3(1,5), r03(21)
            ! write(*,*) "zero_m_2(8), zero_m_2(8)
            ! write(*,*) "rd3(10,5), r03(30)

            r04(1) = r04(1) + (zero_m_4(1)*aqx2*aqx2 - zero_m_3(9)*6.0_dp*aqx2 + zero_m_2(4)*3.0_dp)*buff(4)
            r04(2) = r04(2) + (zero_m_4(1)*aqx2 - zero_m_3(9)*3.0_dp)*buff(4)*aqxy
            r04(3) = r04(3) - (zero_m_4(2)*aqx2 - zero_m_3(10)*3.0_dp)*buff(4)*aqx
            r04(4) = r04(4) + (zero_m_4(1)*aqx2*acy2 - zero_m_3(9)*(aqx2 + acy2) + zero_m_2(4))*buff(4)
            r04(5) = r04(5) - (zero_m_4(2)*aqx2 - zero_m_3(10))*buff(4)*acy
            r04(6) = r04(6) + (zero_m_4(3)*aqx2 - zero_m_3(9)*aqx2 - zero_m_3(11) + zero_m_2(4))*buff(4)
            r04(7) = r04(7) + (zero_m_4(1)*acy2 - zero_m_3(9)*3.0_dp)*buff(4)*aqxy
            r04(8) = r04(8) - (zero_m_4(2)*acy2 - zero_m_3(10))*buff(4)*aqx
            r04(9) = r04(9) + (zero_m_4(3) - zero_m_3(9))*buff(4)*aqxy
            r04(10) = r04(10) - (zero_m_4(4) - zero_m_3(10)*3.0_dp)*buff(4)*aqx
            r04(11) = r04(11) + (zero_m_4(1)*acy2*acy2 - zero_m_3(9)*6.0_dp*acy2 + zero_m_2(4)*3.0_dp)*buff(4)
            r04(12) = r04(12) - (zero_m_4(2)*acy2 - zero_m_3(10)*3.0_dp)*buff(4)*acy
            r04(13) = r04(13) + (zero_m_4(3)*acy2 - zero_m_3(9)*acy2 - zero_m_3(11) + zero_m_2(4))*buff(4)
            r04(14) = r04(14) - (zero_m_4(4) - zero_m_3(10)*3.0_dp)*buff(4)*acy
            r04(15) = r04(15) + (zero_m_4(5) - zero_m_3(11)*6.0_dp + zero_m_2(4)*3.0_dp)*buff(4)
            ! do i = 1, 15
            !    write(*,*) "rd4", r04(i)
            ! enddo
          end do !end ket loop
          qx = rcd*sing
          qz = rcd*cosg

          eri_value(1) = r00(2) + r02(7) + r02(25) + qx*(2*r01(4) + qx*(r00(4) + r02(31)) + 2*r03(21)) + r04(1)
          eri_value(2) = r00(2) + r02(8) + r02(25) + r04(4)
          eri_value(3) = r00(2) + r02(9) + r02(25) + qz*(2*r01(6) + qz*(r00(4) + r02(31)) + 2*r03(23)) + r04(6)
          eri_value(4) = r02(10) + qx*(r01(5) + r03(22)) + r04(2)
          eri_value(5) = r02(11) + qz*(r01(4) + qx*(r00(4) + r02(31)) + r03(21)) + qx*(r01(6) + r03(23)) + r04(3)
          eri_value(6) = r02(12) + qz*(r01(5) + r03(22)) + r04(5)
          eri_value(7) = r02(28) + qx*(qx*r02(34) + 2*r03(22)) + r04(2)
          eri_value(8) = r02(28) + r04(7)
          eri_value(9) = r02(28) + qz*(qz*r02(34) + 2*r03(25)) + r04(9)
          eri_value(10) = qx*r03(24) + r04(4)
          eri_value(11) = qz*(qx*r02(34) + r03(22)) + qx*r03(25) + r04(5)
          eri_value(12) = qz*r03(24) + r04(8)
          eri_value(13) = -r01(7) + r02(29) - r03(11) + qx*(-2*r02(13) + qx*(-r01(10) + r02(35)) + 2*r03(23)) + r04(3)
          eri_value(14) = -r01(7) + r02(29) - r03(14) + r04(8)
          eri_value(15) = -r01(7) + r02(29) - r03(16) + qz*(-2*r02(17) + qz*(-r01(10) + r02(35)) + 2*r03(26)) + r04(10)
          eri_value(16) = -r03(12) + qx*(-r02(16) + r03(25)) + r04(5)
          eri_value(17) = -r03(13) + qz*(-r02(13) + qx*(-r01(10) + r02(35)) + r03(23)) + qx*(-r02(17) + r03(26)) + r04(6)
          eri_value(18) = -r03(15) + qz*(-r02(16) + r03(25)) + r04(9)
          eri_value(19) = r02(28) + qx*(qx*r02(34) + 2*r03(22)) + r04(2)
          eri_value(20) = r02(28) + r04(7)
          eri_value(21) = r02(28) + qz*(qz*r02(34) + 2*r03(25)) + r04(9)
          eri_value(22) = qx*r03(24) + r04(4)
          eri_value(23) = qz*(qx*r02(34) + r03(22)) + qx*r03(25) + r04(5)
          eri_value(24) = qz*r03(24) + r04(8)
          eri_value(25) = r00(2) + r02(7) + r02(26) + qx*(2*r01(4) + qx*(r00(4) + r02(32)) + 2*r03(24)) + r04(4)
          eri_value(26) = r00(2) + r02(8) + r02(26) + r04(11)
          eri_value(27) = r00(2) + r02(9) + r02(26) + qz*(2*r01(6) + qz*(r00(4) + r02(32)) + 2*r03(28)) + r04(13)
          eri_value(28) = r02(10) + qx*(r01(5) + r03(27)) + r04(7)
          eri_value(29) = r02(11) + qz*(r01(4) + qx*(r00(4) + r02(32)) + r03(24)) + qx*(r01(6) + r03(28)) + r04(8)
          eri_value(30) = r02(12) + qz*(r01(5) + r03(27)) + r04(12)
          eri_value(31) = -r01(8) + r02(30) - r03(12) + qx*(-2*r02(16) + qx*(-r01(11) + r02(36)) + 2*r03(25)) + r04(5)
          eri_value(32) = -r01(8) + r02(30) - r03(17) + r04(12)
          eri_value(33) = -r01(8) + r02(30) - r03(19) + qz*(-2*r02(18) + qz*(-r01(11) + r02(36)) + 2*r03(29)) + r04(14)
          eri_value(34) = -r03(14) + qx*(-r02(14) + r03(28)) + r04(8)
          eri_value(35) = -r03(15) + qz*(-r02(16) + qx*(-r01(11) + r02(36)) + r03(25)) + qx*(-r02(18) + r03(29)) + r04(9)
          eri_value(36) = -r03(18) + qz*(-r02(14) + r03(28)) + r04(13)
          eri_value(37) = -r01(13) + r02(29) - r03(1) + qx*(-2*r02(19) + qx*(-r01(16) + r02(35)) + 2*r03(23)) + r04(3)
          eri_value(38) = -r01(13) + r02(29) - r03(4) + r04(8)
          eri_value(39) = -r01(13) + r02(29) - r03(6) + qz*(-2*r02(23) + qz*(-r01(16) + r02(35)) + 2*r03(26)) + r04(10)
          eri_value(40) = -r03(2) + qx*(-r02(22) + r03(25)) + r04(5)
          eri_value(41) = -r03(3) + qz*(-r02(19) + qx*(-r01(16) + r02(35)) + r03(23)) + qx*(-r02(23) + r03(26)) + r04(6)
          eri_value(42) = -r03(5) + qz*(-r02(22) + r03(25)) + r04(9)
          eri_value(43) = -r01(14) + r02(30) - r03(2) + qx*(-2*r02(22) + qx*(-r01(17) + r02(36)) + 2*r03(25)) + r04(5)
          eri_value(44) = -r01(14) + r02(30) - r03(7) + r04(12)
          eri_value(45) = -r01(14) + r02(30) - r03(9) + qz*(-2*r02(24) + qz*(-r01(17) + r02(36)) + 2*r03(29)) + r04(14)
          eri_value(46) = -r03(4) + qx*(-r02(20) + r03(28)) + r04(8)
          eri_value(47) = -r03(5) + qz*(-r02(22) + qx*(-r01(17) + r02(36)) + r03(25)) + qx*(-r02(24) + r03(29)) + r04(9)
          eri_value(48) = -r03(8) + qz*(-r02(20) + r03(28)) + r04(13)
  eri_value(49) = r00(1) + r00(2) - r01(9) - r01(15) + r02(1) + r02(7) + r02(27) - r03(3) - r03(13) + qx*(2*r01(1) + 2*r01(4) - 2*r02(17) - 2*r02(23) + qx*(r00(3) + r00(4) - r01(12) - r01(18) + r02(33)) + 2*r03(26)) + r04(6)
          eri_value(50) = r00(1) + r00(2) - r01(9) - r01(15) + r02(2) + r02(8) + r02(27) - r03(8) - r03(18) + r04(13)
  eri_value(51) = r00(1) + r00(2) - r01(9) - r01(15) + r02(3) + r02(9) + r02(27) - r03(10) - r03(20) + qz*(2*r01(3) + 2*r01(6) - 2*r02(15) - 2*r02(21) + qz*(r00(3) + r00(4) - r01(12) - r01(18) + r02(33)) + 2*r03(30)) + r04(15)
          eri_value(52) = r02(4) + r02(10) - r03(5) - r03(15) + qx*(r01(2) + r01(5) - r02(18) - r02(24) + r03(29)) + r04(9)
  eri_value(53) = r02(5) + r02(11) - r03(6) - r03(16) + qz*(r01(1) + r01(4) - r02(17) - r02(23) + qx*(r00(3) + r00(4) - r01(12) - r01(18) + r02(33)) + r03(26)) + qx*(r01(3) + r01(6) - r02(15) - r02(21) + r03(30)) + r04(10)
          eri_value(54) = r02(6) + r02(12) - r03(9) - r03(19) + qz*(r01(2) + r01(5) - r02(18) - r02(24) + r03(29)) + r04(14)
          ! write(*,*) eri 16", eri_value(16)
          ! ! eri_value(17) = -r03(3) + qz*(-r02(13) + qx*(-r01(10) + r02(35)) + r03(13)) + qx*(-r02(17) + r03(16)) + r04(6)
          ! write(*,*) eri 17", eri_value(17)

          ! eri_value(50) = r00(1) + r00(2) - r01(9) - r01(15) + r02(2) + r02(8) + r02(27) - r03(8) - r03(18) + r04(13)
          ! write(*,*) "r4", r04(13)
          ! write(*,*) "r3", r03(8)
          ! write(*,*) "r3", r03(18)
          ! write(*,*) "r2", r02(2)
          ! write(*,*) "r2", r02(8)
          ! write(*,*) "c3", r02(2) + r02(8) - r03(8) - r03(18) + r04(13)

          ! write(*,*) "r2", r02(27)
          ! write(*,*) "r1", r01(9) !problem
          ! write(*,*) "r1", r01(15) !problem
          ! write(*,*) "r0", r00(1)
          ! write(*,*) "r0", r00(2)
          ! write(*,*) "c1", r00(1) + r00(2) - r01(9) - r01(15) + r02(27)
          ! do i =1,54
          !     write(*,*) "eri", eri_value(i)
          ! enddo !ok

          buff(10) = (buff(8)*buff(3) - buff(9)*buff(2))
          buff(11) = (buff(9)*buff(1) - buff(7)*buff(3))
          buff(12) = (buff(7)*buff(2) - buff(8)*buff(1))

          trans(1) = eri_value(1)
          trans(2) = eri_value(19)
          trans(3) = eri_value(37)
          eri_value(1) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(19) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(37) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(2)
          trans(2) = eri_value(20)
          trans(3) = eri_value(38)
          eri_value(2) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(20) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(38) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(3)
          trans(2) = eri_value(21)
          trans(3) = eri_value(39)
          eri_value(3) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(21) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(39) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(4)
          trans(2) = eri_value(22)
          trans(3) = eri_value(40)
          eri_value(4) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(22) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(40) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(5)
          trans(2) = eri_value(23)
          trans(3) = eri_value(41)
          eri_value(5) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(23) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(41) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(6)
          trans(2) = eri_value(24)
          trans(3) = eri_value(42)
          eri_value(6) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(24) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(42) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(7)
          trans(2) = eri_value(25)
          trans(3) = eri_value(43)
          eri_value(7) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(25) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(43) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(8)
          trans(2) = eri_value(26)
          trans(3) = eri_value(44)
          eri_value(8) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(26) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(44) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(9)
          trans(2) = eri_value(27)
          trans(3) = eri_value(45)
          eri_value(9) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(27) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(45) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(10)
          trans(2) = eri_value(28)
          trans(3) = eri_value(46)
          eri_value(10) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(28) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(46) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(11)
          trans(2) = eri_value(29)
          trans(3) = eri_value(47)
          eri_value(11) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(29) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(47) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(12)
          trans(2) = eri_value(30)
          trans(3) = eri_value(48)
          eri_value(12) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(30) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(48) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(13)
          trans(2) = eri_value(31)
          trans(3) = eri_value(49)
          eri_value(13) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(31) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(49) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(14)
          trans(2) = eri_value(32)
          trans(3) = eri_value(50)
          eri_value(14) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(32) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(50) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(15)
          trans(2) = eri_value(33)
          trans(3) = eri_value(51)
          eri_value(15) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(33) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(51) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(16)
          trans(2) = eri_value(34)
          trans(3) = eri_value(52)
          eri_value(16) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(34) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(52) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(17)
          trans(2) = eri_value(35)
          trans(3) = eri_value(53)
          eri_value(17) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(35) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(53) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(18)
          trans(2) = eri_value(36)
          trans(3) = eri_value(54)
          eri_value(18) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(36) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(54) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)

          trans(1) = eri_value(1)
          trans(2) = eri_value(7)
          trans(3) = eri_value(13)
          eri_value(1) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(7) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(13) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(2)
          trans(2) = eri_value(8)
          trans(3) = eri_value(14)
          eri_value(2) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(8) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(14) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(3)
          trans(2) = eri_value(9)
          trans(3) = eri_value(15)
          eri_value(3) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(9) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(15) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(4)
          trans(2) = eri_value(10)
          trans(3) = eri_value(16)
          eri_value(4) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(10) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(16) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(5)
          trans(2) = eri_value(11)
          trans(3) = eri_value(17)
          eri_value(5) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(11) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(17) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(6)
          trans(2) = eri_value(12)
          trans(3) = eri_value(18)
          eri_value(6) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(12) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(18) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(19)
          trans(2) = eri_value(25)
          trans(3) = eri_value(31)
          eri_value(19) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(25) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(31) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(20)
          trans(2) = eri_value(26)
          trans(3) = eri_value(32)
          eri_value(20) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(26) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(32) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(21)
          trans(2) = eri_value(27)
          trans(3) = eri_value(33)
          eri_value(21) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(27) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(33) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(22)
          trans(2) = eri_value(28)
          trans(3) = eri_value(34)
          eri_value(22) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(28) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(34) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(23)
          trans(2) = eri_value(29)
          trans(3) = eri_value(35)
          eri_value(23) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(29) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(35) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(24)
          trans(2) = eri_value(30)
          trans(3) = eri_value(36)
          eri_value(24) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(30) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(36) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(37)
          trans(2) = eri_value(43)
          trans(3) = eri_value(49)
          eri_value(37) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(43) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(49) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(38)
          trans(2) = eri_value(44)
          trans(3) = eri_value(50)
          eri_value(38) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(44) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(50) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(39)
          trans(2) = eri_value(45)
          trans(3) = eri_value(51)
          eri_value(39) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(45) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(51) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(40)
          trans(2) = eri_value(46)
          trans(3) = eri_value(52)
          eri_value(40) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(46) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(52) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(41)
          trans(2) = eri_value(47)
          trans(3) = eri_value(53)
          eri_value(41) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(47) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(53) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)
          trans(1) = eri_value(42)
          trans(2) = eri_value(48)
          trans(3) = eri_value(54)
          eri_value(42) = buff(10)*trans(1) + buff(7)*trans(2) + buff(1)*trans(3)
          eri_value(48) = buff(11)*trans(1) + buff(8)*trans(2) + buff(2)*trans(3)
          eri_value(54) = buff(12)*trans(1) + buff(9)*trans(2) + buff(3)*trans(3)

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

          ! do i =1,54
          !     write(*,*) "eri", eri_value(i)
          ! enddo !ok

          mini = 1  !1
          maxi = 3  !1
          minj = 1  !0
          maxj = 3  !2
          mink = 1
          maxk = 1
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
            ip = (i - 1)*18 + 1

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
              ijp = (j - 1)*6 + ip !gpople index 16
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
                ijkp = (k - 1)*0 + ijp

                do l = minl, maxl2
                  ijklp = (l - 1)*1 + ijkp

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
      ! write(*,*) "time in 1102 kernel", kernel_only2-kernel_only1
    end do !tiles
    kernel_full2 = omp_get_wtime()
    ! write(*,*) "time in 1102 full", kernel_full2-kernel_full1
    ! do i = 1, size_of_matrix
    !         if(fock(i).ne.0.0_dp) write(*,*) "fa final", i, fock(i)
    !         ! if(fock(i).eq.12.73268741050258) write(*,*) "fa final", i, fock(i); i= 318003
    ! end do
    deallocate (n11bra)
    deallocate (xint11bra)
    deallocate (n02ket)
    deallocate (xint02ket)
  end subroutine int1102
end submodule
