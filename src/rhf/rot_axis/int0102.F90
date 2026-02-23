submodule(rot_axis_kernels) int0102_impl
contains
  module subroutine int0102(sp_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sp_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n01bra(:), n02ket(:)
    real(dp), allocatable :: xint01bra(:), xint02ket(:)
    integer(kind=int64) :: nspbra, nsdket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutspbra, scutsdket
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d01p, d02p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) ::  zero_m_0(1), zero_m_1(4), zero_m_2(6), zero_m_3(4), boys0, boys1, boys2, boys3
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21
    !r-integrals
    real(dp) :: y03
    real(dp) :: r00(2), r01(9), r02(12), r03(10)
    !eri_value
    real(dp) :: eri_value(18), trans(6), qx, qz
    !digestion
    logical :: kandl, same
  integer(kind=int64) :: ii1,kk1,nij,maxl2,jj1,i2,j2,ijp,nkl,itmp,ijklp,ll2,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ll1,ijkp,ij_index,ijk_index,ijkl_index
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: nchunksize_int64
    integer(kind=int64) :: istart, iend, itile, ntile
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp
    real(dp) :: kernel_full1, kernel_full2, kernel_only1, kernel_only2, first_screen1, first_screen2
    integer :: shp_thresh
    basd = 1
    minl = 1
    maxl = 1
    basc = 0
    mink = 1
    maxk = 1
    basb = 6
    minj = 1
    maxj = 6
    basa = 0
    mini = 1
    maxi = 1
    allocate (n01bra(res%n_s_shl*res%n_p_shl))
    allocate (xint01bra(res%n_s_shl*res%n_p_shl))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))
    !start screening
    first_screen1 = omp_get_wtime()
    scutspbra = cutoff_schwarz/maxval(sd_pair%xints)
    nspbra = 0
    do ij = 1, res%n_s_shl*res%n_p_shl
    if (sp_pair%xints(ij) .ge. scutspbra) then
      nspbra = nspbra + 1
      xint01bra(nspbra) = sp_pair%xints(ij)
      n01bra(nspbra) = ij
    end if
    end do
    scutsdket = cutoff_schwarz/maxval(sp_pair%xints)
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
    if ((nspbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = nspbra*nsdket
    ntile = int(nspbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nspbra*nsdket
      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, sd_pair, sp_pair) &
 !$omp shared(nquart_start, nquart_end, xint01bra, xint02ket) &
 !$omp shared(nsdket, n01bra, n02ket) &
 !$omp private(shp_thresh, itmp, tx21, test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d01p, t_alpha, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, r00, r01, r02, r03, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, t_inverse, t_new ) &
 !$omp private(qx, qz, eri_value, kandl, same, maxl2, ijk_index, nij, nkl, ijkp, ijkl_index, ijklp) &
 !$omp private(trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsdket + 1
        kl_tmp = mod(iquart - 1, nsdket) + 1

        test = xint01bra(ij_tmp)*xint02ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then
          ij = n01bra(ij_tmp)
          kl = n02ket(kl_tmp)

          ish_tmp = (ij - 1)/res%n_p_shl + 1
          jsh_tmp = mod(ij - 1, res%n_p_shl) + 1
          ksh_tmp = (kl - 1)/res%n_d_shl + 1
          lsh_tmp = mod(kl - 1, res%n_d_shl) + 1

          ii = res%i_s_shl(ish_tmp)
          jj = res%i_p_shl(jsh_tmp)
          kk = res%i_s_shl(ksh_tmp)
          ll = res%i_d_shl(lsh_tmp)

          ish = ii
          jsh = jj
          ksh = kk
          lsh = ll
          ! write(*,*) "ii", ii
          ! write(*,*) "jj", jj
          ! write(*,*) "kk", kk
          ! write(*,*) "ll", ll

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
          !-----------acy2 stuff check -------
          acy2 = acy*acy
          r00 = 0.0_dp
          r01 = 0.0_dp
          r02 = 0.0_dp
          r03 = 0.0_dp
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
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              bra_loop = bra_loop + 1
              !get bra shell pair info
              shp_thresh = sp_pair%ismlp(sp_pair%pair_loc(ij) + bra_loop) + sd_pair%ismlp(sd_pair%pair_loc(kl) + ket_loop)
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
                t_new = t*1.5711584096264371d+01 !f_increment(4)
                t_int = nint(t_new)
                fmt = boys_grid_zero((3*451*5) + (t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((3*451*5) + (t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((3*451*5) + (t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((3*451*5) + (t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((3*451*5) + (t_int*5) + 1)
                !use extrapolation for exp(-T)
                t_new = t*m_increment !rxinc
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)
                boys2 = ((t + t)*fmt + expt)*0.2000000000000000d00 !*rmr(3)
                boys1 = ((t + t)*boys2 + expt)*0.33333333333333333d00 !*rmr(2)
                boys0 = ((t + t)*boys1 + expt)!*rmr(1)
                boys0 = boys0*sqrt(expon_abcd_inverse)*d01p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d01p*rho
                boys2 = boys2*sqrt(expon_abcd_inverse)*d01p*rho*rho
                boys3 = fmt*sqrt(expon_abcd_inverse)*d01p*rho*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d01p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
                boys3 = boys2*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho)
              end if
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
              zero_m_3(1) = zero_m_3(1) + boys3*tx21
              zero_m_3(2) = zero_m_3(2) + boys3*tx21*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*tx21*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*tx21*pqr*pqs
            end do
            ! write(*,*) "zero_m_0(1) = ", zero_m_0(1)
            ! write(*,*) "zero_m_1(1) = ", zero_m_1(1)
            ! write(*,*) "zero_m_1(2) = ", zero_m_1(2)
            ! write(*,*) "zero_m_1(3) = ", zero_m_1(3)
            ! write(*,*) "zero_m_1(4) = ", zero_m_1(4)
            ! write(*,*) "zero_m_2(1) = ", zero_m_2(1)
            ! write(*,*) "zero_m_2(2) = ", zero_m_2(2)
            ! write(*,*) "zero_m_2(3) = ", zero_m_2(3)
            ! write(*,*) "zero_m_2(4) = ", zero_m_2(4)
            ! write(*,*) "zero_m_2(5) = ", zero_m_2(5)
            ! write(*,*) "zero_m_2(6) = ", zero_m_2(6)
            ! write(*,*) "zero_m_3(1) = ", zero_m_3(1)
            ! write(*,*) "zero_m_3(2) = ", zero_m_3(2)
            ! write(*,*) "zero_m_3(3) = ", zero_m_3(3)
            ! write(*,*) "zero_m_3(4) = ", zero_m_3(4)
            !r integrals here
            r00(1) = r00(1) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq
            r00(2) = r00(2) + zero_m_0(1)*y03*y03*sq
            ! write(*,*) "r00(1) = ", r00(1)
            ! write(*,*) "r00(2) = ", r00(2)

            r01(1) = r01(1) + zero_m_1(1)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03
            r01(2) = r01(2) + zero_m_1(1)*acy*t_inverse_expon_cd*0.5_dp*sq*y03
            r01(3) = r01(3) - zero_m_1(3)*t_inverse_expon_cd*0.5_dp*sq*y03
            r01(4) = r01(4) - zero_m_1(2)*aqx*t_inverse_expon_cd*0.5_dp*sq
            r01(5) = r01(5) - zero_m_1(2)*acy*t_inverse_expon_cd*0.5_dp*sq
            r01(6) = r01(6) + zero_m_1(4)*t_inverse_expon_cd*0.5_dp*sq
            r01(7) = r01(7) - zero_m_1(2)*aqx*y03*y03*sq
            r01(8) = r01(8) - zero_m_1(2)*acy*y03*y03*sq
            r01(9) = r01(9) + zero_m_1(4)*y03*y03*sq

            ! write(*,*) "r01(1) = ", r01(1)
            ! write(*,*) "r01(2) = ", r01(2)
            ! write(*,*) "r01(3) = ", r01(3)
            ! write(*,*) "r01(4) = ", r01(4)
            ! write(*,*) "r01(5) = ", r01(5)
            ! write(*,*) "r01(6) = ", r01(6)
            ! write(*,*) "r01(7) = ", r01(7)
            ! write(*,*) "r01(8) = ", r01(8)
            ! write(*,*) "r01(9) = ", r01(9)

            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(3) = r02(3) + (zero_m_2(5) - zero_m_1(1))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(4) = r02(4) + zero_m_2(1)*aqxy*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(5) = r02(5) - zero_m_2(3)*aqx*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(6) = r02(6) - zero_m_2(3)*acy*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(7) = r02(7) - (zero_m_2(2)*aqx2 - zero_m_1(2))*t_inverse_expon_cd*0.5_dp*sq*y03
            r02(8) = r02(8) - (zero_m_2(2)*acy2 - zero_m_1(2))*t_inverse_expon_cd*0.5_dp*sq*y03
            r02(9) = r02(9) - (zero_m_2(6) - zero_m_1(2))*t_inverse_expon_cd*0.5_dp*sq*y03
            r02(10) = r02(10) - zero_m_2(2)*aqxy*t_inverse_expon_cd*0.5_dp*sq*y03
            r02(11) = r02(11) + zero_m_2(4)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03
            r02(12) = r02(12) + zero_m_2(4)*acy*t_inverse_expon_cd*0.5_dp*sq*y03

            ! write(*,*) "sq = ", sq
            ! write(*,*) "x43 = ", t_inverse_expon_cd
            ! write(*,*) "y03 = ", y03

            ! write(*,*) "zero_m_1(2) = ", zero_m_1(2)
            ! write(*,*) "fqd2(1,3) = ", zero_m_2(2)
            ! write(*,*) "aqx2 = ", aqx2
            ! write(*,*) "bfr work = ", (zero_m_2(2)*aqx2 - zero_m_1(2))
            ! write(*,*) "work(3) = ", t_inverse_expon_cd *0.5_dp*sq*y03
            ! write(*,*) "r02(1) = ", r02(1)
            ! write(*,*) "r02(2) = ", r02(2)
            ! write(*,*) "r02(3) = ", r02(3)
            ! write(*,*) "r02(4) = ", r02(4)
            ! write(*,*) "r02(5) = ", r02(5)
            ! write(*,*) "r02(6) = ", r02(6)
            ! write(*,*) "r02(7) = ", r02(7)
            ! write(*,*) "r02(8) = ", r02(8)
            ! write(*,*) "r02(9) = ", r02(9)
            ! write(*,*) "r02(10) = ", r02(10)
            ! write(*,*) "r02(11) = ", r02(11)
            ! write(*,*) "r02(12) = ", r02(12)
            r03(1) = r03(1) - (zero_m_3(1)*aqx2 - zero_m_2(2)*3.0_dp)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*aqx
            r03(2) = r03(2) - (zero_m_3(1)*aqx2 - zero_m_2(2))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*acy
            r03(3) = r03(3) + (zero_m_3(2)*aqx2 - zero_m_2(4))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r03(4) = r03(4) - (zero_m_3(1)*acy2 - zero_m_2(2))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*aqx
            r03(5) = r03(5) + zero_m_3(2)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*aqxy
            r03(6) = r03(6) - (zero_m_3(3) - zero_m_2(2))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*aqx
            r03(7) = r03(7) - (zero_m_3(1)*acy2 - zero_m_2(2)*3.0_dp)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*acy
            r03(8) = r03(8) + (zero_m_3(2)*acy2 - zero_m_2(4))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r03(9) = r03(9) - (zero_m_3(3) - zero_m_2(2))*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp*acy
            r03(10) = r03(10) + (zero_m_3(4) - zero_m_2(4)*3.0_dp)*t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp

            ! write(*,*) "r03(1) = ", r03(1)
            ! write(*,*) "r03(2) = ", r03(2)
            ! write(*,*) "r03(3) = ", r03(3)
            ! write(*,*) "r03(4) = ", r03(4)
            ! write(*,*) "r03(5) = ", r03(5)
            ! write(*,*) "r03(6) = ", r03(6)
            ! write(*,*) "r03(7) = ", r03(7)
            ! write(*,*) "r03(8) = ", r03(8)
            ! write(*,*) "r03(9) = ", r03(9)
            ! write(*,*) "r03(10) = ", r03(10)
          end do !end ket loop
          qx = rcd*sing
          qz = rcd*cosg
          eri_value(1) = -r01(4) + qx*(-(qx*r01(7)) - 2*r02(7)) - r03(1)
          eri_value(2) = -r01(4) - r03(4)
          eri_value(3) = -r01(4) + qz*(-(qz*r01(7)) - 2*r02(11)) - r03(6)
          eri_value(4) = -(qx*r02(10)) - r03(2)
          eri_value(5) = qz*(-(qx*r01(7)) - r02(7)) - qx*r02(11) - r03(3)
          eri_value(6) = -(qz*r02(10)) - r03(5)
          eri_value(7) = -r01(5) + qx*(-(qx*r01(8)) - 2*r02(10)) - r03(2)
          eri_value(8) = -r01(5) - r03(7)
          eri_value(9) = -r01(5) + qz*(-(qz*r01(8)) - 2*r02(12)) - r03(9)
          eri_value(10) = -(qx*r02(8)) - r03(4)
          eri_value(11) = qz*(-(qx*r01(8)) - r02(10)) - qx*r02(12) - r03(5)
          eri_value(12) = -(qz*r02(8)) - r03(8)
          eri_value(13) = r00(1) - r01(6) + r02(1) + qx*(2*r01(1) + qx*(r00(2) - r01(9)) - 2*r02(11)) - r03(3)
          eri_value(14) = r00(1) - r01(6) + r02(2) - r03(8)
          eri_value(15) = r00(1) - r01(6) + r02(3) + qz*(2*r01(3) + qz*(r00(2) - r01(9)) - 2*r02(9)) - r03(10)
          eri_value(16) = r02(4) + qx*(r01(2) - r02(12)) - r03(5)
          eri_value(17) = qx*(r01(3) - r02(9)) + r02(5) + qz*(r01(1) + qx*(r00(2) - r01(9)) - r02(11)) - r03(6)
          eri_value(18) = r02(6) + qz*(r01(2) - r02(12)) - r03(9)
          ! do i = 1, 18
          !     write(*,*) "eri_value= ", eri_value(i)
          ! enddo
          buff(10) = buff(8)*buff(3) - buff(9)*buff(2)
          buff(11) = buff(9)*buff(1) - buff(7)*buff(3)
          buff(12) = buff(7)*buff(2) - buff(8)*buff(1)

          trans(1) = eri_value(1)
          trans(2) = eri_value(7)
          trans(3) = eri_value(13)
          eri_value(1) = buff(10)*trans(1) + buff(1)*trans(3) + buff(7)*trans(2)
          eri_value(7) = buff(11)*trans(1) + buff(2)*trans(3) + buff(8)*trans(2)
          eri_value(13) = buff(12)*trans(1) + buff(3)*trans(3) + buff(9)*trans(2)
          trans(1) = eri_value(2)
          trans(2) = eri_value(8)
          trans(3) = eri_value(14)
          eri_value(2) = buff(1)*trans(3) + buff(10)*trans(1) + buff(7)*trans(2)
          eri_value(8) = buff(2)*trans(3) + buff(11)*trans(1) + buff(8)*trans(2)
          eri_value(14) = buff(3)*trans(3) + buff(12)*trans(1) + buff(9)*trans(2)

          trans(1) = eri_value(3)
          trans(2) = eri_value(9)
          trans(3) = eri_value(15)
          eri_value(3) = buff(1)*trans(3) + buff(10)*trans(1) + buff(7)*trans(2)
          eri_value(9) = buff(2)*trans(3) + buff(11)*trans(1) + buff(8)*trans(2)
          eri_value(15) = buff(3)*trans(3) + buff(12)*trans(1) + buff(9)*trans(2)
          trans(1) = eri_value(4)
          trans(2) = eri_value(10)
          trans(3) = eri_value(16)
          eri_value(4) = buff(7)*trans(2) + buff(1)*trans(3) + buff(10)*trans(1)
          eri_value(10) = buff(8)*trans(2) + buff(2)*trans(3) + buff(11)*trans(1)
          eri_value(16) = buff(9)*trans(2) + buff(3)*trans(3) + buff(12)*trans(1)
          trans(1) = eri_value(5)
          trans(2) = eri_value(11)
          trans(3) = eri_value(17)
          eri_value(5) = buff(7)*trans(2) + buff(1)*trans(3) + buff(10)*trans(1)
          eri_value(11) = buff(8)*trans(2) + buff(2)*trans(3) + buff(11)*trans(1)
          eri_value(17) = buff(9)*trans(2) + buff(3)*trans(3) + buff(12)*trans(1)
          trans(1) = eri_value(6)
          trans(2) = eri_value(12)
          trans(3) = eri_value(18)
          eri_value(6) = buff(7)*trans(2) + buff(1)*trans(3) + buff(10)*trans(1)
          eri_value(12) = buff(8)*trans(2) + buff(2)*trans(3) + buff(11)*trans(1)
          eri_value(18) = buff(9)*trans(2) + buff(3)*trans(3) + buff(12)*trans(1)

          trans(1) = eri_value(1)
          trans(2) = eri_value(2)
          trans(3) = eri_value(3)
          trans(4) = eri_value(4)
          trans(5) = eri_value(5)
          trans(6) = eri_value(6)
  eri_value(1) = buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(2) = buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(3) = buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(4) = buff(10)*buff(11)*SQRT3*trans(1) + buff(7)*buff(8)*SQRT3*trans(2) + buff(1)*buff(2)*SQRT3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*SQRT3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*SQRT3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*SQRT3*trans(6)
  eri_value(5) = buff(10)*buff(12)*SQRT3*trans(1) + buff(7)*buff(9)*SQRT3*trans(2) + buff(1)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*SQRT3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*SQRT3*trans(6)
  eri_value(6) = buff(11)*buff(12)*SQRT3*trans(1) + buff(8)*buff(9)*SQRT3*trans(2) + buff(2)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*SQRT3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*SQRT3*trans(6)

          trans(1) = eri_value(7)
          trans(2) = eri_value(8)
          trans(3) = eri_value(9)
          trans(4) = eri_value(10)
          trans(5) = eri_value(11)
          trans(6) = eri_value(12)
  eri_value(7) = buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(8) = buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(9) = buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(10) = buff(10)*buff(11)*SQRT3*trans(1) + buff(7)*buff(8)*SQRT3*trans(2) + buff(1)*buff(2)*SQRT3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*SQRT3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*SQRT3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*SQRT3*trans(6)
  eri_value(11) = buff(10)*buff(12)*SQRT3*trans(1) + buff(7)*buff(9)*SQRT3*trans(2) + buff(1)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*SQRT3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*SQRT3*trans(6)
  eri_value(12) = buff(11)*buff(12)*SQRT3*trans(1) + buff(8)*buff(9)*SQRT3*trans(2) + buff(2)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*SQRT3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*SQRT3*trans(6)

          trans(1) = eri_value(13)
          trans(2) = eri_value(14)
          trans(3) = eri_value(15)
          trans(4) = eri_value(16)
          trans(5) = eri_value(17)
          trans(6) = eri_value(18)
  eri_value(13) = buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(14) = buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(15) = buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(16) = buff(10)*buff(11)*SQRT3*trans(1) + buff(7)*buff(8)*SQRT3*trans(2) + buff(1)*buff(2)*SQRT3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*SQRT3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*SQRT3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*SQRT3*trans(6)
  eri_value(17) = buff(10)*buff(12)*SQRT3*trans(1) + buff(7)*buff(9)*SQRT3*trans(2) + buff(1)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*SQRT3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*SQRT3*trans(6)
  eri_value(18) = buff(11)*buff(12)*SQRT3*trans(1) + buff(8)*buff(9)*SQRT3*trans(2) + buff(2)*buff(3)*SQRT3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*SQRT3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*SQRT3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*SQRT3*trans(6)

          ! do i = 1, 18
          !   write(*,*) "eri_value= ", eri_value(i)
          ! enddo
          same = ish .eq. ksh .and. jsh .eq. lsh

          ii1 = res%atom_loc(ish)
          kk1 = res%atom_loc(ksh)
          maxl2 = 6

          nij = 0
          i = 1
          k = 1

          ijk_index = 1
          do j = 1, 3
            nij = nij + 1
            jj1 = j + res%atom_loc(jsh) - 1
            i2 = ii1
            j2 = jj1
            if (ii1 .lt. jj1) then ! sort <ij|
              i2 = jj1
              j2 = ii1
            end if

            ijkl_index = ijk_index
            ijk_index = ijk_index + 6
            nkl = nij

            if (same) then ! account for non-unique permutations
              itmp = min(maxl2 + 1, nkl)
              if (itmp .eq. 0) cycle
              maxl2 = itmp + 1
              nkl = nkl - itmp
            end if

            do l = 1, maxl2

              buff(1) = eri_value(ijkl_index)
              ijkl_index = ijkl_index + 1

              if (abs(buff(1)) .lt. 5.0d-11) cycle ! goto 300

              ll1 = l + res%atom_loc(lsh) - 1
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
        end if !screening
      end do !main loop
      !$omp end target teams distribute parallel do
      kernel_only2 = omp_get_wtime()
      ! write(*,*) "time in 0102 kernel", kernel_only2-kernel_only1
    end do !tiles
    kernel_full2 = omp_get_wtime()
    ! write(*,*) "time in 0102 full", kernel_full2-kernel_full1
    deallocate (n01bra, xint01bra, n02ket, xint02ket)
    ! do i = 1, size_of_matrix
    !         if(fock(i).ne.0.0_dp) write(*,*) "fa final", fock(i)
    ! end do
  end subroutine int0102
end submodule
