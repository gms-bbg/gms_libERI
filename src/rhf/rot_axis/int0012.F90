submodule(rot_axis_kernels) int0012_impl
contains
  module subroutine int0012(ss_pair, pd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, pd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n00bra(:), n12ket(:)
    real(dp), allocatable :: xint00bra(:), xint12ket(:)
    integer(kind=int64) :: nssbra, npdket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ijkl_collapsed
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutssbra, scutpdket
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d00p, d12p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) :: zero_m_0(1), zero_m_1(2), zero_m_2(3), zero_m_3(4), boys0, boys1, boys2, boys3
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21
    !r-integrals
    real(dp) :: r00(3), r01(9), r02(12), r03(10)
    real(dp) :: y03, eri_value(18), qx, qz, trans(6)
    !digestion
    logical :: kandl, same
   integer(kind=int64) :: ii1,kk1,nij,maxl2,jj1,i2,j2,ijp,nkl,itmp,ijklp,ll2,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ll1,ijkp,ijk_index,ijkl_index
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
    basc = 6
    mink = 1
    maxk = 6
    basb = 0
    minj = 1
    maxj = 1
    basa = 0
    mini = 1
    maxi = 1
    allocate (n00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (n12ket(res%n_p_shl*res%n_d_shl))
    allocate (xint12ket(res%n_p_shl*res%n_d_shl))
    !start screening
    first_screen1 = omp_get_wtime()
    scutssbra = cutoff_schwarz/maxval(pd_pair%xints)
    nssbra = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
    if (ss_pair%xints(ij) .ge. scutssbra) then
      nssbra = nssbra + 1
      xint00bra(nssbra) = ss_pair%xints(ij)
      n00bra(nssbra) = ij
    end if
    end do
    scutpdket = cutoff_schwarz/maxval(ss_pair%xints)
    npdket = 0
    do kl = 1, res%n_p_shl*res%n_d_shl
    if (pd_pair%xints(kl) .ge. scutpdket) then
      npdket = npdket + 1
      xint12ket(npdket) = pd_pair%xints(kl)
      n12ket(npdket) = kl
    end if
    end do
    first_screen2 = omp_get_wtime()

    nchunksize_int64 = 375000000
    kernel_full1 = omp_get_wtime()
    if ((nssbra*npdket) .le. nchunksize_int64) nchunksize_int64 = nssbra*npdket
    ntile = int(nssbra*npdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nssbra*npdket
      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, pd_pair, ss_pair) &
 !$omp shared(nquart_start, nquart_end, xint00bra, xint12ket) &
 !$omp shared(npdket, n00bra, n12ket) &
 !$omp private(shp_thresh, tx21, test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d00p, t_alpha, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, r00, r01, r02, r03, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, t_inverse, t_new ) &
 !$omp private(qx, qz, eri_value, kandl, same, maxl2, ijk_index, nij, nkl, ijkp, ijkl_index, ijklp) &
 !$omp private(trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/npdket + 1
        kl_tmp = mod(iquart - 1, npdket) + 1

        test = xint00bra(ij_tmp)*xint12ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then
          ij = n00bra(ij_tmp)
          kl = n12ket(kl_tmp)
          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_d_shl + 1
          lsh_tmp = mod(kl - 1, res%n_d_shl) + 1

          ii = res%i_s_shl(ish_tmp)
          jj = res%i_s_shl(jsh_tmp)
          kk = res%i_p_shl(ksh_tmp)
          ll = res%i_d_shl(lsh_tmp)

          ish = ii
          jsh = jj
          ksh = kk
          lsh = ll
          ! if(ii.eq.8.and.jj.eq.7.and.kk.eq.9.and.ll.eq.6) then
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
            t_expon_cd = pd_pair%t_expon_ab(pd_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_sp
            t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
            y03 = t_inverse_expon_cd*pd_pair%expon_a(pd_pair%pair_loc(kl) + ket_loop) !y03
            sq = pd_pair%sq(pd_pair%pair_loc(kl) + ket_loop)
            cq = pd_pair%t_alpha(pd_pair%pair_loc(kl) + ket_loop)
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
              shp_thresh = ss_pair%ismlp(ss_pair%pair_loc(ij) + bra_loop) + pd_pair%ismlp(pd_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              t_expon_ab = ss_pair%t_expon_ab(ss_pair%pair_loc(ij) + bra_loop) !tx12
              t_alpha = ss_pair%t_alpha(ss_pair%pair_loc(ij) + bra_loop) !ty02
              t_beta = ss_pair%t_beta(ss_pair%pair_loc(ij) + bra_loop) !ty01
              tx21 = ss_pair%t_inverse_expon_ab(ss_pair%pair_loc(ij) + bra_loop) !tx21
              d00p = ss_pair%d_coeff(ss_pair%pair_loc(ij) + bra_loop) !d00p
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
                !    write(*,*) "fmt", fmt
                !use extrapolation for exp(-T)
                t_new = t*m_increment !rxinc
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)
                boys2 = ((t + t)*fmt + expt)*0.2000000000000000d00 !*rmr(3)
                !    write(*,*) "boys2", boys2
                boys1 = ((t + t)*boys2 + expt)*0.333333333333333333d00 !*rmr(2)
                !    write(*,*) "boys1", boys1
                boys0 = ((t + t)*boys1 + expt)!*rmr(1)
                !    write(*,*) "boys0", boys0
                boys0 = boys0*sqrt(expon_abcd_inverse)*d00p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d00p*rho
                boys2 = boys2*sqrt(expon_abcd_inverse)*d00p*rho*rho
                boys3 = fmt*sqrt(expon_abcd_inverse)*d00p*rho*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d00p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
                boys3 = boys2*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho)
              end if
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0
              zero_m_1(1) = zero_m_1(1) + boys1
              zero_m_1(2) = zero_m_1(2) + boys1*pqr
              zero_m_2(1) = zero_m_2(1) + boys2
              zero_m_2(2) = zero_m_2(2) + boys2*pqr
              zero_m_2(3) = zero_m_2(3) + boys2*pqs
              zero_m_3(1) = zero_m_3(1) + boys3
              zero_m_3(2) = zero_m_3(2) + boys3*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*pqr*pqs
            end do
            ! write(*,*) "FQD0(1) = ", zero_m_0(1)
            ! write(*,*) "FQD1(1) = ", zero_m_1(1)
            ! write(*,*) "FQD1(2) = ", zero_m_1(2)
            ! write(*,*) "FQD2(1) = ", zero_m_2(1)
            ! write(*,*) "FQD2(2) = ", zero_m_2(2)
            ! write(*,*) "FQD2(3) = ", zero_m_2(3)
            ! write(*,*) "FQD3(1) = ", zero_m_3(1)
            ! write(*,*) "FQD3(2) = ", zero_m_3(2)
            ! write(*,*) "FQD3(3) = ", zero_m_3(3)
            ! write(*,*) "FQD3(4) = ", zero_m_3(4)
            !r integrals here
            r00(1) = r00(1) - zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq*y03
            r00(2) = r00(2) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq*(1 - y03)
            r00(3) = r00(3) + zero_m_0(1)*y03*y03*(1 - y03)*sq
            ! write(*,*) "r00(1) = ", r00(1)
            ! write(*,*) "r00(2) = ", r00(2)
            ! write(*,*) "r00(3) = ", r00(3)
            r01(1) = r01(1) - zero_m_1(1)*aqx*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq
            r01(2) = r01(2) - zero_m_1(1)*acy*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq
            r01(3) = r01(3) + zero_m_1(2)*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq
            r01(4) = r01(4) - zero_m_1(1)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r01(5) = r01(5) - zero_m_1(1)*acy*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r01(6) = r01(6) + zero_m_1(2)*t_inverse_expon_cd*0.5_dp*sq*y03*y03
            r01(7) = r01(7) + zero_m_1(1)*aqx*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            r01(8) = r01(8) + zero_m_1(1)*acy*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            r01(9) = r01(9) - zero_m_1(2)*t_inverse_expon_cd*0.5_dp*sq*y03*(1 - y03)
            ! write(*,*) "r01(1) = ", r01(1)
            ! write(*,*) "r01(2) = ", r01(2)
            ! write(*,*) "r01(3) = ", r01(3)
            ! write(*,*) "r01(4) = ", r01(4)
            ! write(*,*) "r01(5) = ", r01(5)
            ! write(*,*) "r01(6) = ", r01(6)
            ! write(*,*) "r01(7) = ", r01(7)
            ! write(*,*) "r01(8) = ", r01(8)
            ! write(*,*) "r01(9) = ", r01(9)
            r02(1) = r02(1) - (zero_m_2(1)*aqx2 - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03
            r02(2) = r02(2) - (zero_m_2(1)*acy2 - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03
            r02(3) = r02(3) - (zero_m_2(3) - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03
            r02(4) = r02(4) - zero_m_2(1)*aqxy*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03
            r02(5) = r02(5) + zero_m_2(2)*aqx*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03
            r02(6) = r02(6) + zero_m_2(2)*acy*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*y03

            r02(7) = r02(7) + (zero_m_2(1)*aqx2 - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)
            r02(8) = r02(8) + (zero_m_2(1)*acy2 - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)
            r02(9) = r02(9) + (zero_m_2(3) - zero_m_1(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)
            r02(10) = r02(10) + zero_m_2(1)*aqxy*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)
            r02(11) = r02(11) - zero_m_2(2)*aqx*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)
            r02(12) = r02(12) - zero_m_2(2)*acy*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*sq*(1 - y03)

            ! r02(1,2)= r02(1,2)+fqd2(1)*aqx2-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(1,3)= r02(1,3)+fqd2(1)*aqx2-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(2,2)= r02(2,2)+fqd2(1)*acy2-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(2,3)= r02(2,3)+fqd2(1)*acy2-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(3,2)= r02(3,2)+fqd2(3)-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(3,3)= r02(3,3)+fqd2(3)-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(4,2)= r02(4,2)+fqd2(1)*aqxy*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(4,3)= r02(4,3)+fqd2(1)*aqxy*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(5,2)= r02(5,2)+-fqd2(2)*aqx*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(5,3)= r02(5,3)+-fqd2(2)*aqx*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(6,2)= r02(6,2)+-fqd2(2)*acy*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(6,3)= r02(6,3)+-fqd2(2)*acy*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04

            ! r02(1,2)= r02(1,2)+fqd2(1)*aqx2-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(2,2)= r02(2,2)+fqd2(1)*acy2-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(3,2)= r02(3,2)+fqd2(3)-fqd1(1)*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(4,2)= r02(4,2)+fqd2(1)*aqxy*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(5,2)= r02(5,2)+-fqd2(2)*aqx*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03
            ! r02(6,2)= r02(6,2)+-fqd2(2)*acy*-x43 *0.5_dp*x43 *0.5_dp*sq(2)*y03

            ! r02(1,3)= r02(1,3)+fqd2(1)*aqx2-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(2,3)= r02(2,3)+fqd2(1)*acy2-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(3,3)= r02(3,3)+fqd2(3)-fqd1(1)*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(4,3)= r02(4,3)+fqd2(1)*aqxy*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(5,3)= r02(5,3)+-fqd2(2)*aqx*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04
            ! r02(6,3)= r02(6,3)+-fqd2(2)*acy*x43 *0.5_dp*x43 *0.5_dp*sq(2)*y04

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
 r03(1) = r03(1) - (zero_m_3(1)*aqx2 - zero_m_2(1)*3.0_dp)*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*aqx
     r03(2) = r03(2) - (zero_m_3(1)*aqx2 - zero_m_2(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*acy
         r03(3) = r03(3) + (zero_m_3(2)*aqx2 - zero_m_2(2))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq
     r03(4) = r03(4) - (zero_m_3(1)*acy2 - zero_m_2(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*aqx
            r03(5) = r03(5) + zero_m_3(2)*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*aqxy
          r03(6) = r03(6) - (zero_m_3(3) - zero_m_2(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*aqx
 r03(7) = r03(7) - (zero_m_3(1)*acy2 - zero_m_2(1)*3.0_dp)*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*acy
         r03(8) = r03(8) + (zero_m_3(2)*acy2 - zero_m_2(2))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq
          r03(9) = r03(9) - (zero_m_3(3) - zero_m_2(1))*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq*acy
        r03(10) = r03(10) + (zero_m_3(4) - zero_m_2(2)*3.0_dp)*t_inverse_expon_cd*t_inverse_expon_cd*0.25_dp*t_inverse_expon_cd*0.5_dp*sq
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

 eri_value(1) = qx*qx*qx*r00(3) + 3.0_dp*r01(1) + qx*qx*(r01(4) + 2.0_dp*r01(7)) + qx*(2.0_dp*r00(1) + r00(2) + 2.0_dp*r02(1) + r02(7)) + r03(1)
          eri_value(2) = r01(1) + qx*(r00(2) + r02(8)) + r03(4)
          eri_value(3) = qx*qz*qz*r00(3) + r01(1) + qz*qz*r01(4) + 2*qx*qz*r01(9) + qx*(r00(2) + r02(9)) + 2*qz*r02(5) + r03(6)
          eri_value(4) = r01(2) + qx*qx*r01(8) + qx*(r02(4) + r02(10)) + r03(2)
  eri_value(5) = qx*qx*qz*r00(3) + qx*qz*(r01(4) + r01(7)) + r01(3) + qx*qx*r01(9) + qz*(r00(1) + r02(1)) + qx*(r02(5) + r02(11)) + r03(3)
          eri_value(6) = qx*qz*r01(8) + qz*r02(4) + qx*r02(12) + r03(5)
          eri_value(7) = r01(2) + qx*qx*r01(5) + 2*qx*r02(4) + r03(2)
          eri_value(8) = 3.0_dp*r01(2) + r03(7)
          eri_value(9) = r01(2) + qz*qz*r01(5) + 2*qz*r02(6) + r03(9)
          eri_value(10) = r01(1) + qx*(r00(1) + r02(2)) + r03(4)
          eri_value(11) = qx*qz*r01(5) + qz*r02(4) + qx*r02(6) + r03(5)
          eri_value(12) = r01(3) + qz*(r00(1) + r02(2)) + r03(8)
          eri_value(13) = qx*qx*qz*r00(3) + 2*qx*qz*r01(7) + r01(3) + qx*qx*r01(6) + qz*(r00(2) + r02(7)) + 2*qx*r02(5) + r03(3)
          eri_value(14) = r01(3) + qz*(r00(2) + r02(8)) + r03(8)
  eri_value(15) = qz*qz*qz*r00(3) + 3.0_dp*r01(3) + qz*qz*(r01(6) + 2.0_dp*r01(9)) + qz*(2.0_dp*r00(1) + r00(2) + 2.0_dp*r02(3) + r02(9)) + r03(10)
          eri_value(16) = qx*qz*r01(8) + qz*r02(10) + qx*r02(6) + r03(5)
  eri_value(17) = qx*qz*qz*r00(3) + r01(1) + qz*qz*r01(7) + qx*qz*(r01(6) + r01(9)) + qx*(r00(1) + r02(3)) + qz*(r02(5) + r02(11)) + r03(6)
          eri_value(18) = r01(2) + qz*qz*r01(8) + qz*(r02(6) + r02(12)) + r03(9)
          ! do i=1,18
          !     write(*,*) "eri_value first= ", eri_value(i)
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
          ! write(*,*) "trans(1) = ", trans(1)
          ! write(*,*) "trans(2) = ", trans(2)
          ! write(*,*) "trans(3) = ", trans(3)
          eri_value(2) = buff(1)*trans(3) + buff(10)*trans(1) + buff(7)*trans(2)
          eri_value(8) = buff(2)*trans(3) + buff(11)*trans(1) + buff(8)*trans(2)
          eri_value(14) = buff(3)*trans(3) + buff(12)*trans(1) + buff(9)*trans(2)
          ! write(*,*) "eri_value(1) = ", eri_value(2)
          ! write(*,*) "eri_value(7) = ", eri_value(8)
          ! write(*,*) "eri_value(13) = ", eri_value(14)
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

          ! do i=1,18
          !     write(*,*) "eri_value= ", eri_value(i)
          ! enddo
          ! endif

          kandl = (ksh .eq. lsh)
          same = (ish .eq. ksh .and. jsh .eq. lsh)

          ii1 = res%atom_loc(ish)
          jj1 = res%atom_loc(jsh)
          ! write(*,*) "ii1", ii1
          ! write(*,*) "jj1", jj1

          maxl2 = 6
          IJK_INDEX = 1
          i2 = ii1
          j2 = jj1
          if (ii1 .lt. jj1) then ! sort <ij|
            i2 = jj1
            j2 = ii1
          end if
          NIJ = 0
          NKL = NIJ
          do k = 1, 3
            if (kandl) maxl2 = k

            kk1 = k + res%atom_loc(ksh) - 1

            IJKL_INDEX = IJK_INDEX
            IJK_INDEX = IJK_INDEX + 6

            do l = 1, maxl2

              buff(1) = eri_value(IJKL_INDEX)

              IJKL_INDEX = IJKL_INDEX + 1
              !    write(*,*) "IJKL_INDEX", IJKL_INDEX

              if (abs(buff(1)) .lt. 5.0d-11) cycle ! goto 300
              !    write(*,*) "val", buff(1)

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

    end do !tiles

    deallocate (n00bra, xint00bra, n12ket, xint12ket)

  end subroutine int0012
end submodule
