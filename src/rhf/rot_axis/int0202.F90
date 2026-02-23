submodule(rot_axis_kernels) int0202_impl
contains
  module subroutine int0202(sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n02bra(:)
    real(dp), allocatable :: xint02bra(:)
    integer(kind=int64) :: nsdbra
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutsdbra
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d02p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21, y03
    real(dp) :: zero_m_0(2), zero_m_1(7), zero_m_2(12), zero_m_3(8), zero_m_4(5), boys0, boys1, boys2, boys3, boys4
    !r-integrals
    real(dp) :: r00(4), r01(12), r02(30), r03(20), r04(15), xmdt
    !eri_value
    real(dp) :: eri_value(36), trans(6), qx, qz
    !digestion
    logical :: kandl, iandj, same
  integer(kind=int64) :: ii1,kk1,nij,maxl2,jj1,i2,j2,ijp,nkl,itmp,ijklp,ll2,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ll1,ijkp,ij_index,ijk_index,ijkl_index
    integer(kind=int64) :: maxj2, loci, locj, lock, locl, ip
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: istart, iend, itile, ntile, ijkl_collapsed
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp
    integer(kind=int64) :: niijj, iijj_tmp, kkll_tmp, iijj, kkll, iijjkkll
    integer :: shp_thresh
    mini = 1  !0
    maxi = 1  !2
    minj = 1  !0
    maxj = 6  !2
    mink = 1
    maxk = 1
    minl = 1
    maxl = 6
    allocate (n02bra(res%n_s_shl*res%n_d_shl))
    allocate (xint02bra(res%n_s_shl*res%n_d_shl))

    !start screening
    scutsdbra = cutoff_schwarz/maxval(sd_pair%xints)
    nsdbra = 0
    do ij = 1, res%n_s_shl*res%n_d_shl
    if (sd_pair%xints(ij) .ge. scutsdbra) then
      nsdbra = nsdbra + 1
      xint02bra(nsdbra) = sd_pair%xints(ij)
      n02bra(nsdbra) = ij
    end if
    end do

    !--multi-GPU--work
    nchunk = nsdbra/res%n_size
    nquart_start = nchunk*res%n_rank + 1
    nquart_end = nquart_start + nchunk - 1
    if (res%n_rank .EQ. res%n_size - 1) nquart_end = nsdbra
 !$omp target teams distribute parallel do collapse(2) default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, sd_pair) &
 !$omp shared(nquart_start, nquart_end, xint02bra) &
 !$omp shared(nsdbra, n02bra) &
 !$omp shared(mini, maxi, minj, maxj, mink, maxk, minl, maxl) &
 !$omp private(shp_thresh, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d02p, tx21, t_beta, t_alpha, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, zero_m_4, r00, r01, r02, r03, r04, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, boys4, t_inverse, t_new ) &
 !$omp private(xmdt, qx, qz, eri_value, same, maxl2, ijk_index, nij, nkl, ijkl_index, ijklp) &
 !$omp private(iandj, kandl, loci, locj, lock, locl, ip, ijp, itmp, test, maxj2,ijkp) &
 !$omp private(trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
    do iijj = nquart_start, nquart_end
      do kkll = 1, nsdbra

        if (kkll .gt. iijj) cycle !needed to add this
        test = xint02bra(iijj)*xint02bra(kkll)
        if (test .gt. cutoff_schwarz) then

          ij = n02bra(iijj)
          kl = n02bra(kkll)

          ish_tmp = (ij - 1)/res%n_d_shl + 1
          jsh_tmp = mod(ij - 1, res%n_d_shl) + 1
          ksh_tmp = (kl - 1)/res%n_d_shl + 1
          lsh_tmp = mod(kl - 1, res%n_d_shl) + 1

          ish = res%i_s_shl(ish_tmp)
          jsh = res%i_d_shl(jsh_tmp)
          ksh = res%i_s_shl(ksh_tmp)
          lsh = res%i_d_shl(lsh_tmp)

          ! if(ii.eq.3.and.jj.eq.6.and.kk.eq.1.and.ll.eq.6) then
          ! write(*,*) "---"
          ! write(*,*) "ii", ii
          ! write(*,*) "jj", jj
          ! write(*,*) "kk", kk
          ! write(*,*) "ll", ll
          ! write(*,*) "---"
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
  acx = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* (buff(8)*buff(3) - buff(9)*buff(2))+(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* (buff(9)*buff(1) - buff(7)*buff(3)) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* (buff(7)*buff(2) - buff(8)*buff(1)) !acx
  acy = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* buff(7) +(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* buff(8) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* buff(9) !acy
  acz = (res%coord_sh(ksh,1) - res%coord_sh(ish,1))* buff(1) +(res%coord_sh(ksh,2) - res%coord_sh(ish,2))* buff(2) +(res%coord_sh(ksh,3) - res%coord_sh(ish,3))* buff(3) !acz
          !-----------acy2 stuff check -------
          acy2 = acy*acy !acy2
          r00 = 0.0_dp
          r01 = 0.0_dp
          ! r02 = 0.0_dp
          ! r03 = 0.0_dp
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
          r04 = 0.0_dp
          ket_loop = 0
          do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1
            t_expon_cd = sd_pair%t_expon_ab(sd_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_sd
            t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
            y03 = t_inverse_expon_cd*sd_pair%expon_a(sd_pair%pair_loc(kl) + ket_loop) !y03
            sq = sd_pair%sq(sd_pair%pair_loc(kl) + ket_loop)
            cq = sd_pair%t_alpha(sd_pair%pair_loc(kl) + ket_loop) !cq
            cqx = cq*sing !cqx
            cqz = cq*cosg !cqz
            aqx = acx + cqx !aqx
            aqx2 = aqx*aqx !aqx2
            aqxy = aqx*acy !aqxy
            aqz = acz + cqz !aqz
            qps = aqx2 + acy2 !qps
            zero_m_0 = 0.0_dp
            zero_m_1 = 0.0_dp
            zero_m_2 = 0.0_dp
            zero_m_3 = 0.0_dp
            zero_m_4 = 0.0_dp
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              bra_loop = bra_loop + 1
              !get bra shell pair info
              shp_thresh = sd_pair%ismlp(sd_pair%pair_loc(ij) + bra_loop) + sd_pair%ismlp(sd_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              t_expon_ab = sd_pair%t_expon_ab(sd_pair%pair_loc(ij) + bra_loop) !tx12
              t_alpha = sd_pair%t_alpha(sd_pair%pair_loc(ij) + bra_loop) !ty02
              t_beta = sd_pair%t_beta(sd_pair%pair_loc(ij) + bra_loop) !ty01
              tx21 = sd_pair%t_inverse_expon_ab(sd_pair%pair_loc(ij) + bra_loop) !tx21
              d02p = sd_pair%d_coeff(sd_pair%pair_loc(ij) + bra_loop) !d02p
              expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd) !x41
              pqr = t_alpha - aqz
              pqs = pqr*pqr
              rho = t_expon_ab*t_expon_cd*expon_abcd_inverse
              t = (pqs + qps)*rho
              rho = rho + rho
              if (t .le. t_max) then
                !evaluate boys function with recursion
                t_new = t*1.5365936651378012d+01 !f_increment(5)
                t_int = nint(t_new)
                fmt = boys_grid_zero((4*451*5) + (t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((4*451*5) + (t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((4*451*5) + (t_int*5) + 1)
                !use extrapolation for exp(-T)
                t_new = t*m_increment !rxinc
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)
                boys3 = ((t + t)*fmt + expt)*0.142857142857142857d00 !*rmr(4)
                boys2 = ((t + t)*boys3 + expt)*0.2000000000000000d00 !*rmr(3)
                boys1 = ((t + t)*boys2 + expt)*0.333333333333333333d00 !*rmr(2)
                boys0 = ((t + t)*boys1 + expt)!*rmr(1)
                boys0 = boys0*sqrt(expon_abcd_inverse)*d02p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d02p*rho
                boys2 = boys2*sqrt(expon_abcd_inverse)*d02p*rho*rho
                boys3 = boys3*sqrt(expon_abcd_inverse)*d02p*rho*rho*rho
                boys4 = fmt*sqrt(expon_abcd_inverse)*d02p*rho*rho*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d02p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
                boys3 = boys2*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho)
                boys4 = boys3*(t_inverse*0.5_dp*rho + t_inverse*rho + t_inverse*rho + t_inverse*rho)
              end if
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0*tx21
              zero_m_0(2) = zero_m_0(2) + boys0*t_beta*t_beta !ok

              zero_m_1(1) = zero_m_1(1) + boys1*tx21
              zero_m_1(2) = zero_m_1(2) + boys1*tx21*pqr
              zero_m_1(3) = zero_m_1(3) + boys1*t_beta*t_beta
              zero_m_1(4) = zero_m_1(4) + boys1*t_beta*t_beta*pqr
              zero_m_1(5) = zero_m_1(5) + boys1*tx21*t_beta
              zero_m_1(6) = zero_m_1(6) + boys1*tx21*t_beta*pqr
              zero_m_1(7) = zero_m_1(7) + boys1*tx21*tx21 !ok

              zero_m_2(1) = zero_m_2(1) + boys2*tx21
              zero_m_2(2) = zero_m_2(2) + boys2*tx21*pqr
              zero_m_2(3) = zero_m_2(3) + boys2*tx21*pqs
              zero_m_2(4) = zero_m_2(4) + boys2*t_beta*t_beta
              zero_m_2(5) = zero_m_2(5) + boys2*t_beta*t_beta*pqr
              zero_m_2(6) = zero_m_2(6) + boys2*t_beta*t_beta*pqs
              zero_m_2(7) = zero_m_2(7) + boys2*tx21*t_beta
              zero_m_2(8) = zero_m_2(8) + boys2*tx21*t_beta*pqr
              zero_m_2(9) = zero_m_2(9) + boys2*tx21*t_beta*pqs
              zero_m_2(10) = zero_m_2(10) + boys2*tx21*tx21
              zero_m_2(11) = zero_m_2(11) + boys2*tx21*tx21*pqr
              zero_m_2(12) = zero_m_2(12) + boys2*tx21*tx21*pqs

              zero_m_3(1) = zero_m_3(1) + boys3*tx21*t_beta
              zero_m_3(2) = zero_m_3(2) + boys3*tx21*t_beta*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*tx21*t_beta*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*tx21*t_beta*pqr*pqs
              zero_m_3(5) = zero_m_3(5) + boys3*tx21*tx21
              zero_m_3(6) = zero_m_3(6) + boys3*tx21*tx21*pqr
              zero_m_3(7) = zero_m_3(7) + boys3*tx21*tx21*pqs
              zero_m_3(8) = zero_m_3(8) + boys3*tx21*tx21*pqr*pqs

              zero_m_4(1) = zero_m_4(1) + boys4*tx21*tx21
              zero_m_4(2) = zero_m_4(2) + boys4*tx21*tx21*pqr
              zero_m_4(3) = zero_m_4(3) + boys4*tx21*tx21*pqs
              zero_m_4(4) = zero_m_4(4) + boys4*tx21*tx21*pqr*pqs
              zero_m_4(5) = zero_m_4(5) + boys4*tx21*tx21*pqr*pqs*pqr
            end do
            ! write(*,*) "boys4", boys4
            ! do i = 1, 2
            !   write(*,*) "zero_m_0", zero_m_0(i)
            ! enddo

            ! do i = 1, 7
            !   write(*,*) "zero_m_1", zero_m_1(i)
            ! enddo

            ! do i = 1, 12
            !   write(*,*) "zero_m_2", zero_m_2(i)
            ! enddo

            ! do i = 1, 8
            !   write(*,*) "zero_m_3", zero_m_3(i)
            ! enddo

            ! do i = 1, 5
            !   write(*,*) "zero_m_4", zero_m_4(i)
            ! enddo
            !r integrals here
            r00(1) = r00(1) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq
            r00(2) = r00(2) + zero_m_0(1)*y03*y03*sq
            r00(3) = r00(3) + zero_m_0(2)*t_inverse_expon_cd*0.5_dp*sq
            r00(4) = r00(4) + zero_m_0(2)*y03*y03*sq
            ! do i = 1, 4
            !    write(*,*) "r00", r00(i)
            ! enddo
            xmdt = t_inverse_expon_cd*0.5_dp*sq*y03
            r01(1) = r01(1) + zero_m_1(1)*xmdt*aqx
            r01(2) = r01(2) + zero_m_1(1)*xmdt*acy
            r01(3) = r01(3) - zero_m_1(2)*xmdt
            r01(4) = r01(4) + zero_m_1(3)*xmdt*aqx
            r01(5) = r01(5) + zero_m_1(3)*xmdt*acy
            r01(6) = r01(6) - zero_m_1(4)*xmdt
            xmdt = t_inverse_expon_cd*0.5_dp*sq
            r01(7) = r01(7) - zero_m_1(5)*aqx*xmdt
            r01(8) = r01(8) - zero_m_1(5)*acy*xmdt
            r01(9) = r01(9) + zero_m_1(6)*xmdt
            xmdt = y03*y03*sq
            r01(10) = r01(10) - zero_m_1(5)*aqx*xmdt
            r01(11) = r01(11) - zero_m_1(5)*acy*xmdt
            r01(12) = r01(12) + zero_m_1(6)*xmdt
            ! do i = 1, 12
            ! write(*,*) "r01", r01(i)
            ! enddo
            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*xmdt
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*xmdt
            r02(3) = r02(3) + (zero_m_2(3) - zero_m_1(1))*xmdt
            r02(4) = r02(4) + zero_m_2(1)*aqxy*xmdt
            r02(5) = r02(5) - zero_m_2(2)*aqx*xmdt
            r02(6) = r02(6) - zero_m_2(2)*acy*xmdt
            r02(7) = r02(7) + (zero_m_2(4)*aqx2 - zero_m_1(3))*xmdt
            r02(8) = r02(8) + (zero_m_2(4)*acy2 - zero_m_1(3))*xmdt
            r02(9) = r02(9) + (zero_m_2(6) - zero_m_1(3))*xmdt
            r02(10) = r02(10) + zero_m_2(4)*aqxy*xmdt
            r02(11) = r02(11) - zero_m_2(5)*aqx*xmdt
            r02(12) = r02(12) - zero_m_2(5)*acy*xmdt
            xmdt = t_inverse_expon_cd*0.5_dp*sq*y03
            r02(13) = r02(13) - (zero_m_2(7)*aqx2 - zero_m_1(5))*xmdt
            r02(14) = r02(14) - (zero_m_2(7)*acy2 - zero_m_1(5))*xmdt
            r02(15) = r02(15) - (zero_m_2(9) - zero_m_1(5))*xmdt
            r02(16) = r02(16) - zero_m_2(7)*aqxy*xmdt
            r02(17) = r02(17) + zero_m_2(8)*aqx*xmdt
            r02(18) = r02(18) + zero_m_2(8)*acy*xmdt
            xmdt = t_inverse_expon_cd*0.5_dp*sq
            r02(19) = r02(19) + (zero_m_2(10)*aqx2 - zero_m_1(7))*xmdt
            r02(20) = r02(20) + (zero_m_2(10)*acy2 - zero_m_1(7))*xmdt
            r02(21) = r02(21) + (zero_m_2(12) - zero_m_1(7))*xmdt
            r02(22) = r02(22) + zero_m_2(10)*aqxy*xmdt
            r02(23) = r02(23) - zero_m_2(11)*aqx*xmdt
            r02(24) = r02(24) - zero_m_2(11)*acy*xmdt
            xmdt = y03*y03*sq
            r02(25) = r02(25) + (zero_m_2(10)*aqx2 - zero_m_1(7))*xmdt
            r02(26) = r02(26) + (zero_m_2(10)*acy2 - zero_m_1(7))*xmdt
            r02(27) = r02(27) + (zero_m_2(12) - zero_m_1(7))*xmdt
            r02(28) = r02(28) + zero_m_2(10)*aqxy*xmdt
            r02(29) = r02(29) - zero_m_2(11)*aqx*xmdt
            r02(30) = r02(30) - zero_m_2(11)*acy*xmdt
            ! do i = 1, 30
            !   write(*,*) "r02", r02(i)
            ! enddo
            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r03(1) = r03(1) - (zero_m_3(1)*aqx2 - zero_m_2(7)*3.0_dp)*xmdt*aqx
            r03(2) = r03(2) - (zero_m_3(1)*aqx2 - zero_m_2(7))*xmdt*acy
            r03(3) = r03(3) + (zero_m_3(2)*aqx2 - zero_m_2(8))*xmdt
            r03(4) = r03(4) - (zero_m_3(1)*acy2 - zero_m_2(7))*xmdt*aqx
            r03(5) = r03(5) + zero_m_3(2)*xmdt*aqxy
            r03(6) = r03(6) - (zero_m_3(3) - zero_m_2(7))*xmdt*aqx
            r03(7) = r03(7) - (zero_m_3(1)*acy2 - zero_m_2(7)*3.0_dp)*xmdt*acy
            r03(8) = r03(8) + (zero_m_3(2)*acy2 - zero_m_2(8))*xmdt
            r03(9) = r03(9) - (zero_m_3(3) - zero_m_2(7))*xmdt*acy
            r03(10) = r03(10) + (zero_m_3(4) - zero_m_2(8)*3.0_dp)*xmdt
            xmdt = t_inverse_expon_cd*0.5_dp*sq*y03
            r03(11) = r03(11) + (zero_m_3(5)*aqx2 - zero_m_2(10)*3.0_dp)*xmdt*aqx
            r03(12) = r03(12) + (zero_m_3(5)*aqx2 - zero_m_2(10))*xmdt*acy
            r03(13) = r03(13) - (zero_m_3(6)*aqx2 - zero_m_2(11))*xmdt
            r03(14) = r03(14) + (zero_m_3(5)*acy2 - zero_m_2(10))*xmdt*aqx
            r03(15) = r03(15) - zero_m_3(6)*xmdt*aqxy
            r03(16) = r03(16) + (zero_m_3(7) - zero_m_2(10))*xmdt*aqx
            r03(17) = r03(17) + (zero_m_3(5)*acy2 - zero_m_2(10)*3.0_dp)*xmdt*acy
            r03(18) = r03(18) - (zero_m_3(6)*acy2 - zero_m_2(11))*xmdt
            r03(19) = r03(19) + (zero_m_3(7) - zero_m_2(10))*xmdt*acy
            r03(20) = r03(20) - (zero_m_3(8) - zero_m_2(11)*3.0_dp)*xmdt
            ! do i = 1, 20
            !   write(*,*) "r03", r03(i)
            ! enddo
            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r04(1) = r04(1) + (zero_m_4(1)*aqx2*aqx2 - zero_m_3(5)*6.0_dp*aqx2 + zero_m_2(10)*3.0_dp)*xmdt
            r04(2) = r04(2) + (zero_m_4(1)*aqx2 - zero_m_3(5)*3.0_dp)*xmdt*aqxy
            r04(3) = r04(3) - (zero_m_4(2)*aqx2 - zero_m_3(6)*3.0_dp)*xmdt*aqx
            r04(4) = r04(4) + (zero_m_4(1)*aqx2*acy2 - zero_m_3(5)*(aqx2 + acy2) + zero_m_2(10))*xmdt
            r04(5) = r04(5) - (zero_m_4(2)*aqx2 - zero_m_3(6))*xmdt*acy
            r04(6) = r04(6) + (zero_m_4(3)*aqx2 - zero_m_3(5)*aqx2 - zero_m_3(7) + zero_m_2(10))*xmdt
            r04(7) = r04(7) + (zero_m_4(1)*acy2 - zero_m_3(5)*3.0_dp)*xmdt*aqxy
            r04(8) = r04(8) - (zero_m_4(2)*acy2 - zero_m_3(6))*xmdt*aqx
            r04(9) = r04(9) + (zero_m_4(3) - zero_m_3(5))*xmdt*aqxy
            r04(10) = r04(10) - (zero_m_4(4) - zero_m_3(6)*3.0_dp)*xmdt*aqx
            r04(11) = r04(11) + (zero_m_4(1)*acy2*acy2 - zero_m_3(5)*6.0_dp*acy2 + zero_m_2(10)*3.0_dp)*xmdt
            r04(12) = r04(12) - (zero_m_4(2)*acy2 - zero_m_3(6)*3.0_dp)*xmdt*acy
            r04(13) = r04(13) + (zero_m_4(3)*acy2 - zero_m_3(5)*acy2 - zero_m_3(7) + zero_m_2(10))*xmdt
            r04(14) = r04(14) - (zero_m_4(4) - zero_m_3(6)*3.0_dp)*xmdt*acy
            r04(15) = r04(15) + (zero_m_4(5) - zero_m_3(7)*6.0_dp + zero_m_2(10)*3.0_dp)*xmdt
            ! do i = 1, 15
            !   write(*,*) "r04", r04(i)
            ! enddo
          end do !end ket loop
          qx = rcd*sing
          qz = rcd*cosg

          eri_value(1) = r02(1) + r02(19) + r00(1) + qx*(2*r01(1) + 2*r03(11) + qx*(r02(25) + r00(2))) + r04(1)
          eri_value(2) = r02(2) + r02(19) + r00(1) + r04(4)
          eri_value(3) = r02(3) + r02(19) + r00(1) + qz*(2*r01(3) + 2*r03(13) + qz*(r02(25) + r00(2))) + r04(6)
          eri_value(4) = r02(4) + qx*(r01(2) + r03(12)) + r04(2)
          eri_value(5) = r02(5) + qx*(r01(3) + r03(13)) + qz*(r01(1) + r03(11) + qx*(r02(25) + r00(2))) + r04(3)
          eri_value(6) = r02(6) + qz*(r01(2) + r03(12)) + r04(5)
          eri_value(7) = r02(1) + r02(20) + r00(1) + qx*(2*r01(1) + 2*r03(14) + qx*(r02(26) + r00(2))) + r04(4)
          eri_value(8) = r02(2) + r02(20) + r00(1) + r04(11)
          eri_value(9) = r02(3) + r02(20) + r00(1) + qz*(2*r01(3) + 2*r03(18) + qz*(r02(26) + r00(2))) + r04(13)
          eri_value(10) = r02(4) + qx*(r01(2) + r03(17)) + r04(7)
          eri_value(11) = r02(5) + qx*(r01(3) + r03(18)) + qz*(r01(1) + r03(14) + qx*(r02(26) + r00(2))) + r04(8)
          eri_value(12) = r02(6) + qz*(r01(2) + r03(17)) + r04(12)
  eri_value(13) = -(2.0_dp*r01(9)) + r02(1) + r02(7) + r02(21) - 2.0_dp*r03(3) + r00(1) + r00(3) + qx*(2*r01(1) + 2*r01(4) - 2*2.0_dp*r02(17) + 2*r03(16) + qx*(-(2.0_dp*r01(12)) + r02(27) + r00(2) + r00(4))) + r04(6)
          eri_value(14) = -(2.0_dp*r01(9)) + r02(2) + r02(8) + r02(21) - 2.0_dp*r03(8) + r00(1) + r00(3) + r04(13)
  eri_value(15) = -(2.0_dp*r01(9)) + r02(3) + r02(9) + r02(21) - 2.0_dp*r03(10) + r00(1) + r00(3) + qz*(2*r01(3) + 2*r01(6) - 2*2.0_dp*r02(15) + 2*r03(20) + qz*(-(2.0_dp*r01(12)) + r02(27) + r00(2) + r00(4))) + r04(15)
          eri_value(16) = r02(4) + r02(10) - 2.0_dp*r03(5) + qx*(r01(2) + r01(5) - 2.0_dp*r02(18) + r03(19)) + r04(9)
  eri_value(17) = r02(5) + r02(11) - 2.0_dp*r03(6) + qx*(r01(3) + r01(6) - 2.0_dp*r02(15) + r03(20)) + qz*(r01(1) + r01(4) - 2.0_dp*r02(17) + r03(16) + qx*(-(2.0_dp*r01(12)) + r02(27) + r00(2) + r00(4))) + r04(10)
          eri_value(18) = r02(6) + r02(12) - 2.0_dp*r03(9) + qz*(r01(2) + r01(5) - 2.0_dp*r02(18) + r03(19)) + r04(14)
          eri_value(19) = r02(22) + qx*(qx*r02(28) + 2*r03(12)) + r04(2)
          eri_value(20) = r02(22) + r04(7)
          eri_value(21) = r02(22) + qz*(qz*r02(28) + 2*r03(15)) + r04(9)
          eri_value(22) = qx*r03(14) + r04(4)
          eri_value(23) = qz*(qx*r02(28) + r03(12)) + qx*r03(15) + r04(5)
          eri_value(24) = qz*r03(14) + r04(8)
          eri_value(25) = -r01(7) + r02(23) - r03(1) + qx*(-2*r02(13) + qx*(-r01(10) + r02(29)) + 2*r03(13)) + r04(3)
          eri_value(26) = -r01(7) + r02(23) - r03(4) + r04(8)
          eri_value(27) = -r01(7) + r02(23) - r03(6) + qz*(-2*r02(17) + qz*(-r01(10) + r02(29)) + 2*r03(16)) + r04(10)
          eri_value(28) = -r03(2) + qx*(-r02(16) + r03(15)) + r04(5)
          eri_value(29) = -r03(3) + qz*(-r02(13) + qx*(-r01(10) + r02(29)) + r03(13)) + qx*(-r02(17) + r03(16)) + r04(6)
          eri_value(30) = -r03(5) + qz*(-r02(16) + r03(15)) + r04(9)
          eri_value(31) = -r01(8) + r02(24) - r03(2) + qx*(-2*r02(16) + qx*(-r01(11) + r02(30)) + 2*r03(15)) + r04(5)
          eri_value(32) = -r01(8) + r02(24) - r03(7) + r04(12)
          eri_value(33) = -r01(8) + r02(24) - r03(9) + qz*(-2*r02(18) + qz*(-r01(11) + r02(30)) + 2*r03(19)) + r04(14)
          eri_value(34) = -r03(4) + qx*(-r02(14) + r03(18)) + r04(8)
          eri_value(35) = -r03(5) + qz*(-r02(16) + qx*(-r01(11) + r02(30)) + r03(15)) + qx*(-r02(18) + r03(19)) + r04(9)
          eri_value(36) = -r03(8) + qz*(-r02(14) + r03(18)) + r04(13)

          ! do i = 1, 36
          !   write(*,*) "eri_value", eri_value(i)
          ! enddo
          buff(10) = (buff(8)*buff(3) - buff(9)*buff(2))
          buff(11) = (buff(9)*buff(1) - buff(7)*buff(3))
          buff(12) = (buff(7)*buff(2) - buff(8)*buff(1))

          trans(1) = eri_value(1)
          trans(2) = eri_value(7)
          trans(3) = eri_value(13)
          trans(4) = eri_value(19)
          trans(5) = eri_value(25)
          trans(6) = eri_value(31)
  eri_value(1)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(7)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(13)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(19)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(25)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(31)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(2)
          trans(2) = eri_value(8)
          trans(3) = eri_value(14)
          trans(4) = eri_value(20)
          trans(5) = eri_value(26)
          trans(6) = eri_value(32)
  eri_value(2)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(8)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(14)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(20)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(26)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(32)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
          trans(1) = eri_value(3)
          trans(2) = eri_value(9)
          trans(3) = eri_value(15)
          trans(4) = eri_value(21)
          trans(5) = eri_value(27)
          trans(6) = eri_value(33)
  eri_value(3)=buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(9)=buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(15)=buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(21)=buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(27)=buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(33)=buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)
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

          ! do i = 1, 36
          !   write(*,*) "eri_value", eri_value(i)
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
          do i = mini, maxi !1/2,4
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
                ijkp = (k - 1)*0 + ijp !gpople index 4

                do l = minl, maxl2
                  ijklp = (l - 1)*1 + ijkp !gpople index 1

                  buff(1) = eri_value(ijklp)

                  if (abs(buff(1)) .lt. 5.0d-11) cycle

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

      end do
    end do !main loop
    !$omp end target teams distribute parallel do

    ! do i = 1, size_of_matrix
    !         if(fock(i).ne.0.0_dp) write(*,*) "fa final", fock(i)
    !         ! if(fock(i).eq.12.73268741050258) write(*,*) "fa final", i, fock(i); i= 318003
    ! end do
  end subroutine int0202
end submodule
