submodule(rot_axis_kernels) int0111_impl
contains
  module subroutine int0111(sp_pair, pp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sp_pair, pp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n01bra(:), n11ket(:)
    real(dp), allocatable :: xint01bra(:), xint11ket(:)
    integer(kind=int64) :: nspbra, nppket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutspbra, scutppket
    real(dp) :: r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) :: d01p, d11p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) ::  zero_m_0(1), zero_m_1(4), zero_m_2(6), zero_m_3(4), boys0, boys1, boys2, boys3
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    real(dp) :: tx21, y03
    !r_integrals
    real(dp) :: r00(2), r01(12), r02(18), r03(10)
    real(dp) :: xmdt
    !eri_values
    real(dp) :: eri_value(27), qx, qz, t1, t2, t3
    !digestion
    logical :: same, iandj, kandl
    integer(kind=int64) :: loci, locj, lock, locl, nij, maxj2, ii1, ip, ijp, ijkp, ijklp, maxl2
    integer(kind=int64) :: jj1, i2, j2, nkl, kk1, itmp, ll1, k2, l2, ii2, jj2, kk2, ll2, ik, il, jk, jl
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: nchunksize_int64
    logical(kind=int64), allocatable :: lscreen(:)
    integer(kind=int64) :: istart, iend, itile, ntile, ijkl_collapsed
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp
    real(dp) :: ghondo(90)
    real(dp) :: kernel_full1, kernel_full2, kernel_only1, kernel_only2, first_screen1, first_screen2
    integer :: shp_thresh
    real(dp) :: five, ten, t12, nchunksize_dp
    basd = 1
    minl = 1
    maxl = 1
    basc = 3
    mink = 1
    maxk = 3
    basb = 9
    basa = 0
    mini = 1
    maxi = 1
    allocate (n01bra(res%n_s_shl*res%n_p_shl))
    allocate (xint01bra(res%n_s_shl*res%n_p_shl))
    allocate (n11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    !start screening
    first_screen1 = omp_get_wtime()
    scutspbra = cutoff_schwarz/maxval(pp_pair%xints)
    nspbra = 0
    do ij = 1, res%n_s_shl*res%n_p_shl
    if (sp_pair%xints(ij) .ge. scutspbra) then
      nspbra = nspbra + 1
      xint01bra(nspbra) = sp_pair%xints(ij)
      n01bra(nspbra) = ij
    end if
    end do
    scutppket = cutoff_schwarz/maxval(sp_pair%xints)
    nppket = 0
    do kl = 1, res%n_p_shl*(res%n_p_shl + 1)/2
    if (pp_pair%xints(kl) .ge. scutppket) then
      nppket = nppket + 1
      xint11ket(nppket) = pp_pair%xints(kl)
      n11ket(nppket) = kl
    end if
    end do
    first_screen2 = omp_get_wtime()

    !Total available memory in 32 GB V100s is 32000000000 bytes
    !Leaves us with 4000000000 bytes to account for 8-byte integrals
    ten = 10d00
    nchunksize_dp = (375000000d00*ten)
    nchunksize_int64 = int(nchunksize_dp, int64) !for 32 GB V100s, 8-byte integers

    if ((nspbra*nppket) .le. nchunksize_int64) nchunksize_int64 = nspbra*nppket
    ntile = int(nspbra*nppket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nspbra*nppket

      !--multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, pp_pair, sp_pair) &
 !$omp shared(nquart_start, nquart_end, xint01bra, xint11ket) &
 !$omp shared(nppket, n01bra, n11ket) &
 !$omp private(shp_thresh,test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d01p, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt, t_alpha, tx21) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, zero_m_3, r00, r01, r02, r03, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, boys3, t_inverse, t_new ) &
 !$omp private(xmdt, qx, qz, eri_value, ghondo, iandj, kandl, same, maxl2, nij, nkl, ip, ijp, ijkp, ijklp) &
 !$omp private(mini, minj, mink, minl, loci, locj, lock, locl, maxi, maxj, maxk, maxl, maxj2, itmp) &
 !$omp private(t1, t2, t3, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nppket + 1
        kl_tmp = mod(iquart - 1, nppket) + 1

        test = xint01bra(ij_tmp)*xint11ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then
          ij = n01bra(ij_tmp)
          kl = n11ket(kl_tmp)
          ish_tmp = (ij - 1)/res%n_p_shl + 1
          jsh_tmp = mod(ij - 1, res%n_p_shl) + 1
          ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
          lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2
          ii = res%i_s_shl(ish_tmp)
          jj = res%i_p_shl(jsh_tmp)
          kk = res%i_p_shl(ksh_tmp)
          ll = res%i_p_shl(lsh_tmp)

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

          acy2 = acy*acy
          r00 = 0.0_dp
          r01 = 0.0_dp

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

          r03 = 0.0_dp
          ket_loop = 0
          do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1
            t_expon_cd = pp_pair%t_expon_ab(pp_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_pp
            t_inverse_expon_cd = 1.0_dp/t_expon_cd ! x43
            y03 = t_inverse_expon_cd*pp_pair%expon_a(pp_pair%pair_loc(kl) + ket_loop) !y03
            sq = pp_pair%sq(pp_pair%pair_loc(kl) + ket_loop)
            cq = pp_pair%t_alpha(pp_pair%pair_loc(kl) + ket_loop)
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
              shp_thresh = sp_pair%ismlp(sp_pair%pair_loc(ij) + bra_loop) + pp_pair%ismlp(pp_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              !get bra shell pair info
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
              zero_m_1(2) = zero_m_1(2) + boys1*t_beta*pqr
              zero_m_1(3) = zero_m_1(3) + boys1*tx21
              zero_m_1(4) = zero_m_1(4) + boys1*tx21*pqr

              zero_m_2(1) = zero_m_2(1) + boys2*t_beta
              zero_m_2(2) = zero_m_2(2) + boys2*t_beta*pqr
              zero_m_2(3) = zero_m_2(3) + boys2*t_beta*pqs
              zero_m_2(4) = zero_m_2(4) + boys2*tx21
              zero_m_2(5) = zero_m_2(5) + boys2*tx21*pqr
              zero_m_2(6) = zero_m_2(6) + boys2*tx21*pqs

              zero_m_3(1) = zero_m_3(1) + boys3*tx21
              zero_m_3(2) = zero_m_3(2) + boys3*tx21*pqr
              zero_m_3(3) = zero_m_3(3) + boys3*tx21*pqs
              zero_m_3(4) = zero_m_3(4) + boys3*tx21*pqr*pqs
            end do
            !r integrals here
            xmdt = t_inverse_expon_cd*0.5_dp*sq
            r00(1) = r00(1) - zero_m_0(1)*y03*sq*(1 - y03)
            r00(2) = r00(2) + zero_m_0(1)*xmdt

            r01(1) = r01(1) + zero_m_1(1)*aqx*xmdt*y03
            r01(2) = r01(2) + zero_m_1(1)*acy*xmdt*y03
            r01(3) = r01(3) - zero_m_1(2)*xmdt*y03
            r01(4) = r01(4) - zero_m_1(1)*aqx*xmdt*(1 - y03)
            r01(5) = r01(5) - zero_m_1(1)*acy*xmdt*(1 - y03)
            r01(6) = r01(6) + zero_m_1(2)*xmdt*(1 - y03)
            r01(7) = r01(7) + zero_m_1(3)*aqx*y03*sq*(1 - y03)
            r01(8) = r01(8) + zero_m_1(3)*acy*y03*sq*(1 - y03)
            r01(9) = r01(9) - zero_m_1(4)*y03*sq*(1 - y03)
            r01(10) = r01(10) - zero_m_1(3)*aqx*xmdt
            r01(11) = r01(11) - zero_m_1(3)*acy*xmdt
            r01(12) = r01(12) + zero_m_1(4)*xmdt

            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*xmdt
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*xmdt
            r02(3) = r02(3) + (zero_m_2(3) - zero_m_1(1))*xmdt
            r02(4) = r02(4) + zero_m_2(1)*aqxy*xmdt
            r02(5) = r02(5) - zero_m_2(2)*aqx*xmdt
            r02(6) = r02(6) - zero_m_2(2)*acy*xmdt

            xmdt = t_inverse_expon_cd*sq*0.5_dp
            r02(7) = r02(7) - (zero_m_2(4)*aqx2 - zero_m_1(3))*y03*xmdt
            r02(8) = r02(8) - (zero_m_2(4)*acy2 - zero_m_1(3))*y03*xmdt
            r02(9) = r02(9) - (zero_m_2(6) - zero_m_1(3))*y03*xmdt
            r02(10) = r02(10) - zero_m_2(4)*aqxy*y03*xmdt
            r02(11) = r02(11) + zero_m_2(5)*aqx*y03*xmdt
            r02(12) = r02(12) + zero_m_2(5)*acy*y03*xmdt

            r02(13) = r02(13) + (zero_m_2(4)*aqx2 - zero_m_1(3))*xmdt*(1 - y03)
            r02(14) = r02(14) + (zero_m_2(4)*acy2 - zero_m_1(3))*xmdt*(1 - y03)
            r02(15) = r02(15) + (zero_m_2(6) - zero_m_1(3))*xmdt*(1 - y03)
            r02(16) = r02(16) + zero_m_2(4)*aqxy*xmdt*(1 - y03)
            r02(17) = r02(17) - zero_m_2(5)*aqx*xmdt*(1 - y03)
            r02(18) = r02(18) - zero_m_2(5)*acy*xmdt*(1 - y03)

            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r03(1) = r03(1) - (zero_m_3(1)*aqx2 - zero_m_2(4)*3.0_dp)*xmdt*aqx
            r03(2) = r03(2) - (zero_m_3(1)*aqx2 - zero_m_2(4))*xmdt*acy
            r03(3) = r03(3) + (zero_m_3(2)*aqx2 - zero_m_2(5))*xmdt
            r03(4) = r03(4) - (zero_m_3(1)*acy2 - zero_m_2(4))*xmdt*aqx
            r03(5) = r03(5) + zero_m_3(2)*xmdt*aqxy
            r03(6) = r03(6) - (zero_m_3(3) - zero_m_2(4))*xmdt*aqx
            r03(7) = r03(7) - (zero_m_3(1)*acy2 - zero_m_2(4)*3.0_dp)*xmdt*acy
            r03(8) = r03(8) + (zero_m_3(2)*acy2 - zero_m_2(5))*xmdt
            r03(9) = r03(9) - (zero_m_3(3) - zero_m_2(4))*xmdt*acy
            r03(10) = r03(10) + (zero_m_3(4) - zero_m_2(5)*3.0_dp)*xmdt
          end do !end ket loop

          qx = rcd*sing
          qz = rcd*cosg

          eri_value(1) = -r01(10) + qx*(-(qx*r01(7)) - r02(7) - r02(13)) - r03(1)
          eri_value(2) = -(qx*r02(16)) - r03(2)
          eri_value(3) = qz*(-(qx*r01(7)) - r02(7)) - qx*r02(17) - r03(3)
          eri_value(4) = -(qx*r02(10)) - r03(2)

          eri_value(5) = -r01(10) - r03(4)
          eri_value(6) = -(qz*r02(10)) - r03(5)
          eri_value(7) = -(qz*r02(13)) + qx*(-(qz*r01(7)) - r02(11)) - r03(3)
          eri_value(8) = -(qz*r02(16)) - r03(5)
          eri_value(9) = -r01(10) + qz*(-(qz*r01(7)) - r02(11) - r02(17)) - r03(6)

          eri_value(10) = -r01(11) + qx*(-(qx*r01(8)) - r02(10) - r02(16)) - r03(2)
          eri_value(11) = -(qx*r02(14)) - r03(4)
          eri_value(12) = qz*(-(qx*r01(8)) - r02(10)) - qx*r02(18) - r03(5)
          eri_value(13) = -(qx*r02(8)) - r03(4)
          eri_value(14) = -r01(11) - r03(7)
          eri_value(15) = -(qz*r02(8)) - r03(8)
          eri_value(16) = -(qz*r02(16)) + qx*(-(qz*r01(8)) - r02(12)) - r03(5)
          eri_value(17) = -(qz*r02(14)) - r03(8)
          eri_value(18) = -r01(11) + qz*(-(qz*r01(8)) - r02(12) - r02(18)) - r03(9)

          eri_value(19) = r00(2) - r01(12) + r02(1) + qx*(r01(1) + r01(4) + qx*(r00(1) - r01(9)) - r02(11) - r02(17)) - r03(3)
          eri_value(20) = r02(4) + qx*(r01(5) - r02(18)) - r03(5)
          eri_value(21) = qx*(r01(6) - r02(15)) + r02(5) + qz*(r01(1) + qx*(r00(1) - r01(9)) - r02(11)) - r03(6)
          eri_value(22) = r02(4) + qx*(r01(2) - r02(12)) - r03(5)
          eri_value(23) = r00(2) - r01(12) + r02(2) - r03(8)
          eri_value(24) = r02(6) + qz*(r01(2) - r02(12)) - r03(9)
          eri_value(25) = qx*(r01(3) + qz*(r00(1) - r01(9)) - r02(9)) + r02(5) + qz*(r01(4) - r02(17)) - r03(6)
          eri_value(26) = r02(6) + qz*(r01(5) - r02(18)) - r03(9)
          eri_value(27) = r00(2) - r01(12) + r02(3) + qz*(r01(3) + r01(6) + qz*(r00(1) - r01(9)) - r02(9) - r02(15)) - r03(10)

          buff(10) = (buff(8)*buff(3) - buff(9)*buff(2))
          buff(11) = (buff(9)*buff(1) - buff(7)*buff(3))
          buff(12) = (buff(7)*buff(2) - buff(8)*buff(1))
          !----first----
          t1 = eri_value(1) !1,1,1,0
          t2 = eri_value(10) !1,1,2,0
          t3 = eri_value(19) !1,1,3,0
          eri_value(1) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(10) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(19) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(2)  !2,1,1,0
          t2 = eri_value(11) !2,1,2,0
          t3 = eri_value(20) !2,1,3,0
          eri_value(2) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(11) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(20) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(3)  !3,1,1,0
          t2 = eri_value(12) !3,1,1,0
          t3 = eri_value(21) !3,1,1,0
          eri_value(3) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(12) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(21) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(4) !1,2,1,0
          t2 = eri_value(13)
          t3 = eri_value(22)
          eri_value(4) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(13) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(22) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(5)
          t2 = eri_value(14)
          t3 = eri_value(23)
          eri_value(5) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(14) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(23) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(6)
          t2 = eri_value(15)
          t3 = eri_value(24)
          eri_value(6) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(15) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(24) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(7)
          t2 = eri_value(16)
          t3 = eri_value(25)
          eri_value(7) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(16) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(25) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(8)
          t2 = eri_value(17)
          t3 = eri_value(26)
          eri_value(8) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(17) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(26) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(9)
          t2 = eri_value(18)
          t3 = eri_value(27)
          eri_value(9) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(18) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(27) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          !-----second------
          t1 = eri_value(1) !F(1,1,1,0)
          t2 = eri_value(4) !F(1,2,1,0)
          t3 = eri_value(7) !F(1,3,1,0)
          eri_value(1) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(4) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(7) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(2)  !F(2,1,1,0)
          t2 = eri_value(5) !F(2,2,1,0)
          t3 = eri_value(8) !F(2,3,1,0)
          eri_value(2) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(5) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(8) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(3) !F(3,1,1,0)
          t2 = eri_value(6)  !F(3,2,1,0)
          t3 = eri_value(9)  !F(3,3,1,0)
          eri_value(3) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(6) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(9) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(10) !F(1,1,2,0)
          t2 = eri_value(13) !F(1,2,2,0)
          t3 = eri_value(16) !F(1,3,2,0)
          eri_value(10) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(13) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(16) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(11) !F(2,1,2,0)
          t2 = eri_value(14) !F(2,2,2,0)
          t3 = eri_value(17) !F(2,3,2,0)
          eri_value(11) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(14) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(17) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(12)
          t2 = eri_value(15)
          t3 = eri_value(18)
          eri_value(12) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(15) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(18) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(19) !F(1,1,3,0)
          t2 = eri_value(22) !F(1,2,3,0)
          t3 = eri_value(25) !F(1,3,3,0)
          eri_value(19) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(22) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(25) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(20)  !F(2,1,3,0)
          t2 = eri_value(23) !F(2,2,3,0)
          t3 = eri_value(26) !F(2,3,3,0)
          eri_value(20) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(23) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(26) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(21) !F(3,1,3,0)
          t2 = eri_value(24)  !F(3,2,3,0)
          t3 = eri_value(27)  !F(3,3,3,0)
          eri_value(21) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(24) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(27) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          !-----third------
          t1 = eri_value(1) !F(1,1,1,0)
          t2 = eri_value(2) !F(2,1,1,0)
          t3 = eri_value(3) !F(3,1,1,0)
          eri_value(1) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(2) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(3) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(4)  !F(1,2,1,0)
          t2 = eri_value(5) !F(2,2,1,0)
          t3 = eri_value(6) !F(3,2,1,0)
          eri_value(4) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(5) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(6) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(7) !F(1,3,1,0)
          t2 = eri_value(8)  !F(1,3,1,0)
          t3 = eri_value(9)  !F(1,3,1,0)
          eri_value(7) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(8) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(9) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(10) !F(1,1,2,0)
          t2 = eri_value(11) !F(2,1,2,0)
          t3 = eri_value(12) !F(3,1,2,0)
          eri_value(10) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(11) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(12) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(13)  !F(1,2,2,0)
          t2 = eri_value(14) !F(2,2,2,0)
          t3 = eri_value(15) !F(3,2,2,0)
          eri_value(13) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(14) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(15) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(16) !F(1,3,2,0)
          t2 = eri_value(17)  !F(1,3,2,0)
          t3 = eri_value(18)  !F(1,3,2,0)
          eri_value(16) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(17) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(18) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          t1 = eri_value(19) !F(1,1,3,0)
          t2 = eri_value(20) !F(2,1,3,0)
          t3 = eri_value(21) !F(3,1,3,0)
          eri_value(19) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(20) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(21) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(22)  !F(1,2,3,0)
          t2 = eri_value(23) !F(2,2,3,0)
          t3 = eri_value(24) !F(3,2,3,0)
          eri_value(22) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(23) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(24) = t1*buff(12) + t2*buff(9) + t3*buff(3)
          t1 = eri_value(25) !F(1,3,3,0)
          t2 = eri_value(26)  !F(2,3,3,0)
          t3 = eri_value(27)  !F(3,3,3,0)
          eri_value(25) = t1*buff(10) + t2*buff(7) + t3*buff(1)
          eri_value(26) = t1*buff(11) + t2*buff(8) + t3*buff(2)
          eri_value(27) = t1*buff(12) + t2*buff(9) + t3*buff(3)

          ghondo(1) = eri_value(1)
          ghondo(2) = eri_value(2)
          ghondo(3) = eri_value(3)
          ghondo(4:6) = 0.0_dp
          ghondo(7) = eri_value(4)
          ghondo(8) = eri_value(5)
          ghondo(9) = eri_value(6)
          ghondo(10:12) = 0.0_dp
          ghondo(13) = eri_value(7)
          ghondo(14) = eri_value(8)
          ghondo(15) = eri_value(9)

          ghondo(16) = 0.0_dp
          ghondo(17) = 0.0_dp
          ghondo(18) = 0.0_dp
          ghondo(19) = 0.0_dp
          ghondo(20) = 0.0_dp
          ghondo(21) = 0.0_dp
          ghondo(22) = 0.0_dp
          ghondo(23) = 0.0_dp
          ghondo(24) = 0.0_dp
          ghondo(25) = 0.0_dp
          ghondo(26) = 0.0_dp
          ghondo(27) = 0.0_dp
          ghondo(28) = 0.0_dp
          ghondo(29) = 0.0_dp
          ghondo(30) = 0.0_dp
          ghondo(31) = 0.0_dp
          ghondo(32) = 0.0_dp
          ghondo(33) = 0.0_dp
          ghondo(34) = 0.0_dp
          ghondo(35) = 0.0_dp
          ghondo(36) = 0.0_dp
          ghondo(37) = eri_value(10)
          ghondo(38) = eri_value(11)
          ghondo(39) = eri_value(12)
          ghondo(40:42) = 0.0_dp
          ghondo(43) = eri_value(13)
          ghondo(44) = eri_value(14)
          ghondo(45) = eri_value(15)
          ghondo(46:48) = 0.0_dp
          ghondo(49) = eri_value(16)
          ghondo(50) = eri_value(17)
          ghondo(51) = eri_value(18)

          ghondo(52) = 0.0_dp
          ghondo(53) = 0.0_dp
          ghondo(54) = 0.0_dp
          ghondo(55) = 0.0_dp
          ghondo(56) = 0.0_dp
          ghondo(57) = 0.0_dp
          ghondo(58) = 0.0_dp
          ghondo(59) = 0.0_dp
          ghondo(60) = 0.0_dp
          ghondo(61) = 0.0_dp
          ghondo(62) = 0.0_dp
          ghondo(63) = 0.0_dp
          ghondo(64) = 0.0_dp
          ghondo(65) = 0.0_dp
          ghondo(66) = 0.0_dp
          ghondo(67) = 0.0_dp
          ghondo(68) = 0.0_dp
          ghondo(69) = 0.0_dp
          ghondo(70) = 0.0_dp
          ghondo(71) = 0.0_dp
          ghondo(72) = 0.0_dp
          ghondo(73) = eri_value(19)
          ghondo(74) = eri_value(20)
          ghondo(75) = eri_value(21)
          ghondo(76:78) = 0.0_dp
          ghondo(79) = eri_value(22)
          ghondo(80) = eri_value(23)
          ghondo(81) = eri_value(24)
          ghondo(82:84) = 0.0_dp
          ghondo(85) = eri_value(25)
          ghondo(86) = eri_value(26)
          ghondo(87) = eri_value(27)

          same = (ish .eq. ksh .and. jsh .eq. lsh)
          iandj = (ish .eq. jsh)
          kandl = (ksh .eq. lsh)

          mini = 1
          minj = 1
          mink = 1
          minl = 1
          loci = res%atom_loc(ish) - mini
          locj = res%atom_loc(jsh) - minj
          lock = res%atom_loc(ksh) - mink
          locl = res%atom_loc(lsh) - minl

          maxi = 1
          maxj = 3
          maxk = 3
          maxl = 3

          nij = 0
          maxj2 = maxj
          do i = mini, maxi !1/2,4
            if (iandj) maxj2 = i
            ii1 = i + loci
            ip = (i - 1)*216 + 1 !gpople index 64

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
                  ijklp = (l - 1) + ijkp !gpople index 1

                  buff(1) = ghondo(ijklp)

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

                end do
              end do
            end do
          end do

        end if
      end do !main loop
      !$omp end target teams distribute parallel do
      kernel_only2 = omp_get_wtime()
    end do !tiles
    kernel_full2 = omp_get_wtime()

    deallocate (n01bra)
    deallocate (xint01bra)
    deallocate (n11ket)
    deallocate (xint11ket)

  end subroutine int0111
end submodule
