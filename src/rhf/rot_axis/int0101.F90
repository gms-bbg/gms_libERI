submodule(rot_axis_kernels) int0101_impl
contains
  module subroutine int0101(sp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n01bra(:)
    real(dp), allocatable :: xint01bra(:)
    integer(kind=int64) :: nspbra
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ijij, klkl
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: r12, r34, buff(9), cosg, sing, rcd, rab, acx, acy, acz, tmp, tim1, tim2
    real(dp) :: scutspbra
    real(dp) :: d01p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) ::  zero_m_0(1), zero_m_1(4), zero_m_2(3), boys0, boys1, boys2
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    !r-integrals
    real(dp) :: tx21, y03, xmdt
    real(dp) :: r00(1), r01(6), r02(6)
    real(dp) :: eri_value(9), qx, qz, t1, t2, t3
    !digestion
    integer(kind=int64) :: ij_index, ijk_index, ijkl_index
    real(dp) :: ghondo(75) !for gamess, for now, contains L-shells
    !digestion
    logical :: same, iandj, kandl
    integer(kind=int64) :: loci, locj, lock, locl, nij, maxj2, ii1, ip, ijp, ijkp, ijklp, maxl2
    integer(kind=int64) :: jj1, i2, j2, nkl, kk1, itmp, ll1, k2, l2, ii2, jj2, kk2, ll2, ik, il, jk, jl
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: istart, iend, itile, ntile, ijkl_collapsed
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp, nquart, iijj, kkll
    integer :: new_integer

    mini = 1
    minj = 1
    mink = 1
    minl = 1

    maxi = 1
    maxj = 3
    maxk = 1
    maxl = 3
    allocate (n01bra(res%n_s_shl*res%n_p_shl))
    allocate (xint01bra(res%n_s_shl*res%n_p_shl))

    !start screening
    tim1 = omp_get_wtime()
    scutspbra = cutoff_schwarz/maxval(sp_pair%xints)
    nspbra = 0
    do ij = 1, res%n_s_shl*res%n_p_shl
    if (sp_pair%xints(ij) .ge. scutspbra) then
      nspbra = nspbra + 1
      xint01bra(nspbra) = sp_pair%xints(ij)
      n01bra(nspbra) = ij
    end if
    end do
    tim2 = omp_get_wtime()

    !--multi-GPU--work
    nchunk = nspbra/res%n_size
    nquart_start = nchunk*res%n_rank + 1
    nquart_end = nquart_start + nchunk - 1
    if (res%n_rank .EQ. res%n_size - 1) nquart_end = nspbra

 !$omp target teams distribute parallel do collapse(2) default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, sp_pair) &
 !$omp shared(mini, minj, mink, minl, maxi, maxj, maxk, maxl) &
 !$omp shared(nquart_start, nquart_end) &
 !$omp shared(xint01bra, nspbra, n01bra) &
 !$omp private(test,new_integer,ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d01p, tx21, t_beta, t_alpha, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, r00, r01, r02, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, t_inverse, t_new) &
 !$omp private(xmdt, qx, qz, eri_value, ghondo, same, maxl2, ijk_index, nij, nkl, ijkl_index, ijklp) &
 !$omp private(t1, t2, t3, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l) &
 !$omp private(iandj, kandl, loci, locj, lock, locl, maxj2, ip, ijkp) &
 !$omp private(ij_index, ijp, itmp)
    do iijj = nquart_start, nquart_end
      do kkll = 1, nspbra
        if (kkll .gt. iijj) cycle !needed to add this
        test = xint01bra(iijj)*xint01bra(kkll)
        if (test .lt. cutoff_schwarz) cycle
        ij = n01bra(iijj)
        kl = n01bra(kkll)

        ish_tmp = (ij - 1)/res%n_p_shl + 1
        jsh_tmp = mod(ij - 1, res%n_p_shl) + 1
        ksh_tmp = (kl - 1)/res%n_p_shl + 1
        lsh_tmp = mod(kl - 1, res%n_p_shl) + 1

        ish = res%i_s_shl(ish_tmp)
        jsh = res%i_p_shl(jsh_tmp)
        ksh = res%i_s_shl(ksh_tmp)
        lsh = res%i_p_shl(lsh_tmp)

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
        acy2 = acy*acy !acy2
        r00 = 0.0_dp
        r01 = 0.0_dp
        r02 = 0.0_dp
        ket_loop = 0
        do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
          ket_loop = ket_loop + 1
          t_expon_cd = sp_pair%t_expon_ab(sp_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_sp
          t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
          y03 = t_inverse_expon_cd*sp_pair%expon_a(sp_pair%pair_loc(kl) + ket_loop) !y03
          sq = sp_pair%sq(sp_pair%pair_loc(kl) + ket_loop)
          cq = sp_pair%t_alpha(sp_pair%pair_loc(kl) + ket_loop) !cq
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
          fmt = 0.0_dp
          bra_loop = 0
          do i = 1, res%contr_num(ish)*res%contr_num(jsh)
            bra_loop = bra_loop + 1
            !get bra shell pair info
            new_integer = sp_pair%ismlp(sp_pair%pair_loc(ij) + bra_loop) + sp_pair%ismlp(sp_pair%pair_loc(kl) + ket_loop)
            if (new_integer .ge. 2) cycle
            t_expon_ab = sp_pair%t_expon_ab(sp_pair%pair_loc(ij) + bra_loop) !tx12
            t_beta = sp_pair%t_beta(sp_pair%pair_loc(ij) + bra_loop) !ty01
            t_alpha = sp_pair%t_alpha(sp_pair%pair_loc(ij) + bra_loop) !ty02
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
              t_new = t*1.6109850374431357d+01 !f_increment(3)
              t_int = nint(t_new)
              fmt = boys_grid_zero((2*451*5) + (t_int*5) + 5)*t_new
              fmt = (fmt + boys_grid_zero((2*451*5) + (t_int*5) + 4))*t_new
              fmt = (fmt + boys_grid_zero((2*451*5) + (t_int*5) + 3))*t_new
              fmt = (fmt + boys_grid_zero((2*451*5) + (t_int*5) + 2))*t_new
              fmt = fmt + boys_grid_zero((2*451*5) + (t_int*5) + 1)
              !use extrapolation for exp(-T)
              t_new = t*m_increment !rxinc
              t_int = nint(t_new)
              expt = exponent_grid((t_int*5) + 5)*t_new
              expt = (expt + exponent_grid((t_int*5) + 4))*t_new
              expt = (expt + exponent_grid((t_int*5) + 3))*t_new
              expt = (expt + exponent_grid((t_int*5) + 2))*t_new
              expt = expt + exponent_grid((t_int*5) + 1)
              boys1 = ((t + t)*fmt + expt)*0.33333333333333333d00 !*rmr(2)
              boys0 = ((t + t)*boys1 + expt)!*rmr(1)
              boys0 = boys0*sqrt(expon_abcd_inverse)*d01p
              boys1 = boys1*sqrt(expon_abcd_inverse)*d01p*rho
              boys2 = fmt*sqrt(expon_abcd_inverse)*d01p*rho*rho
            else
              !when t is large, use special formula
              t_inverse = 1.0_dp/t
              boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d01p
              boys1 = boys0*t_inverse*0.5_dp*rho
              boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
            end if
            !zero_m (ss|ss) fundamental integrals generation
            zero_m_0(1) = zero_m_0(1) + boys0*t_beta

            zero_m_1(1) = zero_m_1(1) + boys1*tx21
            zero_m_1(2) = zero_m_1(2) + boys1*tx21*pqr
            zero_m_1(3) = zero_m_1(3) + boys1*t_beta
            zero_m_1(4) = zero_m_1(4) + boys1*t_beta*pqr

            zero_m_2(1) = zero_m_2(1) + boys2*tx21
            zero_m_2(2) = zero_m_2(2) + boys2*tx21*pqr
            zero_m_2(3) = zero_m_2(3) + boys2*tx21*pqs
          end do
          !r integrals here
          r00(1) = r00(1) - zero_m_0(1)*y03*sq

          xmdt = t_inverse_expon_cd*0.5_dp*sq
          r01(1) = r01(1) - zero_m_1(3)*aqx*xmdt
          r01(2) = r01(2) - zero_m_1(3)*acy*xmdt
          r01(3) = r01(3) + zero_m_1(4)*xmdt
          r01(4) = r01(4) + zero_m_1(1)*aqx*y03*sq
          r01(5) = r01(5) + zero_m_1(1)*acy*y03*sq
          r01(6) = r01(6) - zero_m_1(2)*y03*sq

          r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*xmdt
          r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*xmdt
          r02(3) = r02(3) + (zero_m_2(3) - zero_m_1(1))*xmdt
          r02(4) = r02(4) + zero_m_2(1)*aqxy*xmdt
          r02(5) = r02(5) - zero_m_2(2)*aqx*xmdt
          r02(6) = r02(6) - zero_m_2(2)*acy*xmdt

        end do !end ket loop
        qx = rcd*sing
        qz = rcd*cosg
        eri_value(1) = -r02(1) - r01(4)*qx
        eri_value(2) = -r02(4)
        eri_value(3) = -r02(5) - r01(4)*qz
        eri_value(4) = -r02(4) - r01(5)*qx
        eri_value(5) = -r02(2)
        eri_value(6) = -r02(6) - r01(5)*qz
        eri_value(7) = -r02(5) - r01(6)*qx + r01(1) + r00(1)*qx
        eri_value(8) = -r02(6) + r01(2)
        eri_value(9) = -r02(3) - r01(6)*qz + r01(3) + r00(1)*qz

        t1 = eri_value(1)
        t2 = eri_value(4)
        t3 = eri_value(7)
        eri_value(1) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(4) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(7) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)
        t1 = eri_value(2)
        t2 = eri_value(5)
        t3 = eri_value(8)
        eri_value(2) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(5) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(8) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)
        t1 = eri_value(3)
        t2 = eri_value(6)
        t3 = eri_value(9)
        eri_value(3) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(6) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(9) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)

        t1 = eri_value(1)
        t2 = eri_value(2)
        t3 = eri_value(3)
        eri_value(1) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(2) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(3) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)
        t1 = eri_value(4)
        t2 = eri_value(5)
        t3 = eri_value(6)
        eri_value(4) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(5) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(6) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)
        t1 = eri_value(7)
        t2 = eri_value(8)
        t3 = eri_value(9)
        eri_value(7) = t1*(buff(8)*buff(3) - buff(9)*buff(2)) + t2*buff(7) + t3*buff(1)
        eri_value(8) = t1*(buff(9)*buff(1) - buff(7)*buff(3)) + t2*buff(8) + t3*buff(2)
        eri_value(9) = t1*(buff(7)*buff(2) - buff(8)*buff(1)) + t2*buff(9) + t3*buff(3)

        ghondo(1) = eri_value(1)
        ghondo(2) = eri_value(2)
        ghondo(3) = eri_value(3)
        ! ghondo(4:36) = 0.0 unrolling for DC code
        ghondo(4) = 0.0
        ghondo(5) = 0.0
        ghondo(6) = 0.0
        ghondo(7) = 0.0
        ghondo(8) = 0.0
        ghondo(9) = 0.0
        ghondo(10) = 0.0
        ghondo(11) = 0.0
        ghondo(12) = 0.0
        ghondo(13) = 0.0
        ghondo(14) = 0.0
        ghondo(15) = 0.0
        ghondo(16) = 0.0
        ghondo(17) = 0.0
        ghondo(18) = 0.0
        ghondo(19) = 0.0
        ghondo(20) = 0.0
        ghondo(21) = 0.0
        ghondo(22) = 0.0
        ghondo(23) = 0.0
        ghondo(24) = 0.0
        ghondo(25) = 0.0
        ghondo(26) = 0.0
        ghondo(27) = 0.0
        ghondo(28) = 0.0
        ghondo(29) = 0.0
        ghondo(30) = 0.0
        ghondo(31) = 0.0
        ghondo(32) = 0.0
        ghondo(33) = 0.0
        ghondo(34) = 0.0
        ghondo(35) = 0.0
        ghondo(36) = 0.0
        ghondo(37) = eri_value(4)
        ghondo(38) = eri_value(5)
        ghondo(39) = eri_value(6)
        ! ghondo(40:72) = 0.0 unrolling for DC code
        ghondo(40) = 0.0
        ghondo(41) = 0.0
        ghondo(42) = 0.0
        ghondo(43) = 0.0
        ghondo(44) = 0.0
        ghondo(45) = 0.0
        ghondo(46) = 0.0
        ghondo(47) = 0.0
        ghondo(48) = 0.0
        ghondo(49) = 0.0
        ghondo(50) = 0.0
        ghondo(51) = 0.0
        ghondo(52) = 0.0
        ghondo(53) = 0.0
        ghondo(54) = 0.0
        ghondo(55) = 0.0
        ghondo(56) = 0.0
        ghondo(57) = 0.0
        ghondo(58) = 0.0
        ghondo(59) = 0.0
        ghondo(60) = 0.0
        ghondo(61) = 0.0
        ghondo(62) = 0.0
        ghondo(63) = 0.0
        ghondo(64) = 0.0
        ghondo(65) = 0.0
        ghondo(66) = 0.0
        ghondo(67) = 0.0
        ghondo(68) = 0.0
        ghondo(69) = 0.0
        ghondo(70) = 0.0
        ghondo(71) = 0.0
        ghondo(72) = 0.0
        ghondo(73) = eri_value(7)
        ghondo(74) = eri_value(8)
        ghondo(75) = eri_value(9)

        same = (ish .eq. ksh .and. jsh .eq. lsh)
        iandj = (ish .eq. jsh)
        kandl = (ksh .eq. lsh)

        loci = res%atom_loc(ish) - mini
        locj = res%atom_loc(jsh) - minj
        lock = res%atom_loc(ksh) - mink
        locl = res%atom_loc(lsh) - minl

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

      end do !main loop
    end do
    !$omp end target teams distribute parallel do

    deallocate (n01bra)
    deallocate (xint01bra)
  end subroutine int0101
end submodule
