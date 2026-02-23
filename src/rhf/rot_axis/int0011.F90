submodule(rot_axis_kernels) int0011_impl
contains
  module subroutine int0011(ss_pair, pp_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, pp_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n00bra(:), n11ket(:)
    real(dp), allocatable :: xint00bra(:), xint11ket(:)
    integer(kind=int64) :: nssbra, nppket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop, bra, ket
    real(dp) :: scutssbra, scutppket
    real(dp) :: r12, r34, buff(9), cosg, sing, rcd, rab, acx, acy, acz, tmp, tim1, tim2
    real(dp) :: d00p, d11p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) :: sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    real(dp) ::  zero_m_0(1), zero_m_1(2), zero_m_2(3), boys0, boys1, boys2
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    !r-integrals
    real(dp) :: r00(2), r01(6), r02(6)
    real(dp) :: xmdt, y03
    !eri_value
    real(dp) :: eri_value(9), val, qx, qz, t1, t2, t3
    !digestion
    logical :: kandl, same
   integer(kind=int64) :: ii1,kk1,nij,maxl2,jj1,i2,j2,ijp,nkl,itmp,ijklp,ll2,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ll1,ijkp,ijk_index,ijkl_index
    real(dp) :: ghondo(16) !for gamess, for now, contains L-shells
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunk, nquart_start, nquart_end
    integer(kind=int64) :: nchunksize_int64
    integer(kind=int64) :: istart, iend, itile, ntile, ijkl_collapsed
    integer(kind=int64) :: istart_tmp, iend_tmp, nchunksize_tmp
    real(dp) :: kernel_full1, kernel_full2, kernel_only1, kernel_only2, first_screen1, first_screen2
    integer :: new_integer
    real(dp) :: five, ten, t12, nchunksize_dp
    basd = 3
    minl = 1
    maxl = 3
    basc = 3
    mink = 1
    maxk = 3
    basb = 0
    minj = 1
    maxj = 1
    basa = 0
    mini = 1
    maxi = 1
    allocate (n00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (n11ket(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11ket(res%n_p_shl*(res%n_p_shl + 1)/2))

    !start screening
    first_screen1 = omp_get_wtime()
    tim1 = omp_get_wtime()
    scutssbra = cutoff_schwarz/maxval(pp_pair%xints)
    nssbra = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
    if (ss_pair%xints(ij) .ge. scutssbra) then
      nssbra = nssbra + 1
      xint00bra(nssbra) = ss_pair%xints(ij)
      n00bra(nssbra) = ij
    end if
    end do
    scutppket = cutoff_schwarz/maxval(ss_pair%xints)
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

    if ((nssbra*nppket) .le. nchunksize_int64) nchunksize_int64 = nssbra*nppket
    ntile = int(nssbra*nppket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nssbra*nppket

      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, pp_pair, ss_pair) &
 !$omp shared(nquart_start, nquart_end, xint00bra, xint11ket) &
 !$omp shared(nppket, n00bra, n11ket) &
 !$omp private(new_integer,test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d00p, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, r00, r01, r02, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, t_inverse, t_new ) &
 !$omp private(xmdt, qx, qz, eri_value, ghondo, kandl, same, maxl2, ijk_index, nij, nkl, ijkp, ijkl_index, ijklp) &
 !$omp private(t1, t2, t3, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nppket + 1
        kl_tmp = mod(iquart - 1, nppket) + 1

        test = xint00bra(ij_tmp)*xint11ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then

          ij = n00bra(ij_tmp)
          kl = n11ket(kl_tmp)
          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
          lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

          ii = res%i_s_shl(ish_tmp)
          jj = res%i_s_shl(jsh_tmp)
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
          r02 = 0.0_dp
          ket_loop = 0
          do ket = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1
            t_expon_cd = pp_pair%t_expon_ab(pp_pair%pair_loc(kl) + ket_loop) ! buff(3) = expon_cd_pp
            t_inverse_expon_cd = 1.0_dp/t_expon_cd !added new, x43
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
            fmt = 0.0_dp
            bra_loop = 0
            do bra = 1, res%contr_num(ish)*res%contr_num(jsh)
              bra_loop = bra_loop + 1
              !get bra shell pair info
              new_integer = ss_pair%ismlp(ss_pair%pair_loc(ij) + bra_loop) + pp_pair%ismlp(pp_pair%pair_loc(kl) + ket_loop)
              if (new_integer .ge. 2) cycle
              t_expon_ab = ss_pair%t_expon_ab(ss_pair%pair_loc(ij) + bra_loop) !tx12
              t_beta = ss_pair%t_alpha(ss_pair%pair_loc(ij) + bra_loop) !ty02
              d00p = ss_pair%d_coeff(ss_pair%pair_loc(ij) + bra_loop) !d00p
              expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd) !x41
              pqr = t_beta - aqz
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
                boys0 = boys0*sqrt(expon_abcd_inverse)*d00p
                boys1 = boys1*sqrt(expon_abcd_inverse)*d00p*rho
                boys2 = fmt*sqrt(expon_abcd_inverse)*d00p*rho*rho
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d00p
                boys1 = boys0*t_inverse*0.5_dp*rho
                boys2 = boys1*(t_inverse*0.5_dp*rho + t_inverse*rho)
              end if
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0
              zero_m_1(1) = zero_m_1(1) + boys1
              zero_m_1(2) = zero_m_1(2) + boys1*pqr
              zero_m_2(1) = zero_m_2(1) + boys2
              zero_m_2(2) = zero_m_2(2) + boys2*pqr
              zero_m_2(3) = zero_m_2(3) + boys2*pqs
            end do
            !r integrals here
            xmdt = t_inverse_expon_cd*0.5_dp*sq
            r00(1) = r00(1) + zero_m_0(1)*xmdt
            r00(2) = r00(2) - zero_m_0(1)*y03*sq*(1 - y03)

            r01(1) = r01(1) + zero_m_1(1)*aqx*xmdt*y03
            r01(2) = r01(2) + zero_m_1(1)*acy*xmdt*y03
            r01(3) = r01(3) - zero_m_1(2)*xmdt*y03
            r01(4) = r01(4) - zero_m_1(1)*aqx*xmdt*(1 - y03)
            r01(5) = r01(5) - zero_m_1(1)*acy*xmdt*(1 - y03)
            r01(6) = r01(6) + zero_m_1(2)*xmdt*(1 - y03)

            xmdt = t_inverse_expon_cd*sq*t_inverse_expon_cd*0.25_dp
            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*xmdt
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*xmdt
            r02(3) = r02(3) + (zero_m_2(3) - zero_m_1(1))*xmdt
            r02(4) = r02(4) + zero_m_2(1)*aqxy*xmdt
            r02(5) = r02(5) - zero_m_2(2)*aqx*xmdt
            r02(6) = r02(6) - zero_m_2(2)*acy*xmdt
          end do !end ket loop

          qx = rcd*sing
          qz = rcd*cosg

          eri_value(1) = r02(1) + r00(1) + (+r01(1) + r01(4) + r00(2)*qx)*qx
          eri_value(2) = r02(4) + r01(5)*qx
          eri_value(3) = r02(5) + r01(6)*qx + (+r01(1) + r00(2)*qx)*qz
          eri_value(4) = r02(4) + r01(2)*qx
          eri_value(5) = r02(2) + r00(1)
          eri_value(6) = r02(6) + r01(2)*qz
          eri_value(7) = r02(5) + r01(4)*qz + (+r01(3) + r00(2)*qz)*qx
          eri_value(8) = r02(6) + r01(5)*qz
          eri_value(9) = r02(3) + r00(1) + (+r01(3) + r01(6) + r00(2)*qz)*qz

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

          ! write(*,*) "eri_value(1) = ", eri_value(1)
          ! write(*,*) "eri_value(2) = ", eri_value(2)
          ! write(*,*) "eri_value(3) = ", eri_value(3)
          ! write(*,*) "eri_value(4) = ", eri_value(4)
          ! write(*,*) "eri_value(5) = ", eri_value(5)
          ! write(*,*) "eri_value(6) = ", eri_value(6)
          ! write(*,*) "eri_value(7) = ", eri_value(7)
          ! write(*,*) "eri_value(8) = ", eri_value(8)
          ! write(*,*) "eri_value(9) = ", eri_value(9)

          ghondo(1) = eri_value(1)
          ghondo(2) = eri_value(2)
          ghondo(3) = eri_value(3)
          ghondo(4) = 0.0_dp
          ghondo(5) = 0.0_dp
          ghondo(6) = 0.0_dp
          ghondo(7) = eri_value(4)
          ghondo(8) = eri_value(5)
          ghondo(9) = eri_value(6)
          ghondo(10) = 0.0_dp
          ghondo(11) = 0.0_dp
          ghondo(12) = 0.0_dp
          ghondo(13) = eri_value(7)
          ghondo(14) = eri_value(8)
          ghondo(15) = eri_value(9)

          kandl = (ksh .eq. lsh)
          same = (ish .eq. ksh .and. jsh .eq. lsh)

          ii1 = res%atom_loc(ish)
          jj1 = res%atom_loc(jsh)
          ! write(*,*) "ii1", ii1
          ! write(*,*) "jj1", jj1

          maxl2 = 3
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

              buff(1) = ghondo(IJKL_INDEX)

              IJKL_INDEX = IJKL_INDEX + 1
              !    write(*,*) "IJKL_INDEX", IJKL_INDEX

              if (abs(buff(1)) .lt. cutoff_schwarz) cycle ! goto 300
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

        end if
      end do !main loop
      !$omp end target teams distribute parallel do
      kernel_only2 = omp_get_wtime()
    end do !tiles

    deallocate (n00bra)
    deallocate (xint00bra)
    deallocate (n11ket)
    deallocate (xint11ket)
  end subroutine int0011
end submodule
