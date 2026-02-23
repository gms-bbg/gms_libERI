submodule(rot_axis_kernels) int0002_impl
contains
  module subroutine int0002(ss_pair, sd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, sd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n00bra(:), n02ket(:)
    real(dp), allocatable :: xint00bra(:), xint02ket(:)
    integer(kind=int64) :: nssbra, nsdket
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk, minl, maxl
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, ijkl_collapsed
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end
    integer(kind=int64) ::  ish, jsh, ksh, lsh, i, j, k, l, bra_loop, ket_loop
    real(dp) :: scutssbra, scutsdket
    real(dp) ::  r12, r34, buff(12), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) ::  d00p, d02p, t_expon_ab, t_inverse_expon_ab, t_inverse_expon_cd, t_alpha, t_beta, t_alpha_cd, t_beta_cd
    real(dp) ::  sq, cq, cqx, cqz, aqx, aqx2, aqxy, aqz, qps, acy2
    !for boys function
    integer(kind=int64) :: t_int
    ! integer, intent(in) :: res%num_bas, res%ia(res%num_bas)
    real(dp) ::  zero_m_0(1), zero_m_1(2), zero_m_2(3), boys0, boys1, boys2
    real(dp) :: expon_abcd_inverse, t_expon_cd, pqr, pqs, fmt, t_new, t_inverse, rho, t, expt
    !r-integrals
    real(dp) :: r00(2), r01(3), r02(6), eri_value(6), trans(6)
    real(dp) :: xmdt, xmd2, y03
    !digestion
    integer(kind=int64) :: ii1, jj1, kk1, i2, j2, ll1, k2, l2, ii2, jj2, kk2, ik, il, jk, jl
    !multi-GPU
    real(dp) :: test
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
    basb = 0
    minj = 1
    maxj = 1
    basa = 0
    mini = 1
    maxi = 1
    allocate (n00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (n02ket(res%n_s_shl*res%n_d_shl))
    allocate (xint02ket(res%n_s_shl*res%n_d_shl))
    !start screening
    first_screen1 = omp_get_wtime()
    scutssbra = cutoff_schwarz/maxval(sd_pair%xints)
    nssbra = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
    if (ss_pair%xints(ij) .ge. scutssbra) then
      nssbra = nssbra + 1
      xint00bra(nssbra) = ss_pair%xints(ij)
      n00bra(nssbra) = ij
    end if
    end do
    scutsdket = cutoff_schwarz/maxval(ss_pair%xints)
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
    if ((nssbra*nsdket) .le. nchunksize_int64) nchunksize_int64 = nssbra*nsdket
    ntile = int(nssbra*nsdket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nssbra*nsdket
      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend
      kernel_only1 = omp_get_wtime()

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, boys_grid_zero, exponent_grid, sd_pair, ss_pair) &
 !$omp shared(nquart_start, nquart_end, xint00bra, xint02ket) &
 !$omp shared(nsdket, n00bra, n02ket) &
 !$omp private(shp_thresh, xmdt, xmd2,test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
 !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, cq, d00p, t_beta, aqz, t_inverse_expon_cd, expon_abcd_inverse) &
 !$omp private(bra_loop, ket_loop, t_expon_cd, t_expon_ab, y03, cqx, cqz, aqx, aqx2, aqxy, qps, sq, fmt, expt ) &
 !$omp private(zero_m_0, zero_m_1, zero_m_2, r00, r01, r02, pqr, pqs, rho, t, t_int, boys0, boys1, boys2, t_inverse, t_new ) &
 !$omp private(eri_value,trans, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsdket + 1
        kl_tmp = mod(iquart - 1, nsdket) + 1

        test = xint00bra(ij_tmp)*xint02ket(kl_tmp)
        if (test .gt. cutoff_schwarz) then
          ij = n00bra(ij_tmp)
          kl = n02ket(kl_tmp)
          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_d_shl + 1
          lsh_tmp = mod(kl - 1, res%n_d_shl) + 1

          ii = res%i_s_shl(ish_tmp)
          jj = res%i_s_shl(jsh_tmp)
          kk = res%i_s_shl(ksh_tmp)
          ll = res%i_d_shl(lsh_tmp)

          ish = ii
          jsh = jj
          ksh = kk
          lsh = ll
          ! write(*,*) "ish", ish
          ! write(*,*) "jsh", jsh
          ! write(*,*) "ksh", ksh
          ! write(*,*) "lsh", lsh
          ! if(ish.eq.11.and.jsh.eq.1.and.ksh.eq.7.and.lsh.eq.6) then
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
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              bra_loop = bra_loop + 1
              !get bra shell pair info
              shp_thresh = ss_pair%ismlp(ss_pair%pair_loc(ij) + bra_loop) + sd_pair%ismlp(sd_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
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
            ! write(*,*) "zero_m_0(1)", zero_m_0(1)
            ! write(*,*) "zero_m_1(1)", zero_m_1(1)
            ! write(*,*) "zero_m_1(2)", zero_m_1(2)
            ! write(*,*) "zero_m_2(1)", zero_m_2(1)
            ! write(*,*) "zero_m_2(2)", zero_m_2(2)
            ! write(*,*) "zero_m_2(3)", zero_m_2(3)

            xmdt = t_inverse_expon_cd*t_inverse_expon_cd*sq*0.25_dp
            xmd2 = -t_inverse_expon_cd*0.5_dp*sq*y03

            r00(1) = r00(1) + zero_m_0(1)*t_inverse_expon_cd*0.5_dp*sq
            r00(2) = r00(2) + zero_m_0(1)*y03*y03*sq
            ! write(*,*) "aqx", aqx
            ! write(*,*) "acy", acy
            r01(1) = r01(1) - zero_m_1(1)*aqx*xmd2
            r01(2) = r01(2) - zero_m_1(1)*acy*xmd2
            r01(3) = r01(3) + zero_m_1(2)*xmd2

            r02(1) = r02(1) + (zero_m_2(1)*aqx2 - zero_m_1(1))*xmdt
            r02(2) = r02(2) + (zero_m_2(1)*acy2 - zero_m_1(1))*xmdt
            r02(3) = r02(3) + (zero_m_2(3) - zero_m_1(1))*xmdt
            r02(4) = r02(4) + zero_m_2(1)*aqxy*xmdt
            r02(5) = r02(5) - zero_m_2(2)*aqx*xmdt
            r02(6) = r02(6) - zero_m_2(2)*acy*xmdt

            ! write(*,*) "r00(1)=",r00(1)
            ! write(*,*) "r00(2)=",r00(2)
            ! write(*,*) "r01(1)=",r01(1)
            ! write(*,*) "r01(2)=",r01(2)
            ! write(*,*) "r01(3)=",r01(3)
            ! write(*,*) "r02(1)=",r02(1)
            ! write(*,*) "r02(2)=",r02(2)
            ! write(*,*) "r02(3)=",r02(3)
            ! write(*,*) "r02(4)=",r02(4)
            ! write(*,*) "r02(5)=",r02(5)
            ! write(*,*) "r02(6)=",r02(6)

            !r integrals here
          end do !end ket loop
          ! write(*,*) "r00(1)=",r00(1)
          ! write(*,*) "r00(2)=",r00(2)
          ! write(*,*) "r01(1)=",r01(1)
          ! write(*,*) "r01(2)=",r01(2)
          ! write(*,*) "r01(3)=",r01(3)
          ! write(*,*) "r02(1)=",r02(1)
          ! write(*,*) "r02(2)=",r02(2)
          ! write(*,*) "r02(3)=",r02(3)
          ! write(*,*) "r02(4)=",r02(4)
          ! write(*,*) "r02(5)=",r02(5)
          ! write(*,*) "r02(6)=",r02(6)
          ! write(*,*) "rcd=",rcd
          ! write(*,*) "sing=",sing
          ! write(*,*) "cosg=",cosg
          ! eri_value(1) =  +r02(1)+r00(1)+(+r01(1)+r01(1)+r00(2)*rcd*sing)*rcd*sing
          ! eri_value(2) =  +r02(2)+r00(1)
          ! eri_value(3) =  +r02(3)+r00(1)+(+r01(3)+r01(3)+r00(2)*rcd*cosg)*rcd*cosg
          ! eri_value(4) =  +r02(4)+r01(2)*rcd*sing
          ! eri_value(5) =  +r02(5)+r01(3)*rcd*sing+(+r01(1)+r00(2)*rcd*sing)*rcd*cosg
          ! eri_value(6) =  +r02(6)+r01(2)*rcd*cosg

          eri_value(1) = r00(1) + rcd*sing*(rcd*sing*r00(2) + 2*r01(1)) + r02(1)
          eri_value(2) = r00(1) + r02(2)
          eri_value(3) = r00(1) + rcd*cosg*(rcd*cosg*r00(2) + 2*r01(3)) + r02(3)
          eri_value(4) = rcd*sing*r01(2) + r02(4)
          eri_value(5) = rcd*cosg*(rcd*sing*r00(2) + r01(1)) + rcd*sing*r01(3) + r02(5)
          eri_value(6) = rcd*cosg*r01(2) + r02(6)

          ! write(*,*) "eri_value(1) first=",eri_value(1)
          ! write(*,*) "eri_value(2) first=",eri_value(2)
          ! write(*,*) "eri_value(3) first=",eri_value(3)
          ! write(*,*) "eri_value(4) first=",eri_value(4)
          ! write(*,*) "eri_value(5) first=",eri_value(5)
          ! write(*,*) "eri_value(6) first=",eri_value(6)

          buff(10) = buff(8)*buff(3) - buff(9)*buff(2)
          buff(11) = buff(9)*buff(1) - buff(7)*buff(3)
          buff(12) = buff(7)*buff(2) - buff(8)*buff(1)
          trans(1) = eri_value(1)
          trans(2) = eri_value(2)
          trans(3) = eri_value(3)
          trans(4) = eri_value(4)
          trans(5) = eri_value(5)
          trans(6) = eri_value(6)

  eri_value(1) = buff(10)*buff(10)*trans(1) + buff(7)*buff(7)*trans(2) + buff(1)*buff(1)*trans(3) + buff(10)*buff(7)*2.0_dp*trans(4) + buff(1)*buff(10)*2.0_dp*trans(5) + buff(1)*buff(7)*2.0_dp*trans(6)
  eri_value(2) = buff(11)*buff(11)*trans(1) + buff(8)*buff(8)*trans(2) + buff(2)*buff(2)*trans(3) + buff(11)*buff(8)*2.0_dp*trans(4) + buff(11)*buff(2)*2.0_dp*trans(5) + buff(2)*buff(8)*2.0_dp*trans(6)
  eri_value(3) = buff(12)*buff(12)*trans(1) + buff(9)*buff(9)*trans(2) + buff(3)*buff(3)*trans(3) + buff(12)*buff(9)*2.0_dp*trans(4) + buff(12)*buff(3)*2.0_dp*trans(5) + buff(3)*buff(9)*2.0_dp*trans(6)
  eri_value(4) = buff(10)*buff(11)*sqrt3*trans(1) + buff(7)*buff(8)*sqrt3*trans(2) + buff(1)*buff(2)*sqrt3*trans(3) + (buff(11)*buff(7) + buff(10)*buff(8))*sqrt3*trans(4) + (buff(1)*buff(11) + buff(10)*buff(2))*sqrt3*trans(5) + (buff(2)*buff(7) + buff(1)*buff(8))*sqrt3*trans(6)
  eri_value(5) = buff(10)*buff(12)*sqrt3*trans(1) + buff(7)*buff(9)*sqrt3*trans(2) + buff(1)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(7) + buff(10)*buff(9))*sqrt3*trans(4) + (buff(1)*buff(12) + buff(10)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(7) + buff(1)*buff(9))*sqrt3*trans(6)
  eri_value(6) = buff(11)*buff(12)*sqrt3*trans(1) + buff(8)*buff(9)*sqrt3*trans(2) + buff(2)*buff(3)*sqrt3*trans(3) + (buff(12)*buff(8) + buff(11)*buff(9))*sqrt3*trans(4) + (buff(12)*buff(2) + buff(11)*buff(3))*sqrt3*trans(5) + (buff(3)*buff(8) + buff(2)*buff(9))*sqrt3*trans(6)

          ! write(*,*) "eri_value(1)=",eri_value(1)
          ! write(*,*) "eri_value(2)=",eri_value(2)
          ! write(*,*) "eri_value(3)=",eri_value(3)
          ! write(*,*) "eri_value(4)=",eri_value(4)
          ! write(*,*) "eri_value(5)=",eri_value(5)
          ! write(*,*) "eri_value(6)=",eri_value(6)
          ! endif
          ii1 = res%atom_loc(ish)
          jj1 = res%atom_loc(jsh)
          kk1 = res%atom_loc(ksh)
          ! write(*,*) "ii1", ii1
          ! write(*,*) "jj1", jj1
          ! write(*,*) "kk1", kk1

          i2 = ii1
          j2 = jj1
          if (ii1 .lt. jj1) then ! sort <ij|
            i2 = jj1
            j2 = ii1
          end if

          do l = 1, 6

            buff(1) = eri_value(l)
            !  write(*,*) "buff(1)", buff(1)

            if (abs(buff(1)) .lt. 5.0d-11) cycle ! goto 300
            !  write(*,*) " val", buff(1)
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
            !  write(*,*) " ii ", ii
            !  write(*,*) " ii2 ", ii2
            jj2 = res%ia(jj)
            kk2 = res%ia(kk)
            ij = ii2 + jj
            !  write(*,*) " ij ", ij
            ik = ii2 + kk
            il = ii2 + ll
            jk = jj2 + kk
            jl = jj2 + ll
            kl = kk2 + ll
            if (jj .lt. kk) jk = kk2 + jj
            if (jj .lt. ll) jl = res%ia(ll) + jj

            !   account for identical permutations.

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
            buff(9) = -buff(2)*density(il) !1.5s until here, 4s total with atomics
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
        end if !screening
      end do !main loop
      !$omp end target teams distribute parallel do

    end do !tiles

    deallocate (n00bra, xint00bra, n02ket, xint02ket)
  end subroutine int0002
end submodule
