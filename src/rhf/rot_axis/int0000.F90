submodule(rot_axis_kernels) int0000_impl
contains
  module subroutine int0000(ss_pair, density, fock, res)

    use liberi_types, only: int32, int64, dp, shell_pair_t, eri_resources_t
    use parameters, only: cutoff_schwarz, sqrt_pi_div_2
    ! this is terrible practice, FIX
    use boys, only: t_max, boys_grid_zero
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    !necessary variables for the integral class
    integer(kind=int64), allocatable :: n00bra(:)
    real(dp), allocatable :: xint00bra(:)
    integer(kind=int64) :: nssbra
    integer(kind=int64) :: mini, maxi, minj, maxj, mink, maxk
    integer(kind=int64) :: basa, basb, basc, basd
    integer(kind=int64) :: ij, kl, ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, k, l, bra_loop, ket_loop
    real(dp) ::  r12, r34, buff(9), cosg, sing, rcd, rab, acx, acy, acz, tmp
    real(dp) ::  scutssbra
    real(dp) ::  d00p, t_expon_ab 
    real(dp) ::  aqz, qps, acy2
    !for boys function
    integer(kind=int32) :: t_int
    real(dp) :: zero_m_0(1), boys0
    real(dp) :: expon_abcd_inverse, pqr, pqs, fmt, t_new, t_inverse, rho, t
    !r-integrals
    !digestion
    integer(kind=int64) :: ii1, jj1, kk1, ll1, i2, j2, k2, ii2, jj2, kk2, ik, il, jk, jl
    !multi-GPU
    real(dp) :: test
    integer(kind=int64) :: nchunksize_int64, nchunk, nquart_start, nquart_end
    integer(kind=int64) :: istart, iend, itile, ntile
   integer(kind=int64) :: l2, iii,jjj,kkk,lll,iijj, kkll
    integer(kind=int64) :: shp_thresh
   real(dp) :: ten, nchunksize_dp

    basd = 0
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

    !start screening
    scutssbra = cutoff_schwarz/maxval(ss_pair%xints)
    nssbra = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
    if (ss_pair%xints(ij) .ge. scutssbra) then
      nssbra = nssbra + 1
      xint00bra(nssbra) = ss_pair%xints(ij)
      n00bra(nssbra) = ij
    end if
    end do

    !this part is extremely important, because compilers don't always
    !handle longlong integers properly
    !this is only necessary for (ss|ss) kernels for now, because other kernels
    !will have much smaller numbers of quartets
    !of course, as GPUs get more memory, we will need to adjust these
    !Total available memory in 32 GB V100s is 32000000000 bytes
    !Leaves us with 4000000000 bytes to account for 8-byte integrals
    !This is (ss|ss) type kernel, so we need to take the square root
    ten = 10d00
    nchunksize_dp = sqrt(375000000d00*ten)
    nchunksize_int64 = int(nchunksize_dp, int64) !for 32 GB V100s, 8-byte integers

    if ((nssbra) .le. nchunksize_int64) nchunksize_int64 = nssbra
    ntile = int(nssbra/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = nssbra

      !--multi-GPU--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .EQ. res%n_size - 1) nquart_end = iend

      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(res, density, fock, boys_grid_zero) &
      !$omp shared(ss_pair, nssbra, n00bra, nquart_start, nquart_end) &
      !$omp shared(xint00bra) &
      !$omp private(iijj,kkll, shp_thresh, test, ij, kl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp, ii, jj, kk, ll, ish, jsh, ksh, lsh, ijkl, ij_tmp, kl_tmp ) &
      !$omp private(r12, r34, buff, rab, rcd, tmp, sing, cosg, acx, acy, acz, acy2, d00p, aqz, expon_abcd_inverse) &
      !$omp private(bra_loop, ket_loop, t_expon_ab, qps, fmt ) &
      !$omp private(zero_m_0, pqr, pqs, rho, t, t_int, boys0, t_inverse, t_new ) & !fundamental integrals and r-integrals
      !$omp private(iii, jjj, kkk, lll, jj1, kk1, i2, j2, ll1, k2, l2, ik, il, jk, jl, ii1, jj2, kk2, ii2, k,i,l)
      do iijj = nquart_start, nquart_end
        do kkll = 1, nssbra

          if (kkll .gt. iijj) cycle
          test = xint00bra(iijj)*xint00bra(kkll)
          if (test .lt. cutoff_schwarz) cycle

          ij = n00bra(iijj)
          kl = n00bra(kkll)

          ish_tmp = int((1.0 + sqrt(1.0 + 8.0*(real(ij - 1,dp))))/2.0,int64)
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = int((1.0 + sqrt(1.0 + 8.0*(real(kl - 1,dp))))/2.0,int64)
          lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

          ish = res%i_s_shl(ish_tmp)
          jsh = res%i_s_shl(jsh_tmp)
          ksh = res%i_s_shl(ksh_tmp)
          lsh = res%i_s_shl(lsh_tmp)

          ! get distances r12 and r34
          r12 = 0.0_dp 
          r12 = r12 + (res%coord_sh(jsh, 1) - res%coord_sh(ish, 1))*(res%coord_sh(jsh, 1) - res%coord_sh(ish, 1))
          r12 = r12 + (res%coord_sh(jsh, 2) - res%coord_sh(ish, 2))*(res%coord_sh(jsh, 2) - res%coord_sh(ish, 2))
          r12 = r12 + (res%coord_sh(jsh, 3) - res%coord_sh(ish, 3))*(res%coord_sh(jsh, 3) - res%coord_sh(ish, 3))

          r34 = 0.0_dp
          r34 = r34 + (res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1))*(res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1))
          r34 = r34 + (res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2))*(res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2))
          r34 = r34 + (res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3))*(res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3))

          ! rotate axis
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
          if (abs(cosg) .le. 0.9_dp .or. sing .ge. 1.0e-12_dp) then
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
          buff(1) = 0.0_dp !r-integral
          ket_loop = 0
          do k = 1, res%contr_num(ksh)*res%contr_num(lsh)
            ket_loop = ket_loop + 1

            aqz = acz + ss_pair%t_alpha(ss_pair%pair_loc(kl) + ket_loop)*cosg !aqz
            qps = (acx + ss_pair%t_alpha(ss_pair%pair_loc(kl) + ket_loop)*sing)*(acx + ss_pair%t_alpha(ss_pair%pair_loc(kl) + ket_loop)*sing) + acy2 !qps
            zero_m_0 = 0.0_dp
            fmt = 0.0_dp
            bra_loop = 0
            do i = 1, res%contr_num(ish)*res%contr_num(jsh)
              !get bra shell pair info
              bra_loop = bra_loop + 1
              shp_thresh = ss_pair%ismlp(ss_pair%pair_loc(ij) + bra_loop) + ss_pair%ismlp(ss_pair%pair_loc(kl) + ket_loop)
              if (shp_thresh .ge. 2) cycle
              d00p = ss_pair%d_coeff(ss_pair%pair_loc(ij) + bra_loop) !d00p

              expon_abcd_inverse = 1.0_dp/(ss_pair%t_expon_ab(ss_pair%pair_loc(ij) + bra_loop) + ss_pair%t_expon_ab(ss_pair%pair_loc(kl) + ket_loop)) !x41
              pqr = ss_pair%t_alpha(ss_pair%pair_loc(ij) + bra_loop) - aqz
              pqs = pqr*pqr
              rho = ss_pair%t_expon_ab(ss_pair%pair_loc(ij) + bra_loop)*ss_pair%t_expon_ab(ss_pair%pair_loc(kl) + ket_loop)*expon_abcd_inverse
              t = (pqs + qps)*rho
              rho = rho + rho

              if (t .le. t_max) then
                !evaluate boys function with recursion
                t_new = t*1.7140805288547284D+01 !f_increment(2)
                t_int = nint(t_new)
                fmt = boys_grid_zero((t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((t_int*5) + 1)
                boys0 = fmt*sqrt(expon_abcd_inverse)*d00p
              else
                !when t is large, use special formula
                t_inverse = 1.0_dp/t
                boys0 = sqrt(t_inverse*expon_abcd_inverse)*sqrt_pi_div_2*d00p
              end if
              !zero_m (ss|ss) fundamental integrals generation
              zero_m_0(1) = zero_m_0(1) + boys0
            end do
            !r integrals here
            buff(1) = buff(1) + zero_m_0(1)*ss_pair%sq(ss_pair%pair_loc(kl) + ket_loop)
          end do !end ket loop
          ii1 = res%atom_loc(ish)
          jj1 = res%atom_loc(jsh)
          kk1 = res%atom_loc(ksh)
          ll1 = res%atom_loc(lsh)

          if (abs(buff(1)) .lt. cutoff_schwarz) cycle

          i2 = ii1
          j2 = jj1
          k2 = kk1
          l2 = ll1
          if (ii1 .lt. jj1) then ! sort <ij|
            i2 = jj1
            j2 = ii1
          end if
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
          elseif (ii .eq. kk .and. jj .lt. ll) then ! sort <ij|il>
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
      !$omp end target teams distribute parallel do

    end do !tiles

    deallocate (n00bra)
    deallocate (xint00bra)
  end subroutine int0000
end submodule
