! The total angular momentum of this class is:           4
! The algorithm chosen is: PHR a.k.a. ERIC
! Writing an ERIC kernel
submodule(eric_kernels) int0031_impl
contains
  module subroutine int0031(ss_pair, pf_pair, density, fock, res)

    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, pf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

! Variables for the class
    integer(kind=int64), allocatable :: n00bra(:), n13ket(:)
    real(dp), allocatable :: xint00bra(:), xint13ket(:)
    integer(kind=int64) :: nssbra, npfket
    real(dp) :: scutssbra, scutpfket, test
    integer(kind=int64) :: indxi, maxi, indxj, maxj, indxk, maxk, indxl, maxl
    integer(kind=int64) :: leni, lenk, idim, ioff
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, iquart, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, start, end, ijtop, kltop
    integer(kind=int64) :: ish, jsh, ksh, lsh, i, j, k, l, n, iii, t_int, bra_loop, ket_loop
    real(dp) :: t_expon_ab, t_inverse_expon_ab, t_expon_a, t_expon_b
    real(dp) :: t_expon_cd, t_inverse_expon_cd, t_expon_c, t_expon_d
    real(dp) :: ccfket, slket, xkl, ykl, zkl, rxket
    real(dp) :: ccfbra, slbra, xij, yij, zij, rxbra
    real(dp) :: buff(9), scale_factor1(9), scale_factor2(9), cnf(11), fac1, fac2
    real(dp) :: inverse_expon_abcd, ccfbraket, rx, ry, rz, rsq, rho, tt, t2, ftf, xin
    real(dp) :: t_new, fmt, expt
    real(dp) :: ft(5), phi(178), c_factor(35), rxyz(3)
    real(dp) :: work2(100, 1)
    real(dp) :: eri_value(30), angl(20)
    real(dp) :: ai, aij, aijk, aijkl
    integer(kind=int64) :: iord(20)
    integer(kind=int64) :: kstride, lstride
    integer(kind=int64) :: ijk, lo, jc, ir, io, jo, ko
    integer(kind=int64) :: ii1, kk1, nij, maxl2, jj1, j2, ijp, nkl, ijklp
    integer(kind=int64) :: l2, ii2, jj2, kk2, ik, il, jk, jl, ll1, ijkp
    integer(kind=int64) :: maxj2, loci, locj, lock, locl, ip, i2, k2
    integer(kind=int64) :: nchunksize_k10, istart, iend, itile, ntile

    data iord/1, &
      2, 3, 4, &
      5, 7, 10, 6, 8, 9, &
      11, 14, 20, 12, 15, 13, 17, 18, 19, 16/

    data angl/1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, sqrt3, 1.0_dp, sqrt3, sqrt3, 1.0_dp, &
      1.0_dp, sqrt5, sqrt5, 1.0_dp, sqrt5, sqrt15, sqrt5, sqrt5, sqrt5, 1.0_dp/

    allocate (n00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (xint00bra(res%n_s_shl*(res%n_s_shl + 1)/2))
    allocate (n13ket(res%n_p_shl*res%n_f_shl))
    allocate (xint13ket(res%n_p_shl*res%n_f_shl))

! Start screening

    scutssbra = cutoff_schwarz/maxval(ss_pair%xints)
    nssbra = 0
    do ij = 1, res%n_s_shl*(res%n_s_shl + 1)/2
      if (ss_pair%xints(ij) .ge. scutssbra) then
        nssbra = nssbra + 1
        xint00bra(nssbra) = ss_pair%xints(ij)
        n00bra(nssbra) = ij
      end if
    end do

    scutpfket = cutoff_schwarz/maxval(pf_pair%xints)
    npfket = 0
    do ij = 1, res%n_p_shl*res%n_f_shl
      if (pf_pair%xints(ij) .ge. scutpfket) then
        npfket = npfket + 1
        xint13ket(npfket) = pf_pair%xints(ij)
        n13ket(npfket) = ij
      end if
    end do

    nchunksize_k10 = 375000000

    if ((nssbra*npfket) .le. nchunksize_k10) nchunksize_k10 = nssbra*npfket
    ntile = int(nssbra*npfket/nchunksize_k10)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_k10 + 1
      iend = itile*nchunksize_k10
      if (itile .eq. ntile) iend = nssbra*npfket

! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

! Mappings to GPU

!$omp target teams distribute parallel do default(none) &
!$omp shared(res, density, fock, nquart_start, nquart_end, nssbra, xint00bra, n00bra, xint13ket, n13ket, npfket, pf_pair, ss_pair) &
!$omp shared(boys_grid_zero, exponent_grid, iord, angl) &
!$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
!$omp private(ish,jsh,ksh,lsh,ijtop,kltop,eri_value) &
!$omp private(ket_loop,phi,n,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
!$omp private(ccfket,slket,xkl,ykl,zkl,rxket,bra_loop,i) &
!$omp private(t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
!$omp private(ccfbra,slbra,xij,yij,zij,rxbra) &
!$omp private(inverse_expon_abcd,ccfbraket,rx,ry,rz,rsq,rho,tt) &
!$omp private(t_new,t_int,fmt,expt,ft,t2,ftf,xin) &
!$omp private(rxyz,c_factor,fac1,scale_factor1,fac2,scale_factor2,j,iii) &
!$omp private(work2,cnf) &
!$omp private(indxi,maxi,indxj,maxj,indxk,maxk,indxl,maxl,idim,ioff,leni,lenk,ijkl) &
!$omp private(io,ai,jo,aij,jc,ko,aijk,l,lo,aijkl,ir) &
!$omp private(loci,locj,lock,locl,nij,kstride,lstride) &
!$omp private(maxl2,ip,ii1,ijp,jj1,i2,j2,ijkp,nkl,kk1) &
!$omp private(ijklp,buff,ll1,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ii,jj,kk,ll)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/npfket + 1
        kl_tmp = mod(iquart - 1, npfket) + 1

        test = xint00bra(ij_tmp)*xint13ket(kl_tmp)

        if (test .gt. cutoff_schwarz) then

          ij = n00bra(ij_tmp)
          kl = n13ket(kl_tmp)

          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_f_shl + 1
          lsh_tmp = mod(kl - 1, res%n_f_shl) + 1

          ish = res%i_s_shl(ish_tmp)
          jsh = res%i_s_shl(jsh_tmp)
          ksh = res%i_f_shl(lsh_tmp)
          lsh = res%i_p_shl(ksh_tmp)

          ijtop = res%contr_num(ish)*res%contr_num(jsh)

          kltop = res%contr_num(ksh)*res%contr_num(lsh)

          eri_value = 0.0_dp
          ket_loop = 0

          ! 0.0_dp out phi elements to contract over kl
          ! PHI = "pre-Hermite integrals"

          phi(48) = 0.0_dp
          phi(49) = 0.0_dp
          phi(50) = 0.0_dp
          phi(51) = 0.0_dp
          phi(52) = 0.0_dp
          phi(56) = 0.0_dp
          phi(57) = 0.0_dp
          phi(58) = 0.0_dp
          phi(59) = 0.0_dp
          phi(60) = 0.0_dp
          phi(61) = 0.0_dp
          phi(62) = 0.0_dp
          phi(63) = 0.0_dp
          phi(64) = 0.0_dp
          phi(65) = 0.0_dp
          phi(66) = 0.0_dp
          phi(67) = 0.0_dp
          phi(75) = 0.0_dp
          phi(76) = 0.0_dp
          phi(77) = 0.0_dp
          phi(78) = 0.0_dp
          phi(79) = 0.0_dp
          phi(80) = 0.0_dp
          phi(81) = 0.0_dp
          phi(82) = 0.0_dp
          phi(83) = 0.0_dp
          phi(84) = 0.0_dp
          phi(85) = 0.0_dp
          phi(86) = 0.0_dp
          phi(87) = 0.0_dp
          phi(88) = 0.0_dp
          phi(89) = 0.0_dp
          phi(90) = 0.0_dp
          phi(91) = 0.0_dp
          phi(92) = 0.0_dp
          phi(93) = 0.0_dp
          phi(94) = 0.0_dp
          phi(95) = 0.0_dp
          phi(109) = 0.0_dp
          phi(110) = 0.0_dp
          phi(111) = 0.0_dp
          phi(112) = 0.0_dp
          phi(113) = 0.0_dp
          phi(114) = 0.0_dp
          phi(115) = 0.0_dp
          phi(116) = 0.0_dp
          phi(117) = 0.0_dp
          phi(118) = 0.0_dp
          phi(119) = 0.0_dp
          phi(120) = 0.0_dp
          phi(121) = 0.0_dp
          phi(122) = 0.0_dp
          phi(123) = 0.0_dp
          phi(124) = 0.0_dp
          phi(125) = 0.0_dp
          phi(126) = 0.0_dp
          phi(127) = 0.0_dp
          phi(128) = 0.0_dp
          phi(129) = 0.0_dp
          phi(130) = 0.0_dp
          phi(131) = 0.0_dp
          phi(132) = 0.0_dp
          phi(133) = 0.0_dp
          phi(134) = 0.0_dp
          phi(157) = 0.0_dp
          phi(158) = 0.0_dp
          phi(159) = 0.0_dp
          phi(160) = 0.0_dp
          phi(161) = 0.0_dp
          phi(162) = 0.0_dp
          phi(163) = 0.0_dp
          phi(164) = 0.0_dp
          phi(165) = 0.0_dp
          phi(166) = 0.0_dp
          phi(167) = 0.0_dp
          phi(168) = 0.0_dp
          phi(169) = 0.0_dp
          phi(170) = 0.0_dp
          phi(171) = 0.0_dp
          phi(172) = 0.0_dp
          phi(173) = 0.0_dp
          phi(174) = 0.0_dp
          phi(175) = 0.0_dp
          phi(176) = 0.0_dp
          phi(177) = 0.0_dp
          phi(178) = 0.0_dp

! Begin looping over kl primitives

          do k = 1, kltop

            t_expon_cd = pf_pair%t_expon_ab(pf_pair%pair_loc(kl) + 1 + ket_loop) ! exp_c + exp_d
            t_expon_c = pf_pair%expon_b(pf_pair%pair_loc(kl) + 1 + ket_loop)
            t_expon_d = pf_pair%expon_a(pf_pair%pair_loc(kl) + 1 + ket_loop)

            t_inverse_expon_cd = 1.0_dp/t_expon_cd
            ccfket = pf_pair%sq(pf_pair%pair_loc(kl) + 1 + ket_loop)*sqrt2_pi_5_4
            slket = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_cd)

            xkl = ((t_expon_c*res%coord_sh(ksh, 1)) + (t_expon_d*res%coord_sh(lsh, 1)))*t_inverse_expon_cd
            ykl = ((t_expon_c*res%coord_sh(ksh, 2)) + (t_expon_d*res%coord_sh(lsh, 2)))*t_inverse_expon_cd
            zkl = ((t_expon_c*res%coord_sh(ksh, 3)) + (t_expon_d*res%coord_sh(lsh, 3)))*t_inverse_expon_cd

            rxket = pf_pair%t_inverse_expon_ab(pf_pair%pair_loc(kl) + 1 + ket_loop)! inverse_expon_ab*0.5_dp

            ket_loop = ket_loop + 1
            bra_loop = 0

            ! 0.0_dp out phi elements to contract over ij

            phi(47) = 0.0_dp
            phi(53) = 0.0_dp
            phi(54) = 0.0_dp
            phi(55) = 0.0_dp
            phi(68) = 0.0_dp
            phi(69) = 0.0_dp
            phi(70) = 0.0_dp
            phi(71) = 0.0_dp
            phi(72) = 0.0_dp
            phi(73) = 0.0_dp
            phi(74) = 0.0_dp
            phi(96) = 0.0_dp
            phi(97) = 0.0_dp
            phi(98) = 0.0_dp
            phi(99) = 0.0_dp
            phi(100) = 0.0_dp
            phi(101) = 0.0_dp
            phi(102) = 0.0_dp
            phi(103) = 0.0_dp
            phi(104) = 0.0_dp
            phi(105) = 0.0_dp
            phi(106) = 0.0_dp
            phi(107) = 0.0_dp
            phi(108) = 0.0_dp
            phi(135) = 0.0_dp
            phi(136) = 0.0_dp
            phi(137) = 0.0_dp
            phi(138) = 0.0_dp
            phi(139) = 0.0_dp
            phi(140) = 0.0_dp
            phi(141) = 0.0_dp
            phi(142) = 0.0_dp
            phi(143) = 0.0_dp
            phi(144) = 0.0_dp
            phi(145) = 0.0_dp
            phi(146) = 0.0_dp
            phi(147) = 0.0_dp
            phi(148) = 0.0_dp
            phi(149) = 0.0_dp
            phi(150) = 0.0_dp
            phi(151) = 0.0_dp
            phi(152) = 0.0_dp
            phi(153) = 0.0_dp
            phi(154) = 0.0_dp
            phi(155) = 0.0_dp
            phi(156) = 0.0_dp

! Begin looping over ij primitives

            do i = 1, ijtop

              t_expon_ab = ss_pair%t_expon_ab(ss_pair%pair_loc(ij) + 1 + bra_loop) ! exp_c + exp_d
              t_expon_a = ss_pair%expon_a(ss_pair%pair_loc(ij) + 1 + bra_loop)
              t_expon_b = ss_pair%expon_b(ss_pair%pair_loc(ij) + 1 + bra_loop)

              t_inverse_expon_ab = 1.0_dp/t_expon_ab
              ccfbra = ss_pair%sq(ss_pair%pair_loc(ij) + 1 + bra_loop)*sqrt2_pi_5_4
              slbra = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_ab)

              xij = ((t_expon_a*res%coord_sh(ish, 1)) + (t_expon_b*res%coord_sh(jsh, 1)))*t_inverse_expon_ab
              yij = ((t_expon_a*res%coord_sh(ish, 2)) + (t_expon_b*res%coord_sh(jsh, 2)))*t_inverse_expon_ab
              zij = ((t_expon_a*res%coord_sh(ish, 3)) + (t_expon_b*res%coord_sh(jsh, 3)))*t_inverse_expon_ab

              rxbra = ss_pair%t_inverse_expon_ab(ss_pair%pair_loc(ij) + 1 + bra_loop)! inverse_expon_ab*0.5_dp

              bra_loop = bra_loop + 1

! Get 4-index data
              inverse_expon_abcd = 1.0_dp/(t_expon_ab + t_expon_cd)
              ccfbraket = ccfbra*ccfket
              rx = xkl - xij
              ry = ykl - yij
              rz = zkl - zij
              rsq = (rx*rx) + (ry*ry) + (rz*rz)
              rho = t_expon_cd*t_expon_ab*inverse_expon_abcd
              tt = rsq*rho

              n = 4
              if (tt .le. t_max) then

                ! Boys function evaluation with Chebyshev interpolation

                t_new = tt*0.1536593665137801D+02
                t_int = nint(t_new)
                fmt = boys_grid_zero((n*451*5) + (t_int*5) + 5)*t_new
                fmt = (fmt + boys_grid_zero((n*451*5) + (t_int*5) + 4))*t_new
                fmt = (fmt + boys_grid_zero((n*451*5) + (t_int*5) + 3))*t_new
                fmt = (fmt + boys_grid_zero((n*451*5) + (t_int*5) + 2))*t_new
                fmt = fmt + boys_grid_zero((n*451*5) + (t_int*5) + 1)
                t_new = tt*2.768915858120726D+01
                t_int = nint(t_new)
                expt = exponent_grid((t_int*5) + 5)*t_new
                expt = (expt + exponent_grid((t_int*5) + 4))*t_new
                expt = (expt + exponent_grid((t_int*5) + 3))*t_new
                expt = (expt + exponent_grid((t_int*5) + 2))*t_new
                expt = expt + exponent_grid((t_int*5) + 1)

                ! Generate necessary f_m

                ft(5) = fmt
                t2 = tt + tt
                ft(4) = (t2*ft(5) + expt)*0.142857142857143D+00
                ft(3) = (t2*ft(4) + expt)*0.200000000000000D+00
                ft(2) = (t2*ft(3) + expt)*0.333333333333333D+00
                ft(1) = (t2*ft(2) + expt)*0.100000000000000D+01
                !
                rho = rho + rho
                ftf = ccfbraket*sqrt(inverse_expon_abcd)
                ft(1) = ft(1)*ftf
                ftf = ftf*rho
                ft(2) = ft(2)*ftf
                ftf = ftf*rho
                ft(3) = ft(3)*ftf
                ftf = ftf*rho
                ft(4) = ft(4)*ftf
                ftf = ftf*rho
                ft(5) = ft(5)*ftf

              else

                ! Large T formula

                xin = 1.0_dp/rsq
                ftf = ccfbraket*sqrt(xin)*slbra*slket
                ft(1) = ftf
                ftf = ftf*xin
                ft(2) = ftf*0.100000000000000D+01
                ftf = ftf*xin
                ft(3) = ftf*0.300000000000000D+01
                ftf = ftf*xin
                ft(4) = ftf*0.150000000000000D+02
                ftf = ftf*xin
                ft(5) = ftf*0.105000000000000D+03
              end if

! ft(1 to L+1) now contains the [0]^0, [0]^(1), ... [0]^(L)

! Generate Pre-Hermite Integrals (PHI)

              rxyz(1) = rx
              rxyz(2) = ry
              rxyz(3) = rz
              c_factor(1) = 1.0_dp

! Cartesian factors for px,py,pz
              c_factor(2) = rx
              c_factor(3) = ry
              c_factor(4) = rz
              ! Cartesian factors for (in order):
              ! dxx,dxy,dyy,dxz,dyz,dzz
              c_factor(5) = c_factor(2)*rxyz(1)
              c_factor(6) = c_factor(3)*rxyz(1)
              c_factor(7) = c_factor(3)*rxyz(2)
              c_factor(8) = c_factor(4)*rxyz(1)
              c_factor(9) = c_factor(4)*rxyz(2)
              c_factor(10) = c_factor(4)*rxyz(3)
              ! Cartesian factors for (in order):
              ! fxxx,fxxy,fxyy,fyyy,fxxz,fxyz,fyyz,fzzx,fzzy,fzzz
              c_factor(11) = c_factor(5)*rxyz(1)
              c_factor(12) = c_factor(6)*rxyz(1)
              c_factor(13) = c_factor(7)*rxyz(1)
              c_factor(14) = c_factor(7)*rxyz(2)
              c_factor(15) = c_factor(8)*rxyz(1)
              c_factor(16) = c_factor(9)*rxyz(1)
              c_factor(17) = c_factor(9)*rxyz(2)
              c_factor(18) = c_factor(10)*rxyz(1)
              c_factor(19) = c_factor(10)*rxyz(2)
              c_factor(20) = c_factor(10)*rxyz(3)
! Still need to figure these out :P
              c_factor(21) = c_factor(11)*rxyz(1)
              c_factor(22) = c_factor(12)*rxyz(1)
              c_factor(23) = c_factor(13)*rxyz(1)
              c_factor(24) = c_factor(14)*rxyz(1)
              c_factor(25) = c_factor(14)*rxyz(2)
              c_factor(26) = c_factor(15)*rxyz(1)
              c_factor(27) = c_factor(16)*rxyz(1)
              c_factor(28) = c_factor(17)*rxyz(1)
              c_factor(29) = c_factor(17)*rxyz(2)
              c_factor(30) = c_factor(18)*rxyz(1)
              c_factor(31) = c_factor(19)*rxyz(1)
              c_factor(32) = c_factor(19)*rxyz(2)
              c_factor(33) = c_factor(20)*rxyz(1)
              c_factor(34) = c_factor(20)*rxyz(2)
              c_factor(35) = c_factor(20)*rxyz(3)
              ! Generate pre-Hermite integrals with cartesian factors
              ! and Boys function
              phi(1) = ft(1)
              phi(2) = c_factor(2)*ft(2)
              phi(3) = c_factor(3)*ft(2)
              phi(4) = c_factor(4)*ft(2)
              phi(5) = c_factor(1)*ft(2)
              phi(6) = c_factor(5)*ft(3)
              phi(7) = c_factor(6)*ft(3)
              phi(8) = c_factor(7)*ft(3)
              phi(9) = c_factor(8)*ft(3)
              phi(10) = c_factor(9)*ft(3)
              phi(11) = c_factor(10)*ft(3)
              phi(12) = c_factor(2)*ft(3)
              phi(13) = c_factor(3)*ft(3)
              phi(14) = c_factor(4)*ft(3)
              phi(15) = c_factor(11)*ft(4)
              phi(16) = c_factor(12)*ft(4)
              phi(17) = c_factor(13)*ft(4)
              phi(18) = c_factor(14)*ft(4)
              phi(19) = c_factor(15)*ft(4)
              phi(20) = c_factor(16)*ft(4)
              phi(21) = c_factor(17)*ft(4)
              phi(22) = c_factor(18)*ft(4)
              phi(23) = c_factor(19)*ft(4)
              phi(24) = c_factor(20)*ft(4)
              phi(25) = c_factor(1)*ft(3)
              phi(26) = c_factor(5)*ft(4)
              phi(27) = c_factor(6)*ft(4)
              phi(28) = c_factor(7)*ft(4)
              phi(29) = c_factor(8)*ft(4)
              phi(30) = c_factor(9)*ft(4)
              phi(31) = c_factor(10)*ft(4)
              phi(32) = c_factor(21)*ft(5)
              phi(33) = c_factor(22)*ft(5)
              phi(34) = c_factor(23)*ft(5)
              phi(35) = c_factor(24)*ft(5)
              phi(36) = c_factor(25)*ft(5)
              phi(37) = c_factor(26)*ft(5)
              phi(38) = c_factor(27)*ft(5)
              phi(39) = c_factor(28)*ft(5)
              phi(40) = c_factor(29)*ft(5)
              phi(41) = c_factor(30)*ft(5)
              phi(42) = c_factor(31)*ft(5)
              phi(43) = c_factor(32)*ft(5)
              phi(44) = c_factor(33)*ft(5)
              phi(45) = c_factor(34)*ft(5)
              phi(46) = c_factor(35)*ft(5)

              ! Begin contraction & scaling over ij
              ! i shell scaling factors
              fac1 = 1.0_dp
              scale_factor1(1) = 1.0_dp
              ! j shell scaling factors
              fac2 = 1.0_dp
              scale_factor2(1) = 1.0_dp

              ! ij contraction
              phi(47) = phi(47) + phi(1)

              phi(53) = phi(53) + phi(2)
              phi(54) = phi(54) + phi(3)
              phi(55) = phi(55) + phi(4)

              phi(68) = phi(68) + phi(5)
              phi(69) = phi(69) + phi(6)
              phi(70) = phi(70) + phi(7)
              phi(71) = phi(71) + phi(8)
              phi(72) = phi(72) + phi(9)
              phi(73) = phi(73) + phi(10)
              phi(74) = phi(74) + phi(11)

              phi(96) = phi(96) + phi(12)
              phi(97) = phi(97) + phi(13)
              phi(98) = phi(98) + phi(14)
              phi(99) = phi(99) + phi(15)
              phi(100) = phi(100) + phi(16)
              phi(101) = phi(101) + phi(17)
              phi(102) = phi(102) + phi(18)
              phi(103) = phi(103) + phi(19)
              phi(104) = phi(104) + phi(20)
              phi(105) = phi(105) + phi(21)
              phi(106) = phi(106) + phi(22)
              phi(107) = phi(107) + phi(23)
              phi(108) = phi(108) + phi(24)

              phi(135) = phi(135) + phi(25)
              phi(136) = phi(136) + phi(26)
              phi(137) = phi(137) + phi(27)
              phi(138) = phi(138) + phi(28)
              phi(139) = phi(139) + phi(29)
              phi(140) = phi(140) + phi(30)
              phi(141) = phi(141) + phi(31)
              phi(142) = phi(142) + phi(32)
              phi(143) = phi(143) + phi(33)
              phi(144) = phi(144) + phi(34)
              phi(145) = phi(145) + phi(35)
              phi(146) = phi(146) + phi(36)
              phi(147) = phi(147) + phi(37)
              phi(148) = phi(148) + phi(38)
              phi(149) = phi(149) + phi(39)
              phi(150) = phi(150) + phi(40)
              phi(151) = phi(151) + phi(41)
              phi(152) = phi(152) + phi(42)
              phi(153) = phi(153) + phi(43)
              phi(154) = phi(154) + phi(44)
              phi(155) = phi(155) + phi(45)
              phi(156) = phi(156) + phi(46)

            end do ! i

            ! Begin contraction & scaling over kl
            ! k shell scaling factors
            fac1 = 1.0_dp
            scale_factor1(1) = 1.0_dp
            fac1 = fac1*rxket
            scale_factor1(2) = fac1 ! (0.5*inv_exp_cd)** 1
            fac1 = fac1*rxket
            scale_factor1(3) = fac1 ! (0.5*inv_exp_cd)** 2
            fac1 = fac1*rxket
            scale_factor1(4) = fac1 ! (0.5*inv_exp_cd)** 3
            fac1 = fac1*rxket
            scale_factor1(5) = fac1 ! (0.5*inv_exp_cd)** 4
            ! l shell scaling factors
            fac2 = 1.0_dp
            scale_factor2(1) = 1.0_dp
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(2) = fac2 ! (2*t_expon_d)** 1
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(3) = fac2 ! (2*t_expon_d)** 2
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(4) = fac2 ! (2*t_expon_d)** 3
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(5) = fac2 ! (2*t_expon_d)** 4

            ! kl contraction
            phi(48) = phi(48) + scale_factor1(3)*phi(47)
            phi(49) = phi(49) + scale_factor2(2)*scale_factor1(3)*phi(47)
            phi(50) = phi(50) + scale_factor2(3)*scale_factor1(4)*phi(47)
            phi(51) = phi(51) + scale_factor2(4)*scale_factor1(4)*phi(47)
            phi(52) = phi(52) + scale_factor2(5)*scale_factor1(5)*phi(47)

            phi(56) = phi(56) + scale_factor1(3)*phi(53)
            phi(57) = phi(57) + scale_factor1(3)*phi(54)
            phi(58) = phi(58) + scale_factor1(3)*phi(55)
            phi(59) = phi(59) + scale_factor2(2)*scale_factor1(4)*phi(53)
            phi(60) = phi(60) + scale_factor2(2)*scale_factor1(4)*phi(54)
            phi(61) = phi(61) + scale_factor2(2)*scale_factor1(4)*phi(55)
            phi(62) = phi(62) + scale_factor2(3)*scale_factor1(4)*phi(53)
            phi(63) = phi(63) + scale_factor2(3)*scale_factor1(4)*phi(54)
            phi(64) = phi(64) + scale_factor2(3)*scale_factor1(4)*phi(55)
            phi(65) = phi(65) + scale_factor2(4)*scale_factor1(5)*phi(53)
            phi(66) = phi(66) + scale_factor2(4)*scale_factor1(5)*phi(54)
            phi(67) = phi(67) + scale_factor2(4)*scale_factor1(5)*phi(55)

            phi(75) = phi(75) + scale_factor1(4)*phi(68)
            phi(76) = phi(76) + scale_factor1(4)*phi(69)
            phi(77) = phi(77) + scale_factor1(4)*phi(70)
            phi(78) = phi(78) + scale_factor1(4)*phi(71)
            phi(79) = phi(79) + scale_factor1(4)*phi(72)
            phi(80) = phi(80) + scale_factor1(4)*phi(73)
            phi(81) = phi(81) + scale_factor1(4)*phi(74)
            phi(82) = phi(82) + scale_factor2(2)*scale_factor1(4)*phi(68)
            phi(83) = phi(83) + scale_factor2(2)*scale_factor1(4)*phi(69)
            phi(84) = phi(84) + scale_factor2(2)*scale_factor1(4)*phi(70)
            phi(85) = phi(85) + scale_factor2(2)*scale_factor1(4)*phi(71)
            phi(86) = phi(86) + scale_factor2(2)*scale_factor1(4)*phi(72)
            phi(87) = phi(87) + scale_factor2(2)*scale_factor1(4)*phi(73)
            phi(88) = phi(88) + scale_factor2(2)*scale_factor1(4)*phi(74)
            phi(89) = phi(89) + scale_factor2(3)*scale_factor1(5)*phi(68)
            phi(90) = phi(90) + scale_factor2(3)*scale_factor1(5)*phi(69)
            phi(91) = phi(91) + scale_factor2(3)*scale_factor1(5)*phi(70)
            phi(92) = phi(92) + scale_factor2(3)*scale_factor1(5)*phi(71)
            phi(93) = phi(93) + scale_factor2(3)*scale_factor1(5)*phi(72)
            phi(94) = phi(94) + scale_factor2(3)*scale_factor1(5)*phi(73)
            phi(95) = phi(95) + scale_factor2(3)*scale_factor1(5)*phi(74)

            phi(109) = phi(109) + scale_factor1(4)*phi(96)
            phi(110) = phi(110) + scale_factor1(4)*phi(97)
            phi(111) = phi(111) + scale_factor1(4)*phi(98)
            phi(112) = phi(112) + scale_factor1(4)*phi(99)
            phi(113) = phi(113) + scale_factor1(4)*phi(100)
            phi(114) = phi(114) + scale_factor1(4)*phi(101)
            phi(115) = phi(115) + scale_factor1(4)*phi(102)
            phi(116) = phi(116) + scale_factor1(4)*phi(103)
            phi(117) = phi(117) + scale_factor1(4)*phi(104)
            phi(118) = phi(118) + scale_factor1(4)*phi(105)
            phi(119) = phi(119) + scale_factor1(4)*phi(106)
            phi(120) = phi(120) + scale_factor1(4)*phi(107)
            phi(121) = phi(121) + scale_factor1(4)*phi(108)
            phi(122) = phi(122) + scale_factor2(2)*scale_factor1(5)*phi(96)
            phi(123) = phi(123) + scale_factor2(2)*scale_factor1(5)*phi(97)
            phi(124) = phi(124) + scale_factor2(2)*scale_factor1(5)*phi(98)
            phi(125) = phi(125) + scale_factor2(2)*scale_factor1(5)*phi(99)
            phi(126) = phi(126) + scale_factor2(2)*scale_factor1(5)*phi(100)
            phi(127) = phi(127) + scale_factor2(2)*scale_factor1(5)*phi(101)
            phi(128) = phi(128) + scale_factor2(2)*scale_factor1(5)*phi(102)
            phi(129) = phi(129) + scale_factor2(2)*scale_factor1(5)*phi(103)
            phi(130) = phi(130) + scale_factor2(2)*scale_factor1(5)*phi(104)
            phi(131) = phi(131) + scale_factor2(2)*scale_factor1(5)*phi(105)
            phi(132) = phi(132) + scale_factor2(2)*scale_factor1(5)*phi(106)
            phi(133) = phi(133) + scale_factor2(2)*scale_factor1(5)*phi(107)
            phi(134) = phi(134) + scale_factor2(2)*scale_factor1(5)*phi(108)

            phi(157) = phi(157) + scale_factor1(5)*phi(135)
            phi(158) = phi(158) + scale_factor1(5)*phi(136)
            phi(159) = phi(159) + scale_factor1(5)*phi(137)
            phi(160) = phi(160) + scale_factor1(5)*phi(138)
            phi(161) = phi(161) + scale_factor1(5)*phi(139)
            phi(162) = phi(162) + scale_factor1(5)*phi(140)
            phi(163) = phi(163) + scale_factor1(5)*phi(141)
            phi(164) = phi(164) + scale_factor1(5)*phi(142)
            phi(165) = phi(165) + scale_factor1(5)*phi(143)
            phi(166) = phi(166) + scale_factor1(5)*phi(144)
            phi(167) = phi(167) + scale_factor1(5)*phi(145)
            phi(168) = phi(168) + scale_factor1(5)*phi(146)
            phi(169) = phi(169) + scale_factor1(5)*phi(147)
            phi(170) = phi(170) + scale_factor1(5)*phi(148)
            phi(171) = phi(171) + scale_factor1(5)*phi(149)
            phi(172) = phi(172) + scale_factor1(5)*phi(150)
            phi(173) = phi(173) + scale_factor1(5)*phi(151)
            phi(174) = phi(174) + scale_factor1(5)*phi(152)
            phi(175) = phi(175) + scale_factor1(5)*phi(153)
            phi(176) = phi(176) + scale_factor1(5)*phi(154)
            phi(177) = phi(177) + scale_factor1(5)*phi(155)
            phi(178) = phi(178) + scale_factor1(5)*phi(156)

          end do ! k

!                   ******************************
!                   *                            *
!                   *   Post-contraction phase   *
!                   *                            *
!                   ******************************

          phi(76) = phi(76) - phi(75)
          phi(78) = phi(78) - phi(75)
          phi(81) = phi(81) - phi(75)

          phi(83) = phi(83) - phi(82)
          phi(85) = phi(85) - phi(82)
          phi(88) = phi(88) - phi(82)

          phi(90) = phi(90) - phi(89)
          phi(92) = phi(92) - phi(89)
          phi(95) = phi(95) - phi(89)

          phi(112) = phi(112) - phi(109)
          phi(113) = phi(113) - phi(110)
          phi(115) = phi(115) - phi(110)
          phi(116) = phi(116) - phi(111)
          phi(118) = phi(118) - phi(111)
          phi(121) = phi(121) - phi(111)
          phi(112) = phi(112) - (2.0D+00)*phi(109)
          phi(114) = phi(114) - phi(109)
          phi(115) = phi(115) - (2.0D+00)*phi(110)
          phi(119) = phi(119) - phi(109)
          phi(120) = phi(120) - phi(110)
          phi(121) = phi(121) - (2.0D+00)*phi(111)

          phi(125) = phi(125) - phi(122)
          phi(126) = phi(126) - phi(123)
          phi(128) = phi(128) - phi(123)
          phi(129) = phi(129) - phi(124)
          phi(131) = phi(131) - phi(124)
          phi(134) = phi(134) - phi(124)
          phi(125) = phi(125) - (2.0D+00)*phi(122)
          phi(127) = phi(127) - phi(122)
          phi(128) = phi(128) - (2.0D+00)*phi(123)
          phi(132) = phi(132) - phi(122)
          phi(133) = phi(133) - phi(123)
          phi(134) = phi(134) - (2.0D+00)*phi(124)

          phi(164) = phi(164) - phi(158)
          phi(165) = phi(165) - phi(159)
          phi(166) = phi(166) - phi(160)
          phi(168) = phi(168) - phi(160)
          phi(169) = phi(169) - phi(161)
          phi(170) = phi(170) - phi(162)
          phi(172) = phi(172) - phi(162)
          phi(173) = phi(173) - phi(163)
          phi(175) = phi(175) - phi(163)
          phi(178) = phi(178) - phi(163)
          phi(164) = phi(164) - (2.0D+00)*phi(158)
          phi(165) = phi(165) - (2.0D+00)*phi(159)
          phi(167) = phi(167) - phi(159)
          phi(168) = phi(168) - (2.0D+00)*phi(160)
          phi(169) = phi(169) - (2.0D+00)*phi(161)
          phi(171) = phi(171) - phi(161)
          phi(172) = phi(172) - (2.0D+00)*phi(162)
          phi(176) = phi(176) - phi(161)
          phi(177) = phi(177) - phi(162)
          phi(178) = phi(178) - (2.0D+00)*phi(163)
          phi(158) = phi(158) - phi(157)
          phi(160) = phi(160) - phi(157)
          phi(163) = phi(163) - phi(157)
          phi(164) = phi(164) - (3.0D+00)*phi(158)
          phi(166) = phi(166) - phi(158)
          phi(167) = phi(167) - (2.0D+00)*phi(159)
          phi(168) = phi(168) - (3.0D+00)*phi(160)
          phi(173) = phi(173) - phi(158)
          phi(174) = phi(174) - phi(159)
          phi(175) = phi(175) - phi(160)
          phi(176) = phi(176) - (2.0D+00)*phi(161)
          phi(177) = phi(177) - (2.0D+00)*phi(162)
          phi(178) = phi(178) - (3.0D+00)*phi(163)

!           *************************************
!           *                                   *
!           * -- Angular momentum conversion -- *
!           *   1-center (r) to 2-center (p|q)  *
!           *  Bra contracted transfer equation *
!           *   Horizontal recurrence relation  *
!           *                                   *
!           *************************************

          work2(1, 1) = phi(51)
          work2(2, 1) = phi(49)
          work2(3, 1) = -phi(62)
          work2(4, 1) = -phi(63)
          work2(5, 1) = -phi(64)
          work2(6, 1) = -phi(56)
          work2(7, 1) = -phi(57)
          work2(8, 1) = -phi(58)
          work2(9, 1) = phi(83)
          work2(10, 1) = phi(84)
          work2(11, 1) = phi(85)
          work2(12, 1) = phi(86)
          work2(13, 1) = phi(87)
          work2(14, 1) = phi(88)
          work2(15, 1) = -phi(112)
          work2(16, 1) = -phi(113)
          work2(17, 1) = -phi(114)
          work2(18, 1) = -phi(115)
          work2(19, 1) = -phi(116)
          work2(20, 1) = -phi(117)
          work2(21, 1) = -phi(118)
          work2(22, 1) = -phi(119)
          work2(23, 1) = -phi(120)
          work2(24, 1) = -phi(121)
          work2(25, 1) = phi(52)
          work2(26, 1) = phi(50)
          work2(27, 1) = phi(48)
          work2(28, 1) = -phi(65)
          work2(29, 1) = -phi(66)
          work2(30, 1) = -phi(67)
          work2(31, 1) = -phi(59)
          work2(32, 1) = -phi(60)
          work2(33, 1) = -phi(61)
          work2(34, 1) = phi(90)
          work2(35, 1) = phi(91)
          work2(36, 1) = phi(92)
          work2(37, 1) = phi(93)
          work2(38, 1) = phi(94)
          work2(39, 1) = phi(95)
          work2(40, 1) = phi(76)
          work2(41, 1) = phi(77)
          work2(42, 1) = phi(78)
          work2(43, 1) = phi(79)
          work2(44, 1) = phi(80)
          work2(45, 1) = phi(81)
          work2(46, 1) = -phi(125)
          work2(47, 1) = -phi(126)
          work2(48, 1) = -phi(127)
          work2(49, 1) = -phi(128)
          work2(50, 1) = -phi(129)
          work2(51, 1) = -phi(130)
          work2(52, 1) = -phi(131)
          work2(53, 1) = -phi(132)
          work2(54, 1) = -phi(133)
          work2(55, 1) = -phi(134)
          work2(56, 1) = phi(164)
          work2(57, 1) = phi(165)
          work2(58, 1) = phi(166)
          work2(59, 1) = phi(167)
          work2(60, 1) = phi(168)
          work2(61, 1) = phi(169)
          work2(62, 1) = phi(170)
          work2(63, 1) = phi(171)
          work2(64, 1) = phi(172)
          work2(65, 1) = phi(173)
          work2(66, 1) = phi(174)
          work2(67, 1) = phi(175)
          work2(68, 1) = phi(176)
          work2(69, 1) = phi(177)
          work2(70, 1) = phi(178)

!           *************************************
!           *                                   *
!           *  ket contracted transfer equation *
!           *   Horizontal recurrence relation  *
!           *                                   *
!           *************************************

          cnf(1) = res%coord_sh(lsh, 1) - res%coord_sh(ksh, 1)
          cnf(2) = res%coord_sh(lsh, 2) - res%coord_sh(ksh, 2)
          cnf(3) = res%coord_sh(lsh, 3) - res%coord_sh(ksh, 3)
          cnf(4) = 1.0D+00
          cnf(5) = 2.0D+00
          cnf(6) = 3.0D+00
          cnf(7) = 4.0D+00
          cnf(8) = 5.0D+00
          cnf(9) = 6.0D+00
          cnf(10) = 7.0D+00
          cnf(11) = 8.0D+00

          work2(15, 1) = work2(15, 1) + cnf(1)*work2(9, 1)
          work2(15, 1) = work2(15, 1) + cnf(5)*work2(6, 1)
          work2(16, 1) = work2(16, 1) + cnf(1)*work2(10, 1)
          work2(16, 1) = work2(16, 1) + cnf(4)*work2(7, 1)
          work2(17, 1) = work2(17, 1) + cnf(1)*work2(11, 1)
          work2(18, 1) = work2(18, 1) + cnf(2)*work2(11, 1)
          work2(18, 1) = work2(18, 1) + cnf(5)*work2(7, 1)
          work2(19, 1) = work2(19, 1) + cnf(1)*work2(12, 1)
          work2(19, 1) = work2(19, 1) + cnf(4)*work2(8, 1)
          work2(20, 1) = work2(20, 1) + cnf(1)*work2(13, 1)
          work2(21, 1) = work2(21, 1) + cnf(2)*work2(13, 1)
          work2(21, 1) = work2(21, 1) + cnf(4)*work2(8, 1)
          work2(22, 1) = work2(22, 1) + cnf(1)*work2(14, 1)
          work2(23, 1) = work2(23, 1) + cnf(2)*work2(14, 1)
          work2(24, 1) = work2(24, 1) + cnf(3)*work2(14, 1)
          work2(24, 1) = work2(24, 1) + cnf(5)*work2(8, 1)
          work2(9, 1) = work2(9, 1) + cnf(1)*work2(3, 1)
          work2(9, 1) = work2(9, 1) + cnf(4)*work2(2, 1)
          work2(10, 1) = work2(10, 1) + cnf(1)*work2(4, 1)
          work2(11, 1) = work2(11, 1) + cnf(2)*work2(4, 1)
          work2(11, 1) = work2(11, 1) + cnf(4)*work2(2, 1)
          work2(12, 1) = work2(12, 1) + cnf(1)*work2(5, 1)
          work2(13, 1) = work2(13, 1) + cnf(2)*work2(5, 1)
          work2(14, 1) = work2(14, 1) + cnf(3)*work2(5, 1)
          work2(14, 1) = work2(14, 1) + cnf(4)*work2(2, 1)
          work2(3, 1) = work2(3, 1) + cnf(1)*work2(1, 1)
          work2(4, 1) = work2(4, 1) + cnf(2)*work2(1, 1)
          work2(5, 1) = work2(5, 1) + cnf(3)*work2(1, 1)
          work2(6, 1) = work2(6, 1) + cnf(1)*work2(2, 1)
          work2(7, 1) = work2(7, 1) + cnf(2)*work2(2, 1)
          work2(8, 1) = work2(8, 1) + cnf(3)*work2(2, 1)
          work2(56, 1) = work2(56, 1) + cnf(1)*work2(46, 1)
          work2(56, 1) = work2(56, 1) + cnf(6)*work2(40, 1)
          work2(57, 1) = work2(57, 1) + cnf(1)*work2(47, 1)
          work2(57, 1) = work2(57, 1) + cnf(5)*work2(41, 1)
          work2(58, 1) = work2(58, 1) + cnf(1)*work2(48, 1)
          work2(58, 1) = work2(58, 1) + cnf(4)*work2(42, 1)
          work2(59, 1) = work2(59, 1) + cnf(1)*work2(49, 1)
          work2(60, 1) = work2(60, 1) + cnf(2)*work2(49, 1)
          work2(60, 1) = work2(60, 1) + cnf(6)*work2(42, 1)
          work2(61, 1) = work2(61, 1) + cnf(1)*work2(50, 1)
          work2(61, 1) = work2(61, 1) + cnf(5)*work2(43, 1)
          work2(62, 1) = work2(62, 1) + cnf(1)*work2(51, 1)
          work2(62, 1) = work2(62, 1) + cnf(4)*work2(44, 1)
          work2(63, 1) = work2(63, 1) + cnf(1)*work2(52, 1)
          work2(64, 1) = work2(64, 1) + cnf(2)*work2(52, 1)
          work2(64, 1) = work2(64, 1) + cnf(5)*work2(44, 1)
          work2(65, 1) = work2(65, 1) + cnf(1)*work2(53, 1)
          work2(65, 1) = work2(65, 1) + cnf(4)*work2(45, 1)
          work2(66, 1) = work2(66, 1) + cnf(1)*work2(54, 1)
          work2(67, 1) = work2(67, 1) + cnf(2)*work2(54, 1)
          work2(67, 1) = work2(67, 1) + cnf(4)*work2(45, 1)
          work2(68, 1) = work2(68, 1) + cnf(1)*work2(55, 1)
          work2(69, 1) = work2(69, 1) + cnf(2)*work2(55, 1)
          work2(70, 1) = work2(70, 1) + cnf(3)*work2(55, 1)
          work2(70, 1) = work2(70, 1) + cnf(6)*work2(45, 1)
          work2(46, 1) = work2(46, 1) + cnf(1)*work2(34, 1)
          work2(46, 1) = work2(46, 1) + cnf(5)*work2(31, 1)
          work2(47, 1) = work2(47, 1) + cnf(1)*work2(35, 1)
          work2(47, 1) = work2(47, 1) + cnf(4)*work2(32, 1)
          work2(48, 1) = work2(48, 1) + cnf(1)*work2(36, 1)
          work2(49, 1) = work2(49, 1) + cnf(2)*work2(36, 1)
          work2(49, 1) = work2(49, 1) + cnf(5)*work2(32, 1)
          work2(50, 1) = work2(50, 1) + cnf(1)*work2(37, 1)
          work2(50, 1) = work2(50, 1) + cnf(4)*work2(33, 1)
          work2(51, 1) = work2(51, 1) + cnf(1)*work2(38, 1)
          work2(52, 1) = work2(52, 1) + cnf(2)*work2(38, 1)
          work2(52, 1) = work2(52, 1) + cnf(4)*work2(33, 1)
          work2(53, 1) = work2(53, 1) + cnf(1)*work2(39, 1)
          work2(54, 1) = work2(54, 1) + cnf(2)*work2(39, 1)
          work2(55, 1) = work2(55, 1) + cnf(3)*work2(39, 1)
          work2(55, 1) = work2(55, 1) + cnf(5)*work2(33, 1)
          work2(34, 1) = work2(34, 1) + cnf(1)*work2(28, 1)
          work2(34, 1) = work2(34, 1) + cnf(4)*work2(26, 1)
          work2(35, 1) = work2(35, 1) + cnf(1)*work2(29, 1)
          work2(36, 1) = work2(36, 1) + cnf(2)*work2(29, 1)
          work2(36, 1) = work2(36, 1) + cnf(4)*work2(26, 1)
          work2(37, 1) = work2(37, 1) + cnf(1)*work2(30, 1)
          work2(38, 1) = work2(38, 1) + cnf(2)*work2(30, 1)
          work2(39, 1) = work2(39, 1) + cnf(3)*work2(30, 1)
          work2(39, 1) = work2(39, 1) + cnf(4)*work2(26, 1)
          work2(40, 1) = work2(40, 1) + cnf(1)*work2(31, 1)
          work2(40, 1) = work2(40, 1) + cnf(4)*work2(27, 1)
          work2(41, 1) = work2(41, 1) + cnf(1)*work2(32, 1)
          work2(42, 1) = work2(42, 1) + cnf(2)*work2(32, 1)
          work2(42, 1) = work2(42, 1) + cnf(4)*work2(27, 1)
          work2(43, 1) = work2(43, 1) + cnf(1)*work2(33, 1)
          work2(44, 1) = work2(44, 1) + cnf(2)*work2(33, 1)
          work2(45, 1) = work2(45, 1) + cnf(3)*work2(33, 1)
          work2(45, 1) = work2(45, 1) + cnf(4)*work2(27, 1)
          work2(28, 1) = work2(28, 1) + cnf(1)*work2(25, 1)
          work2(29, 1) = work2(29, 1) + cnf(2)*work2(25, 1)
          work2(30, 1) = work2(30, 1) + cnf(3)*work2(25, 1)
          work2(31, 1) = work2(31, 1) + cnf(1)*work2(26, 1)
          work2(32, 1) = work2(32, 1) + cnf(2)*work2(26, 1)
          work2(33, 1) = work2(33, 1) + cnf(3)*work2(26, 1)
          work2(15, 1) = work2(15, 1) + cnf(1)*work2(9, 1)
          work2(15, 1) = work2(15, 1) + cnf(4)*work2(6, 1)
          work2(16, 1) = work2(16, 1) + cnf(1)*work2(10, 1)
          work2(17, 1) = work2(17, 1) + cnf(2)*work2(10, 1)
          work2(17, 1) = work2(17, 1) + cnf(4)*work2(6, 1)
          work2(18, 1) = work2(18, 1) + cnf(2)*work2(11, 1)
          work2(18, 1) = work2(18, 1) + cnf(4)*work2(7, 1)
          work2(19, 1) = work2(19, 1) + cnf(1)*work2(12, 1)
          work2(20, 1) = work2(20, 1) + cnf(2)*work2(12, 1)
          work2(21, 1) = work2(21, 1) + cnf(2)*work2(13, 1)
          work2(22, 1) = work2(22, 1) + cnf(3)*work2(12, 1)
          work2(22, 1) = work2(22, 1) + cnf(4)*work2(6, 1)
          work2(23, 1) = work2(23, 1) + cnf(3)*work2(13, 1)
          work2(23, 1) = work2(23, 1) + cnf(4)*work2(7, 1)
          work2(24, 1) = work2(24, 1) + cnf(3)*work2(14, 1)
          work2(24, 1) = work2(24, 1) + cnf(4)*work2(8, 1)
          work2(9, 1) = work2(9, 1) + cnf(1)*work2(3, 1)
          work2(10, 1) = work2(10, 1) + cnf(2)*work2(3, 1)
          work2(11, 1) = work2(11, 1) + cnf(2)*work2(4, 1)
          work2(12, 1) = work2(12, 1) + cnf(3)*work2(3, 1)
          work2(13, 1) = work2(13, 1) + cnf(3)*work2(4, 1)
          work2(14, 1) = work2(14, 1) + cnf(3)*work2(5, 1)
          work2(15, 1) = work2(15, 1) + cnf(1)*work2(9, 1)
          work2(16, 1) = work2(16, 1) + cnf(2)*work2(9, 1)
          work2(17, 1) = work2(17, 1) + cnf(2)*work2(10, 1)
          work2(18, 1) = work2(18, 1) + cnf(2)*work2(11, 1)
          work2(19, 1) = work2(19, 1) + cnf(3)*work2(9, 1)
          work2(20, 1) = work2(20, 1) + cnf(3)*work2(10, 1)
          work2(21, 1) = work2(21, 1) + cnf(3)*work2(11, 1)
          work2(22, 1) = work2(22, 1) + cnf(3)*work2(12, 1)
          work2(23, 1) = work2(23, 1) + cnf(3)*work2(13, 1)
          work2(24, 1) = work2(24, 1) + cnf(3)*work2(14, 1)
          work2(56, 1) = work2(56, 1) + cnf(1)*work2(46, 1)
          work2(56, 1) = work2(56, 1) + cnf(5)*work2(40, 1)
          work2(57, 1) = work2(57, 1) + cnf(1)*work2(47, 1)
          work2(57, 1) = work2(57, 1) + cnf(4)*work2(41, 1)
          work2(58, 1) = work2(58, 1) + cnf(1)*work2(48, 1)
          work2(59, 1) = work2(59, 1) + cnf(2)*work2(48, 1)
          work2(59, 1) = work2(59, 1) + cnf(5)*work2(41, 1)
          work2(60, 1) = work2(60, 1) + cnf(2)*work2(49, 1)
          work2(60, 1) = work2(60, 1) + cnf(5)*work2(42, 1)
          work2(61, 1) = work2(61, 1) + cnf(1)*work2(50, 1)
          work2(61, 1) = work2(61, 1) + cnf(4)*work2(43, 1)
          work2(62, 1) = work2(62, 1) + cnf(1)*work2(51, 1)
          work2(63, 1) = work2(63, 1) + cnf(2)*work2(51, 1)
          work2(63, 1) = work2(63, 1) + cnf(4)*work2(43, 1)
          work2(64, 1) = work2(64, 1) + cnf(2)*work2(52, 1)
          work2(64, 1) = work2(64, 1) + cnf(4)*work2(44, 1)
          work2(65, 1) = work2(65, 1) + cnf(1)*work2(53, 1)
          work2(66, 1) = work2(66, 1) + cnf(2)*work2(53, 1)
          work2(67, 1) = work2(67, 1) + cnf(2)*work2(54, 1)
          work2(68, 1) = work2(68, 1) + cnf(3)*work2(53, 1)
          work2(68, 1) = work2(68, 1) + cnf(5)*work2(43, 1)
          work2(69, 1) = work2(69, 1) + cnf(3)*work2(54, 1)
          work2(69, 1) = work2(69, 1) + cnf(5)*work2(44, 1)
          work2(70, 1) = work2(70, 1) + cnf(3)*work2(55, 1)
          work2(70, 1) = work2(70, 1) + cnf(5)*work2(45, 1)
          work2(46, 1) = work2(46, 1) + cnf(1)*work2(34, 1)
          work2(46, 1) = work2(46, 1) + cnf(4)*work2(31, 1)
          work2(47, 1) = work2(47, 1) + cnf(1)*work2(35, 1)
          work2(48, 1) = work2(48, 1) + cnf(2)*work2(35, 1)
          work2(48, 1) = work2(48, 1) + cnf(4)*work2(31, 1)
          work2(49, 1) = work2(49, 1) + cnf(2)*work2(36, 1)
          work2(49, 1) = work2(49, 1) + cnf(4)*work2(32, 1)
          work2(50, 1) = work2(50, 1) + cnf(1)*work2(37, 1)
          work2(51, 1) = work2(51, 1) + cnf(2)*work2(37, 1)
          work2(52, 1) = work2(52, 1) + cnf(2)*work2(38, 1)
          work2(53, 1) = work2(53, 1) + cnf(3)*work2(37, 1)
          work2(53, 1) = work2(53, 1) + cnf(4)*work2(31, 1)
          work2(54, 1) = work2(54, 1) + cnf(3)*work2(38, 1)
          work2(54, 1) = work2(54, 1) + cnf(4)*work2(32, 1)
          work2(55, 1) = work2(55, 1) + cnf(3)*work2(39, 1)
          work2(55, 1) = work2(55, 1) + cnf(4)*work2(33, 1)
          work2(34, 1) = work2(34, 1) + cnf(1)*work2(28, 1)
          work2(35, 1) = work2(35, 1) + cnf(2)*work2(28, 1)
          work2(36, 1) = work2(36, 1) + cnf(2)*work2(29, 1)
          work2(37, 1) = work2(37, 1) + cnf(3)*work2(28, 1)
          work2(38, 1) = work2(38, 1) + cnf(3)*work2(29, 1)
          work2(39, 1) = work2(39, 1) + cnf(3)*work2(30, 1)
          work2(40, 1) = work2(40, 1) + cnf(1)*work2(31, 1)
          work2(41, 1) = work2(41, 1) + cnf(2)*work2(31, 1)
          work2(42, 1) = work2(42, 1) + cnf(2)*work2(32, 1)
          work2(43, 1) = work2(43, 1) + cnf(3)*work2(31, 1)
          work2(44, 1) = work2(44, 1) + cnf(3)*work2(32, 1)
          work2(45, 1) = work2(45, 1) + cnf(3)*work2(33, 1)
          work2(56, 1) = work2(56, 1) + cnf(1)*work2(46, 1)
          work2(56, 1) = work2(56, 1) + cnf(4)*work2(40, 1)
          work2(57, 1) = work2(57, 1) + cnf(1)*work2(47, 1)
          work2(58, 1) = work2(58, 1) + cnf(2)*work2(47, 1)
          work2(58, 1) = work2(58, 1) + cnf(4)*work2(40, 1)
          work2(59, 1) = work2(59, 1) + cnf(2)*work2(48, 1)
          work2(59, 1) = work2(59, 1) + cnf(4)*work2(41, 1)
          work2(60, 1) = work2(60, 1) + cnf(2)*work2(49, 1)
          work2(60, 1) = work2(60, 1) + cnf(4)*work2(42, 1)
          work2(61, 1) = work2(61, 1) + cnf(1)*work2(50, 1)
          work2(62, 1) = work2(62, 1) + cnf(2)*work2(50, 1)
          work2(63, 1) = work2(63, 1) + cnf(2)*work2(51, 1)
          work2(64, 1) = work2(64, 1) + cnf(2)*work2(52, 1)
          work2(65, 1) = work2(65, 1) + cnf(3)*work2(50, 1)
          work2(65, 1) = work2(65, 1) + cnf(4)*work2(40, 1)
          work2(66, 1) = work2(66, 1) + cnf(3)*work2(51, 1)
          work2(66, 1) = work2(66, 1) + cnf(4)*work2(41, 1)
          work2(67, 1) = work2(67, 1) + cnf(3)*work2(52, 1)
          work2(67, 1) = work2(67, 1) + cnf(4)*work2(42, 1)
          work2(68, 1) = work2(68, 1) + cnf(3)*work2(53, 1)
          work2(68, 1) = work2(68, 1) + cnf(4)*work2(43, 1)
          work2(69, 1) = work2(69, 1) + cnf(3)*work2(54, 1)
          work2(69, 1) = work2(69, 1) + cnf(4)*work2(44, 1)
          work2(70, 1) = work2(70, 1) + cnf(3)*work2(55, 1)
          work2(70, 1) = work2(70, 1) + cnf(4)*work2(45, 1)
          work2(46, 1) = work2(46, 1) + cnf(1)*work2(34, 1)
          work2(47, 1) = work2(47, 1) + cnf(2)*work2(34, 1)
          work2(48, 1) = work2(48, 1) + cnf(2)*work2(35, 1)
          work2(49, 1) = work2(49, 1) + cnf(2)*work2(36, 1)
          work2(50, 1) = work2(50, 1) + cnf(3)*work2(34, 1)
          work2(51, 1) = work2(51, 1) + cnf(3)*work2(35, 1)
          work2(52, 1) = work2(52, 1) + cnf(3)*work2(36, 1)
          work2(53, 1) = work2(53, 1) + cnf(3)*work2(37, 1)
          work2(54, 1) = work2(54, 1) + cnf(3)*work2(38, 1)
          work2(55, 1) = work2(55, 1) + cnf(3)*work2(39, 1)
          work2(56, 1) = work2(56, 1) + cnf(1)*work2(46, 1)
          work2(57, 1) = work2(57, 1) + cnf(2)*work2(46, 1)
          work2(58, 1) = work2(58, 1) + cnf(2)*work2(47, 1)
          work2(59, 1) = work2(59, 1) + cnf(2)*work2(48, 1)
          work2(60, 1) = work2(60, 1) + cnf(2)*work2(49, 1)
          work2(61, 1) = work2(61, 1) + cnf(3)*work2(46, 1)
          work2(62, 1) = work2(62, 1) + cnf(3)*work2(47, 1)
          work2(63, 1) = work2(63, 1) + cnf(3)*work2(48, 1)
          work2(64, 1) = work2(64, 1) + cnf(3)*work2(49, 1)
          work2(65, 1) = work2(65, 1) + cnf(3)*work2(50, 1)
          work2(66, 1) = work2(66, 1) + cnf(3)*work2(51, 1)
          work2(67, 1) = work2(67, 1) + cnf(3)*work2(52, 1)
          work2(68, 1) = work2(68, 1) + cnf(3)*work2(53, 1)
          work2(69, 1) = work2(69, 1) + cnf(3)*work2(54, 1)
          work2(70, 1) = work2(70, 1) + cnf(3)*work2(55, 1)
          work2(71, 1) = work2(56, 1) - cnf(1)*work2(15, 1)
          work2(72, 1) = work2(57, 1) - cnf(1)*work2(16, 1)
          work2(73, 1) = work2(58, 1) - cnf(1)*work2(17, 1)
          work2(74, 1) = work2(59, 1) - cnf(1)*work2(18, 1)
          work2(75, 1) = work2(61, 1) - cnf(1)*work2(19, 1)
          work2(76, 1) = work2(62, 1) - cnf(1)*work2(20, 1)
          work2(77, 1) = work2(63, 1) - cnf(1)*work2(21, 1)
          work2(78, 1) = work2(65, 1) - cnf(1)*work2(22, 1)
          work2(79, 1) = work2(66, 1) - cnf(1)*work2(23, 1)
          work2(80, 1) = work2(68, 1) - cnf(1)*work2(24, 1)
          work2(81, 1) = work2(57, 1) - cnf(2)*work2(15, 1)
          work2(82, 1) = work2(58, 1) - cnf(2)*work2(16, 1)
          work2(83, 1) = work2(59, 1) - cnf(2)*work2(17, 1)
          work2(84, 1) = work2(60, 1) - cnf(2)*work2(18, 1)
          work2(85, 1) = work2(62, 1) - cnf(2)*work2(19, 1)
          work2(86, 1) = work2(63, 1) - cnf(2)*work2(20, 1)
          work2(87, 1) = work2(64, 1) - cnf(2)*work2(21, 1)
          work2(88, 1) = work2(66, 1) - cnf(2)*work2(22, 1)
          work2(89, 1) = work2(67, 1) - cnf(2)*work2(23, 1)
          work2(90, 1) = work2(69, 1) - cnf(2)*work2(24, 1)
          work2(91, 1) = work2(61, 1) - cnf(3)*work2(15, 1)
          work2(92, 1) = work2(62, 1) - cnf(3)*work2(16, 1)
          work2(93, 1) = work2(63, 1) - cnf(3)*work2(17, 1)
          work2(94, 1) = work2(64, 1) - cnf(3)*work2(18, 1)
          work2(95, 1) = work2(65, 1) - cnf(3)*work2(19, 1)
          work2(96, 1) = work2(66, 1) - cnf(3)*work2(20, 1)
          work2(97, 1) = work2(67, 1) - cnf(3)*work2(21, 1)
          work2(98, 1) = work2(68, 1) - cnf(3)*work2(22, 1)
          work2(99, 1) = work2(69, 1) - cnf(3)*work2(23, 1)
          work2(100, 1) = work2(70, 1) - cnf(3)*work2(24, 1)

! ************************
! *                      *
! * Form final integrals *
! *                      *
! ************************

          eri_value(1) = work2(71, 1)
          eri_value(2) = work2(81, 1)
          eri_value(3) = work2(91, 1)
          eri_value(4) = work2(74, 1)
          eri_value(5) = work2(84, 1)
          eri_value(6) = work2(94, 1)
          eri_value(7) = work2(80, 1)
          eri_value(8) = work2(90, 1)
          eri_value(9) = work2(100, 1)
          eri_value(10) = work2(72, 1)*sqrt5
          eri_value(11) = work2(82, 1)*sqrt5
          eri_value(12) = work2(92, 1)*sqrt5
          eri_value(13) = work2(75, 1)*sqrt5
          eri_value(14) = work2(85, 1)*sqrt5
          eri_value(15) = work2(95, 1)*sqrt5
          eri_value(16) = work2(73, 1)*sqrt5
          eri_value(17) = work2(83, 1)*sqrt5
          eri_value(18) = work2(93, 1)*sqrt5
          eri_value(19) = work2(77, 1)*sqrt5
          eri_value(20) = work2(87, 1)*sqrt5
          eri_value(21) = work2(97, 1)*sqrt5
          eri_value(22) = work2(78, 1)*sqrt5
          eri_value(23) = work2(88, 1)*sqrt5
          eri_value(24) = work2(98, 1)*sqrt5
          eri_value(25) = work2(79, 1)*sqrt5
          eri_value(26) = work2(89, 1)*sqrt5
          eri_value(27) = work2(99, 1)*sqrt5
          eri_value(28) = work2(76, 1)*sqrt15
          eri_value(29) = work2(86, 1)*sqrt15
          eri_value(30) = work2(96, 1)*sqrt15

          maxk = 10
          maxl = 3

! ******************************
! *                            *
! * Digestion into Fock matrix *
! *                            *
! ******************************

          loci = res%atom_loc(ish) - 1
          locj = res%atom_loc(jsh) - 1
          lock = res%atom_loc(ksh) - 1
          locl = res%atom_loc(lsh) - 1

          nij = 0

          lstride = 1
          kstride = 3

          ip = 1
          i = 1
          ii1 = i + loci
          ijp = ip

          j = 1
          nij = nij + 1
          maxl2 = maxl
          jj1 = j + locj
          i2 = ii1
          j2 = jj1
          if (ii1 .lt. jj1) then ! sort <ij|
            i2 = jj1
            j2 = ii1
          end if

          ijkp = ijp
          nkl = nij

          do k = 1, maxk
            kk1 = k + lock
            kk1 = k + lock

            ijklp = ijkp
            ijkp = ijkp + kstride

            do l = 1, maxl
              buff(1) = eri_value(ijklp)
              ijklp = ijklp + lstride
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
              !   Account for identical permutations
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
          !
!

        end if ! test.gt.cutoff_schwarz
      end do ! iquart
!$omp end target teams distribute parallel do

    end do ! itile

  end subroutine int0031
end submodule
