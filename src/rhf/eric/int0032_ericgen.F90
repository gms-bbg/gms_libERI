! The total angular momentum of this class is:           5
! The algorithm chosen is: PHR a.k.a. ERIC
! Writing an ERIC kernel
submodule(eric_kernels) int0032_impl
contains
  module subroutine int0032(ss_pair, df_pair, density, fock, res)

    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, df_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

! Variables for the class
    integer(kind=int64), allocatable :: n00bra(:), n23ket(:)
    real(dp), allocatable :: xint00bra(:), xint23ket(:)
    integer(kind=int64) :: nssbra, ndfket
    real(dp) :: scutssbra, scutdfket, test
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
    real(dp) :: ft(9), phi(354), c_factor(56), rxyz(3)
    real(dp) :: work2(285, 1)
    real(dp) :: eri_value(60), angl(20)
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
    allocate (n23ket(res%n_d_shl*res%n_f_shl))
    allocate (xint23ket(res%n_d_shl*res%n_f_shl))

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

    scutdfket = cutoff_schwarz/maxval(df_pair%xints)
    ndfket = 0
    do ij = 1, res%n_d_shl*res%n_f_shl
      if (df_pair%xints(ij) .ge. scutdfket) then
        ndfket = ndfket + 1
        xint23ket(ndfket) = df_pair%xints(ij)
        n23ket(ndfket) = ij
      end if
    end do

    nchunksize_k10 = 375000000

    if ((nssbra*ndfket) .le. nchunksize_k10) nchunksize_k10 = nssbra*ndfket
    ntile = int(nssbra*ndfket/nchunksize_k10)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_k10 + 1
      iend = itile*nchunksize_k10
      if (itile .eq. ntile) iend = nssbra*ndfket

! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

! Mappings to GPU

!$omp target teams distribute parallel do default(none) &
!$omp shared(res, density, fock, nquart_start, nquart_end, nssbra, xint00bra, n00bra, xint23ket, n23ket, ndfket, df_pair, ss_pair) &
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

        ij_tmp = (iquart - 1)/ndfket + 1
        kl_tmp = mod(iquart - 1, ndfket) + 1

        test = xint00bra(ij_tmp)*xint23ket(kl_tmp)

        if (test .gt. cutoff_schwarz) then

          ij = n00bra(ij_tmp)
          kl = n23ket(kl_tmp)

          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_f_shl + 1
          lsh_tmp = mod(kl - 1, res%n_f_shl) + 1

          ish = res%i_s_shl(ish_tmp)
          jsh = res%i_s_shl(jsh_tmp)
          ksh = res%i_f_shl(lsh_tmp)
          lsh = res%i_d_shl(ksh_tmp)

          ijtop = res%contr_num(ish)*res%contr_num(jsh)

          kltop = res%contr_num(ksh)*res%contr_num(lsh)

          eri_value = 0.0_dp
          ket_loop = 0

          ! 0.0_dp out phi elements to contract over kl
          ! PHI = "pre-Hermite integrals"

          phi(82) = 0.0_dp
          phi(83) = 0.0_dp
          phi(84) = 0.0_dp
          phi(85) = 0.0_dp
          phi(86) = 0.0_dp
          phi(87) = 0.0_dp
          phi(88) = 0.0_dp
          phi(89) = 0.0_dp
          phi(93) = 0.0_dp
          phi(94) = 0.0_dp
          phi(95) = 0.0_dp
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
          phi(109) = 0.0_dp
          phi(110) = 0.0_dp
          phi(111) = 0.0_dp
          phi(112) = 0.0_dp
          phi(113) = 0.0_dp
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
          phi(179) = 0.0_dp
          phi(180) = 0.0_dp
          phi(181) = 0.0_dp
          phi(182) = 0.0_dp
          phi(183) = 0.0_dp
          phi(184) = 0.0_dp
          phi(185) = 0.0_dp
          phi(186) = 0.0_dp
          phi(187) = 0.0_dp
          phi(188) = 0.0_dp
          phi(189) = 0.0_dp
          phi(190) = 0.0_dp
          phi(191) = 0.0_dp
          phi(192) = 0.0_dp
          phi(193) = 0.0_dp
          phi(194) = 0.0_dp
          phi(195) = 0.0_dp
          phi(196) = 0.0_dp
          phi(197) = 0.0_dp
          phi(198) = 0.0_dp
          phi(199) = 0.0_dp
          phi(200) = 0.0_dp
          phi(201) = 0.0_dp
          phi(202) = 0.0_dp
          phi(203) = 0.0_dp
          phi(204) = 0.0_dp
          phi(205) = 0.0_dp
          phi(206) = 0.0_dp
          phi(207) = 0.0_dp
          phi(208) = 0.0_dp
          phi(209) = 0.0_dp
          phi(210) = 0.0_dp
          phi(211) = 0.0_dp
          phi(212) = 0.0_dp
          phi(213) = 0.0_dp
          phi(214) = 0.0_dp
          phi(215) = 0.0_dp
          phi(216) = 0.0_dp
          phi(217) = 0.0_dp
          phi(218) = 0.0_dp
          phi(219) = 0.0_dp
          phi(220) = 0.0_dp
          phi(243) = 0.0_dp
          phi(244) = 0.0_dp
          phi(245) = 0.0_dp
          phi(246) = 0.0_dp
          phi(247) = 0.0_dp
          phi(248) = 0.0_dp
          phi(249) = 0.0_dp
          phi(250) = 0.0_dp
          phi(251) = 0.0_dp
          phi(252) = 0.0_dp
          phi(253) = 0.0_dp
          phi(254) = 0.0_dp
          phi(255) = 0.0_dp
          phi(256) = 0.0_dp
          phi(257) = 0.0_dp
          phi(258) = 0.0_dp
          phi(259) = 0.0_dp
          phi(260) = 0.0_dp
          phi(261) = 0.0_dp
          phi(262) = 0.0_dp
          phi(263) = 0.0_dp
          phi(264) = 0.0_dp
          phi(265) = 0.0_dp
          phi(266) = 0.0_dp
          phi(267) = 0.0_dp
          phi(268) = 0.0_dp
          phi(269) = 0.0_dp
          phi(270) = 0.0_dp
          phi(271) = 0.0_dp
          phi(272) = 0.0_dp
          phi(273) = 0.0_dp
          phi(274) = 0.0_dp
          phi(275) = 0.0_dp
          phi(276) = 0.0_dp
          phi(277) = 0.0_dp
          phi(278) = 0.0_dp
          phi(279) = 0.0_dp
          phi(280) = 0.0_dp
          phi(281) = 0.0_dp
          phi(282) = 0.0_dp
          phi(283) = 0.0_dp
          phi(284) = 0.0_dp
          phi(285) = 0.0_dp
          phi(286) = 0.0_dp
          phi(321) = 0.0_dp
          phi(322) = 0.0_dp
          phi(323) = 0.0_dp
          phi(324) = 0.0_dp
          phi(325) = 0.0_dp
          phi(326) = 0.0_dp
          phi(327) = 0.0_dp
          phi(328) = 0.0_dp
          phi(329) = 0.0_dp
          phi(330) = 0.0_dp
          phi(331) = 0.0_dp
          phi(332) = 0.0_dp
          phi(333) = 0.0_dp
          phi(334) = 0.0_dp
          phi(335) = 0.0_dp
          phi(336) = 0.0_dp
          phi(337) = 0.0_dp
          phi(338) = 0.0_dp
          phi(339) = 0.0_dp
          phi(340) = 0.0_dp
          phi(341) = 0.0_dp
          phi(342) = 0.0_dp
          phi(343) = 0.0_dp
          phi(344) = 0.0_dp
          phi(345) = 0.0_dp
          phi(346) = 0.0_dp
          phi(347) = 0.0_dp
          phi(348) = 0.0_dp
          phi(349) = 0.0_dp
          phi(350) = 0.0_dp
          phi(351) = 0.0_dp
          phi(352) = 0.0_dp
          phi(353) = 0.0_dp
          phi(354) = 0.0_dp

! Begin looping over kl primitives

          do k = 1, kltop

            t_expon_cd = df_pair%t_expon_ab(df_pair%pair_loc(kl) + 1 + ket_loop) ! exp_c + exp_d
            t_expon_c = df_pair%expon_b(df_pair%pair_loc(kl) + 1 + ket_loop)
            t_expon_d = df_pair%expon_a(df_pair%pair_loc(kl) + 1 + ket_loop)

            t_inverse_expon_cd = 1.0_dp/t_expon_cd
            ccfket = df_pair%sq(df_pair%pair_loc(kl) + 1 + ket_loop)*sqrt2_pi_5_4
            slket = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_cd)

            xkl = ((t_expon_c*res%coord_sh(ksh, 1)) + (t_expon_d*res%coord_sh(lsh, 1)))*t_inverse_expon_cd
            ykl = ((t_expon_c*res%coord_sh(ksh, 2)) + (t_expon_d*res%coord_sh(lsh, 2)))*t_inverse_expon_cd
            zkl = ((t_expon_c*res%coord_sh(ksh, 3)) + (t_expon_d*res%coord_sh(lsh, 3)))*t_inverse_expon_cd

            rxket = df_pair%t_inverse_expon_ab(df_pair%pair_loc(kl) + 1 + ket_loop)! inverse_expon_ab*0.5_dp

            ket_loop = ket_loop + 1
            bra_loop = 0

            ! 0.0_dp out phi elements to contract over ij

            phi(81) = 0.0_dp
            phi(90) = 0.0_dp
            phi(91) = 0.0_dp
            phi(92) = 0.0_dp
            phi(114) = 0.0_dp
            phi(115) = 0.0_dp
            phi(116) = 0.0_dp
            phi(117) = 0.0_dp
            phi(118) = 0.0_dp
            phi(119) = 0.0_dp
            phi(120) = 0.0_dp
            phi(156) = 0.0_dp
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
            phi(221) = 0.0_dp
            phi(222) = 0.0_dp
            phi(223) = 0.0_dp
            phi(224) = 0.0_dp
            phi(225) = 0.0_dp
            phi(226) = 0.0_dp
            phi(227) = 0.0_dp
            phi(228) = 0.0_dp
            phi(229) = 0.0_dp
            phi(230) = 0.0_dp
            phi(231) = 0.0_dp
            phi(232) = 0.0_dp
            phi(233) = 0.0_dp
            phi(234) = 0.0_dp
            phi(235) = 0.0_dp
            phi(236) = 0.0_dp
            phi(237) = 0.0_dp
            phi(238) = 0.0_dp
            phi(239) = 0.0_dp
            phi(240) = 0.0_dp
            phi(241) = 0.0_dp
            phi(242) = 0.0_dp
            phi(287) = 0.0_dp
            phi(288) = 0.0_dp
            phi(289) = 0.0_dp
            phi(290) = 0.0_dp
            phi(291) = 0.0_dp
            phi(292) = 0.0_dp
            phi(293) = 0.0_dp
            phi(294) = 0.0_dp
            phi(295) = 0.0_dp
            phi(296) = 0.0_dp
            phi(297) = 0.0_dp
            phi(298) = 0.0_dp
            phi(299) = 0.0_dp
            phi(300) = 0.0_dp
            phi(301) = 0.0_dp
            phi(302) = 0.0_dp
            phi(303) = 0.0_dp
            phi(304) = 0.0_dp
            phi(305) = 0.0_dp
            phi(306) = 0.0_dp
            phi(307) = 0.0_dp
            phi(308) = 0.0_dp
            phi(309) = 0.0_dp
            phi(310) = 0.0_dp
            phi(311) = 0.0_dp
            phi(312) = 0.0_dp
            phi(313) = 0.0_dp
            phi(314) = 0.0_dp
            phi(315) = 0.0_dp
            phi(316) = 0.0_dp
            phi(317) = 0.0_dp
            phi(318) = 0.0_dp
            phi(319) = 0.0_dp
            phi(320) = 0.0_dp

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

              n = 5
              if (tt .le. t_max) then

                ! Boys function evaluation with Chebyshev interpolation

                t_new = tt*0.1432309939656232D+02
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

                ft(9) = fmt
                t2 = tt + tt
                ft(8) = (t2*ft(9) + expt)*0.666666666666667D-01
                ft(7) = (t2*ft(8) + expt)*0.769230769230769D-01
                ft(6) = (t2*ft(7) + expt)*0.909090909090909D-01
                ft(5) = (t2*ft(6) + expt)*0.111111111111111D+00
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
                ftf = ftf*rho
                ft(6) = ft(6)*ftf

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
                ftf = ftf*xin
                ft(6) = ftf*0.945000000000000D+03
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
! Still need to figure these out :P
              c_factor(36) = c_factor(21)*rxyz(1)
              c_factor(37) = c_factor(22)*rxyz(1)
              c_factor(38) = c_factor(23)*rxyz(1)
              c_factor(39) = c_factor(24)*rxyz(1)
              c_factor(40) = c_factor(25)*rxyz(1)
              c_factor(41) = c_factor(25)*rxyz(2)
              c_factor(42) = c_factor(26)*rxyz(1)
              c_factor(43) = c_factor(27)*rxyz(1)
              c_factor(44) = c_factor(28)*rxyz(1)
              c_factor(45) = c_factor(29)*rxyz(1)
              c_factor(46) = c_factor(29)*rxyz(2)
              c_factor(47) = c_factor(30)*rxyz(1)
              c_factor(48) = c_factor(31)*rxyz(1)
              c_factor(49) = c_factor(32)*rxyz(1)
              c_factor(50) = c_factor(32)*rxyz(2)
              c_factor(51) = c_factor(33)*rxyz(1)
              c_factor(52) = c_factor(34)*rxyz(1)
              c_factor(53) = c_factor(34)*rxyz(2)
              c_factor(54) = c_factor(35)*rxyz(1)
              c_factor(55) = c_factor(35)*rxyz(2)
              c_factor(56) = c_factor(35)*rxyz(3)
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
              phi(47) = c_factor(2)*ft(4)
              phi(48) = c_factor(3)*ft(4)
              phi(49) = c_factor(4)*ft(4)
              phi(50) = c_factor(11)*ft(5)
              phi(51) = c_factor(12)*ft(5)
              phi(52) = c_factor(13)*ft(5)
              phi(53) = c_factor(14)*ft(5)
              phi(54) = c_factor(15)*ft(5)
              phi(55) = c_factor(16)*ft(5)
              phi(56) = c_factor(17)*ft(5)
              phi(57) = c_factor(18)*ft(5)
              phi(58) = c_factor(19)*ft(5)
              phi(59) = c_factor(20)*ft(5)
              phi(60) = c_factor(36)*ft(6)
              phi(61) = c_factor(37)*ft(6)
              phi(62) = c_factor(38)*ft(6)
              phi(63) = c_factor(39)*ft(6)
              phi(64) = c_factor(40)*ft(6)
              phi(65) = c_factor(41)*ft(6)
              phi(66) = c_factor(42)*ft(6)
              phi(67) = c_factor(43)*ft(6)
              phi(68) = c_factor(44)*ft(6)
              phi(69) = c_factor(45)*ft(6)
              phi(70) = c_factor(46)*ft(6)
              phi(71) = c_factor(47)*ft(6)
              phi(72) = c_factor(48)*ft(6)
              phi(73) = c_factor(49)*ft(6)
              phi(74) = c_factor(50)*ft(6)
              phi(75) = c_factor(51)*ft(6)
              phi(76) = c_factor(52)*ft(6)
              phi(77) = c_factor(53)*ft(6)
              phi(78) = c_factor(54)*ft(6)
              phi(79) = c_factor(55)*ft(6)
              phi(80) = c_factor(56)*ft(6)

              ! Begin contraction & scaling over ij
              ! i shell scaling factors
              fac1 = 1.0_dp
              scale_factor1(1) = 1.0_dp
              ! j shell scaling factors
              fac2 = 1.0_dp
              scale_factor2(1) = 1.0_dp

              ! ij contraction
              phi(81) = phi(81) + phi(1)

              phi(90) = phi(90) + phi(2)
              phi(91) = phi(91) + phi(3)
              phi(92) = phi(92) + phi(4)

              phi(114) = phi(114) + phi(5)
              phi(115) = phi(115) + phi(6)
              phi(116) = phi(116) + phi(7)
              phi(117) = phi(117) + phi(8)
              phi(118) = phi(118) + phi(9)
              phi(119) = phi(119) + phi(10)
              phi(120) = phi(120) + phi(11)

              phi(156) = phi(156) + phi(12)
              phi(157) = phi(157) + phi(13)
              phi(158) = phi(158) + phi(14)
              phi(159) = phi(159) + phi(15)
              phi(160) = phi(160) + phi(16)
              phi(161) = phi(161) + phi(17)
              phi(162) = phi(162) + phi(18)
              phi(163) = phi(163) + phi(19)
              phi(164) = phi(164) + phi(20)
              phi(165) = phi(165) + phi(21)
              phi(166) = phi(166) + phi(22)
              phi(167) = phi(167) + phi(23)
              phi(168) = phi(168) + phi(24)

              phi(221) = phi(221) + phi(25)
              phi(222) = phi(222) + phi(26)
              phi(223) = phi(223) + phi(27)
              phi(224) = phi(224) + phi(28)
              phi(225) = phi(225) + phi(29)
              phi(226) = phi(226) + phi(30)
              phi(227) = phi(227) + phi(31)
              phi(228) = phi(228) + phi(32)
              phi(229) = phi(229) + phi(33)
              phi(230) = phi(230) + phi(34)
              phi(231) = phi(231) + phi(35)
              phi(232) = phi(232) + phi(36)
              phi(233) = phi(233) + phi(37)
              phi(234) = phi(234) + phi(38)
              phi(235) = phi(235) + phi(39)
              phi(236) = phi(236) + phi(40)
              phi(237) = phi(237) + phi(41)
              phi(238) = phi(238) + phi(42)
              phi(239) = phi(239) + phi(43)
              phi(240) = phi(240) + phi(44)
              phi(241) = phi(241) + phi(45)
              phi(242) = phi(242) + phi(46)

              phi(287) = phi(287) + phi(47)
              phi(288) = phi(288) + phi(48)
              phi(289) = phi(289) + phi(49)
              phi(290) = phi(290) + phi(50)
              phi(291) = phi(291) + phi(51)
              phi(292) = phi(292) + phi(52)
              phi(293) = phi(293) + phi(53)
              phi(294) = phi(294) + phi(54)
              phi(295) = phi(295) + phi(55)
              phi(296) = phi(296) + phi(56)
              phi(297) = phi(297) + phi(57)
              phi(298) = phi(298) + phi(58)
              phi(299) = phi(299) + phi(59)
              phi(300) = phi(300) + phi(60)
              phi(301) = phi(301) + phi(61)
              phi(302) = phi(302) + phi(62)
              phi(303) = phi(303) + phi(63)
              phi(304) = phi(304) + phi(64)
              phi(305) = phi(305) + phi(65)
              phi(306) = phi(306) + phi(66)
              phi(307) = phi(307) + phi(67)
              phi(308) = phi(308) + phi(68)
              phi(309) = phi(309) + phi(69)
              phi(310) = phi(310) + phi(70)
              phi(311) = phi(311) + phi(71)
              phi(312) = phi(312) + phi(72)
              phi(313) = phi(313) + phi(73)
              phi(314) = phi(314) + phi(74)
              phi(315) = phi(315) + phi(75)
              phi(316) = phi(316) + phi(76)
              phi(317) = phi(317) + phi(77)
              phi(318) = phi(318) + phi(78)
              phi(319) = phi(319) + phi(79)
              phi(320) = phi(320) + phi(80)

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
            fac1 = fac1*rxket
            scale_factor1(6) = fac1 ! (0.5*inv_exp_cd)** 5
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
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(6) = fac2 ! (2*t_expon_d)** 5

            ! kl contraction
            phi(82) = phi(82) + scale_factor1(3)*phi(81)
            phi(83) = phi(83) + scale_factor2(2)*scale_factor1(3)*phi(81)
            phi(84) = phi(84) + scale_factor2(2)*scale_factor1(4)*phi(81)
            phi(85) = phi(85) + scale_factor2(3)*scale_factor1(4)*phi(81)
            phi(86) = phi(86) + scale_factor2(4)*scale_factor1(4)*phi(81)
            phi(87) = phi(87) + scale_factor2(4)*scale_factor1(5)*phi(81)
            phi(88) = phi(88) + scale_factor2(5)*scale_factor1(5)*phi(81)
            phi(89) = phi(89) + scale_factor2(6)*scale_factor1(6)*phi(81)

            phi(93) = phi(93) + scale_factor1(3)*phi(90)
            phi(94) = phi(94) + scale_factor1(3)*phi(91)
            phi(95) = phi(95) + scale_factor1(3)*phi(92)
            phi(96) = phi(96) + scale_factor1(4)*phi(90)
            phi(97) = phi(97) + scale_factor1(4)*phi(91)
            phi(98) = phi(98) + scale_factor1(4)*phi(92)
            phi(99) = phi(99) + scale_factor2(2)*scale_factor1(4)*phi(90)
            phi(100) = phi(100) + scale_factor2(2)*scale_factor1(4)*phi(91)
            phi(101) = phi(101) + scale_factor2(2)*scale_factor1(4)*phi(92)
            phi(102) = phi(102) + scale_factor2(3)*scale_factor1(4)*phi(90)
            phi(103) = phi(103) + scale_factor2(3)*scale_factor1(4)*phi(91)
            phi(104) = phi(104) + scale_factor2(3)*scale_factor1(4)*phi(92)
            phi(105) = phi(105) + scale_factor2(3)*scale_factor1(5)*phi(90)
            phi(106) = phi(106) + scale_factor2(3)*scale_factor1(5)*phi(91)
            phi(107) = phi(107) + scale_factor2(3)*scale_factor1(5)*phi(92)
            phi(108) = phi(108) + scale_factor2(4)*scale_factor1(5)*phi(90)
            phi(109) = phi(109) + scale_factor2(4)*scale_factor1(5)*phi(91)
            phi(110) = phi(110) + scale_factor2(4)*scale_factor1(5)*phi(92)
            phi(111) = phi(111) + scale_factor2(5)*scale_factor1(6)*phi(90)
            phi(112) = phi(112) + scale_factor2(5)*scale_factor1(6)*phi(91)
            phi(113) = phi(113) + scale_factor2(5)*scale_factor1(6)*phi(92)

            phi(121) = phi(121) + scale_factor1(4)*phi(114)
            phi(122) = phi(122) + scale_factor1(4)*phi(115)
            phi(123) = phi(123) + scale_factor1(4)*phi(116)
            phi(124) = phi(124) + scale_factor1(4)*phi(117)
            phi(125) = phi(125) + scale_factor1(4)*phi(118)
            phi(126) = phi(126) + scale_factor1(4)*phi(119)
            phi(127) = phi(127) + scale_factor1(4)*phi(120)
            phi(128) = phi(128) + scale_factor2(2)*scale_factor1(4)*phi(114)
            phi(129) = phi(129) + scale_factor2(2)*scale_factor1(4)*phi(115)
            phi(130) = phi(130) + scale_factor2(2)*scale_factor1(4)*phi(116)
            phi(131) = phi(131) + scale_factor2(2)*scale_factor1(4)*phi(117)
            phi(132) = phi(132) + scale_factor2(2)*scale_factor1(4)*phi(118)
            phi(133) = phi(133) + scale_factor2(2)*scale_factor1(4)*phi(119)
            phi(134) = phi(134) + scale_factor2(2)*scale_factor1(4)*phi(120)
            phi(135) = phi(135) + scale_factor2(2)*scale_factor1(5)*phi(114)
            phi(136) = phi(136) + scale_factor2(2)*scale_factor1(5)*phi(115)
            phi(137) = phi(137) + scale_factor2(2)*scale_factor1(5)*phi(116)
            phi(138) = phi(138) + scale_factor2(2)*scale_factor1(5)*phi(117)
            phi(139) = phi(139) + scale_factor2(2)*scale_factor1(5)*phi(118)
            phi(140) = phi(140) + scale_factor2(2)*scale_factor1(5)*phi(119)
            phi(141) = phi(141) + scale_factor2(2)*scale_factor1(5)*phi(120)
            phi(142) = phi(142) + scale_factor2(3)*scale_factor1(5)*phi(114)
            phi(143) = phi(143) + scale_factor2(3)*scale_factor1(5)*phi(115)
            phi(144) = phi(144) + scale_factor2(3)*scale_factor1(5)*phi(116)
            phi(145) = phi(145) + scale_factor2(3)*scale_factor1(5)*phi(117)
            phi(146) = phi(146) + scale_factor2(3)*scale_factor1(5)*phi(118)
            phi(147) = phi(147) + scale_factor2(3)*scale_factor1(5)*phi(119)
            phi(148) = phi(148) + scale_factor2(3)*scale_factor1(5)*phi(120)
            phi(149) = phi(149) + scale_factor2(4)*scale_factor1(6)*phi(114)
            phi(150) = phi(150) + scale_factor2(4)*scale_factor1(6)*phi(115)
            phi(151) = phi(151) + scale_factor2(4)*scale_factor1(6)*phi(116)
            phi(152) = phi(152) + scale_factor2(4)*scale_factor1(6)*phi(117)
            phi(153) = phi(153) + scale_factor2(4)*scale_factor1(6)*phi(118)
            phi(154) = phi(154) + scale_factor2(4)*scale_factor1(6)*phi(119)
            phi(155) = phi(155) + scale_factor2(4)*scale_factor1(6)*phi(120)

            phi(169) = phi(169) + scale_factor1(4)*phi(156)
            phi(170) = phi(170) + scale_factor1(4)*phi(157)
            phi(171) = phi(171) + scale_factor1(4)*phi(158)
            phi(172) = phi(172) + scale_factor1(4)*phi(159)
            phi(173) = phi(173) + scale_factor1(4)*phi(160)
            phi(174) = phi(174) + scale_factor1(4)*phi(161)
            phi(175) = phi(175) + scale_factor1(4)*phi(162)
            phi(176) = phi(176) + scale_factor1(4)*phi(163)
            phi(177) = phi(177) + scale_factor1(4)*phi(164)
            phi(178) = phi(178) + scale_factor1(4)*phi(165)
            phi(179) = phi(179) + scale_factor1(4)*phi(166)
            phi(180) = phi(180) + scale_factor1(4)*phi(167)
            phi(181) = phi(181) + scale_factor1(4)*phi(168)
            phi(182) = phi(182) + scale_factor1(5)*phi(156)
            phi(183) = phi(183) + scale_factor1(5)*phi(157)
            phi(184) = phi(184) + scale_factor1(5)*phi(158)
            phi(185) = phi(185) + scale_factor1(5)*phi(159)
            phi(186) = phi(186) + scale_factor1(5)*phi(160)
            phi(187) = phi(187) + scale_factor1(5)*phi(161)
            phi(188) = phi(188) + scale_factor1(5)*phi(162)
            phi(189) = phi(189) + scale_factor1(5)*phi(163)
            phi(190) = phi(190) + scale_factor1(5)*phi(164)
            phi(191) = phi(191) + scale_factor1(5)*phi(165)
            phi(192) = phi(192) + scale_factor1(5)*phi(166)
            phi(193) = phi(193) + scale_factor1(5)*phi(167)
            phi(194) = phi(194) + scale_factor1(5)*phi(168)
            phi(195) = phi(195) + scale_factor2(2)*scale_factor1(5)*phi(156)
            phi(196) = phi(196) + scale_factor2(2)*scale_factor1(5)*phi(157)
            phi(197) = phi(197) + scale_factor2(2)*scale_factor1(5)*phi(158)
            phi(198) = phi(198) + scale_factor2(2)*scale_factor1(5)*phi(159)
            phi(199) = phi(199) + scale_factor2(2)*scale_factor1(5)*phi(160)
            phi(200) = phi(200) + scale_factor2(2)*scale_factor1(5)*phi(161)
            phi(201) = phi(201) + scale_factor2(2)*scale_factor1(5)*phi(162)
            phi(202) = phi(202) + scale_factor2(2)*scale_factor1(5)*phi(163)
            phi(203) = phi(203) + scale_factor2(2)*scale_factor1(5)*phi(164)
            phi(204) = phi(204) + scale_factor2(2)*scale_factor1(5)*phi(165)
            phi(205) = phi(205) + scale_factor2(2)*scale_factor1(5)*phi(166)
            phi(206) = phi(206) + scale_factor2(2)*scale_factor1(5)*phi(167)
            phi(207) = phi(207) + scale_factor2(2)*scale_factor1(5)*phi(168)
            phi(208) = phi(208) + scale_factor2(3)*scale_factor1(6)*phi(156)
            phi(209) = phi(209) + scale_factor2(3)*scale_factor1(6)*phi(157)
            phi(210) = phi(210) + scale_factor2(3)*scale_factor1(6)*phi(158)
            phi(211) = phi(211) + scale_factor2(3)*scale_factor1(6)*phi(159)
            phi(212) = phi(212) + scale_factor2(3)*scale_factor1(6)*phi(160)
            phi(213) = phi(213) + scale_factor2(3)*scale_factor1(6)*phi(161)
            phi(214) = phi(214) + scale_factor2(3)*scale_factor1(6)*phi(162)
            phi(215) = phi(215) + scale_factor2(3)*scale_factor1(6)*phi(163)
            phi(216) = phi(216) + scale_factor2(3)*scale_factor1(6)*phi(164)
            phi(217) = phi(217) + scale_factor2(3)*scale_factor1(6)*phi(165)
            phi(218) = phi(218) + scale_factor2(3)*scale_factor1(6)*phi(166)
            phi(219) = phi(219) + scale_factor2(3)*scale_factor1(6)*phi(167)
            phi(220) = phi(220) + scale_factor2(3)*scale_factor1(6)*phi(168)

            phi(243) = phi(243) + scale_factor1(5)*phi(221)
            phi(244) = phi(244) + scale_factor1(5)*phi(222)
            phi(245) = phi(245) + scale_factor1(5)*phi(223)
            phi(246) = phi(246) + scale_factor1(5)*phi(224)
            phi(247) = phi(247) + scale_factor1(5)*phi(225)
            phi(248) = phi(248) + scale_factor1(5)*phi(226)
            phi(249) = phi(249) + scale_factor1(5)*phi(227)
            phi(250) = phi(250) + scale_factor1(5)*phi(228)
            phi(251) = phi(251) + scale_factor1(5)*phi(229)
            phi(252) = phi(252) + scale_factor1(5)*phi(230)
            phi(253) = phi(253) + scale_factor1(5)*phi(231)
            phi(254) = phi(254) + scale_factor1(5)*phi(232)
            phi(255) = phi(255) + scale_factor1(5)*phi(233)
            phi(256) = phi(256) + scale_factor1(5)*phi(234)
            phi(257) = phi(257) + scale_factor1(5)*phi(235)
            phi(258) = phi(258) + scale_factor1(5)*phi(236)
            phi(259) = phi(259) + scale_factor1(5)*phi(237)
            phi(260) = phi(260) + scale_factor1(5)*phi(238)
            phi(261) = phi(261) + scale_factor1(5)*phi(239)
            phi(262) = phi(262) + scale_factor1(5)*phi(240)
            phi(263) = phi(263) + scale_factor1(5)*phi(241)
            phi(264) = phi(264) + scale_factor1(5)*phi(242)
            phi(265) = phi(265) + scale_factor2(2)*scale_factor1(6)*phi(221)
            phi(266) = phi(266) + scale_factor2(2)*scale_factor1(6)*phi(222)
            phi(267) = phi(267) + scale_factor2(2)*scale_factor1(6)*phi(223)
            phi(268) = phi(268) + scale_factor2(2)*scale_factor1(6)*phi(224)
            phi(269) = phi(269) + scale_factor2(2)*scale_factor1(6)*phi(225)
            phi(270) = phi(270) + scale_factor2(2)*scale_factor1(6)*phi(226)
            phi(271) = phi(271) + scale_factor2(2)*scale_factor1(6)*phi(227)
            phi(272) = phi(272) + scale_factor2(2)*scale_factor1(6)*phi(228)
            phi(273) = phi(273) + scale_factor2(2)*scale_factor1(6)*phi(229)
            phi(274) = phi(274) + scale_factor2(2)*scale_factor1(6)*phi(230)
            phi(275) = phi(275) + scale_factor2(2)*scale_factor1(6)*phi(231)
            phi(276) = phi(276) + scale_factor2(2)*scale_factor1(6)*phi(232)
            phi(277) = phi(277) + scale_factor2(2)*scale_factor1(6)*phi(233)
            phi(278) = phi(278) + scale_factor2(2)*scale_factor1(6)*phi(234)
            phi(279) = phi(279) + scale_factor2(2)*scale_factor1(6)*phi(235)
            phi(280) = phi(280) + scale_factor2(2)*scale_factor1(6)*phi(236)
            phi(281) = phi(281) + scale_factor2(2)*scale_factor1(6)*phi(237)
            phi(282) = phi(282) + scale_factor2(2)*scale_factor1(6)*phi(238)
            phi(283) = phi(283) + scale_factor2(2)*scale_factor1(6)*phi(239)
            phi(284) = phi(284) + scale_factor2(2)*scale_factor1(6)*phi(240)
            phi(285) = phi(285) + scale_factor2(2)*scale_factor1(6)*phi(241)
            phi(286) = phi(286) + scale_factor2(2)*scale_factor1(6)*phi(242)

            phi(321) = phi(321) + scale_factor1(6)*phi(287)
            phi(322) = phi(322) + scale_factor1(6)*phi(288)
            phi(323) = phi(323) + scale_factor1(6)*phi(289)
            phi(324) = phi(324) + scale_factor1(6)*phi(290)
            phi(325) = phi(325) + scale_factor1(6)*phi(291)
            phi(326) = phi(326) + scale_factor1(6)*phi(292)
            phi(327) = phi(327) + scale_factor1(6)*phi(293)
            phi(328) = phi(328) + scale_factor1(6)*phi(294)
            phi(329) = phi(329) + scale_factor1(6)*phi(295)
            phi(330) = phi(330) + scale_factor1(6)*phi(296)
            phi(331) = phi(331) + scale_factor1(6)*phi(297)
            phi(332) = phi(332) + scale_factor1(6)*phi(298)
            phi(333) = phi(333) + scale_factor1(6)*phi(299)
            phi(334) = phi(334) + scale_factor1(6)*phi(300)
            phi(335) = phi(335) + scale_factor1(6)*phi(301)
            phi(336) = phi(336) + scale_factor1(6)*phi(302)
            phi(337) = phi(337) + scale_factor1(6)*phi(303)
            phi(338) = phi(338) + scale_factor1(6)*phi(304)
            phi(339) = phi(339) + scale_factor1(6)*phi(305)
            phi(340) = phi(340) + scale_factor1(6)*phi(306)
            phi(341) = phi(341) + scale_factor1(6)*phi(307)
            phi(342) = phi(342) + scale_factor1(6)*phi(308)
            phi(343) = phi(343) + scale_factor1(6)*phi(309)
            phi(344) = phi(344) + scale_factor1(6)*phi(310)
            phi(345) = phi(345) + scale_factor1(6)*phi(311)
            phi(346) = phi(346) + scale_factor1(6)*phi(312)
            phi(347) = phi(347) + scale_factor1(6)*phi(313)
            phi(348) = phi(348) + scale_factor1(6)*phi(314)
            phi(349) = phi(349) + scale_factor1(6)*phi(315)
            phi(350) = phi(350) + scale_factor1(6)*phi(316)
            phi(351) = phi(351) + scale_factor1(6)*phi(317)
            phi(352) = phi(352) + scale_factor1(6)*phi(318)
            phi(353) = phi(353) + scale_factor1(6)*phi(319)
            phi(354) = phi(354) + scale_factor1(6)*phi(320)

          end do ! k

!                   ******************************
!                   *                            *
!                   *   Post-contraction phase   *
!                   *                            *
!                   ******************************

          phi(122) = phi(122) - phi(121)
          phi(124) = phi(124) - phi(121)
          phi(127) = phi(127) - phi(121)

          phi(129) = phi(129) - phi(128)
          phi(131) = phi(131) - phi(128)
          phi(134) = phi(134) - phi(128)

          phi(136) = phi(136) - phi(135)
          phi(138) = phi(138) - phi(135)
          phi(141) = phi(141) - phi(135)

          phi(143) = phi(143) - phi(142)
          phi(145) = phi(145) - phi(142)
          phi(148) = phi(148) - phi(142)

          phi(150) = phi(150) - phi(149)
          phi(152) = phi(152) - phi(149)
          phi(155) = phi(155) - phi(149)

          phi(172) = phi(172) - phi(169)
          phi(173) = phi(173) - phi(170)
          phi(175) = phi(175) - phi(170)
          phi(176) = phi(176) - phi(171)
          phi(178) = phi(178) - phi(171)
          phi(181) = phi(181) - phi(171)
          phi(172) = phi(172) - (2.0D+00)*phi(169)
          phi(174) = phi(174) - phi(169)
          phi(175) = phi(175) - (2.0D+00)*phi(170)
          phi(179) = phi(179) - phi(169)
          phi(180) = phi(180) - phi(170)
          phi(181) = phi(181) - (2.0D+00)*phi(171)

          phi(185) = phi(185) - phi(182)
          phi(186) = phi(186) - phi(183)
          phi(188) = phi(188) - phi(183)
          phi(189) = phi(189) - phi(184)
          phi(191) = phi(191) - phi(184)
          phi(194) = phi(194) - phi(184)
          phi(185) = phi(185) - (2.0D+00)*phi(182)
          phi(187) = phi(187) - phi(182)
          phi(188) = phi(188) - (2.0D+00)*phi(183)
          phi(192) = phi(192) - phi(182)
          phi(193) = phi(193) - phi(183)
          phi(194) = phi(194) - (2.0D+00)*phi(184)

          phi(198) = phi(198) - phi(195)
          phi(199) = phi(199) - phi(196)
          phi(201) = phi(201) - phi(196)
          phi(202) = phi(202) - phi(197)
          phi(204) = phi(204) - phi(197)
          phi(207) = phi(207) - phi(197)
          phi(198) = phi(198) - (2.0D+00)*phi(195)
          phi(200) = phi(200) - phi(195)
          phi(201) = phi(201) - (2.0D+00)*phi(196)
          phi(205) = phi(205) - phi(195)
          phi(206) = phi(206) - phi(196)
          phi(207) = phi(207) - (2.0D+00)*phi(197)

          phi(211) = phi(211) - phi(208)
          phi(212) = phi(212) - phi(209)
          phi(214) = phi(214) - phi(209)
          phi(215) = phi(215) - phi(210)
          phi(217) = phi(217) - phi(210)
          phi(220) = phi(220) - phi(210)
          phi(211) = phi(211) - (2.0D+00)*phi(208)
          phi(213) = phi(213) - phi(208)
          phi(214) = phi(214) - (2.0D+00)*phi(209)
          phi(218) = phi(218) - phi(208)
          phi(219) = phi(219) - phi(209)
          phi(220) = phi(220) - (2.0D+00)*phi(210)

          phi(250) = phi(250) - phi(244)
          phi(251) = phi(251) - phi(245)
          phi(252) = phi(252) - phi(246)
          phi(254) = phi(254) - phi(246)
          phi(255) = phi(255) - phi(247)
          phi(256) = phi(256) - phi(248)
          phi(258) = phi(258) - phi(248)
          phi(259) = phi(259) - phi(249)
          phi(261) = phi(261) - phi(249)
          phi(264) = phi(264) - phi(249)
          phi(250) = phi(250) - (2.0D+00)*phi(244)
          phi(251) = phi(251) - (2.0D+00)*phi(245)
          phi(253) = phi(253) - phi(245)
          phi(254) = phi(254) - (2.0D+00)*phi(246)
          phi(255) = phi(255) - (2.0D+00)*phi(247)
          phi(257) = phi(257) - phi(247)
          phi(258) = phi(258) - (2.0D+00)*phi(248)
          phi(262) = phi(262) - phi(247)
          phi(263) = phi(263) - phi(248)
          phi(264) = phi(264) - (2.0D+00)*phi(249)
          phi(244) = phi(244) - phi(243)
          phi(246) = phi(246) - phi(243)
          phi(249) = phi(249) - phi(243)
          phi(250) = phi(250) - (3.0D+00)*phi(244)
          phi(252) = phi(252) - phi(244)
          phi(253) = phi(253) - (2.0D+00)*phi(245)
          phi(254) = phi(254) - (3.0D+00)*phi(246)
          phi(259) = phi(259) - phi(244)
          phi(260) = phi(260) - phi(245)
          phi(261) = phi(261) - phi(246)
          phi(262) = phi(262) - (2.0D+00)*phi(247)
          phi(263) = phi(263) - (2.0D+00)*phi(248)
          phi(264) = phi(264) - (3.0D+00)*phi(249)

          phi(272) = phi(272) - phi(266)
          phi(273) = phi(273) - phi(267)
          phi(274) = phi(274) - phi(268)
          phi(276) = phi(276) - phi(268)
          phi(277) = phi(277) - phi(269)
          phi(278) = phi(278) - phi(270)
          phi(280) = phi(280) - phi(270)
          phi(281) = phi(281) - phi(271)
          phi(283) = phi(283) - phi(271)
          phi(286) = phi(286) - phi(271)
          phi(272) = phi(272) - (2.0D+00)*phi(266)
          phi(273) = phi(273) - (2.0D+00)*phi(267)
          phi(275) = phi(275) - phi(267)
          phi(276) = phi(276) - (2.0D+00)*phi(268)
          phi(277) = phi(277) - (2.0D+00)*phi(269)
          phi(279) = phi(279) - phi(269)
          phi(280) = phi(280) - (2.0D+00)*phi(270)
          phi(284) = phi(284) - phi(269)
          phi(285) = phi(285) - phi(270)
          phi(286) = phi(286) - (2.0D+00)*phi(271)
          phi(266) = phi(266) - phi(265)
          phi(268) = phi(268) - phi(265)
          phi(271) = phi(271) - phi(265)
          phi(272) = phi(272) - (3.0D+00)*phi(266)
          phi(274) = phi(274) - phi(266)
          phi(275) = phi(275) - (2.0D+00)*phi(267)
          phi(276) = phi(276) - (3.0D+00)*phi(268)
          phi(281) = phi(281) - phi(266)
          phi(282) = phi(282) - phi(267)
          phi(283) = phi(283) - phi(268)
          phi(284) = phi(284) - (2.0D+00)*phi(269)
          phi(285) = phi(285) - (2.0D+00)*phi(270)
          phi(286) = phi(286) - (3.0D+00)*phi(271)

          phi(334) = phi(334) - phi(324)
          phi(335) = phi(335) - phi(325)
          phi(336) = phi(336) - phi(326)
          phi(337) = phi(337) - phi(327)
          phi(339) = phi(339) - phi(327)
          phi(340) = phi(340) - phi(328)
          phi(341) = phi(341) - phi(329)
          phi(342) = phi(342) - phi(330)
          phi(344) = phi(344) - phi(330)
          phi(345) = phi(345) - phi(331)
          phi(346) = phi(346) - phi(332)
          phi(348) = phi(348) - phi(332)
          phi(349) = phi(349) - phi(333)
          phi(351) = phi(351) - phi(333)
          phi(354) = phi(354) - phi(333)
          phi(334) = phi(334) - (2.0D+00)*phi(324)
          phi(335) = phi(335) - (2.0D+00)*phi(325)
          phi(336) = phi(336) - (2.0D+00)*phi(326)
          phi(338) = phi(338) - phi(326)
          phi(339) = phi(339) - (2.0D+00)*phi(327)
          phi(340) = phi(340) - (2.0D+00)*phi(328)
          phi(341) = phi(341) - (2.0D+00)*phi(329)
          phi(343) = phi(343) - phi(329)
          phi(344) = phi(344) - (2.0D+00)*phi(330)
          phi(345) = phi(345) - (2.0D+00)*phi(331)
          phi(347) = phi(347) - phi(331)
          phi(348) = phi(348) - (2.0D+00)*phi(332)
          phi(352) = phi(352) - phi(331)
          phi(353) = phi(353) - phi(332)
          phi(354) = phi(354) - (2.0D+00)*phi(333)
          phi(324) = phi(324) - phi(321)
          phi(325) = phi(325) - phi(322)
          phi(327) = phi(327) - phi(322)
          phi(328) = phi(328) - phi(323)
          phi(330) = phi(330) - phi(323)
          phi(333) = phi(333) - phi(323)
          phi(334) = phi(334) - (3.0D+00)*phi(324)
          phi(335) = phi(335) - (3.0D+00)*phi(325)
          phi(337) = phi(337) - phi(325)
          phi(338) = phi(338) - (2.0D+00)*phi(326)
          phi(339) = phi(339) - (3.0D+00)*phi(327)
          phi(340) = phi(340) - (3.0D+00)*phi(328)
          phi(342) = phi(342) - phi(328)
          phi(343) = phi(343) - (2.0D+00)*phi(329)
          phi(344) = phi(344) - (3.0D+00)*phi(330)
          phi(349) = phi(349) - phi(328)
          phi(350) = phi(350) - phi(329)
          phi(351) = phi(351) - phi(330)
          phi(352) = phi(352) - (2.0D+00)*phi(331)
          phi(353) = phi(353) - (2.0D+00)*phi(332)
          phi(354) = phi(354) - (3.0D+00)*phi(333)
          phi(324) = phi(324) - (2.0D+00)*phi(321)
          phi(326) = phi(326) - phi(321)
          phi(327) = phi(327) - (2.0D+00)*phi(322)
          phi(331) = phi(331) - phi(321)
          phi(332) = phi(332) - phi(322)
          phi(333) = phi(333) - (2.0D+00)*phi(323)
          phi(334) = phi(334) - (4.0D+00)*phi(324)
          phi(336) = phi(336) - phi(324)
          phi(337) = phi(337) - (2.0D+00)*phi(325)
          phi(338) = phi(338) - (3.0D+00)*phi(326)
          phi(339) = phi(339) - (4.0D+00)*phi(327)
          phi(345) = phi(345) - phi(324)
          phi(346) = phi(346) - phi(325)
          phi(347) = phi(347) - phi(326)
          phi(348) = phi(348) - phi(327)
          phi(349) = phi(349) - (2.0D+00)*phi(328)
          phi(350) = phi(350) - (2.0D+00)*phi(329)
          phi(351) = phi(351) - (2.0D+00)*phi(330)
          phi(352) = phi(352) - (3.0D+00)*phi(331)
          phi(353) = phi(353) - (3.0D+00)*phi(332)
          phi(354) = phi(354) - (4.0D+00)*phi(333)

!           *************************************
!           *                                   *
!           * -- Angular momentum conversion -- *
!           *   1-center (r) to 2-center (p|q)  *
!           *  Bra contracted transfer equation *
!           *   Horizontal recurrence relation  *
!           *                                   *
!           *************************************

          work2(1, 1) = phi(86)
          work2(2, 1) = phi(83)
          work2(3, 1) = -phi(102)
          work2(4, 1) = -phi(103)
          work2(5, 1) = -phi(104)
          work2(6, 1) = -phi(93)
          work2(7, 1) = -phi(94)
          work2(8, 1) = -phi(95)
          work2(9, 1) = phi(129)
          work2(10, 1) = phi(130)
          work2(11, 1) = phi(131)
          work2(12, 1) = phi(132)
          work2(13, 1) = phi(133)
          work2(14, 1) = phi(134)
          work2(15, 1) = -phi(172)
          work2(16, 1) = -phi(173)
          work2(17, 1) = -phi(174)
          work2(18, 1) = -phi(175)
          work2(19, 1) = -phi(176)
          work2(20, 1) = -phi(177)
          work2(21, 1) = -phi(178)
          work2(22, 1) = -phi(179)
          work2(23, 1) = -phi(180)
          work2(24, 1) = -phi(181)
          work2(25, 1) = phi(88)
          work2(26, 1) = phi(85)
          work2(27, 1) = phi(82)
          work2(28, 1) = -phi(108)
          work2(29, 1) = -phi(109)
          work2(30, 1) = -phi(110)
          work2(31, 1) = -phi(99)
          work2(32, 1) = -phi(100)
          work2(33, 1) = -phi(101)
          work2(34, 1) = phi(143)
          work2(35, 1) = phi(144)
          work2(36, 1) = phi(145)
          work2(37, 1) = phi(146)
          work2(38, 1) = phi(147)
          work2(39, 1) = phi(148)
          work2(40, 1) = phi(122)
          work2(41, 1) = phi(123)
          work2(42, 1) = phi(124)
          work2(43, 1) = phi(125)
          work2(44, 1) = phi(126)
          work2(45, 1) = phi(127)
          work2(46, 1) = -phi(198)
          work2(47, 1) = -phi(199)
          work2(48, 1) = -phi(200)
          work2(49, 1) = -phi(201)
          work2(50, 1) = -phi(202)
          work2(51, 1) = -phi(203)
          work2(52, 1) = -phi(204)
          work2(53, 1) = -phi(205)
          work2(54, 1) = -phi(206)
          work2(55, 1) = -phi(207)
          work2(56, 1) = phi(250)
          work2(57, 1) = phi(251)
          work2(58, 1) = phi(252)
          work2(59, 1) = phi(253)
          work2(60, 1) = phi(254)
          work2(61, 1) = phi(255)
          work2(62, 1) = phi(256)
          work2(63, 1) = phi(257)
          work2(64, 1) = phi(258)
          work2(65, 1) = phi(259)
          work2(66, 1) = phi(260)
          work2(67, 1) = phi(261)
          work2(68, 1) = phi(262)
          work2(69, 1) = phi(263)
          work2(70, 1) = phi(264)
          work2(71, 1) = phi(89)
          work2(72, 1) = phi(87)
          work2(73, 1) = phi(84)
          work2(74, 1) = -phi(111)
          work2(75, 1) = -phi(112)
          work2(76, 1) = -phi(113)
          work2(77, 1) = -phi(105)
          work2(78, 1) = -phi(106)
          work2(79, 1) = -phi(107)
          work2(80, 1) = -phi(96)
          work2(81, 1) = -phi(97)
          work2(82, 1) = -phi(98)
          work2(83, 1) = phi(150)
          work2(84, 1) = phi(151)
          work2(85, 1) = phi(152)
          work2(86, 1) = phi(153)
          work2(87, 1) = phi(154)
          work2(88, 1) = phi(155)
          work2(89, 1) = phi(136)
          work2(90, 1) = phi(137)
          work2(91, 1) = phi(138)
          work2(92, 1) = phi(139)
          work2(93, 1) = phi(140)
          work2(94, 1) = phi(141)
          work2(95, 1) = -phi(211)
          work2(96, 1) = -phi(212)
          work2(97, 1) = -phi(213)
          work2(98, 1) = -phi(214)
          work2(99, 1) = -phi(215)
          work2(100, 1) = -phi(216)
          work2(101, 1) = -phi(217)
          work2(102, 1) = -phi(218)
          work2(103, 1) = -phi(219)
          work2(104, 1) = -phi(220)
          work2(105, 1) = -phi(185)
          work2(106, 1) = -phi(186)
          work2(107, 1) = -phi(187)
          work2(108, 1) = -phi(188)
          work2(109, 1) = -phi(189)
          work2(110, 1) = -phi(190)
          work2(111, 1) = -phi(191)
          work2(112, 1) = -phi(192)
          work2(113, 1) = -phi(193)
          work2(114, 1) = -phi(194)
          work2(115, 1) = phi(272)
          work2(116, 1) = phi(273)
          work2(117, 1) = phi(274)
          work2(118, 1) = phi(275)
          work2(119, 1) = phi(276)
          work2(120, 1) = phi(277)
          work2(121, 1) = phi(278)
          work2(122, 1) = phi(279)
          work2(123, 1) = phi(280)
          work2(124, 1) = phi(281)
          work2(125, 1) = phi(282)
          work2(126, 1) = phi(283)
          work2(127, 1) = phi(284)
          work2(128, 1) = phi(285)
          work2(129, 1) = phi(286)
          work2(130, 1) = -phi(334)
          work2(131, 1) = -phi(335)
          work2(132, 1) = -phi(336)
          work2(133, 1) = -phi(337)
          work2(134, 1) = -phi(338)
          work2(135, 1) = -phi(339)
          work2(136, 1) = -phi(340)
          work2(137, 1) = -phi(341)
          work2(138, 1) = -phi(342)
          work2(139, 1) = -phi(343)
          work2(140, 1) = -phi(344)
          work2(141, 1) = -phi(345)
          work2(142, 1) = -phi(346)
          work2(143, 1) = -phi(347)
          work2(144, 1) = -phi(348)
          work2(145, 1) = -phi(349)
          work2(146, 1) = -phi(350)
          work2(147, 1) = -phi(351)
          work2(148, 1) = -phi(352)
          work2(149, 1) = -phi(353)
          work2(150, 1) = -phi(354)

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
          work2(130, 1) = work2(130, 1) + cnf(1)*work2(115, 1)
          work2(130, 1) = work2(130, 1) + cnf(7)*work2(105, 1)
          work2(131, 1) = work2(131, 1) + cnf(1)*work2(116, 1)
          work2(131, 1) = work2(131, 1) + cnf(6)*work2(106, 1)
          work2(132, 1) = work2(132, 1) + cnf(1)*work2(117, 1)
          work2(132, 1) = work2(132, 1) + cnf(5)*work2(107, 1)
          work2(133, 1) = work2(133, 1) + cnf(1)*work2(118, 1)
          work2(133, 1) = work2(133, 1) + cnf(4)*work2(108, 1)
          work2(134, 1) = work2(134, 1) + cnf(1)*work2(119, 1)
          work2(135, 1) = work2(135, 1) + cnf(2)*work2(119, 1)
          work2(135, 1) = work2(135, 1) + cnf(7)*work2(108, 1)
          work2(136, 1) = work2(136, 1) + cnf(1)*work2(120, 1)
          work2(136, 1) = work2(136, 1) + cnf(6)*work2(109, 1)
          work2(137, 1) = work2(137, 1) + cnf(1)*work2(121, 1)
          work2(137, 1) = work2(137, 1) + cnf(5)*work2(110, 1)
          work2(138, 1) = work2(138, 1) + cnf(1)*work2(122, 1)
          work2(138, 1) = work2(138, 1) + cnf(4)*work2(111, 1)
          work2(139, 1) = work2(139, 1) + cnf(1)*work2(123, 1)
          work2(140, 1) = work2(140, 1) + cnf(2)*work2(123, 1)
          work2(140, 1) = work2(140, 1) + cnf(6)*work2(111, 1)
          work2(141, 1) = work2(141, 1) + cnf(1)*work2(124, 1)
          work2(141, 1) = work2(141, 1) + cnf(5)*work2(112, 1)
          work2(142, 1) = work2(142, 1) + cnf(1)*work2(125, 1)
          work2(142, 1) = work2(142, 1) + cnf(4)*work2(113, 1)
          work2(143, 1) = work2(143, 1) + cnf(1)*work2(126, 1)
          work2(144, 1) = work2(144, 1) + cnf(2)*work2(126, 1)
          work2(144, 1) = work2(144, 1) + cnf(5)*work2(113, 1)
          work2(145, 1) = work2(145, 1) + cnf(1)*work2(127, 1)
          work2(145, 1) = work2(145, 1) + cnf(4)*work2(114, 1)
          work2(146, 1) = work2(146, 1) + cnf(1)*work2(128, 1)
          work2(147, 1) = work2(147, 1) + cnf(2)*work2(128, 1)
          work2(147, 1) = work2(147, 1) + cnf(4)*work2(114, 1)
          work2(148, 1) = work2(148, 1) + cnf(1)*work2(129, 1)
          work2(149, 1) = work2(149, 1) + cnf(2)*work2(129, 1)
          work2(150, 1) = work2(150, 1) + cnf(3)*work2(129, 1)
          work2(150, 1) = work2(150, 1) + cnf(7)*work2(114, 1)
          work2(115, 1) = work2(115, 1) + cnf(1)*work2(95, 1)
          work2(115, 1) = work2(115, 1) + cnf(6)*work2(89, 1)
          work2(116, 1) = work2(116, 1) + cnf(1)*work2(96, 1)
          work2(116, 1) = work2(116, 1) + cnf(5)*work2(90, 1)
          work2(117, 1) = work2(117, 1) + cnf(1)*work2(97, 1)
          work2(117, 1) = work2(117, 1) + cnf(4)*work2(91, 1)
          work2(118, 1) = work2(118, 1) + cnf(1)*work2(98, 1)
          work2(119, 1) = work2(119, 1) + cnf(2)*work2(98, 1)
          work2(119, 1) = work2(119, 1) + cnf(6)*work2(91, 1)
          work2(120, 1) = work2(120, 1) + cnf(1)*work2(99, 1)
          work2(120, 1) = work2(120, 1) + cnf(5)*work2(92, 1)
          work2(121, 1) = work2(121, 1) + cnf(1)*work2(100, 1)
          work2(121, 1) = work2(121, 1) + cnf(4)*work2(93, 1)
          work2(122, 1) = work2(122, 1) + cnf(1)*work2(101, 1)
          work2(123, 1) = work2(123, 1) + cnf(2)*work2(101, 1)
          work2(123, 1) = work2(123, 1) + cnf(5)*work2(93, 1)
          work2(124, 1) = work2(124, 1) + cnf(1)*work2(102, 1)
          work2(124, 1) = work2(124, 1) + cnf(4)*work2(94, 1)
          work2(125, 1) = work2(125, 1) + cnf(1)*work2(103, 1)
          work2(126, 1) = work2(126, 1) + cnf(2)*work2(103, 1)
          work2(126, 1) = work2(126, 1) + cnf(4)*work2(94, 1)
          work2(127, 1) = work2(127, 1) + cnf(1)*work2(104, 1)
          work2(128, 1) = work2(128, 1) + cnf(2)*work2(104, 1)
          work2(129, 1) = work2(129, 1) + cnf(3)*work2(104, 1)
          work2(129, 1) = work2(129, 1) + cnf(6)*work2(94, 1)
          work2(95, 1) = work2(95, 1) + cnf(1)*work2(83, 1)
          work2(95, 1) = work2(95, 1) + cnf(5)*work2(77, 1)
          work2(96, 1) = work2(96, 1) + cnf(1)*work2(84, 1)
          work2(96, 1) = work2(96, 1) + cnf(4)*work2(78, 1)
          work2(97, 1) = work2(97, 1) + cnf(1)*work2(85, 1)
          work2(98, 1) = work2(98, 1) + cnf(2)*work2(85, 1)
          work2(98, 1) = work2(98, 1) + cnf(5)*work2(78, 1)
          work2(99, 1) = work2(99, 1) + cnf(1)*work2(86, 1)
          work2(99, 1) = work2(99, 1) + cnf(4)*work2(79, 1)
          work2(100, 1) = work2(100, 1) + cnf(1)*work2(87, 1)
          work2(101, 1) = work2(101, 1) + cnf(2)*work2(87, 1)
          work2(101, 1) = work2(101, 1) + cnf(4)*work2(79, 1)
          work2(102, 1) = work2(102, 1) + cnf(1)*work2(88, 1)
          work2(103, 1) = work2(103, 1) + cnf(2)*work2(88, 1)
          work2(104, 1) = work2(104, 1) + cnf(3)*work2(88, 1)
          work2(104, 1) = work2(104, 1) + cnf(5)*work2(79, 1)
          work2(105, 1) = work2(105, 1) + cnf(1)*work2(89, 1)
          work2(105, 1) = work2(105, 1) + cnf(5)*work2(80, 1)
          work2(106, 1) = work2(106, 1) + cnf(1)*work2(90, 1)
          work2(106, 1) = work2(106, 1) + cnf(4)*work2(81, 1)
          work2(107, 1) = work2(107, 1) + cnf(1)*work2(91, 1)
          work2(108, 1) = work2(108, 1) + cnf(2)*work2(91, 1)
          work2(108, 1) = work2(108, 1) + cnf(5)*work2(81, 1)
          work2(109, 1) = work2(109, 1) + cnf(1)*work2(92, 1)
          work2(109, 1) = work2(109, 1) + cnf(4)*work2(82, 1)
          work2(110, 1) = work2(110, 1) + cnf(1)*work2(93, 1)
          work2(111, 1) = work2(111, 1) + cnf(2)*work2(93, 1)
          work2(111, 1) = work2(111, 1) + cnf(4)*work2(82, 1)
          work2(112, 1) = work2(112, 1) + cnf(1)*work2(94, 1)
          work2(113, 1) = work2(113, 1) + cnf(2)*work2(94, 1)
          work2(114, 1) = work2(114, 1) + cnf(3)*work2(94, 1)
          work2(114, 1) = work2(114, 1) + cnf(5)*work2(82, 1)
          work2(83, 1) = work2(83, 1) + cnf(1)*work2(74, 1)
          work2(83, 1) = work2(83, 1) + cnf(4)*work2(72, 1)
          work2(84, 1) = work2(84, 1) + cnf(1)*work2(75, 1)
          work2(85, 1) = work2(85, 1) + cnf(2)*work2(75, 1)
          work2(85, 1) = work2(85, 1) + cnf(4)*work2(72, 1)
          work2(86, 1) = work2(86, 1) + cnf(1)*work2(76, 1)
          work2(87, 1) = work2(87, 1) + cnf(2)*work2(76, 1)
          work2(88, 1) = work2(88, 1) + cnf(3)*work2(76, 1)
          work2(88, 1) = work2(88, 1) + cnf(4)*work2(72, 1)
          work2(89, 1) = work2(89, 1) + cnf(1)*work2(77, 1)
          work2(89, 1) = work2(89, 1) + cnf(4)*work2(73, 1)
          work2(90, 1) = work2(90, 1) + cnf(1)*work2(78, 1)
          work2(91, 1) = work2(91, 1) + cnf(2)*work2(78, 1)
          work2(91, 1) = work2(91, 1) + cnf(4)*work2(73, 1)
          work2(92, 1) = work2(92, 1) + cnf(1)*work2(79, 1)
          work2(93, 1) = work2(93, 1) + cnf(2)*work2(79, 1)
          work2(94, 1) = work2(94, 1) + cnf(3)*work2(79, 1)
          work2(94, 1) = work2(94, 1) + cnf(4)*work2(73, 1)
          work2(74, 1) = work2(74, 1) + cnf(1)*work2(71, 1)
          work2(75, 1) = work2(75, 1) + cnf(2)*work2(71, 1)
          work2(76, 1) = work2(76, 1) + cnf(3)*work2(71, 1)
          work2(77, 1) = work2(77, 1) + cnf(1)*work2(72, 1)
          work2(78, 1) = work2(78, 1) + cnf(2)*work2(72, 1)
          work2(79, 1) = work2(79, 1) + cnf(3)*work2(72, 1)
          work2(80, 1) = work2(80, 1) + cnf(1)*work2(73, 1)
          work2(81, 1) = work2(81, 1) + cnf(2)*work2(73, 1)
          work2(82, 1) = work2(82, 1) + cnf(3)*work2(73, 1)
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
          work2(130, 1) = work2(130, 1) + cnf(1)*work2(115, 1)
          work2(130, 1) = work2(130, 1) + cnf(6)*work2(105, 1)
          work2(131, 1) = work2(131, 1) + cnf(1)*work2(116, 1)
          work2(131, 1) = work2(131, 1) + cnf(5)*work2(106, 1)
          work2(132, 1) = work2(132, 1) + cnf(1)*work2(117, 1)
          work2(132, 1) = work2(132, 1) + cnf(4)*work2(107, 1)
          work2(133, 1) = work2(133, 1) + cnf(1)*work2(118, 1)
          work2(134, 1) = work2(134, 1) + cnf(2)*work2(118, 1)
          work2(134, 1) = work2(134, 1) + cnf(6)*work2(107, 1)
          work2(135, 1) = work2(135, 1) + cnf(2)*work2(119, 1)
          work2(135, 1) = work2(135, 1) + cnf(6)*work2(108, 1)
          work2(136, 1) = work2(136, 1) + cnf(1)*work2(120, 1)
          work2(136, 1) = work2(136, 1) + cnf(5)*work2(109, 1)
          work2(137, 1) = work2(137, 1) + cnf(1)*work2(121, 1)
          work2(137, 1) = work2(137, 1) + cnf(4)*work2(110, 1)
          work2(138, 1) = work2(138, 1) + cnf(1)*work2(122, 1)
          work2(139, 1) = work2(139, 1) + cnf(2)*work2(122, 1)
          work2(139, 1) = work2(139, 1) + cnf(5)*work2(110, 1)
          work2(140, 1) = work2(140, 1) + cnf(2)*work2(123, 1)
          work2(140, 1) = work2(140, 1) + cnf(5)*work2(111, 1)
          work2(141, 1) = work2(141, 1) + cnf(1)*work2(124, 1)
          work2(141, 1) = work2(141, 1) + cnf(4)*work2(112, 1)
          work2(142, 1) = work2(142, 1) + cnf(1)*work2(125, 1)
          work2(143, 1) = work2(143, 1) + cnf(2)*work2(125, 1)
          work2(143, 1) = work2(143, 1) + cnf(4)*work2(112, 1)
          work2(144, 1) = work2(144, 1) + cnf(2)*work2(126, 1)
          work2(144, 1) = work2(144, 1) + cnf(4)*work2(113, 1)
          work2(145, 1) = work2(145, 1) + cnf(1)*work2(127, 1)
          work2(146, 1) = work2(146, 1) + cnf(2)*work2(127, 1)
          work2(147, 1) = work2(147, 1) + cnf(2)*work2(128, 1)
          work2(148, 1) = work2(148, 1) + cnf(3)*work2(127, 1)
          work2(148, 1) = work2(148, 1) + cnf(6)*work2(112, 1)
          work2(149, 1) = work2(149, 1) + cnf(3)*work2(128, 1)
          work2(149, 1) = work2(149, 1) + cnf(6)*work2(113, 1)
          work2(150, 1) = work2(150, 1) + cnf(3)*work2(129, 1)
          work2(150, 1) = work2(150, 1) + cnf(6)*work2(114, 1)
          work2(115, 1) = work2(115, 1) + cnf(1)*work2(95, 1)
          work2(115, 1) = work2(115, 1) + cnf(5)*work2(89, 1)
          work2(116, 1) = work2(116, 1) + cnf(1)*work2(96, 1)
          work2(116, 1) = work2(116, 1) + cnf(4)*work2(90, 1)
          work2(117, 1) = work2(117, 1) + cnf(1)*work2(97, 1)
          work2(118, 1) = work2(118, 1) + cnf(2)*work2(97, 1)
          work2(118, 1) = work2(118, 1) + cnf(5)*work2(90, 1)
          work2(119, 1) = work2(119, 1) + cnf(2)*work2(98, 1)
          work2(119, 1) = work2(119, 1) + cnf(5)*work2(91, 1)
          work2(120, 1) = work2(120, 1) + cnf(1)*work2(99, 1)
          work2(120, 1) = work2(120, 1) + cnf(4)*work2(92, 1)
          work2(121, 1) = work2(121, 1) + cnf(1)*work2(100, 1)
          work2(122, 1) = work2(122, 1) + cnf(2)*work2(100, 1)
          work2(122, 1) = work2(122, 1) + cnf(4)*work2(92, 1)
          work2(123, 1) = work2(123, 1) + cnf(2)*work2(101, 1)
          work2(123, 1) = work2(123, 1) + cnf(4)*work2(93, 1)
          work2(124, 1) = work2(124, 1) + cnf(1)*work2(102, 1)
          work2(125, 1) = work2(125, 1) + cnf(2)*work2(102, 1)
          work2(126, 1) = work2(126, 1) + cnf(2)*work2(103, 1)
          work2(127, 1) = work2(127, 1) + cnf(3)*work2(102, 1)
          work2(127, 1) = work2(127, 1) + cnf(5)*work2(92, 1)
          work2(128, 1) = work2(128, 1) + cnf(3)*work2(103, 1)
          work2(128, 1) = work2(128, 1) + cnf(5)*work2(93, 1)
          work2(129, 1) = work2(129, 1) + cnf(3)*work2(104, 1)
          work2(129, 1) = work2(129, 1) + cnf(5)*work2(94, 1)
          work2(95, 1) = work2(95, 1) + cnf(1)*work2(83, 1)
          work2(95, 1) = work2(95, 1) + cnf(4)*work2(77, 1)
          work2(96, 1) = work2(96, 1) + cnf(1)*work2(84, 1)
          work2(97, 1) = work2(97, 1) + cnf(2)*work2(84, 1)
          work2(97, 1) = work2(97, 1) + cnf(4)*work2(77, 1)
          work2(98, 1) = work2(98, 1) + cnf(2)*work2(85, 1)
          work2(98, 1) = work2(98, 1) + cnf(4)*work2(78, 1)
          work2(99, 1) = work2(99, 1) + cnf(1)*work2(86, 1)
          work2(100, 1) = work2(100, 1) + cnf(2)*work2(86, 1)
          work2(101, 1) = work2(101, 1) + cnf(2)*work2(87, 1)
          work2(102, 1) = work2(102, 1) + cnf(3)*work2(86, 1)
          work2(102, 1) = work2(102, 1) + cnf(4)*work2(77, 1)
          work2(103, 1) = work2(103, 1) + cnf(3)*work2(87, 1)
          work2(103, 1) = work2(103, 1) + cnf(4)*work2(78, 1)
          work2(104, 1) = work2(104, 1) + cnf(3)*work2(88, 1)
          work2(104, 1) = work2(104, 1) + cnf(4)*work2(79, 1)
          work2(105, 1) = work2(105, 1) + cnf(1)*work2(89, 1)
          work2(105, 1) = work2(105, 1) + cnf(4)*work2(80, 1)
          work2(106, 1) = work2(106, 1) + cnf(1)*work2(90, 1)
          work2(107, 1) = work2(107, 1) + cnf(2)*work2(90, 1)
          work2(107, 1) = work2(107, 1) + cnf(4)*work2(80, 1)
          work2(108, 1) = work2(108, 1) + cnf(2)*work2(91, 1)
          work2(108, 1) = work2(108, 1) + cnf(4)*work2(81, 1)
          work2(109, 1) = work2(109, 1) + cnf(1)*work2(92, 1)
          work2(110, 1) = work2(110, 1) + cnf(2)*work2(92, 1)
          work2(111, 1) = work2(111, 1) + cnf(2)*work2(93, 1)
          work2(112, 1) = work2(112, 1) + cnf(3)*work2(92, 1)
          work2(112, 1) = work2(112, 1) + cnf(4)*work2(80, 1)
          work2(113, 1) = work2(113, 1) + cnf(3)*work2(93, 1)
          work2(113, 1) = work2(113, 1) + cnf(4)*work2(81, 1)
          work2(114, 1) = work2(114, 1) + cnf(3)*work2(94, 1)
          work2(114, 1) = work2(114, 1) + cnf(4)*work2(82, 1)
          work2(83, 1) = work2(83, 1) + cnf(1)*work2(74, 1)
          work2(84, 1) = work2(84, 1) + cnf(2)*work2(74, 1)
          work2(85, 1) = work2(85, 1) + cnf(2)*work2(75, 1)
          work2(86, 1) = work2(86, 1) + cnf(3)*work2(74, 1)
          work2(87, 1) = work2(87, 1) + cnf(3)*work2(75, 1)
          work2(88, 1) = work2(88, 1) + cnf(3)*work2(76, 1)
          work2(89, 1) = work2(89, 1) + cnf(1)*work2(77, 1)
          work2(90, 1) = work2(90, 1) + cnf(2)*work2(77, 1)
          work2(91, 1) = work2(91, 1) + cnf(2)*work2(78, 1)
          work2(92, 1) = work2(92, 1) + cnf(3)*work2(77, 1)
          work2(93, 1) = work2(93, 1) + cnf(3)*work2(78, 1)
          work2(94, 1) = work2(94, 1) + cnf(3)*work2(79, 1)
          work2(130, 1) = work2(130, 1) + cnf(1)*work2(115, 1)
          work2(130, 1) = work2(130, 1) + cnf(5)*work2(105, 1)
          work2(131, 1) = work2(131, 1) + cnf(1)*work2(116, 1)
          work2(131, 1) = work2(131, 1) + cnf(4)*work2(106, 1)
          work2(132, 1) = work2(132, 1) + cnf(1)*work2(117, 1)
          work2(133, 1) = work2(133, 1) + cnf(2)*work2(117, 1)
          work2(133, 1) = work2(133, 1) + cnf(5)*work2(106, 1)
          work2(134, 1) = work2(134, 1) + cnf(2)*work2(118, 1)
          work2(134, 1) = work2(134, 1) + cnf(5)*work2(107, 1)
          work2(135, 1) = work2(135, 1) + cnf(2)*work2(119, 1)
          work2(135, 1) = work2(135, 1) + cnf(5)*work2(108, 1)
          work2(136, 1) = work2(136, 1) + cnf(1)*work2(120, 1)
          work2(136, 1) = work2(136, 1) + cnf(4)*work2(109, 1)
          work2(137, 1) = work2(137, 1) + cnf(1)*work2(121, 1)
          work2(138, 1) = work2(138, 1) + cnf(2)*work2(121, 1)
          work2(138, 1) = work2(138, 1) + cnf(4)*work2(109, 1)
          work2(139, 1) = work2(139, 1) + cnf(2)*work2(122, 1)
          work2(139, 1) = work2(139, 1) + cnf(4)*work2(110, 1)
          work2(140, 1) = work2(140, 1) + cnf(2)*work2(123, 1)
          work2(140, 1) = work2(140, 1) + cnf(4)*work2(111, 1)
          work2(141, 1) = work2(141, 1) + cnf(1)*work2(124, 1)
          work2(142, 1) = work2(142, 1) + cnf(2)*work2(124, 1)
          work2(143, 1) = work2(143, 1) + cnf(2)*work2(125, 1)
          work2(144, 1) = work2(144, 1) + cnf(2)*work2(126, 1)
          work2(145, 1) = work2(145, 1) + cnf(3)*work2(124, 1)
          work2(145, 1) = work2(145, 1) + cnf(5)*work2(109, 1)
          work2(146, 1) = work2(146, 1) + cnf(3)*work2(125, 1)
          work2(146, 1) = work2(146, 1) + cnf(5)*work2(110, 1)
          work2(147, 1) = work2(147, 1) + cnf(3)*work2(126, 1)
          work2(147, 1) = work2(147, 1) + cnf(5)*work2(111, 1)
          work2(148, 1) = work2(148, 1) + cnf(3)*work2(127, 1)
          work2(148, 1) = work2(148, 1) + cnf(5)*work2(112, 1)
          work2(149, 1) = work2(149, 1) + cnf(3)*work2(128, 1)
          work2(149, 1) = work2(149, 1) + cnf(5)*work2(113, 1)
          work2(150, 1) = work2(150, 1) + cnf(3)*work2(129, 1)
          work2(150, 1) = work2(150, 1) + cnf(5)*work2(114, 1)
          work2(115, 1) = work2(115, 1) + cnf(1)*work2(95, 1)
          work2(115, 1) = work2(115, 1) + cnf(4)*work2(89, 1)
          work2(116, 1) = work2(116, 1) + cnf(1)*work2(96, 1)
          work2(117, 1) = work2(117, 1) + cnf(2)*work2(96, 1)
          work2(117, 1) = work2(117, 1) + cnf(4)*work2(89, 1)
          work2(118, 1) = work2(118, 1) + cnf(2)*work2(97, 1)
          work2(118, 1) = work2(118, 1) + cnf(4)*work2(90, 1)
          work2(119, 1) = work2(119, 1) + cnf(2)*work2(98, 1)
          work2(119, 1) = work2(119, 1) + cnf(4)*work2(91, 1)
          work2(120, 1) = work2(120, 1) + cnf(1)*work2(99, 1)
          work2(121, 1) = work2(121, 1) + cnf(2)*work2(99, 1)
          work2(122, 1) = work2(122, 1) + cnf(2)*work2(100, 1)
          work2(123, 1) = work2(123, 1) + cnf(2)*work2(101, 1)
          work2(124, 1) = work2(124, 1) + cnf(3)*work2(99, 1)
          work2(124, 1) = work2(124, 1) + cnf(4)*work2(89, 1)
          work2(125, 1) = work2(125, 1) + cnf(3)*work2(100, 1)
          work2(125, 1) = work2(125, 1) + cnf(4)*work2(90, 1)
          work2(126, 1) = work2(126, 1) + cnf(3)*work2(101, 1)
          work2(126, 1) = work2(126, 1) + cnf(4)*work2(91, 1)
          work2(127, 1) = work2(127, 1) + cnf(3)*work2(102, 1)
          work2(127, 1) = work2(127, 1) + cnf(4)*work2(92, 1)
          work2(128, 1) = work2(128, 1) + cnf(3)*work2(103, 1)
          work2(128, 1) = work2(128, 1) + cnf(4)*work2(93, 1)
          work2(129, 1) = work2(129, 1) + cnf(3)*work2(104, 1)
          work2(129, 1) = work2(129, 1) + cnf(4)*work2(94, 1)
          work2(95, 1) = work2(95, 1) + cnf(1)*work2(83, 1)
          work2(96, 1) = work2(96, 1) + cnf(2)*work2(83, 1)
          work2(97, 1) = work2(97, 1) + cnf(2)*work2(84, 1)
          work2(98, 1) = work2(98, 1) + cnf(2)*work2(85, 1)
          work2(99, 1) = work2(99, 1) + cnf(3)*work2(83, 1)
          work2(100, 1) = work2(100, 1) + cnf(3)*work2(84, 1)
          work2(101, 1) = work2(101, 1) + cnf(3)*work2(85, 1)
          work2(102, 1) = work2(102, 1) + cnf(3)*work2(86, 1)
          work2(103, 1) = work2(103, 1) + cnf(3)*work2(87, 1)
          work2(104, 1) = work2(104, 1) + cnf(3)*work2(88, 1)
          work2(105, 1) = work2(105, 1) + cnf(1)*work2(89, 1)
          work2(106, 1) = work2(106, 1) + cnf(2)*work2(89, 1)
          work2(107, 1) = work2(107, 1) + cnf(2)*work2(90, 1)
          work2(108, 1) = work2(108, 1) + cnf(2)*work2(91, 1)
          work2(109, 1) = work2(109, 1) + cnf(3)*work2(89, 1)
          work2(110, 1) = work2(110, 1) + cnf(3)*work2(90, 1)
          work2(111, 1) = work2(111, 1) + cnf(3)*work2(91, 1)
          work2(112, 1) = work2(112, 1) + cnf(3)*work2(92, 1)
          work2(113, 1) = work2(113, 1) + cnf(3)*work2(93, 1)
          work2(114, 1) = work2(114, 1) + cnf(3)*work2(94, 1)
          work2(130, 1) = work2(130, 1) + cnf(1)*work2(115, 1)
          work2(130, 1) = work2(130, 1) + cnf(4)*work2(105, 1)
          work2(131, 1) = work2(131, 1) + cnf(1)*work2(116, 1)
          work2(132, 1) = work2(132, 1) + cnf(2)*work2(116, 1)
          work2(132, 1) = work2(132, 1) + cnf(4)*work2(105, 1)
          work2(133, 1) = work2(133, 1) + cnf(2)*work2(117, 1)
          work2(133, 1) = work2(133, 1) + cnf(4)*work2(106, 1)
          work2(134, 1) = work2(134, 1) + cnf(2)*work2(118, 1)
          work2(134, 1) = work2(134, 1) + cnf(4)*work2(107, 1)
          work2(135, 1) = work2(135, 1) + cnf(2)*work2(119, 1)
          work2(135, 1) = work2(135, 1) + cnf(4)*work2(108, 1)
          work2(136, 1) = work2(136, 1) + cnf(1)*work2(120, 1)
          work2(137, 1) = work2(137, 1) + cnf(2)*work2(120, 1)
          work2(138, 1) = work2(138, 1) + cnf(2)*work2(121, 1)
          work2(139, 1) = work2(139, 1) + cnf(2)*work2(122, 1)
          work2(140, 1) = work2(140, 1) + cnf(2)*work2(123, 1)
          work2(141, 1) = work2(141, 1) + cnf(3)*work2(120, 1)
          work2(141, 1) = work2(141, 1) + cnf(4)*work2(105, 1)
          work2(142, 1) = work2(142, 1) + cnf(3)*work2(121, 1)
          work2(142, 1) = work2(142, 1) + cnf(4)*work2(106, 1)
          work2(143, 1) = work2(143, 1) + cnf(3)*work2(122, 1)
          work2(143, 1) = work2(143, 1) + cnf(4)*work2(107, 1)
          work2(144, 1) = work2(144, 1) + cnf(3)*work2(123, 1)
          work2(144, 1) = work2(144, 1) + cnf(4)*work2(108, 1)
          work2(145, 1) = work2(145, 1) + cnf(3)*work2(124, 1)
          work2(145, 1) = work2(145, 1) + cnf(4)*work2(109, 1)
          work2(146, 1) = work2(146, 1) + cnf(3)*work2(125, 1)
          work2(146, 1) = work2(146, 1) + cnf(4)*work2(110, 1)
          work2(147, 1) = work2(147, 1) + cnf(3)*work2(126, 1)
          work2(147, 1) = work2(147, 1) + cnf(4)*work2(111, 1)
          work2(148, 1) = work2(148, 1) + cnf(3)*work2(127, 1)
          work2(148, 1) = work2(148, 1) + cnf(4)*work2(112, 1)
          work2(149, 1) = work2(149, 1) + cnf(3)*work2(128, 1)
          work2(149, 1) = work2(149, 1) + cnf(4)*work2(113, 1)
          work2(150, 1) = work2(150, 1) + cnf(3)*work2(129, 1)
          work2(150, 1) = work2(150, 1) + cnf(4)*work2(114, 1)
          work2(115, 1) = work2(115, 1) + cnf(1)*work2(95, 1)
          work2(116, 1) = work2(116, 1) + cnf(2)*work2(95, 1)
          work2(117, 1) = work2(117, 1) + cnf(2)*work2(96, 1)
          work2(118, 1) = work2(118, 1) + cnf(2)*work2(97, 1)
          work2(119, 1) = work2(119, 1) + cnf(2)*work2(98, 1)
          work2(120, 1) = work2(120, 1) + cnf(3)*work2(95, 1)
          work2(121, 1) = work2(121, 1) + cnf(3)*work2(96, 1)
          work2(122, 1) = work2(122, 1) + cnf(3)*work2(97, 1)
          work2(123, 1) = work2(123, 1) + cnf(3)*work2(98, 1)
          work2(124, 1) = work2(124, 1) + cnf(3)*work2(99, 1)
          work2(125, 1) = work2(125, 1) + cnf(3)*work2(100, 1)
          work2(126, 1) = work2(126, 1) + cnf(3)*work2(101, 1)
          work2(127, 1) = work2(127, 1) + cnf(3)*work2(102, 1)
          work2(128, 1) = work2(128, 1) + cnf(3)*work2(103, 1)
          work2(129, 1) = work2(129, 1) + cnf(3)*work2(104, 1)
          work2(130, 1) = work2(130, 1) + cnf(1)*work2(115, 1)
          work2(131, 1) = work2(131, 1) + cnf(2)*work2(115, 1)
          work2(132, 1) = work2(132, 1) + cnf(2)*work2(116, 1)
          work2(133, 1) = work2(133, 1) + cnf(2)*work2(117, 1)
          work2(134, 1) = work2(134, 1) + cnf(2)*work2(118, 1)
          work2(135, 1) = work2(135, 1) + cnf(2)*work2(119, 1)
          work2(136, 1) = work2(136, 1) + cnf(3)*work2(115, 1)
          work2(137, 1) = work2(137, 1) + cnf(3)*work2(116, 1)
          work2(138, 1) = work2(138, 1) + cnf(3)*work2(117, 1)
          work2(139, 1) = work2(139, 1) + cnf(3)*work2(118, 1)
          work2(140, 1) = work2(140, 1) + cnf(3)*work2(119, 1)
          work2(141, 1) = work2(141, 1) + cnf(3)*work2(120, 1)
          work2(142, 1) = work2(142, 1) + cnf(3)*work2(121, 1)
          work2(143, 1) = work2(143, 1) + cnf(3)*work2(122, 1)
          work2(144, 1) = work2(144, 1) + cnf(3)*work2(123, 1)
          work2(145, 1) = work2(145, 1) + cnf(3)*work2(124, 1)
          work2(146, 1) = work2(146, 1) + cnf(3)*work2(125, 1)
          work2(147, 1) = work2(147, 1) + cnf(3)*work2(126, 1)
          work2(148, 1) = work2(148, 1) + cnf(3)*work2(127, 1)
          work2(149, 1) = work2(149, 1) + cnf(3)*work2(128, 1)
          work2(150, 1) = work2(150, 1) + cnf(3)*work2(129, 1)
          work2(151, 1) = work2(56, 1) - cnf(1)*work2(15, 1)
          work2(152, 1) = work2(57, 1) - cnf(1)*work2(16, 1)
          work2(153, 1) = work2(58, 1) - cnf(1)*work2(17, 1)
          work2(154, 1) = work2(59, 1) - cnf(1)*work2(18, 1)
          work2(155, 1) = work2(61, 1) - cnf(1)*work2(19, 1)
          work2(156, 1) = work2(62, 1) - cnf(1)*work2(20, 1)
          work2(157, 1) = work2(63, 1) - cnf(1)*work2(21, 1)
          work2(158, 1) = work2(65, 1) - cnf(1)*work2(22, 1)
          work2(159, 1) = work2(66, 1) - cnf(1)*work2(23, 1)
          work2(160, 1) = work2(68, 1) - cnf(1)*work2(24, 1)
          work2(161, 1) = work2(130, 1) - cnf(1)*work2(56, 1)
          work2(162, 1) = work2(131, 1) - cnf(1)*work2(57, 1)
          work2(163, 1) = work2(132, 1) - cnf(1)*work2(58, 1)
          work2(164, 1) = work2(133, 1) - cnf(1)*work2(59, 1)
          work2(165, 1) = work2(134, 1) - cnf(1)*work2(60, 1)
          work2(166, 1) = work2(136, 1) - cnf(1)*work2(61, 1)
          work2(167, 1) = work2(137, 1) - cnf(1)*work2(62, 1)
          work2(168, 1) = work2(138, 1) - cnf(1)*work2(63, 1)
          work2(169, 1) = work2(139, 1) - cnf(1)*work2(64, 1)
          work2(170, 1) = work2(141, 1) - cnf(1)*work2(65, 1)
          work2(171, 1) = work2(142, 1) - cnf(1)*work2(66, 1)
          work2(172, 1) = work2(143, 1) - cnf(1)*work2(67, 1)
          work2(173, 1) = work2(145, 1) - cnf(1)*work2(68, 1)
          work2(174, 1) = work2(146, 1) - cnf(1)*work2(69, 1)
          work2(175, 1) = work2(148, 1) - cnf(1)*work2(70, 1)
          work2(176, 1) = work2(57, 1) - cnf(2)*work2(15, 1)
          work2(177, 1) = work2(58, 1) - cnf(2)*work2(16, 1)
          work2(178, 1) = work2(59, 1) - cnf(2)*work2(17, 1)
          work2(179, 1) = work2(60, 1) - cnf(2)*work2(18, 1)
          work2(180, 1) = work2(62, 1) - cnf(2)*work2(19, 1)
          work2(181, 1) = work2(63, 1) - cnf(2)*work2(20, 1)
          work2(182, 1) = work2(64, 1) - cnf(2)*work2(21, 1)
          work2(183, 1) = work2(66, 1) - cnf(2)*work2(22, 1)
          work2(184, 1) = work2(67, 1) - cnf(2)*work2(23, 1)
          work2(185, 1) = work2(69, 1) - cnf(2)*work2(24, 1)
          work2(186, 1) = work2(131, 1) - cnf(2)*work2(56, 1)
          work2(187, 1) = work2(132, 1) - cnf(2)*work2(57, 1)
          work2(188, 1) = work2(133, 1) - cnf(2)*work2(58, 1)
          work2(189, 1) = work2(134, 1) - cnf(2)*work2(59, 1)
          work2(190, 1) = work2(135, 1) - cnf(2)*work2(60, 1)
          work2(191, 1) = work2(137, 1) - cnf(2)*work2(61, 1)
          work2(192, 1) = work2(138, 1) - cnf(2)*work2(62, 1)
          work2(193, 1) = work2(139, 1) - cnf(2)*work2(63, 1)
          work2(194, 1) = work2(140, 1) - cnf(2)*work2(64, 1)
          work2(195, 1) = work2(142, 1) - cnf(2)*work2(65, 1)
          work2(196, 1) = work2(143, 1) - cnf(2)*work2(66, 1)
          work2(197, 1) = work2(144, 1) - cnf(2)*work2(67, 1)
          work2(198, 1) = work2(146, 1) - cnf(2)*work2(68, 1)
          work2(199, 1) = work2(147, 1) - cnf(2)*work2(69, 1)
          work2(200, 1) = work2(149, 1) - cnf(2)*work2(70, 1)
          work2(201, 1) = work2(61, 1) - cnf(3)*work2(15, 1)
          work2(202, 1) = work2(62, 1) - cnf(3)*work2(16, 1)
          work2(203, 1) = work2(63, 1) - cnf(3)*work2(17, 1)
          work2(204, 1) = work2(64, 1) - cnf(3)*work2(18, 1)
          work2(205, 1) = work2(65, 1) - cnf(3)*work2(19, 1)
          work2(206, 1) = work2(66, 1) - cnf(3)*work2(20, 1)
          work2(207, 1) = work2(67, 1) - cnf(3)*work2(21, 1)
          work2(208, 1) = work2(68, 1) - cnf(3)*work2(22, 1)
          work2(209, 1) = work2(69, 1) - cnf(3)*work2(23, 1)
          work2(210, 1) = work2(70, 1) - cnf(3)*work2(24, 1)
          work2(211, 1) = work2(136, 1) - cnf(3)*work2(56, 1)
          work2(212, 1) = work2(137, 1) - cnf(3)*work2(57, 1)
          work2(213, 1) = work2(138, 1) - cnf(3)*work2(58, 1)
          work2(214, 1) = work2(139, 1) - cnf(3)*work2(59, 1)
          work2(215, 1) = work2(140, 1) - cnf(3)*work2(60, 1)
          work2(216, 1) = work2(141, 1) - cnf(3)*work2(61, 1)
          work2(217, 1) = work2(142, 1) - cnf(3)*work2(62, 1)
          work2(218, 1) = work2(143, 1) - cnf(3)*work2(63, 1)
          work2(219, 1) = work2(144, 1) - cnf(3)*work2(64, 1)
          work2(220, 1) = work2(145, 1) - cnf(3)*work2(65, 1)
          work2(221, 1) = work2(146, 1) - cnf(3)*work2(66, 1)
          work2(222, 1) = work2(147, 1) - cnf(3)*work2(67, 1)
          work2(223, 1) = work2(148, 1) - cnf(3)*work2(68, 1)
          work2(224, 1) = work2(149, 1) - cnf(3)*work2(69, 1)
          work2(225, 1) = work2(150, 1) - cnf(3)*work2(70, 1)
          work2(226, 1) = work2(161, 1) - cnf(1)*work2(151, 1)
          work2(227, 1) = work2(162, 1) - cnf(1)*work2(152, 1)
          work2(228, 1) = work2(163, 1) - cnf(1)*work2(153, 1)
          work2(229, 1) = work2(164, 1) - cnf(1)*work2(154, 1)
          work2(230, 1) = work2(166, 1) - cnf(1)*work2(155, 1)
          work2(231, 1) = work2(167, 1) - cnf(1)*work2(156, 1)
          work2(232, 1) = work2(168, 1) - cnf(1)*work2(157, 1)
          work2(233, 1) = work2(170, 1) - cnf(1)*work2(158, 1)
          work2(234, 1) = work2(171, 1) - cnf(1)*work2(159, 1)
          work2(235, 1) = work2(173, 1) - cnf(1)*work2(160, 1)
          work2(236, 1) = work2(186, 1) - cnf(1)*work2(176, 1)
          work2(237, 1) = work2(187, 1) - cnf(1)*work2(177, 1)
          work2(238, 1) = work2(188, 1) - cnf(1)*work2(178, 1)
          work2(239, 1) = work2(189, 1) - cnf(1)*work2(179, 1)
          work2(240, 1) = work2(191, 1) - cnf(1)*work2(180, 1)
          work2(241, 1) = work2(192, 1) - cnf(1)*work2(181, 1)
          work2(242, 1) = work2(193, 1) - cnf(1)*work2(182, 1)
          work2(243, 1) = work2(195, 1) - cnf(1)*work2(183, 1)
          work2(244, 1) = work2(196, 1) - cnf(1)*work2(184, 1)
          work2(245, 1) = work2(198, 1) - cnf(1)*work2(185, 1)
          work2(246, 1) = work2(187, 1) - cnf(2)*work2(176, 1)
          work2(247, 1) = work2(188, 1) - cnf(2)*work2(177, 1)
          work2(248, 1) = work2(189, 1) - cnf(2)*work2(178, 1)
          work2(249, 1) = work2(190, 1) - cnf(2)*work2(179, 1)
          work2(250, 1) = work2(192, 1) - cnf(2)*work2(180, 1)
          work2(251, 1) = work2(193, 1) - cnf(2)*work2(181, 1)
          work2(252, 1) = work2(194, 1) - cnf(2)*work2(182, 1)
          work2(253, 1) = work2(196, 1) - cnf(2)*work2(183, 1)
          work2(254, 1) = work2(197, 1) - cnf(2)*work2(184, 1)
          work2(255, 1) = work2(199, 1) - cnf(2)*work2(185, 1)
          work2(256, 1) = work2(211, 1) - cnf(1)*work2(201, 1)
          work2(257, 1) = work2(212, 1) - cnf(1)*work2(202, 1)
          work2(258, 1) = work2(213, 1) - cnf(1)*work2(203, 1)
          work2(259, 1) = work2(214, 1) - cnf(1)*work2(204, 1)
          work2(260, 1) = work2(216, 1) - cnf(1)*work2(205, 1)
          work2(261, 1) = work2(217, 1) - cnf(1)*work2(206, 1)
          work2(262, 1) = work2(218, 1) - cnf(1)*work2(207, 1)
          work2(263, 1) = work2(220, 1) - cnf(1)*work2(208, 1)
          work2(264, 1) = work2(221, 1) - cnf(1)*work2(209, 1)
          work2(265, 1) = work2(223, 1) - cnf(1)*work2(210, 1)
          work2(266, 1) = work2(212, 1) - cnf(2)*work2(201, 1)
          work2(267, 1) = work2(213, 1) - cnf(2)*work2(202, 1)
          work2(268, 1) = work2(214, 1) - cnf(2)*work2(203, 1)
          work2(269, 1) = work2(215, 1) - cnf(2)*work2(204, 1)
          work2(270, 1) = work2(217, 1) - cnf(2)*work2(205, 1)
          work2(271, 1) = work2(218, 1) - cnf(2)*work2(206, 1)
          work2(272, 1) = work2(219, 1) - cnf(2)*work2(207, 1)
          work2(273, 1) = work2(221, 1) - cnf(2)*work2(208, 1)
          work2(274, 1) = work2(222, 1) - cnf(2)*work2(209, 1)
          work2(275, 1) = work2(224, 1) - cnf(2)*work2(210, 1)
          work2(276, 1) = work2(216, 1) - cnf(3)*work2(201, 1)
          work2(277, 1) = work2(217, 1) - cnf(3)*work2(202, 1)
          work2(278, 1) = work2(218, 1) - cnf(3)*work2(203, 1)
          work2(279, 1) = work2(219, 1) - cnf(3)*work2(204, 1)
          work2(280, 1) = work2(220, 1) - cnf(3)*work2(205, 1)
          work2(281, 1) = work2(221, 1) - cnf(3)*work2(206, 1)
          work2(282, 1) = work2(222, 1) - cnf(3)*work2(207, 1)
          work2(283, 1) = work2(223, 1) - cnf(3)*work2(208, 1)
          work2(284, 1) = work2(224, 1) - cnf(3)*work2(209, 1)
          work2(285, 1) = work2(225, 1) - cnf(3)*work2(210, 1)

! ************************
! *                      *
! * Form final integrals *
! *                      *
! ************************

          eri_value(1) = work2(226, 1)
          eri_value(2) = work2(246, 1)
          eri_value(3) = work2(276, 1)
          eri_value(4) = work2(236, 1)*sqrt3
          eri_value(5) = work2(256, 1)*sqrt3
          eri_value(6) = work2(266, 1)*sqrt3
          eri_value(7) = work2(229, 1)
          eri_value(8) = work2(249, 1)
          eri_value(9) = work2(279, 1)
          eri_value(10) = work2(239, 1)*sqrt3
          eri_value(11) = work2(259, 1)*sqrt3
          eri_value(12) = work2(269, 1)*sqrt3
          eri_value(13) = work2(235, 1)
          eri_value(14) = work2(255, 1)
          eri_value(15) = work2(285, 1)
          eri_value(16) = work2(245, 1)*sqrt3
          eri_value(17) = work2(265, 1)*sqrt3
          eri_value(18) = work2(275, 1)*sqrt3
          eri_value(19) = work2(227, 1)*sqrt5
          eri_value(20) = work2(247, 1)*sqrt5
          eri_value(21) = work2(277, 1)*sqrt5
          eri_value(22) = work2(237, 1)*sqrt5*sqrt3
          eri_value(23) = work2(257, 1)*sqrt5*sqrt3
          eri_value(24) = work2(267, 1)*sqrt5*sqrt3
          eri_value(25) = work2(230, 1)*sqrt5
          eri_value(26) = work2(250, 1)*sqrt5
          eri_value(27) = work2(280, 1)*sqrt5
          eri_value(28) = work2(240, 1)*sqrt5*sqrt3
          eri_value(29) = work2(260, 1)*sqrt5*sqrt3
          eri_value(30) = work2(270, 1)*sqrt5*sqrt3
          eri_value(31) = work2(228, 1)*sqrt5
          eri_value(32) = work2(248, 1)*sqrt5
          eri_value(33) = work2(278, 1)*sqrt5
          eri_value(34) = work2(238, 1)*sqrt5*sqrt3
          eri_value(35) = work2(258, 1)*sqrt5*sqrt3
          eri_value(36) = work2(268, 1)*sqrt5*sqrt3
          eri_value(37) = work2(232, 1)*sqrt5
          eri_value(38) = work2(252, 1)*sqrt5
          eri_value(39) = work2(282, 1)*sqrt5
          eri_value(40) = work2(242, 1)*sqrt5*sqrt3
          eri_value(41) = work2(262, 1)*sqrt5*sqrt3
          eri_value(42) = work2(272, 1)*sqrt5*sqrt3
          eri_value(43) = work2(233, 1)*sqrt5
          eri_value(44) = work2(253, 1)*sqrt5
          eri_value(45) = work2(283, 1)*sqrt5
          eri_value(46) = work2(243, 1)*sqrt5*sqrt3
          eri_value(47) = work2(263, 1)*sqrt5*sqrt3
          eri_value(48) = work2(273, 1)*sqrt5*sqrt3
          eri_value(49) = work2(234, 1)*sqrt5
          eri_value(50) = work2(254, 1)*sqrt5
          eri_value(51) = work2(284, 1)*sqrt5
          eri_value(52) = work2(244, 1)*sqrt5*sqrt3
          eri_value(53) = work2(264, 1)*sqrt5*sqrt3
          eri_value(54) = work2(274, 1)*sqrt5*sqrt3
          eri_value(55) = work2(231, 1)*sqrt15
          eri_value(56) = work2(251, 1)*sqrt15
          eri_value(57) = work2(281, 1)*sqrt15
          eri_value(58) = work2(241, 1)*sqrt15*sqrt3
          eri_value(59) = work2(261, 1)*sqrt15*sqrt3
          eri_value(60) = work2(271, 1)*sqrt15*sqrt3

          maxk = 10
          maxl = 6

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
          kstride = 6

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

  end subroutine int0032
end submodule
