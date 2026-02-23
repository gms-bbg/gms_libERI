! The total angular momentum of this class is:           5
! The algorithm chosen is: PHR a.k.a. ERIC
! Writing an ERIC kernel
submodule(eric_kernels) int1130_impl
contains
  module subroutine int1130(pp_pair, sf_pair, density, fock, res)

    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: pp_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

! Variables for the class
    integer(kind=int64), allocatable :: n11bra(:), n03ket(:)
    real(dp), allocatable :: xint11bra(:), xint03ket(:)
    integer(kind=int64) :: nppbra, nsfket
    real(dp) :: scutppbra, scutsfket, test
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
    real(dp) :: ft(9), phi(585), c_factor(56), rxyz(3)
    real(dp) :: work2(24, 9), work1(24)
    real(dp) :: eri_value(90), angl(20)
    real(dp) :: ai, aij, aijk, aijkl
    integer(kind=int64) :: iord(20)
    integer(kind=int64) :: istride, jstride, kstride
    integer(kind=int64) :: ijk, lo, jc, ir, io, jo, ko
    integer(kind=int64) :: ii1, kk1, nij, maxl2, jj1, j2, ijp, nkl, ijklp
    integer(kind=int64) :: l2, ii2, jj2, kk2, ik, il, jk, jl, ll1, ijkp
    integer(kind=int64) :: maxj2, loci, locj, lock, locl, ip, i2, k2
    integer(kind=int64) :: nchunksize_k10, istart, iend, itile, ntile
    logical :: iandj

    data iord/1, &
      2, 3, 4, &
      5, 7, 10, 6, 8, 9, &
      11, 14, 20, 12, 15, 13, 17, 18, 19, 16/

    data angl/1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, sqrt3, 1.0_dp, sqrt3, sqrt3, 1.0_dp, &
      1.0_dp, sqrt5, sqrt5, 1.0_dp, sqrt5, sqrt15, sqrt5, sqrt5, sqrt5, 1.0_dp/

    allocate (n11bra(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (xint11bra(res%n_p_shl*(res%n_p_shl + 1)/2))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

! Start screening

    scutppbra = cutoff_schwarz/maxval(pp_pair%xints)
    nppbra = 0
    do ij = 1, res%n_p_shl*(res%n_p_shl + 1)/2
      if (pp_pair%xints(ij) .ge. scutppbra) then
        nppbra = nppbra + 1
        xint11bra(nppbra) = pp_pair%xints(ij)
        n11bra(nppbra) = ij
      end if
    end do

    scutsfket = cutoff_schwarz/maxval(sf_pair%xints)
    nsfket = 0
    do ij = 1, res%n_s_shl*res%n_f_shl
      if (sf_pair%xints(ij) .ge. scutsfket) then
        nsfket = nsfket + 1
        xint03ket(nsfket) = sf_pair%xints(ij)
        n03ket(nsfket) = ij
      end if
    end do

    nchunksize_k10 = 375000000

    if ((nppbra*nsfket) .le. nchunksize_k10) nchunksize_k10 = nppbra*nsfket
    ntile = int(nppbra*nsfket/nchunksize_k10)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_k10 + 1
      iend = itile*nchunksize_k10
      if (itile .eq. ntile) iend = nppbra*nsfket

! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

! Mappings to GPU

!$omp target teams distribute parallel do default(none) &
!$omp shared(res, density, fock, nquart_start, nquart_end, nppbra, xint11bra, n11bra, xint03ket, n03ket, nsfket, pp_pair, sf_pair) &
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
!$omp private(work2,cnf,work1) &
!$omp private(indxi,maxi,indxj,maxj,indxk,maxk,indxl,maxl,idim,ioff,leni,lenk,ijkl) &
!$omp private(io,ai,jo,aij,jc,ko,aijk,l,lo,aijkl,ir) &
!$omp private(loci,locj,lock,locl,nij,istride,jstride,kstride) &
!$omp private(iandj,maxj2,ip,ii1,ijp,jj1,i2,j2,ijkp,nkl,kk1) &
!$omp private(ijklp,buff,ll1,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ii,jj,kk,ll)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsfket + 1
        kl_tmp = mod(iquart - 1, nsfket) + 1

        test = xint11bra(ij_tmp)*xint03ket(kl_tmp)

        if (test .gt. cutoff_schwarz) then

          ij = n11bra(ij_tmp)
          kl = n03ket(kl_tmp)

          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_f_shl + 1
          lsh_tmp = mod(kl - 1, res%n_f_shl) + 1

          ish = res%i_p_shl(ish_tmp)
          jsh = res%i_p_shl(jsh_tmp)
          ksh = res%i_f_shl(lsh_tmp)
          lsh = res%i_s_shl(ksh_tmp)

          ijtop = res%contr_num(ish)*res%contr_num(jsh)

          kltop = res%contr_num(ksh)*res%contr_num(lsh)

          eri_value = 0.0_dp
          ket_loop = 0

          ! 0.0_dp out phi elements to contract over kl
          ! PHI = "pre-Hermite integrals"

          phi(82) = 0.0_dp
          phi(83) = 0.0_dp
          phi(85) = 0.0_dp
          phi(86) = 0.0_dp
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
          phi(108) = 0.0_dp
          phi(109) = 0.0_dp
          phi(110) = 0.0_dp
          phi(111) = 0.0_dp
          phi(112) = 0.0_dp
          phi(113) = 0.0_dp
          phi(117) = 0.0_dp
          phi(118) = 0.0_dp
          phi(119) = 0.0_dp
          phi(120) = 0.0_dp
          phi(121) = 0.0_dp
          phi(122) = 0.0_dp
          phi(126) = 0.0_dp
          phi(127) = 0.0_dp
          phi(128) = 0.0_dp
          phi(129) = 0.0_dp
          phi(130) = 0.0_dp
          phi(131) = 0.0_dp
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
          phi(157) = 0.0_dp
          phi(158) = 0.0_dp
          phi(159) = 0.0_dp
          phi(167) = 0.0_dp
          phi(168) = 0.0_dp
          phi(169) = 0.0_dp
          phi(170) = 0.0_dp
          phi(171) = 0.0_dp
          phi(172) = 0.0_dp
          phi(173) = 0.0_dp
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
          phi(223) = 0.0_dp
          phi(224) = 0.0_dp
          phi(225) = 0.0_dp
          phi(226) = 0.0_dp
          phi(227) = 0.0_dp
          phi(228) = 0.0_dp
          phi(229) = 0.0_dp
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
          phi(282) = 0.0_dp
          phi(283) = 0.0_dp
          phi(284) = 0.0_dp
          phi(285) = 0.0_dp
          phi(286) = 0.0_dp
          phi(287) = 0.0_dp
          phi(288) = 0.0_dp
          phi(289) = 0.0_dp
          phi(290) = 0.0_dp
          phi(291) = 0.0_dp
          phi(292) = 0.0_dp
          phi(293) = 0.0_dp
          phi(294) = 0.0_dp
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
          phi(347) = 0.0_dp
          phi(348) = 0.0_dp
          phi(349) = 0.0_dp
          phi(350) = 0.0_dp
          phi(351) = 0.0_dp
          phi(352) = 0.0_dp
          phi(353) = 0.0_dp
          phi(354) = 0.0_dp
          phi(355) = 0.0_dp
          phi(356) = 0.0_dp
          phi(357) = 0.0_dp
          phi(358) = 0.0_dp
          phi(359) = 0.0_dp
          phi(373) = 0.0_dp
          phi(374) = 0.0_dp
          phi(375) = 0.0_dp
          phi(376) = 0.0_dp
          phi(377) = 0.0_dp
          phi(378) = 0.0_dp
          phi(379) = 0.0_dp
          phi(380) = 0.0_dp
          phi(381) = 0.0_dp
          phi(382) = 0.0_dp
          phi(383) = 0.0_dp
          phi(384) = 0.0_dp
          phi(385) = 0.0_dp
          phi(408) = 0.0_dp
          phi(409) = 0.0_dp
          phi(410) = 0.0_dp
          phi(411) = 0.0_dp
          phi(412) = 0.0_dp
          phi(413) = 0.0_dp
          phi(414) = 0.0_dp
          phi(415) = 0.0_dp
          phi(416) = 0.0_dp
          phi(417) = 0.0_dp
          phi(418) = 0.0_dp
          phi(419) = 0.0_dp
          phi(420) = 0.0_dp
          phi(421) = 0.0_dp
          phi(422) = 0.0_dp
          phi(423) = 0.0_dp
          phi(424) = 0.0_dp
          phi(425) = 0.0_dp
          phi(426) = 0.0_dp
          phi(427) = 0.0_dp
          phi(428) = 0.0_dp
          phi(429) = 0.0_dp
          phi(452) = 0.0_dp
          phi(453) = 0.0_dp
          phi(454) = 0.0_dp
          phi(455) = 0.0_dp
          phi(456) = 0.0_dp
          phi(457) = 0.0_dp
          phi(458) = 0.0_dp
          phi(459) = 0.0_dp
          phi(460) = 0.0_dp
          phi(461) = 0.0_dp
          phi(462) = 0.0_dp
          phi(463) = 0.0_dp
          phi(464) = 0.0_dp
          phi(465) = 0.0_dp
          phi(466) = 0.0_dp
          phi(467) = 0.0_dp
          phi(468) = 0.0_dp
          phi(469) = 0.0_dp
          phi(470) = 0.0_dp
          phi(471) = 0.0_dp
          phi(472) = 0.0_dp
          phi(473) = 0.0_dp
          phi(496) = 0.0_dp
          phi(497) = 0.0_dp
          phi(498) = 0.0_dp
          phi(499) = 0.0_dp
          phi(500) = 0.0_dp
          phi(501) = 0.0_dp
          phi(502) = 0.0_dp
          phi(503) = 0.0_dp
          phi(504) = 0.0_dp
          phi(505) = 0.0_dp
          phi(506) = 0.0_dp
          phi(507) = 0.0_dp
          phi(508) = 0.0_dp
          phi(509) = 0.0_dp
          phi(510) = 0.0_dp
          phi(511) = 0.0_dp
          phi(512) = 0.0_dp
          phi(513) = 0.0_dp
          phi(514) = 0.0_dp
          phi(515) = 0.0_dp
          phi(516) = 0.0_dp
          phi(517) = 0.0_dp
          phi(552) = 0.0_dp
          phi(553) = 0.0_dp
          phi(554) = 0.0_dp
          phi(555) = 0.0_dp
          phi(556) = 0.0_dp
          phi(557) = 0.0_dp
          phi(558) = 0.0_dp
          phi(559) = 0.0_dp
          phi(560) = 0.0_dp
          phi(561) = 0.0_dp
          phi(562) = 0.0_dp
          phi(563) = 0.0_dp
          phi(564) = 0.0_dp
          phi(565) = 0.0_dp
          phi(566) = 0.0_dp
          phi(567) = 0.0_dp
          phi(568) = 0.0_dp
          phi(569) = 0.0_dp
          phi(570) = 0.0_dp
          phi(571) = 0.0_dp
          phi(572) = 0.0_dp
          phi(573) = 0.0_dp
          phi(574) = 0.0_dp
          phi(575) = 0.0_dp
          phi(576) = 0.0_dp
          phi(577) = 0.0_dp
          phi(578) = 0.0_dp
          phi(579) = 0.0_dp
          phi(580) = 0.0_dp
          phi(581) = 0.0_dp
          phi(582) = 0.0_dp
          phi(583) = 0.0_dp
          phi(584) = 0.0_dp
          phi(585) = 0.0_dp

! Begin looping over kl primitives

          do k = 1, kltop

            t_expon_cd = sf_pair%t_expon_ab(sf_pair%pair_loc(kl) + 1 + ket_loop) ! exp_c + exp_d
            t_expon_c = sf_pair%expon_b(sf_pair%pair_loc(kl) + 1 + ket_loop)
            t_expon_d = sf_pair%expon_a(sf_pair%pair_loc(kl) + 1 + ket_loop)

            t_inverse_expon_cd = 1.0_dp/t_expon_cd
            ccfket = sf_pair%sq(sf_pair%pair_loc(kl) + 1 + ket_loop)*sqrt2_pi_5_4
            slket = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_cd)

            xkl = ((t_expon_c*res%coord_sh(ksh, 1)) + (t_expon_d*res%coord_sh(lsh, 1)))*t_inverse_expon_cd
            ykl = ((t_expon_c*res%coord_sh(ksh, 2)) + (t_expon_d*res%coord_sh(lsh, 2)))*t_inverse_expon_cd
            zkl = ((t_expon_c*res%coord_sh(ksh, 3)) + (t_expon_d*res%coord_sh(lsh, 3)))*t_inverse_expon_cd

            rxket = sf_pair%t_inverse_expon_ab(sf_pair%pair_loc(kl) + 1 + ket_loop)! inverse_expon_ab*0.5_dp

            ket_loop = ket_loop + 1
            bra_loop = 0

            ! 0.0_dp out phi elements to contract over ij

            phi(81) = 0.0_dp
            phi(84) = 0.0_dp
            phi(87) = 0.0_dp
            phi(90) = 0.0_dp
            phi(91) = 0.0_dp
            phi(92) = 0.0_dp
            phi(105) = 0.0_dp
            phi(106) = 0.0_dp
            phi(107) = 0.0_dp
            phi(114) = 0.0_dp
            phi(115) = 0.0_dp
            phi(116) = 0.0_dp
            phi(123) = 0.0_dp
            phi(124) = 0.0_dp
            phi(125) = 0.0_dp
            phi(132) = 0.0_dp
            phi(133) = 0.0_dp
            phi(134) = 0.0_dp
            phi(135) = 0.0_dp
            phi(136) = 0.0_dp
            phi(137) = 0.0_dp
            phi(138) = 0.0_dp
            phi(160) = 0.0_dp
            phi(161) = 0.0_dp
            phi(162) = 0.0_dp
            phi(163) = 0.0_dp
            phi(164) = 0.0_dp
            phi(165) = 0.0_dp
            phi(166) = 0.0_dp
            phi(174) = 0.0_dp
            phi(175) = 0.0_dp
            phi(176) = 0.0_dp
            phi(177) = 0.0_dp
            phi(178) = 0.0_dp
            phi(179) = 0.0_dp
            phi(180) = 0.0_dp
            phi(195) = 0.0_dp
            phi(196) = 0.0_dp
            phi(197) = 0.0_dp
            phi(198) = 0.0_dp
            phi(199) = 0.0_dp
            phi(200) = 0.0_dp
            phi(201) = 0.0_dp
            phi(216) = 0.0_dp
            phi(217) = 0.0_dp
            phi(218) = 0.0_dp
            phi(219) = 0.0_dp
            phi(220) = 0.0_dp
            phi(221) = 0.0_dp
            phi(222) = 0.0_dp
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
            phi(360) = 0.0_dp
            phi(361) = 0.0_dp
            phi(362) = 0.0_dp
            phi(363) = 0.0_dp
            phi(364) = 0.0_dp
            phi(365) = 0.0_dp
            phi(366) = 0.0_dp
            phi(367) = 0.0_dp
            phi(368) = 0.0_dp
            phi(369) = 0.0_dp
            phi(370) = 0.0_dp
            phi(371) = 0.0_dp
            phi(372) = 0.0_dp
            phi(386) = 0.0_dp
            phi(387) = 0.0_dp
            phi(388) = 0.0_dp
            phi(389) = 0.0_dp
            phi(390) = 0.0_dp
            phi(391) = 0.0_dp
            phi(392) = 0.0_dp
            phi(393) = 0.0_dp
            phi(394) = 0.0_dp
            phi(395) = 0.0_dp
            phi(396) = 0.0_dp
            phi(397) = 0.0_dp
            phi(398) = 0.0_dp
            phi(399) = 0.0_dp
            phi(400) = 0.0_dp
            phi(401) = 0.0_dp
            phi(402) = 0.0_dp
            phi(403) = 0.0_dp
            phi(404) = 0.0_dp
            phi(405) = 0.0_dp
            phi(406) = 0.0_dp
            phi(407) = 0.0_dp
            phi(430) = 0.0_dp
            phi(431) = 0.0_dp
            phi(432) = 0.0_dp
            phi(433) = 0.0_dp
            phi(434) = 0.0_dp
            phi(435) = 0.0_dp
            phi(436) = 0.0_dp
            phi(437) = 0.0_dp
            phi(438) = 0.0_dp
            phi(439) = 0.0_dp
            phi(440) = 0.0_dp
            phi(441) = 0.0_dp
            phi(442) = 0.0_dp
            phi(443) = 0.0_dp
            phi(444) = 0.0_dp
            phi(445) = 0.0_dp
            phi(446) = 0.0_dp
            phi(447) = 0.0_dp
            phi(448) = 0.0_dp
            phi(449) = 0.0_dp
            phi(450) = 0.0_dp
            phi(451) = 0.0_dp
            phi(474) = 0.0_dp
            phi(475) = 0.0_dp
            phi(476) = 0.0_dp
            phi(477) = 0.0_dp
            phi(478) = 0.0_dp
            phi(479) = 0.0_dp
            phi(480) = 0.0_dp
            phi(481) = 0.0_dp
            phi(482) = 0.0_dp
            phi(483) = 0.0_dp
            phi(484) = 0.0_dp
            phi(485) = 0.0_dp
            phi(486) = 0.0_dp
            phi(487) = 0.0_dp
            phi(488) = 0.0_dp
            phi(489) = 0.0_dp
            phi(490) = 0.0_dp
            phi(491) = 0.0_dp
            phi(492) = 0.0_dp
            phi(493) = 0.0_dp
            phi(494) = 0.0_dp
            phi(495) = 0.0_dp
            phi(518) = 0.0_dp
            phi(519) = 0.0_dp
            phi(520) = 0.0_dp
            phi(521) = 0.0_dp
            phi(522) = 0.0_dp
            phi(523) = 0.0_dp
            phi(524) = 0.0_dp
            phi(525) = 0.0_dp
            phi(526) = 0.0_dp
            phi(527) = 0.0_dp
            phi(528) = 0.0_dp
            phi(529) = 0.0_dp
            phi(530) = 0.0_dp
            phi(531) = 0.0_dp
            phi(532) = 0.0_dp
            phi(533) = 0.0_dp
            phi(534) = 0.0_dp
            phi(535) = 0.0_dp
            phi(536) = 0.0_dp
            phi(537) = 0.0_dp
            phi(538) = 0.0_dp
            phi(539) = 0.0_dp
            phi(540) = 0.0_dp
            phi(541) = 0.0_dp
            phi(542) = 0.0_dp
            phi(543) = 0.0_dp
            phi(544) = 0.0_dp
            phi(545) = 0.0_dp
            phi(546) = 0.0_dp
            phi(547) = 0.0_dp
            phi(548) = 0.0_dp
            phi(549) = 0.0_dp
            phi(550) = 0.0_dp
            phi(551) = 0.0_dp

! Begin looping over ij primitives

            do i = 1, ijtop

              t_expon_ab = pp_pair%t_expon_ab(pp_pair%pair_loc(ij) + 1 + bra_loop) ! exp_c + exp_d
              t_expon_a = pp_pair%expon_a(pp_pair%pair_loc(ij) + 1 + bra_loop)
              t_expon_b = pp_pair%expon_b(pp_pair%pair_loc(ij) + 1 + bra_loop)

              t_inverse_expon_ab = 1.0_dp/t_expon_ab
              ccfbra = pp_pair%sq(pp_pair%pair_loc(ij) + 1 + bra_loop)*sqrt2_pi_5_4
              slbra = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_ab)

              xij = ((t_expon_a*res%coord_sh(ish, 1)) + (t_expon_b*res%coord_sh(jsh, 1)))*t_inverse_expon_ab
              yij = ((t_expon_a*res%coord_sh(ish, 2)) + (t_expon_b*res%coord_sh(jsh, 2)))*t_inverse_expon_ab
              zij = ((t_expon_a*res%coord_sh(ish, 3)) + (t_expon_b*res%coord_sh(jsh, 3)))*t_inverse_expon_ab

              rxbra = pp_pair%t_inverse_expon_ab(pp_pair%pair_loc(ij) + 1 + bra_loop)! inverse_expon_ab*0.5_dp

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
              fac1 = fac1*rxbra
              scale_factor1(2) = fac1 ! (0.5*inv_exp_ab)** 1
              fac1 = fac1*rxbra
              scale_factor1(3) = fac1 ! (0.5*inv_exp_ab)** 2
              ! j shell scaling factors
              fac2 = 1.0_dp
              scale_factor2(1) = 1.0_dp
              fac2 = fac2*t_expon_b*2.0_dp
              scale_factor2(2) = fac2 ! (2*t_expon_b)** 1
              fac2 = fac2*t_expon_b*2.0_dp
              scale_factor2(3) = fac2 ! (2*t_expon_b)** 2

              ! ij contraction
              phi(81) = phi(81) + scale_factor1(2)*phi(1)
              phi(84) = phi(84) + scale_factor2(2)*scale_factor1(2)*phi(1)
              phi(87) = phi(87) + scale_factor2(3)*scale_factor1(3)*phi(1)

              phi(90) = phi(90) + scale_factor1(2)*phi(2)
              phi(91) = phi(91) + scale_factor1(2)*phi(3)
              phi(92) = phi(92) + scale_factor1(2)*phi(4)
              phi(105) = phi(105) + scale_factor2(2)*scale_factor1(2)*phi(2)
              phi(106) = phi(106) + scale_factor2(2)*scale_factor1(2)*phi(3)
              phi(107) = phi(107) + scale_factor2(2)*scale_factor1(2)*phi(4)
              phi(114) = phi(114) + scale_factor2(2)*scale_factor1(3)*phi(2)
              phi(115) = phi(115) + scale_factor2(2)*scale_factor1(3)*phi(3)
              phi(116) = phi(116) + scale_factor2(2)*scale_factor1(3)*phi(4)
              phi(123) = phi(123) + scale_factor2(3)*scale_factor1(3)*phi(2)
              phi(124) = phi(124) + scale_factor2(3)*scale_factor1(3)*phi(3)
              phi(125) = phi(125) + scale_factor2(3)*scale_factor1(3)*phi(4)

              phi(132) = phi(132) + scale_factor1(2)*phi(5)
              phi(133) = phi(133) + scale_factor1(2)*phi(6)
              phi(134) = phi(134) + scale_factor1(2)*phi(7)
              phi(135) = phi(135) + scale_factor1(2)*phi(8)
              phi(136) = phi(136) + scale_factor1(2)*phi(9)
              phi(137) = phi(137) + scale_factor1(2)*phi(10)
              phi(138) = phi(138) + scale_factor1(2)*phi(11)
              phi(160) = phi(160) + scale_factor2(2)*scale_factor1(2)*phi(5)
              phi(161) = phi(161) + scale_factor2(2)*scale_factor1(2)*phi(6)
              phi(162) = phi(162) + scale_factor2(2)*scale_factor1(2)*phi(7)
              phi(163) = phi(163) + scale_factor2(2)*scale_factor1(2)*phi(8)
              phi(164) = phi(164) + scale_factor2(2)*scale_factor1(2)*phi(9)
              phi(165) = phi(165) + scale_factor2(2)*scale_factor1(2)*phi(10)
              phi(166) = phi(166) + scale_factor2(2)*scale_factor1(2)*phi(11)
              phi(174) = phi(174) + scale_factor1(3)*phi(5)
              phi(175) = phi(175) + scale_factor1(3)*phi(6)
              phi(176) = phi(176) + scale_factor1(3)*phi(7)
              phi(177) = phi(177) + scale_factor1(3)*phi(8)
              phi(178) = phi(178) + scale_factor1(3)*phi(9)
              phi(179) = phi(179) + scale_factor1(3)*phi(10)
              phi(180) = phi(180) + scale_factor1(3)*phi(11)
              phi(195) = phi(195) + scale_factor2(2)*scale_factor1(3)*phi(5)
              phi(196) = phi(196) + scale_factor2(2)*scale_factor1(3)*phi(6)
              phi(197) = phi(197) + scale_factor2(2)*scale_factor1(3)*phi(7)
              phi(198) = phi(198) + scale_factor2(2)*scale_factor1(3)*phi(8)
              phi(199) = phi(199) + scale_factor2(2)*scale_factor1(3)*phi(9)
              phi(200) = phi(200) + scale_factor2(2)*scale_factor1(3)*phi(10)
              phi(201) = phi(201) + scale_factor2(2)*scale_factor1(3)*phi(11)
              phi(216) = phi(216) + scale_factor2(3)*scale_factor1(3)*phi(5)
              phi(217) = phi(217) + scale_factor2(3)*scale_factor1(3)*phi(6)
              phi(218) = phi(218) + scale_factor2(3)*scale_factor1(3)*phi(7)
              phi(219) = phi(219) + scale_factor2(3)*scale_factor1(3)*phi(8)
              phi(220) = phi(220) + scale_factor2(3)*scale_factor1(3)*phi(9)
              phi(221) = phi(221) + scale_factor2(3)*scale_factor1(3)*phi(10)
              phi(222) = phi(222) + scale_factor2(3)*scale_factor1(3)*phi(11)

              phi(230) = phi(230) + scale_factor1(2)*phi(12)
              phi(231) = phi(231) + scale_factor1(2)*phi(13)
              phi(232) = phi(232) + scale_factor1(2)*phi(14)
              phi(233) = phi(233) + scale_factor1(2)*phi(15)
              phi(234) = phi(234) + scale_factor1(2)*phi(16)
              phi(235) = phi(235) + scale_factor1(2)*phi(17)
              phi(236) = phi(236) + scale_factor1(2)*phi(18)
              phi(237) = phi(237) + scale_factor1(2)*phi(19)
              phi(238) = phi(238) + scale_factor1(2)*phi(20)
              phi(239) = phi(239) + scale_factor1(2)*phi(21)
              phi(240) = phi(240) + scale_factor1(2)*phi(22)
              phi(241) = phi(241) + scale_factor1(2)*phi(23)
              phi(242) = phi(242) + scale_factor1(2)*phi(24)
              phi(269) = phi(269) + scale_factor2(2)*scale_factor1(2)*phi(12)
              phi(270) = phi(270) + scale_factor2(2)*scale_factor1(2)*phi(13)
              phi(271) = phi(271) + scale_factor2(2)*scale_factor1(2)*phi(14)
              phi(272) = phi(272) + scale_factor2(2)*scale_factor1(2)*phi(15)
              phi(273) = phi(273) + scale_factor2(2)*scale_factor1(2)*phi(16)
              phi(274) = phi(274) + scale_factor2(2)*scale_factor1(2)*phi(17)
              phi(275) = phi(275) + scale_factor2(2)*scale_factor1(2)*phi(18)
              phi(276) = phi(276) + scale_factor2(2)*scale_factor1(2)*phi(19)
              phi(277) = phi(277) + scale_factor2(2)*scale_factor1(2)*phi(20)
              phi(278) = phi(278) + scale_factor2(2)*scale_factor1(2)*phi(21)
              phi(279) = phi(279) + scale_factor2(2)*scale_factor1(2)*phi(22)
              phi(280) = phi(280) + scale_factor2(2)*scale_factor1(2)*phi(23)
              phi(281) = phi(281) + scale_factor2(2)*scale_factor1(2)*phi(24)
              phi(295) = phi(295) + scale_factor1(3)*phi(12)
              phi(296) = phi(296) + scale_factor1(3)*phi(13)
              phi(297) = phi(297) + scale_factor1(3)*phi(14)
              phi(298) = phi(298) + scale_factor1(3)*phi(15)
              phi(299) = phi(299) + scale_factor1(3)*phi(16)
              phi(300) = phi(300) + scale_factor1(3)*phi(17)
              phi(301) = phi(301) + scale_factor1(3)*phi(18)
              phi(302) = phi(302) + scale_factor1(3)*phi(19)
              phi(303) = phi(303) + scale_factor1(3)*phi(20)
              phi(304) = phi(304) + scale_factor1(3)*phi(21)
              phi(305) = phi(305) + scale_factor1(3)*phi(22)
              phi(306) = phi(306) + scale_factor1(3)*phi(23)
              phi(307) = phi(307) + scale_factor1(3)*phi(24)
              phi(334) = phi(334) + scale_factor2(2)*scale_factor1(3)*phi(12)
              phi(335) = phi(335) + scale_factor2(2)*scale_factor1(3)*phi(13)
              phi(336) = phi(336) + scale_factor2(2)*scale_factor1(3)*phi(14)
              phi(337) = phi(337) + scale_factor2(2)*scale_factor1(3)*phi(15)
              phi(338) = phi(338) + scale_factor2(2)*scale_factor1(3)*phi(16)
              phi(339) = phi(339) + scale_factor2(2)*scale_factor1(3)*phi(17)
              phi(340) = phi(340) + scale_factor2(2)*scale_factor1(3)*phi(18)
              phi(341) = phi(341) + scale_factor2(2)*scale_factor1(3)*phi(19)
              phi(342) = phi(342) + scale_factor2(2)*scale_factor1(3)*phi(20)
              phi(343) = phi(343) + scale_factor2(2)*scale_factor1(3)*phi(21)
              phi(344) = phi(344) + scale_factor2(2)*scale_factor1(3)*phi(22)
              phi(345) = phi(345) + scale_factor2(2)*scale_factor1(3)*phi(23)
              phi(346) = phi(346) + scale_factor2(2)*scale_factor1(3)*phi(24)
              phi(360) = phi(360) + scale_factor2(3)*scale_factor1(3)*phi(12)
              phi(361) = phi(361) + scale_factor2(3)*scale_factor1(3)*phi(13)
              phi(362) = phi(362) + scale_factor2(3)*scale_factor1(3)*phi(14)
              phi(363) = phi(363) + scale_factor2(3)*scale_factor1(3)*phi(15)
              phi(364) = phi(364) + scale_factor2(3)*scale_factor1(3)*phi(16)
              phi(365) = phi(365) + scale_factor2(3)*scale_factor1(3)*phi(17)
              phi(366) = phi(366) + scale_factor2(3)*scale_factor1(3)*phi(18)
              phi(367) = phi(367) + scale_factor2(3)*scale_factor1(3)*phi(19)
              phi(368) = phi(368) + scale_factor2(3)*scale_factor1(3)*phi(20)
              phi(369) = phi(369) + scale_factor2(3)*scale_factor1(3)*phi(21)
              phi(370) = phi(370) + scale_factor2(3)*scale_factor1(3)*phi(22)
              phi(371) = phi(371) + scale_factor2(3)*scale_factor1(3)*phi(23)
              phi(372) = phi(372) + scale_factor2(3)*scale_factor1(3)*phi(24)

              phi(386) = phi(386) + scale_factor1(2)*phi(25)
              phi(387) = phi(387) + scale_factor1(2)*phi(26)
              phi(388) = phi(388) + scale_factor1(2)*phi(27)
              phi(389) = phi(389) + scale_factor1(2)*phi(28)
              phi(390) = phi(390) + scale_factor1(2)*phi(29)
              phi(391) = phi(391) + scale_factor1(2)*phi(30)
              phi(392) = phi(392) + scale_factor1(2)*phi(31)
              phi(393) = phi(393) + scale_factor1(2)*phi(32)
              phi(394) = phi(394) + scale_factor1(2)*phi(33)
              phi(395) = phi(395) + scale_factor1(2)*phi(34)
              phi(396) = phi(396) + scale_factor1(2)*phi(35)
              phi(397) = phi(397) + scale_factor1(2)*phi(36)
              phi(398) = phi(398) + scale_factor1(2)*phi(37)
              phi(399) = phi(399) + scale_factor1(2)*phi(38)
              phi(400) = phi(400) + scale_factor1(2)*phi(39)
              phi(401) = phi(401) + scale_factor1(2)*phi(40)
              phi(402) = phi(402) + scale_factor1(2)*phi(41)
              phi(403) = phi(403) + scale_factor1(2)*phi(42)
              phi(404) = phi(404) + scale_factor1(2)*phi(43)
              phi(405) = phi(405) + scale_factor1(2)*phi(44)
              phi(406) = phi(406) + scale_factor1(2)*phi(45)
              phi(407) = phi(407) + scale_factor1(2)*phi(46)
              phi(430) = phi(430) + scale_factor1(3)*phi(25)
              phi(431) = phi(431) + scale_factor1(3)*phi(26)
              phi(432) = phi(432) + scale_factor1(3)*phi(27)
              phi(433) = phi(433) + scale_factor1(3)*phi(28)
              phi(434) = phi(434) + scale_factor1(3)*phi(29)
              phi(435) = phi(435) + scale_factor1(3)*phi(30)
              phi(436) = phi(436) + scale_factor1(3)*phi(31)
              phi(437) = phi(437) + scale_factor1(3)*phi(32)
              phi(438) = phi(438) + scale_factor1(3)*phi(33)
              phi(439) = phi(439) + scale_factor1(3)*phi(34)
              phi(440) = phi(440) + scale_factor1(3)*phi(35)
              phi(441) = phi(441) + scale_factor1(3)*phi(36)
              phi(442) = phi(442) + scale_factor1(3)*phi(37)
              phi(443) = phi(443) + scale_factor1(3)*phi(38)
              phi(444) = phi(444) + scale_factor1(3)*phi(39)
              phi(445) = phi(445) + scale_factor1(3)*phi(40)
              phi(446) = phi(446) + scale_factor1(3)*phi(41)
              phi(447) = phi(447) + scale_factor1(3)*phi(42)
              phi(448) = phi(448) + scale_factor1(3)*phi(43)
              phi(449) = phi(449) + scale_factor1(3)*phi(44)
              phi(450) = phi(450) + scale_factor1(3)*phi(45)
              phi(451) = phi(451) + scale_factor1(3)*phi(46)
              phi(474) = phi(474) + scale_factor2(2)*scale_factor1(3)*phi(25)
              phi(475) = phi(475) + scale_factor2(2)*scale_factor1(3)*phi(26)
              phi(476) = phi(476) + scale_factor2(2)*scale_factor1(3)*phi(27)
              phi(477) = phi(477) + scale_factor2(2)*scale_factor1(3)*phi(28)
              phi(478) = phi(478) + scale_factor2(2)*scale_factor1(3)*phi(29)
              phi(479) = phi(479) + scale_factor2(2)*scale_factor1(3)*phi(30)
              phi(480) = phi(480) + scale_factor2(2)*scale_factor1(3)*phi(31)
              phi(481) = phi(481) + scale_factor2(2)*scale_factor1(3)*phi(32)
              phi(482) = phi(482) + scale_factor2(2)*scale_factor1(3)*phi(33)
              phi(483) = phi(483) + scale_factor2(2)*scale_factor1(3)*phi(34)
              phi(484) = phi(484) + scale_factor2(2)*scale_factor1(3)*phi(35)
              phi(485) = phi(485) + scale_factor2(2)*scale_factor1(3)*phi(36)
              phi(486) = phi(486) + scale_factor2(2)*scale_factor1(3)*phi(37)
              phi(487) = phi(487) + scale_factor2(2)*scale_factor1(3)*phi(38)
              phi(488) = phi(488) + scale_factor2(2)*scale_factor1(3)*phi(39)
              phi(489) = phi(489) + scale_factor2(2)*scale_factor1(3)*phi(40)
              phi(490) = phi(490) + scale_factor2(2)*scale_factor1(3)*phi(41)
              phi(491) = phi(491) + scale_factor2(2)*scale_factor1(3)*phi(42)
              phi(492) = phi(492) + scale_factor2(2)*scale_factor1(3)*phi(43)
              phi(493) = phi(493) + scale_factor2(2)*scale_factor1(3)*phi(44)
              phi(494) = phi(494) + scale_factor2(2)*scale_factor1(3)*phi(45)
              phi(495) = phi(495) + scale_factor2(2)*scale_factor1(3)*phi(46)

              phi(518) = phi(518) + scale_factor1(3)*phi(47)
              phi(519) = phi(519) + scale_factor1(3)*phi(48)
              phi(520) = phi(520) + scale_factor1(3)*phi(49)
              phi(521) = phi(521) + scale_factor1(3)*phi(50)
              phi(522) = phi(522) + scale_factor1(3)*phi(51)
              phi(523) = phi(523) + scale_factor1(3)*phi(52)
              phi(524) = phi(524) + scale_factor1(3)*phi(53)
              phi(525) = phi(525) + scale_factor1(3)*phi(54)
              phi(526) = phi(526) + scale_factor1(3)*phi(55)
              phi(527) = phi(527) + scale_factor1(3)*phi(56)
              phi(528) = phi(528) + scale_factor1(3)*phi(57)
              phi(529) = phi(529) + scale_factor1(3)*phi(58)
              phi(530) = phi(530) + scale_factor1(3)*phi(59)
              phi(531) = phi(531) + scale_factor1(3)*phi(60)
              phi(532) = phi(532) + scale_factor1(3)*phi(61)
              phi(533) = phi(533) + scale_factor1(3)*phi(62)
              phi(534) = phi(534) + scale_factor1(3)*phi(63)
              phi(535) = phi(535) + scale_factor1(3)*phi(64)
              phi(536) = phi(536) + scale_factor1(3)*phi(65)
              phi(537) = phi(537) + scale_factor1(3)*phi(66)
              phi(538) = phi(538) + scale_factor1(3)*phi(67)
              phi(539) = phi(539) + scale_factor1(3)*phi(68)
              phi(540) = phi(540) + scale_factor1(3)*phi(69)
              phi(541) = phi(541) + scale_factor1(3)*phi(70)
              phi(542) = phi(542) + scale_factor1(3)*phi(71)
              phi(543) = phi(543) + scale_factor1(3)*phi(72)
              phi(544) = phi(544) + scale_factor1(3)*phi(73)
              phi(545) = phi(545) + scale_factor1(3)*phi(74)
              phi(546) = phi(546) + scale_factor1(3)*phi(75)
              phi(547) = phi(547) + scale_factor1(3)*phi(76)
              phi(548) = phi(548) + scale_factor1(3)*phi(77)
              phi(549) = phi(549) + scale_factor1(3)*phi(78)
              phi(550) = phi(550) + scale_factor1(3)*phi(79)
              phi(551) = phi(551) + scale_factor1(3)*phi(80)

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
            ! l shell scaling factors
            fac2 = 1.0_dp
            scale_factor2(1) = 1.0_dp
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(2) = fac2 ! (2*t_expon_d)** 1
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(3) = fac2 ! (2*t_expon_d)** 2
            fac2 = fac2*t_expon_d*2.0_dp
            scale_factor2(4) = fac2 ! (2*t_expon_d)** 3

            ! kl contraction
            phi(82) = phi(82) + scale_factor2(2)*scale_factor1(3)*phi(81)
            phi(83) = phi(83) + scale_factor2(4)*scale_factor1(4)*phi(81)
            phi(85) = phi(85) + scale_factor2(2)*scale_factor1(3)*phi(84)
            phi(86) = phi(86) + scale_factor2(4)*scale_factor1(4)*phi(84)
            phi(88) = phi(88) + scale_factor2(2)*scale_factor1(3)*phi(87)
            phi(89) = phi(89) + scale_factor2(4)*scale_factor1(4)*phi(87)

            phi(93) = phi(93) + scale_factor1(3)*phi(90)
            phi(94) = phi(94) + scale_factor1(3)*phi(91)
            phi(95) = phi(95) + scale_factor1(3)*phi(92)
            phi(96) = phi(96) + scale_factor2(2)*scale_factor1(3)*phi(90)
            phi(97) = phi(97) + scale_factor2(2)*scale_factor1(3)*phi(91)
            phi(98) = phi(98) + scale_factor2(2)*scale_factor1(3)*phi(92)
            phi(99) = phi(99) + scale_factor2(3)*scale_factor1(4)*phi(90)
            phi(100) = phi(100) + scale_factor2(3)*scale_factor1(4)*phi(91)
            phi(101) = phi(101) + scale_factor2(3)*scale_factor1(4)*phi(92)
            phi(102) = phi(102) + scale_factor2(4)*scale_factor1(4)*phi(90)
            phi(103) = phi(103) + scale_factor2(4)*scale_factor1(4)*phi(91)
            phi(104) = phi(104) + scale_factor2(4)*scale_factor1(4)*phi(92)
            phi(108) = phi(108) + scale_factor1(3)*phi(105)
            phi(109) = phi(109) + scale_factor1(3)*phi(106)
            phi(110) = phi(110) + scale_factor1(3)*phi(107)
            phi(111) = phi(111) + scale_factor2(3)*scale_factor1(4)*phi(105)
            phi(112) = phi(112) + scale_factor2(3)*scale_factor1(4)*phi(106)
            phi(113) = phi(113) + scale_factor2(3)*scale_factor1(4)*phi(107)
            phi(117) = phi(117) + scale_factor2(2)*scale_factor1(3)*phi(114)
            phi(118) = phi(118) + scale_factor2(2)*scale_factor1(3)*phi(115)
            phi(119) = phi(119) + scale_factor2(2)*scale_factor1(3)*phi(116)
            phi(120) = phi(120) + scale_factor2(4)*scale_factor1(4)*phi(114)
            phi(121) = phi(121) + scale_factor2(4)*scale_factor1(4)*phi(115)
            phi(122) = phi(122) + scale_factor2(4)*scale_factor1(4)*phi(116)
            phi(126) = phi(126) + scale_factor1(3)*phi(123)
            phi(127) = phi(127) + scale_factor1(3)*phi(124)
            phi(128) = phi(128) + scale_factor1(3)*phi(125)
            phi(129) = phi(129) + scale_factor2(3)*scale_factor1(4)*phi(123)
            phi(130) = phi(130) + scale_factor2(3)*scale_factor1(4)*phi(124)
            phi(131) = phi(131) + scale_factor2(3)*scale_factor1(4)*phi(125)

            phi(139) = phi(139) + scale_factor1(3)*phi(132)
            phi(140) = phi(140) + scale_factor1(3)*phi(133)
            phi(141) = phi(141) + scale_factor1(3)*phi(134)
            phi(142) = phi(142) + scale_factor1(3)*phi(135)
            phi(143) = phi(143) + scale_factor1(3)*phi(136)
            phi(144) = phi(144) + scale_factor1(3)*phi(137)
            phi(145) = phi(145) + scale_factor1(3)*phi(138)
            phi(146) = phi(146) + scale_factor2(2)*scale_factor1(4)*phi(132)
            phi(147) = phi(147) + scale_factor2(2)*scale_factor1(4)*phi(133)
            phi(148) = phi(148) + scale_factor2(2)*scale_factor1(4)*phi(134)
            phi(149) = phi(149) + scale_factor2(2)*scale_factor1(4)*phi(135)
            phi(150) = phi(150) + scale_factor2(2)*scale_factor1(4)*phi(136)
            phi(151) = phi(151) + scale_factor2(2)*scale_factor1(4)*phi(137)
            phi(152) = phi(152) + scale_factor2(2)*scale_factor1(4)*phi(138)
            phi(153) = phi(153) + scale_factor2(3)*scale_factor1(4)*phi(132)
            phi(154) = phi(154) + scale_factor2(3)*scale_factor1(4)*phi(133)
            phi(155) = phi(155) + scale_factor2(3)*scale_factor1(4)*phi(134)
            phi(156) = phi(156) + scale_factor2(3)*scale_factor1(4)*phi(135)
            phi(157) = phi(157) + scale_factor2(3)*scale_factor1(4)*phi(136)
            phi(158) = phi(158) + scale_factor2(3)*scale_factor1(4)*phi(137)
            phi(159) = phi(159) + scale_factor2(3)*scale_factor1(4)*phi(138)
            phi(167) = phi(167) + scale_factor2(2)*scale_factor1(4)*phi(160)
            phi(168) = phi(168) + scale_factor2(2)*scale_factor1(4)*phi(161)
            phi(169) = phi(169) + scale_factor2(2)*scale_factor1(4)*phi(162)
            phi(170) = phi(170) + scale_factor2(2)*scale_factor1(4)*phi(163)
            phi(171) = phi(171) + scale_factor2(2)*scale_factor1(4)*phi(164)
            phi(172) = phi(172) + scale_factor2(2)*scale_factor1(4)*phi(165)
            phi(173) = phi(173) + scale_factor2(2)*scale_factor1(4)*phi(166)
            phi(181) = phi(181) + scale_factor2(2)*scale_factor1(3)*phi(174)
            phi(182) = phi(182) + scale_factor2(2)*scale_factor1(3)*phi(175)
            phi(183) = phi(183) + scale_factor2(2)*scale_factor1(3)*phi(176)
            phi(184) = phi(184) + scale_factor2(2)*scale_factor1(3)*phi(177)
            phi(185) = phi(185) + scale_factor2(2)*scale_factor1(3)*phi(178)
            phi(186) = phi(186) + scale_factor2(2)*scale_factor1(3)*phi(179)
            phi(187) = phi(187) + scale_factor2(2)*scale_factor1(3)*phi(180)
            phi(188) = phi(188) + scale_factor2(4)*scale_factor1(4)*phi(174)
            phi(189) = phi(189) + scale_factor2(4)*scale_factor1(4)*phi(175)
            phi(190) = phi(190) + scale_factor2(4)*scale_factor1(4)*phi(176)
            phi(191) = phi(191) + scale_factor2(4)*scale_factor1(4)*phi(177)
            phi(192) = phi(192) + scale_factor2(4)*scale_factor1(4)*phi(178)
            phi(193) = phi(193) + scale_factor2(4)*scale_factor1(4)*phi(179)
            phi(194) = phi(194) + scale_factor2(4)*scale_factor1(4)*phi(180)
            phi(202) = phi(202) + scale_factor1(3)*phi(195)
            phi(203) = phi(203) + scale_factor1(3)*phi(196)
            phi(204) = phi(204) + scale_factor1(3)*phi(197)
            phi(205) = phi(205) + scale_factor1(3)*phi(198)
            phi(206) = phi(206) + scale_factor1(3)*phi(199)
            phi(207) = phi(207) + scale_factor1(3)*phi(200)
            phi(208) = phi(208) + scale_factor1(3)*phi(201)
            phi(209) = phi(209) + scale_factor2(3)*scale_factor1(4)*phi(195)
            phi(210) = phi(210) + scale_factor2(3)*scale_factor1(4)*phi(196)
            phi(211) = phi(211) + scale_factor2(3)*scale_factor1(4)*phi(197)
            phi(212) = phi(212) + scale_factor2(3)*scale_factor1(4)*phi(198)
            phi(213) = phi(213) + scale_factor2(3)*scale_factor1(4)*phi(199)
            phi(214) = phi(214) + scale_factor2(3)*scale_factor1(4)*phi(200)
            phi(215) = phi(215) + scale_factor2(3)*scale_factor1(4)*phi(201)
            phi(223) = phi(223) + scale_factor2(2)*scale_factor1(4)*phi(216)
            phi(224) = phi(224) + scale_factor2(2)*scale_factor1(4)*phi(217)
            phi(225) = phi(225) + scale_factor2(2)*scale_factor1(4)*phi(218)
            phi(226) = phi(226) + scale_factor2(2)*scale_factor1(4)*phi(219)
            phi(227) = phi(227) + scale_factor2(2)*scale_factor1(4)*phi(220)
            phi(228) = phi(228) + scale_factor2(2)*scale_factor1(4)*phi(221)
            phi(229) = phi(229) + scale_factor2(2)*scale_factor1(4)*phi(222)

            phi(243) = phi(243) + scale_factor1(4)*phi(230)
            phi(244) = phi(244) + scale_factor1(4)*phi(231)
            phi(245) = phi(245) + scale_factor1(4)*phi(232)
            phi(246) = phi(246) + scale_factor1(4)*phi(233)
            phi(247) = phi(247) + scale_factor1(4)*phi(234)
            phi(248) = phi(248) + scale_factor1(4)*phi(235)
            phi(249) = phi(249) + scale_factor1(4)*phi(236)
            phi(250) = phi(250) + scale_factor1(4)*phi(237)
            phi(251) = phi(251) + scale_factor1(4)*phi(238)
            phi(252) = phi(252) + scale_factor1(4)*phi(239)
            phi(253) = phi(253) + scale_factor1(4)*phi(240)
            phi(254) = phi(254) + scale_factor1(4)*phi(241)
            phi(255) = phi(255) + scale_factor1(4)*phi(242)
            phi(256) = phi(256) + scale_factor2(2)*scale_factor1(4)*phi(230)
            phi(257) = phi(257) + scale_factor2(2)*scale_factor1(4)*phi(231)
            phi(258) = phi(258) + scale_factor2(2)*scale_factor1(4)*phi(232)
            phi(259) = phi(259) + scale_factor2(2)*scale_factor1(4)*phi(233)
            phi(260) = phi(260) + scale_factor2(2)*scale_factor1(4)*phi(234)
            phi(261) = phi(261) + scale_factor2(2)*scale_factor1(4)*phi(235)
            phi(262) = phi(262) + scale_factor2(2)*scale_factor1(4)*phi(236)
            phi(263) = phi(263) + scale_factor2(2)*scale_factor1(4)*phi(237)
            phi(264) = phi(264) + scale_factor2(2)*scale_factor1(4)*phi(238)
            phi(265) = phi(265) + scale_factor2(2)*scale_factor1(4)*phi(239)
            phi(266) = phi(266) + scale_factor2(2)*scale_factor1(4)*phi(240)
            phi(267) = phi(267) + scale_factor2(2)*scale_factor1(4)*phi(241)
            phi(268) = phi(268) + scale_factor2(2)*scale_factor1(4)*phi(242)
            phi(282) = phi(282) + scale_factor1(4)*phi(269)
            phi(283) = phi(283) + scale_factor1(4)*phi(270)
            phi(284) = phi(284) + scale_factor1(4)*phi(271)
            phi(285) = phi(285) + scale_factor1(4)*phi(272)
            phi(286) = phi(286) + scale_factor1(4)*phi(273)
            phi(287) = phi(287) + scale_factor1(4)*phi(274)
            phi(288) = phi(288) + scale_factor1(4)*phi(275)
            phi(289) = phi(289) + scale_factor1(4)*phi(276)
            phi(290) = phi(290) + scale_factor1(4)*phi(277)
            phi(291) = phi(291) + scale_factor1(4)*phi(278)
            phi(292) = phi(292) + scale_factor1(4)*phi(279)
            phi(293) = phi(293) + scale_factor1(4)*phi(280)
            phi(294) = phi(294) + scale_factor1(4)*phi(281)
            phi(308) = phi(308) + scale_factor1(3)*phi(295)
            phi(309) = phi(309) + scale_factor1(3)*phi(296)
            phi(310) = phi(310) + scale_factor1(3)*phi(297)
            phi(311) = phi(311) + scale_factor1(3)*phi(298)
            phi(312) = phi(312) + scale_factor1(3)*phi(299)
            phi(313) = phi(313) + scale_factor1(3)*phi(300)
            phi(314) = phi(314) + scale_factor1(3)*phi(301)
            phi(315) = phi(315) + scale_factor1(3)*phi(302)
            phi(316) = phi(316) + scale_factor1(3)*phi(303)
            phi(317) = phi(317) + scale_factor1(3)*phi(304)
            phi(318) = phi(318) + scale_factor1(3)*phi(305)
            phi(319) = phi(319) + scale_factor1(3)*phi(306)
            phi(320) = phi(320) + scale_factor1(3)*phi(307)
            phi(321) = phi(321) + scale_factor2(3)*scale_factor1(4)*phi(295)
            phi(322) = phi(322) + scale_factor2(3)*scale_factor1(4)*phi(296)
            phi(323) = phi(323) + scale_factor2(3)*scale_factor1(4)*phi(297)
            phi(324) = phi(324) + scale_factor2(3)*scale_factor1(4)*phi(298)
            phi(325) = phi(325) + scale_factor2(3)*scale_factor1(4)*phi(299)
            phi(326) = phi(326) + scale_factor2(3)*scale_factor1(4)*phi(300)
            phi(327) = phi(327) + scale_factor2(3)*scale_factor1(4)*phi(301)
            phi(328) = phi(328) + scale_factor2(3)*scale_factor1(4)*phi(302)
            phi(329) = phi(329) + scale_factor2(3)*scale_factor1(4)*phi(303)
            phi(330) = phi(330) + scale_factor2(3)*scale_factor1(4)*phi(304)
            phi(331) = phi(331) + scale_factor2(3)*scale_factor1(4)*phi(305)
            phi(332) = phi(332) + scale_factor2(3)*scale_factor1(4)*phi(306)
            phi(333) = phi(333) + scale_factor2(3)*scale_factor1(4)*phi(307)
            phi(347) = phi(347) + scale_factor2(2)*scale_factor1(4)*phi(334)
            phi(348) = phi(348) + scale_factor2(2)*scale_factor1(4)*phi(335)
            phi(349) = phi(349) + scale_factor2(2)*scale_factor1(4)*phi(336)
            phi(350) = phi(350) + scale_factor2(2)*scale_factor1(4)*phi(337)
            phi(351) = phi(351) + scale_factor2(2)*scale_factor1(4)*phi(338)
            phi(352) = phi(352) + scale_factor2(2)*scale_factor1(4)*phi(339)
            phi(353) = phi(353) + scale_factor2(2)*scale_factor1(4)*phi(340)
            phi(354) = phi(354) + scale_factor2(2)*scale_factor1(4)*phi(341)
            phi(355) = phi(355) + scale_factor2(2)*scale_factor1(4)*phi(342)
            phi(356) = phi(356) + scale_factor2(2)*scale_factor1(4)*phi(343)
            phi(357) = phi(357) + scale_factor2(2)*scale_factor1(4)*phi(344)
            phi(358) = phi(358) + scale_factor2(2)*scale_factor1(4)*phi(345)
            phi(359) = phi(359) + scale_factor2(2)*scale_factor1(4)*phi(346)
            phi(373) = phi(373) + scale_factor1(4)*phi(360)
            phi(374) = phi(374) + scale_factor1(4)*phi(361)
            phi(375) = phi(375) + scale_factor1(4)*phi(362)
            phi(376) = phi(376) + scale_factor1(4)*phi(363)
            phi(377) = phi(377) + scale_factor1(4)*phi(364)
            phi(378) = phi(378) + scale_factor1(4)*phi(365)
            phi(379) = phi(379) + scale_factor1(4)*phi(366)
            phi(380) = phi(380) + scale_factor1(4)*phi(367)
            phi(381) = phi(381) + scale_factor1(4)*phi(368)
            phi(382) = phi(382) + scale_factor1(4)*phi(369)
            phi(383) = phi(383) + scale_factor1(4)*phi(370)
            phi(384) = phi(384) + scale_factor1(4)*phi(371)
            phi(385) = phi(385) + scale_factor1(4)*phi(372)

            phi(408) = phi(408) + scale_factor1(4)*phi(386)
            phi(409) = phi(409) + scale_factor1(4)*phi(387)
            phi(410) = phi(410) + scale_factor1(4)*phi(388)
            phi(411) = phi(411) + scale_factor1(4)*phi(389)
            phi(412) = phi(412) + scale_factor1(4)*phi(390)
            phi(413) = phi(413) + scale_factor1(4)*phi(391)
            phi(414) = phi(414) + scale_factor1(4)*phi(392)
            phi(415) = phi(415) + scale_factor1(4)*phi(393)
            phi(416) = phi(416) + scale_factor1(4)*phi(394)
            phi(417) = phi(417) + scale_factor1(4)*phi(395)
            phi(418) = phi(418) + scale_factor1(4)*phi(396)
            phi(419) = phi(419) + scale_factor1(4)*phi(397)
            phi(420) = phi(420) + scale_factor1(4)*phi(398)
            phi(421) = phi(421) + scale_factor1(4)*phi(399)
            phi(422) = phi(422) + scale_factor1(4)*phi(400)
            phi(423) = phi(423) + scale_factor1(4)*phi(401)
            phi(424) = phi(424) + scale_factor1(4)*phi(402)
            phi(425) = phi(425) + scale_factor1(4)*phi(403)
            phi(426) = phi(426) + scale_factor1(4)*phi(404)
            phi(427) = phi(427) + scale_factor1(4)*phi(405)
            phi(428) = phi(428) + scale_factor1(4)*phi(406)
            phi(429) = phi(429) + scale_factor1(4)*phi(407)
            phi(452) = phi(452) + scale_factor2(2)*scale_factor1(4)*phi(430)
            phi(453) = phi(453) + scale_factor2(2)*scale_factor1(4)*phi(431)
            phi(454) = phi(454) + scale_factor2(2)*scale_factor1(4)*phi(432)
            phi(455) = phi(455) + scale_factor2(2)*scale_factor1(4)*phi(433)
            phi(456) = phi(456) + scale_factor2(2)*scale_factor1(4)*phi(434)
            phi(457) = phi(457) + scale_factor2(2)*scale_factor1(4)*phi(435)
            phi(458) = phi(458) + scale_factor2(2)*scale_factor1(4)*phi(436)
            phi(459) = phi(459) + scale_factor2(2)*scale_factor1(4)*phi(437)
            phi(460) = phi(460) + scale_factor2(2)*scale_factor1(4)*phi(438)
            phi(461) = phi(461) + scale_factor2(2)*scale_factor1(4)*phi(439)
            phi(462) = phi(462) + scale_factor2(2)*scale_factor1(4)*phi(440)
            phi(463) = phi(463) + scale_factor2(2)*scale_factor1(4)*phi(441)
            phi(464) = phi(464) + scale_factor2(2)*scale_factor1(4)*phi(442)
            phi(465) = phi(465) + scale_factor2(2)*scale_factor1(4)*phi(443)
            phi(466) = phi(466) + scale_factor2(2)*scale_factor1(4)*phi(444)
            phi(467) = phi(467) + scale_factor2(2)*scale_factor1(4)*phi(445)
            phi(468) = phi(468) + scale_factor2(2)*scale_factor1(4)*phi(446)
            phi(469) = phi(469) + scale_factor2(2)*scale_factor1(4)*phi(447)
            phi(470) = phi(470) + scale_factor2(2)*scale_factor1(4)*phi(448)
            phi(471) = phi(471) + scale_factor2(2)*scale_factor1(4)*phi(449)
            phi(472) = phi(472) + scale_factor2(2)*scale_factor1(4)*phi(450)
            phi(473) = phi(473) + scale_factor2(2)*scale_factor1(4)*phi(451)
            phi(496) = phi(496) + scale_factor1(4)*phi(474)
            phi(497) = phi(497) + scale_factor1(4)*phi(475)
            phi(498) = phi(498) + scale_factor1(4)*phi(476)
            phi(499) = phi(499) + scale_factor1(4)*phi(477)
            phi(500) = phi(500) + scale_factor1(4)*phi(478)
            phi(501) = phi(501) + scale_factor1(4)*phi(479)
            phi(502) = phi(502) + scale_factor1(4)*phi(480)
            phi(503) = phi(503) + scale_factor1(4)*phi(481)
            phi(504) = phi(504) + scale_factor1(4)*phi(482)
            phi(505) = phi(505) + scale_factor1(4)*phi(483)
            phi(506) = phi(506) + scale_factor1(4)*phi(484)
            phi(507) = phi(507) + scale_factor1(4)*phi(485)
            phi(508) = phi(508) + scale_factor1(4)*phi(486)
            phi(509) = phi(509) + scale_factor1(4)*phi(487)
            phi(510) = phi(510) + scale_factor1(4)*phi(488)
            phi(511) = phi(511) + scale_factor1(4)*phi(489)
            phi(512) = phi(512) + scale_factor1(4)*phi(490)
            phi(513) = phi(513) + scale_factor1(4)*phi(491)
            phi(514) = phi(514) + scale_factor1(4)*phi(492)
            phi(515) = phi(515) + scale_factor1(4)*phi(493)
            phi(516) = phi(516) + scale_factor1(4)*phi(494)
            phi(517) = phi(517) + scale_factor1(4)*phi(495)

            phi(552) = phi(552) + scale_factor1(4)*phi(518)
            phi(553) = phi(553) + scale_factor1(4)*phi(519)
            phi(554) = phi(554) + scale_factor1(4)*phi(520)
            phi(555) = phi(555) + scale_factor1(4)*phi(521)
            phi(556) = phi(556) + scale_factor1(4)*phi(522)
            phi(557) = phi(557) + scale_factor1(4)*phi(523)
            phi(558) = phi(558) + scale_factor1(4)*phi(524)
            phi(559) = phi(559) + scale_factor1(4)*phi(525)
            phi(560) = phi(560) + scale_factor1(4)*phi(526)
            phi(561) = phi(561) + scale_factor1(4)*phi(527)
            phi(562) = phi(562) + scale_factor1(4)*phi(528)
            phi(563) = phi(563) + scale_factor1(4)*phi(529)
            phi(564) = phi(564) + scale_factor1(4)*phi(530)
            phi(565) = phi(565) + scale_factor1(4)*phi(531)
            phi(566) = phi(566) + scale_factor1(4)*phi(532)
            phi(567) = phi(567) + scale_factor1(4)*phi(533)
            phi(568) = phi(568) + scale_factor1(4)*phi(534)
            phi(569) = phi(569) + scale_factor1(4)*phi(535)
            phi(570) = phi(570) + scale_factor1(4)*phi(536)
            phi(571) = phi(571) + scale_factor1(4)*phi(537)
            phi(572) = phi(572) + scale_factor1(4)*phi(538)
            phi(573) = phi(573) + scale_factor1(4)*phi(539)
            phi(574) = phi(574) + scale_factor1(4)*phi(540)
            phi(575) = phi(575) + scale_factor1(4)*phi(541)
            phi(576) = phi(576) + scale_factor1(4)*phi(542)
            phi(577) = phi(577) + scale_factor1(4)*phi(543)
            phi(578) = phi(578) + scale_factor1(4)*phi(544)
            phi(579) = phi(579) + scale_factor1(4)*phi(545)
            phi(580) = phi(580) + scale_factor1(4)*phi(546)
            phi(581) = phi(581) + scale_factor1(4)*phi(547)
            phi(582) = phi(582) + scale_factor1(4)*phi(548)
            phi(583) = phi(583) + scale_factor1(4)*phi(549)
            phi(584) = phi(584) + scale_factor1(4)*phi(550)
            phi(585) = phi(585) + scale_factor1(4)*phi(551)

          end do ! k

!                   ******************************
!                   *                            *
!                   *   Post-contraction phase   *
!                   *                            *
!                   ******************************

          phi(140) = phi(140) - phi(139)
          phi(142) = phi(142) - phi(139)
          phi(145) = phi(145) - phi(139)

          phi(147) = phi(147) - phi(146)
          phi(149) = phi(149) - phi(146)
          phi(152) = phi(152) - phi(146)

          phi(154) = phi(154) - phi(153)
          phi(156) = phi(156) - phi(153)
          phi(159) = phi(159) - phi(153)

          phi(168) = phi(168) - phi(167)
          phi(170) = phi(170) - phi(167)
          phi(173) = phi(173) - phi(167)

          phi(182) = phi(182) - phi(181)
          phi(184) = phi(184) - phi(181)
          phi(187) = phi(187) - phi(181)

          phi(189) = phi(189) - phi(188)
          phi(191) = phi(191) - phi(188)
          phi(194) = phi(194) - phi(188)

          phi(203) = phi(203) - phi(202)
          phi(205) = phi(205) - phi(202)
          phi(208) = phi(208) - phi(202)

          phi(210) = phi(210) - phi(209)
          phi(212) = phi(212) - phi(209)
          phi(215) = phi(215) - phi(209)

          phi(224) = phi(224) - phi(223)
          phi(226) = phi(226) - phi(223)
          phi(229) = phi(229) - phi(223)

          phi(246) = phi(246) - phi(243)
          phi(247) = phi(247) - phi(244)
          phi(249) = phi(249) - phi(244)
          phi(250) = phi(250) - phi(245)
          phi(252) = phi(252) - phi(245)
          phi(255) = phi(255) - phi(245)
          phi(246) = phi(246) - (2.0D+00)*phi(243)
          phi(248) = phi(248) - phi(243)
          phi(249) = phi(249) - (2.0D+00)*phi(244)
          phi(253) = phi(253) - phi(243)
          phi(254) = phi(254) - phi(244)
          phi(255) = phi(255) - (2.0D+00)*phi(245)

          phi(259) = phi(259) - phi(256)
          phi(260) = phi(260) - phi(257)
          phi(262) = phi(262) - phi(257)
          phi(263) = phi(263) - phi(258)
          phi(265) = phi(265) - phi(258)
          phi(268) = phi(268) - phi(258)
          phi(259) = phi(259) - (2.0D+00)*phi(256)
          phi(261) = phi(261) - phi(256)
          phi(262) = phi(262) - (2.0D+00)*phi(257)
          phi(266) = phi(266) - phi(256)
          phi(267) = phi(267) - phi(257)
          phi(268) = phi(268) - (2.0D+00)*phi(258)

          phi(285) = phi(285) - phi(282)
          phi(286) = phi(286) - phi(283)
          phi(288) = phi(288) - phi(283)
          phi(289) = phi(289) - phi(284)
          phi(291) = phi(291) - phi(284)
          phi(294) = phi(294) - phi(284)
          phi(285) = phi(285) - (2.0D+00)*phi(282)
          phi(287) = phi(287) - phi(282)
          phi(288) = phi(288) - (2.0D+00)*phi(283)
          phi(292) = phi(292) - phi(282)
          phi(293) = phi(293) - phi(283)
          phi(294) = phi(294) - (2.0D+00)*phi(284)

          phi(311) = phi(311) - phi(308)
          phi(312) = phi(312) - phi(309)
          phi(314) = phi(314) - phi(309)
          phi(315) = phi(315) - phi(310)
          phi(317) = phi(317) - phi(310)
          phi(320) = phi(320) - phi(310)
          phi(311) = phi(311) - (2.0D+00)*phi(308)
          phi(313) = phi(313) - phi(308)
          phi(314) = phi(314) - (2.0D+00)*phi(309)
          phi(318) = phi(318) - phi(308)
          phi(319) = phi(319) - phi(309)
          phi(320) = phi(320) - (2.0D+00)*phi(310)

          phi(324) = phi(324) - phi(321)
          phi(325) = phi(325) - phi(322)
          phi(327) = phi(327) - phi(322)
          phi(328) = phi(328) - phi(323)
          phi(330) = phi(330) - phi(323)
          phi(333) = phi(333) - phi(323)
          phi(324) = phi(324) - (2.0D+00)*phi(321)
          phi(326) = phi(326) - phi(321)
          phi(327) = phi(327) - (2.0D+00)*phi(322)
          phi(331) = phi(331) - phi(321)
          phi(332) = phi(332) - phi(322)
          phi(333) = phi(333) - (2.0D+00)*phi(323)

          phi(350) = phi(350) - phi(347)
          phi(351) = phi(351) - phi(348)
          phi(353) = phi(353) - phi(348)
          phi(354) = phi(354) - phi(349)
          phi(356) = phi(356) - phi(349)
          phi(359) = phi(359) - phi(349)
          phi(350) = phi(350) - (2.0D+00)*phi(347)
          phi(352) = phi(352) - phi(347)
          phi(353) = phi(353) - (2.0D+00)*phi(348)
          phi(357) = phi(357) - phi(347)
          phi(358) = phi(358) - phi(348)
          phi(359) = phi(359) - (2.0D+00)*phi(349)

          phi(376) = phi(376) - phi(373)
          phi(377) = phi(377) - phi(374)
          phi(379) = phi(379) - phi(374)
          phi(380) = phi(380) - phi(375)
          phi(382) = phi(382) - phi(375)
          phi(385) = phi(385) - phi(375)
          phi(376) = phi(376) - (2.0D+00)*phi(373)
          phi(378) = phi(378) - phi(373)
          phi(379) = phi(379) - (2.0D+00)*phi(374)
          phi(383) = phi(383) - phi(373)
          phi(384) = phi(384) - phi(374)
          phi(385) = phi(385) - (2.0D+00)*phi(375)

          phi(415) = phi(415) - phi(409)
          phi(416) = phi(416) - phi(410)
          phi(417) = phi(417) - phi(411)
          phi(419) = phi(419) - phi(411)
          phi(420) = phi(420) - phi(412)
          phi(421) = phi(421) - phi(413)
          phi(423) = phi(423) - phi(413)
          phi(424) = phi(424) - phi(414)
          phi(426) = phi(426) - phi(414)
          phi(429) = phi(429) - phi(414)
          phi(415) = phi(415) - (2.0D+00)*phi(409)
          phi(416) = phi(416) - (2.0D+00)*phi(410)
          phi(418) = phi(418) - phi(410)
          phi(419) = phi(419) - (2.0D+00)*phi(411)
          phi(420) = phi(420) - (2.0D+00)*phi(412)
          phi(422) = phi(422) - phi(412)
          phi(423) = phi(423) - (2.0D+00)*phi(413)
          phi(427) = phi(427) - phi(412)
          phi(428) = phi(428) - phi(413)
          phi(429) = phi(429) - (2.0D+00)*phi(414)
          phi(409) = phi(409) - phi(408)
          phi(411) = phi(411) - phi(408)
          phi(414) = phi(414) - phi(408)
          phi(415) = phi(415) - (3.0D+00)*phi(409)
          phi(417) = phi(417) - phi(409)
          phi(418) = phi(418) - (2.0D+00)*phi(410)
          phi(419) = phi(419) - (3.0D+00)*phi(411)
          phi(424) = phi(424) - phi(409)
          phi(425) = phi(425) - phi(410)
          phi(426) = phi(426) - phi(411)
          phi(427) = phi(427) - (2.0D+00)*phi(412)
          phi(428) = phi(428) - (2.0D+00)*phi(413)
          phi(429) = phi(429) - (3.0D+00)*phi(414)

          phi(459) = phi(459) - phi(453)
          phi(460) = phi(460) - phi(454)
          phi(461) = phi(461) - phi(455)
          phi(463) = phi(463) - phi(455)
          phi(464) = phi(464) - phi(456)
          phi(465) = phi(465) - phi(457)
          phi(467) = phi(467) - phi(457)
          phi(468) = phi(468) - phi(458)
          phi(470) = phi(470) - phi(458)
          phi(473) = phi(473) - phi(458)
          phi(459) = phi(459) - (2.0D+00)*phi(453)
          phi(460) = phi(460) - (2.0D+00)*phi(454)
          phi(462) = phi(462) - phi(454)
          phi(463) = phi(463) - (2.0D+00)*phi(455)
          phi(464) = phi(464) - (2.0D+00)*phi(456)
          phi(466) = phi(466) - phi(456)
          phi(467) = phi(467) - (2.0D+00)*phi(457)
          phi(471) = phi(471) - phi(456)
          phi(472) = phi(472) - phi(457)
          phi(473) = phi(473) - (2.0D+00)*phi(458)
          phi(453) = phi(453) - phi(452)
          phi(455) = phi(455) - phi(452)
          phi(458) = phi(458) - phi(452)
          phi(459) = phi(459) - (3.0D+00)*phi(453)
          phi(461) = phi(461) - phi(453)
          phi(462) = phi(462) - (2.0D+00)*phi(454)
          phi(463) = phi(463) - (3.0D+00)*phi(455)
          phi(468) = phi(468) - phi(453)
          phi(469) = phi(469) - phi(454)
          phi(470) = phi(470) - phi(455)
          phi(471) = phi(471) - (2.0D+00)*phi(456)
          phi(472) = phi(472) - (2.0D+00)*phi(457)
          phi(473) = phi(473) - (3.0D+00)*phi(458)

          phi(503) = phi(503) - phi(497)
          phi(504) = phi(504) - phi(498)
          phi(505) = phi(505) - phi(499)
          phi(507) = phi(507) - phi(499)
          phi(508) = phi(508) - phi(500)
          phi(509) = phi(509) - phi(501)
          phi(511) = phi(511) - phi(501)
          phi(512) = phi(512) - phi(502)
          phi(514) = phi(514) - phi(502)
          phi(517) = phi(517) - phi(502)
          phi(503) = phi(503) - (2.0D+00)*phi(497)
          phi(504) = phi(504) - (2.0D+00)*phi(498)
          phi(506) = phi(506) - phi(498)
          phi(507) = phi(507) - (2.0D+00)*phi(499)
          phi(508) = phi(508) - (2.0D+00)*phi(500)
          phi(510) = phi(510) - phi(500)
          phi(511) = phi(511) - (2.0D+00)*phi(501)
          phi(515) = phi(515) - phi(500)
          phi(516) = phi(516) - phi(501)
          phi(517) = phi(517) - (2.0D+00)*phi(502)
          phi(497) = phi(497) - phi(496)
          phi(499) = phi(499) - phi(496)
          phi(502) = phi(502) - phi(496)
          phi(503) = phi(503) - (3.0D+00)*phi(497)
          phi(505) = phi(505) - phi(497)
          phi(506) = phi(506) - (2.0D+00)*phi(498)
          phi(507) = phi(507) - (3.0D+00)*phi(499)
          phi(512) = phi(512) - phi(497)
          phi(513) = phi(513) - phi(498)
          phi(514) = phi(514) - phi(499)
          phi(515) = phi(515) - (2.0D+00)*phi(500)
          phi(516) = phi(516) - (2.0D+00)*phi(501)
          phi(517) = phi(517) - (3.0D+00)*phi(502)

          phi(565) = phi(565) - phi(555)
          phi(566) = phi(566) - phi(556)
          phi(567) = phi(567) - phi(557)
          phi(568) = phi(568) - phi(558)
          phi(570) = phi(570) - phi(558)
          phi(571) = phi(571) - phi(559)
          phi(572) = phi(572) - phi(560)
          phi(573) = phi(573) - phi(561)
          phi(575) = phi(575) - phi(561)
          phi(576) = phi(576) - phi(562)
          phi(577) = phi(577) - phi(563)
          phi(579) = phi(579) - phi(563)
          phi(580) = phi(580) - phi(564)
          phi(582) = phi(582) - phi(564)
          phi(585) = phi(585) - phi(564)
          phi(565) = phi(565) - (2.0D+00)*phi(555)
          phi(566) = phi(566) - (2.0D+00)*phi(556)
          phi(567) = phi(567) - (2.0D+00)*phi(557)
          phi(569) = phi(569) - phi(557)
          phi(570) = phi(570) - (2.0D+00)*phi(558)
          phi(571) = phi(571) - (2.0D+00)*phi(559)
          phi(572) = phi(572) - (2.0D+00)*phi(560)
          phi(574) = phi(574) - phi(560)
          phi(575) = phi(575) - (2.0D+00)*phi(561)
          phi(576) = phi(576) - (2.0D+00)*phi(562)
          phi(578) = phi(578) - phi(562)
          phi(579) = phi(579) - (2.0D+00)*phi(563)
          phi(583) = phi(583) - phi(562)
          phi(584) = phi(584) - phi(563)
          phi(585) = phi(585) - (2.0D+00)*phi(564)
          phi(555) = phi(555) - phi(552)
          phi(556) = phi(556) - phi(553)
          phi(558) = phi(558) - phi(553)
          phi(559) = phi(559) - phi(554)
          phi(561) = phi(561) - phi(554)
          phi(564) = phi(564) - phi(554)
          phi(565) = phi(565) - (3.0D+00)*phi(555)
          phi(566) = phi(566) - (3.0D+00)*phi(556)
          phi(568) = phi(568) - phi(556)
          phi(569) = phi(569) - (2.0D+00)*phi(557)
          phi(570) = phi(570) - (3.0D+00)*phi(558)
          phi(571) = phi(571) - (3.0D+00)*phi(559)
          phi(573) = phi(573) - phi(559)
          phi(574) = phi(574) - (2.0D+00)*phi(560)
          phi(575) = phi(575) - (3.0D+00)*phi(561)
          phi(580) = phi(580) - phi(559)
          phi(581) = phi(581) - phi(560)
          phi(582) = phi(582) - phi(561)
          phi(583) = phi(583) - (2.0D+00)*phi(562)
          phi(584) = phi(584) - (2.0D+00)*phi(563)
          phi(585) = phi(585) - (3.0D+00)*phi(564)
          phi(555) = phi(555) - (2.0D+00)*phi(552)
          phi(557) = phi(557) - phi(552)
          phi(558) = phi(558) - (2.0D+00)*phi(553)
          phi(562) = phi(562) - phi(552)
          phi(563) = phi(563) - phi(553)
          phi(564) = phi(564) - (2.0D+00)*phi(554)
          phi(565) = phi(565) - (4.0D+00)*phi(555)
          phi(567) = phi(567) - phi(555)
          phi(568) = phi(568) - (2.0D+00)*phi(556)
          phi(569) = phi(569) - (3.0D+00)*phi(557)
          phi(570) = phi(570) - (4.0D+00)*phi(558)
          phi(576) = phi(576) - phi(555)
          phi(577) = phi(577) - phi(556)
          phi(578) = phi(578) - phi(557)
          phi(579) = phi(579) - phi(558)
          phi(580) = phi(580) - (2.0D+00)*phi(559)
          phi(581) = phi(581) - (2.0D+00)*phi(560)
          phi(582) = phi(582) - (2.0D+00)*phi(561)
          phi(583) = phi(583) - (3.0D+00)*phi(562)
          phi(584) = phi(584) - (3.0D+00)*phi(563)
          phi(585) = phi(585) - (4.0D+00)*phi(564)

!           *************************************
!           *                                   *
!           * -- Angular momentum conversion -- *
!           *   1-center (r) to 2-center (p|q)  *
!           *  Bra contracted transfer equation *
!           *   Horizontal recurrence relation  *
!           *                                   *
!           *************************************

          cnf(1) = res%coord_sh(jsh, 1) - res%coord_sh(ish, 1)
          cnf(2) = res%coord_sh(jsh, 2) - res%coord_sh(ish, 2)
          cnf(3) = res%coord_sh(jsh, 3) - res%coord_sh(ish, 3)
          cnf(4) = 1.0D+00
          cnf(5) = 2.0D+00
          cnf(6) = 3.0D+00
          cnf(7) = 4.0D+00
          cnf(8) = 5.0D+00
          cnf(9) = 6.0D+00
          cnf(10) = 7.0D+00
          cnf(11) = 8.0D+00

          work1(2) = phi(102) + cnf(1)*phi(86)
          work1(3) = phi(103) + cnf(2)*phi(86)
          work1(4) = phi(104) + cnf(3)*phi(86)
          work1(10) = phi(189) + cnf(1)*phi(120)
          work1(11) = phi(190) + cnf(1)*phi(121)
          work1(12) = phi(191) + cnf(2)*phi(121)
          work1(13) = phi(192) + cnf(1)*phi(122)
          work1(14) = phi(193) + cnf(2)*phi(122)
          work1(15) = phi(194) + cnf(3)*phi(122)
          work1(7) = phi(120) + cnf(1)*phi(89)
          work1(8) = phi(121) + cnf(2)*phi(89)
          work1(9) = phi(122) + cnf(3)*phi(89)
          work1(10) = work1(10) + cnf(4)*phi(83)
          work1(12) = work1(12) + cnf(4)*phi(83)
          work1(15) = work1(15) + cnf(4)*phi(83)
          iii = 1
          ! Generating BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! Bcte done

          work1(2) = phi(96) + cnf(1)*phi(85)
          work1(3) = phi(97) + cnf(2)*phi(85)
          work1(4) = phi(98) + cnf(3)*phi(85)
          work1(10) = phi(182) + cnf(1)*phi(117)
          work1(11) = phi(183) + cnf(1)*phi(118)
          work1(12) = phi(184) + cnf(2)*phi(118)
          work1(13) = phi(185) + cnf(1)*phi(119)
          work1(14) = phi(186) + cnf(2)*phi(119)
          work1(15) = phi(187) + cnf(3)*phi(119)
          work1(7) = phi(117) + cnf(1)*phi(88)
          work1(8) = phi(118) + cnf(2)*phi(88)
          work1(9) = phi(119) + cnf(3)*phi(88)
          work1(10) = work1(10) + cnf(4)*phi(82)
          work1(12) = work1(12) + cnf(4)*phi(82)
          work1(15) = work1(15) + cnf(4)*phi(82)
          iii = 2
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(154) - cnf(1)*phi(111)
          work1(3) = -phi(155) - cnf(2)*phi(111)
          work1(4) = -phi(157) - cnf(3)*phi(111)
          work1(10) = -phi(324) - cnf(1)*phi(210)
          work1(11) = -phi(325) - cnf(1)*phi(211)
          work1(12) = -phi(326) - cnf(2)*phi(211)
          work1(13) = -phi(328) - cnf(1)*phi(213)
          work1(14) = -phi(329) - cnf(2)*phi(213)
          work1(15) = -phi(331) - cnf(3)*phi(213)
          work1(7) = -phi(210) - cnf(1)*phi(129)
          work1(8) = -phi(211) - cnf(2)*phi(129)
          work1(9) = -phi(213) - cnf(3)*phi(129)
          work1(10) = work1(10) - cnf(4)*phi(99)
          work1(12) = work1(12) - cnf(4)*phi(99)
          work1(15) = work1(15) - cnf(4)*phi(99)
          iii = 3
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(155) - cnf(1)*phi(112)
          work1(3) = -phi(156) - cnf(2)*phi(112)
          work1(4) = -phi(158) - cnf(3)*phi(112)
          work1(10) = -phi(325) - cnf(1)*phi(211)
          work1(11) = -phi(326) - cnf(1)*phi(212)
          work1(12) = -phi(327) - cnf(2)*phi(212)
          work1(13) = -phi(329) - cnf(1)*phi(214)
          work1(14) = -phi(330) - cnf(2)*phi(214)
          work1(15) = -phi(332) - cnf(3)*phi(214)
          work1(7) = -phi(211) - cnf(1)*phi(130)
          work1(8) = -phi(212) - cnf(2)*phi(130)
          work1(9) = -phi(214) - cnf(3)*phi(130)
          work1(10) = work1(10) - cnf(4)*phi(100)
          work1(12) = work1(12) - cnf(4)*phi(100)
          work1(15) = work1(15) - cnf(4)*phi(100)
          iii = 4
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(157) - cnf(1)*phi(113)
          work1(3) = -phi(158) - cnf(2)*phi(113)
          work1(4) = -phi(159) - cnf(3)*phi(113)
          work1(10) = -phi(328) - cnf(1)*phi(213)
          work1(11) = -phi(329) - cnf(1)*phi(214)
          work1(12) = -phi(330) - cnf(2)*phi(214)
          work1(13) = -phi(331) - cnf(1)*phi(215)
          work1(14) = -phi(332) - cnf(2)*phi(215)
          work1(15) = -phi(333) - cnf(3)*phi(215)
          work1(7) = -phi(213) - cnf(1)*phi(131)
          work1(8) = -phi(214) - cnf(2)*phi(131)
          work1(9) = -phi(215) - cnf(3)*phi(131)
          work1(10) = work1(10) - cnf(4)*phi(101)
          work1(12) = work1(12) - cnf(4)*phi(101)
          work1(15) = work1(15) - cnf(4)*phi(101)
          iii = 5
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(140) - cnf(1)*phi(108)
          work1(3) = -phi(141) - cnf(2)*phi(108)
          work1(4) = -phi(143) - cnf(3)*phi(108)
          work1(10) = -phi(311) - cnf(1)*phi(203)
          work1(11) = -phi(312) - cnf(1)*phi(204)
          work1(12) = -phi(313) - cnf(2)*phi(204)
          work1(13) = -phi(315) - cnf(1)*phi(206)
          work1(14) = -phi(316) - cnf(2)*phi(206)
          work1(15) = -phi(318) - cnf(3)*phi(206)
          work1(7) = -phi(203) - cnf(1)*phi(126)
          work1(8) = -phi(204) - cnf(2)*phi(126)
          work1(9) = -phi(206) - cnf(3)*phi(126)
          work1(10) = work1(10) - cnf(4)*phi(93)
          work1(12) = work1(12) - cnf(4)*phi(93)
          work1(15) = work1(15) - cnf(4)*phi(93)
          iii = 6
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(141) - cnf(1)*phi(109)
          work1(3) = -phi(142) - cnf(2)*phi(109)
          work1(4) = -phi(144) - cnf(3)*phi(109)
          work1(10) = -phi(312) - cnf(1)*phi(204)
          work1(11) = -phi(313) - cnf(1)*phi(205)
          work1(12) = -phi(314) - cnf(2)*phi(205)
          work1(13) = -phi(316) - cnf(1)*phi(207)
          work1(14) = -phi(317) - cnf(2)*phi(207)
          work1(15) = -phi(319) - cnf(3)*phi(207)
          work1(7) = -phi(204) - cnf(1)*phi(127)
          work1(8) = -phi(205) - cnf(2)*phi(127)
          work1(9) = -phi(207) - cnf(3)*phi(127)
          work1(10) = work1(10) - cnf(4)*phi(94)
          work1(12) = work1(12) - cnf(4)*phi(94)
          work1(15) = work1(15) - cnf(4)*phi(94)
          iii = 7
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(143) - cnf(1)*phi(110)
          work1(3) = -phi(144) - cnf(2)*phi(110)
          work1(4) = -phi(145) - cnf(3)*phi(110)
          work1(10) = -phi(315) - cnf(1)*phi(206)
          work1(11) = -phi(316) - cnf(1)*phi(207)
          work1(12) = -phi(317) - cnf(2)*phi(207)
          work1(13) = -phi(318) - cnf(1)*phi(208)
          work1(14) = -phi(319) - cnf(2)*phi(208)
          work1(15) = -phi(320) - cnf(3)*phi(208)
          work1(7) = -phi(206) - cnf(1)*phi(128)
          work1(8) = -phi(207) - cnf(2)*phi(128)
          work1(9) = -phi(208) - cnf(3)*phi(128)
          work1(10) = work1(10) - cnf(4)*phi(95)
          work1(12) = work1(12) - cnf(4)*phi(95)
          work1(15) = work1(15) - cnf(4)*phi(95)
          iii = 8
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(259) + cnf(1)*phi(168)
          work1(3) = phi(260) + cnf(2)*phi(168)
          work1(4) = phi(263) + cnf(3)*phi(168)
          work1(10) = phi(459) + cnf(1)*phi(350)
          work1(11) = phi(460) + cnf(1)*phi(351)
          work1(12) = phi(461) + cnf(2)*phi(351)
          work1(13) = phi(464) + cnf(1)*phi(354)
          work1(14) = phi(465) + cnf(2)*phi(354)
          work1(15) = phi(468) + cnf(3)*phi(354)
          work1(7) = phi(350) + cnf(1)*phi(224)
          work1(8) = phi(351) + cnf(2)*phi(224)
          work1(9) = phi(354) + cnf(3)*phi(224)
          work1(10) = work1(10) + cnf(4)*phi(147)
          work1(12) = work1(12) + cnf(4)*phi(147)
          work1(15) = work1(15) + cnf(4)*phi(147)
          iii = 9
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(260) + cnf(1)*phi(169)
          work1(3) = phi(261) + cnf(2)*phi(169)
          work1(4) = phi(264) + cnf(3)*phi(169)
          work1(10) = phi(460) + cnf(1)*phi(351)
          work1(11) = phi(461) + cnf(1)*phi(352)
          work1(12) = phi(462) + cnf(2)*phi(352)
          work1(13) = phi(465) + cnf(1)*phi(355)
          work1(14) = phi(466) + cnf(2)*phi(355)
          work1(15) = phi(469) + cnf(3)*phi(355)
          work1(7) = phi(351) + cnf(1)*phi(225)
          work1(8) = phi(352) + cnf(2)*phi(225)
          work1(9) = phi(355) + cnf(3)*phi(225)
          work1(10) = work1(10) + cnf(4)*phi(148)
          work1(12) = work1(12) + cnf(4)*phi(148)
          work1(15) = work1(15) + cnf(4)*phi(148)
          iii = 10
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(261) + cnf(1)*phi(170)
          work1(3) = phi(262) + cnf(2)*phi(170)
          work1(4) = phi(265) + cnf(3)*phi(170)
          work1(10) = phi(461) + cnf(1)*phi(352)
          work1(11) = phi(462) + cnf(1)*phi(353)
          work1(12) = phi(463) + cnf(2)*phi(353)
          work1(13) = phi(466) + cnf(1)*phi(356)
          work1(14) = phi(467) + cnf(2)*phi(356)
          work1(15) = phi(470) + cnf(3)*phi(356)
          work1(7) = phi(352) + cnf(1)*phi(226)
          work1(8) = phi(353) + cnf(2)*phi(226)
          work1(9) = phi(356) + cnf(3)*phi(226)
          work1(10) = work1(10) + cnf(4)*phi(149)
          work1(12) = work1(12) + cnf(4)*phi(149)
          work1(15) = work1(15) + cnf(4)*phi(149)
          iii = 11
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(263) + cnf(1)*phi(171)
          work1(3) = phi(264) + cnf(2)*phi(171)
          work1(4) = phi(266) + cnf(3)*phi(171)
          work1(10) = phi(464) + cnf(1)*phi(354)
          work1(11) = phi(465) + cnf(1)*phi(355)
          work1(12) = phi(466) + cnf(2)*phi(355)
          work1(13) = phi(468) + cnf(1)*phi(357)
          work1(14) = phi(469) + cnf(2)*phi(357)
          work1(15) = phi(471) + cnf(3)*phi(357)
          work1(7) = phi(354) + cnf(1)*phi(227)
          work1(8) = phi(355) + cnf(2)*phi(227)
          work1(9) = phi(357) + cnf(3)*phi(227)
          work1(10) = work1(10) + cnf(4)*phi(150)
          work1(12) = work1(12) + cnf(4)*phi(150)
          work1(15) = work1(15) + cnf(4)*phi(150)
          iii = 12
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(264) + cnf(1)*phi(172)
          work1(3) = phi(265) + cnf(2)*phi(172)
          work1(4) = phi(267) + cnf(3)*phi(172)
          work1(10) = phi(465) + cnf(1)*phi(355)
          work1(11) = phi(466) + cnf(1)*phi(356)
          work1(12) = phi(467) + cnf(2)*phi(356)
          work1(13) = phi(469) + cnf(1)*phi(358)
          work1(14) = phi(470) + cnf(2)*phi(358)
          work1(15) = phi(472) + cnf(3)*phi(358)
          work1(7) = phi(355) + cnf(1)*phi(228)
          work1(8) = phi(356) + cnf(2)*phi(228)
          work1(9) = phi(358) + cnf(3)*phi(228)
          work1(10) = work1(10) + cnf(4)*phi(151)
          work1(12) = work1(12) + cnf(4)*phi(151)
          work1(15) = work1(15) + cnf(4)*phi(151)
          iii = 13
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = phi(266) + cnf(1)*phi(173)
          work1(3) = phi(267) + cnf(2)*phi(173)
          work1(4) = phi(268) + cnf(3)*phi(173)
          work1(10) = phi(468) + cnf(1)*phi(357)
          work1(11) = phi(469) + cnf(1)*phi(358)
          work1(12) = phi(470) + cnf(2)*phi(358)
          work1(13) = phi(471) + cnf(1)*phi(359)
          work1(14) = phi(472) + cnf(2)*phi(359)
          work1(15) = phi(473) + cnf(3)*phi(359)
          work1(7) = phi(357) + cnf(1)*phi(229)
          work1(8) = phi(358) + cnf(2)*phi(229)
          work1(9) = phi(359) + cnf(3)*phi(229)
          work1(10) = work1(10) + cnf(4)*phi(152)
          work1(12) = work1(12) + cnf(4)*phi(152)
          work1(15) = work1(15) + cnf(4)*phi(152)
          iii = 14
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(415) - cnf(1)*phi(285)
          work1(3) = -phi(416) - cnf(2)*phi(285)
          work1(4) = -phi(420) - cnf(3)*phi(285)
          work1(10) = -phi(565) - cnf(1)*phi(503)
          work1(11) = -phi(566) - cnf(1)*phi(504)
          work1(12) = -phi(567) - cnf(2)*phi(504)
          work1(13) = -phi(571) - cnf(1)*phi(508)
          work1(14) = -phi(572) - cnf(2)*phi(508)
          work1(15) = -phi(576) - cnf(3)*phi(508)
          work1(7) = -phi(503) - cnf(1)*phi(376)
          work1(8) = -phi(504) - cnf(2)*phi(376)
          work1(9) = -phi(508) - cnf(3)*phi(376)
          work1(10) = work1(10) - cnf(4)*phi(246)
          work1(12) = work1(12) - cnf(4)*phi(246)
          work1(15) = work1(15) - cnf(4)*phi(246)
          iii = 15
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(416) - cnf(1)*phi(286)
          work1(3) = -phi(417) - cnf(2)*phi(286)
          work1(4) = -phi(421) - cnf(3)*phi(286)
          work1(10) = -phi(566) - cnf(1)*phi(504)
          work1(11) = -phi(567) - cnf(1)*phi(505)
          work1(12) = -phi(568) - cnf(2)*phi(505)
          work1(13) = -phi(572) - cnf(1)*phi(509)
          work1(14) = -phi(573) - cnf(2)*phi(509)
          work1(15) = -phi(577) - cnf(3)*phi(509)
          work1(7) = -phi(504) - cnf(1)*phi(377)
          work1(8) = -phi(505) - cnf(2)*phi(377)
          work1(9) = -phi(509) - cnf(3)*phi(377)
          work1(10) = work1(10) - cnf(4)*phi(247)
          work1(12) = work1(12) - cnf(4)*phi(247)
          work1(15) = work1(15) - cnf(4)*phi(247)
          iii = 16
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(417) - cnf(1)*phi(287)
          work1(3) = -phi(418) - cnf(2)*phi(287)
          work1(4) = -phi(422) - cnf(3)*phi(287)
          work1(10) = -phi(567) - cnf(1)*phi(505)
          work1(11) = -phi(568) - cnf(1)*phi(506)
          work1(12) = -phi(569) - cnf(2)*phi(506)
          work1(13) = -phi(573) - cnf(1)*phi(510)
          work1(14) = -phi(574) - cnf(2)*phi(510)
          work1(15) = -phi(578) - cnf(3)*phi(510)
          work1(7) = -phi(505) - cnf(1)*phi(378)
          work1(8) = -phi(506) - cnf(2)*phi(378)
          work1(9) = -phi(510) - cnf(3)*phi(378)
          work1(10) = work1(10) - cnf(4)*phi(248)
          work1(12) = work1(12) - cnf(4)*phi(248)
          work1(15) = work1(15) - cnf(4)*phi(248)
          iii = 17
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(418) - cnf(1)*phi(288)
          work1(3) = -phi(419) - cnf(2)*phi(288)
          work1(4) = -phi(423) - cnf(3)*phi(288)
          work1(10) = -phi(568) - cnf(1)*phi(506)
          work1(11) = -phi(569) - cnf(1)*phi(507)
          work1(12) = -phi(570) - cnf(2)*phi(507)
          work1(13) = -phi(574) - cnf(1)*phi(511)
          work1(14) = -phi(575) - cnf(2)*phi(511)
          work1(15) = -phi(579) - cnf(3)*phi(511)
          work1(7) = -phi(506) - cnf(1)*phi(379)
          work1(8) = -phi(507) - cnf(2)*phi(379)
          work1(9) = -phi(511) - cnf(3)*phi(379)
          work1(10) = work1(10) - cnf(4)*phi(249)
          work1(12) = work1(12) - cnf(4)*phi(249)
          work1(15) = work1(15) - cnf(4)*phi(249)
          iii = 18
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(420) - cnf(1)*phi(289)
          work1(3) = -phi(421) - cnf(2)*phi(289)
          work1(4) = -phi(424) - cnf(3)*phi(289)
          work1(10) = -phi(571) - cnf(1)*phi(508)
          work1(11) = -phi(572) - cnf(1)*phi(509)
          work1(12) = -phi(573) - cnf(2)*phi(509)
          work1(13) = -phi(576) - cnf(1)*phi(512)
          work1(14) = -phi(577) - cnf(2)*phi(512)
          work1(15) = -phi(580) - cnf(3)*phi(512)
          work1(7) = -phi(508) - cnf(1)*phi(380)
          work1(8) = -phi(509) - cnf(2)*phi(380)
          work1(9) = -phi(512) - cnf(3)*phi(380)
          work1(10) = work1(10) - cnf(4)*phi(250)
          work1(12) = work1(12) - cnf(4)*phi(250)
          work1(15) = work1(15) - cnf(4)*phi(250)
          iii = 19
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(421) - cnf(1)*phi(290)
          work1(3) = -phi(422) - cnf(2)*phi(290)
          work1(4) = -phi(425) - cnf(3)*phi(290)
          work1(10) = -phi(572) - cnf(1)*phi(509)
          work1(11) = -phi(573) - cnf(1)*phi(510)
          work1(12) = -phi(574) - cnf(2)*phi(510)
          work1(13) = -phi(577) - cnf(1)*phi(513)
          work1(14) = -phi(578) - cnf(2)*phi(513)
          work1(15) = -phi(581) - cnf(3)*phi(513)
          work1(7) = -phi(509) - cnf(1)*phi(381)
          work1(8) = -phi(510) - cnf(2)*phi(381)
          work1(9) = -phi(513) - cnf(3)*phi(381)
          work1(10) = work1(10) - cnf(4)*phi(251)
          work1(12) = work1(12) - cnf(4)*phi(251)
          work1(15) = work1(15) - cnf(4)*phi(251)
          iii = 20
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(422) - cnf(1)*phi(291)
          work1(3) = -phi(423) - cnf(2)*phi(291)
          work1(4) = -phi(426) - cnf(3)*phi(291)
          work1(10) = -phi(573) - cnf(1)*phi(510)
          work1(11) = -phi(574) - cnf(1)*phi(511)
          work1(12) = -phi(575) - cnf(2)*phi(511)
          work1(13) = -phi(578) - cnf(1)*phi(514)
          work1(14) = -phi(579) - cnf(2)*phi(514)
          work1(15) = -phi(582) - cnf(3)*phi(514)
          work1(7) = -phi(510) - cnf(1)*phi(382)
          work1(8) = -phi(511) - cnf(2)*phi(382)
          work1(9) = -phi(514) - cnf(3)*phi(382)
          work1(10) = work1(10) - cnf(4)*phi(252)
          work1(12) = work1(12) - cnf(4)*phi(252)
          work1(15) = work1(15) - cnf(4)*phi(252)
          iii = 21
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(424) - cnf(1)*phi(292)
          work1(3) = -phi(425) - cnf(2)*phi(292)
          work1(4) = -phi(427) - cnf(3)*phi(292)
          work1(10) = -phi(576) - cnf(1)*phi(512)
          work1(11) = -phi(577) - cnf(1)*phi(513)
          work1(12) = -phi(578) - cnf(2)*phi(513)
          work1(13) = -phi(580) - cnf(1)*phi(515)
          work1(14) = -phi(581) - cnf(2)*phi(515)
          work1(15) = -phi(583) - cnf(3)*phi(515)
          work1(7) = -phi(512) - cnf(1)*phi(383)
          work1(8) = -phi(513) - cnf(2)*phi(383)
          work1(9) = -phi(515) - cnf(3)*phi(383)
          work1(10) = work1(10) - cnf(4)*phi(253)
          work1(12) = work1(12) - cnf(4)*phi(253)
          work1(15) = work1(15) - cnf(4)*phi(253)
          iii = 22
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(425) - cnf(1)*phi(293)
          work1(3) = -phi(426) - cnf(2)*phi(293)
          work1(4) = -phi(428) - cnf(3)*phi(293)
          work1(10) = -phi(577) - cnf(1)*phi(513)
          work1(11) = -phi(578) - cnf(1)*phi(514)
          work1(12) = -phi(579) - cnf(2)*phi(514)
          work1(13) = -phi(581) - cnf(1)*phi(516)
          work1(14) = -phi(582) - cnf(2)*phi(516)
          work1(15) = -phi(584) - cnf(3)*phi(516)
          work1(7) = -phi(513) - cnf(1)*phi(384)
          work1(8) = -phi(514) - cnf(2)*phi(384)
          work1(9) = -phi(516) - cnf(3)*phi(384)
          work1(10) = work1(10) - cnf(4)*phi(254)
          work1(12) = work1(12) - cnf(4)*phi(254)
          work1(15) = work1(15) - cnf(4)*phi(254)
          iii = 23
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(2) = -phi(427) - cnf(1)*phi(294)
          work1(3) = -phi(428) - cnf(2)*phi(294)
          work1(4) = -phi(429) - cnf(3)*phi(294)
          work1(10) = -phi(580) - cnf(1)*phi(515)
          work1(11) = -phi(581) - cnf(1)*phi(516)
          work1(12) = -phi(582) - cnf(2)*phi(516)
          work1(13) = -phi(583) - cnf(1)*phi(517)
          work1(14) = -phi(584) - cnf(2)*phi(517)
          work1(15) = -phi(585) - cnf(3)*phi(517)
          work1(7) = -phi(515) - cnf(1)*phi(385)
          work1(8) = -phi(516) - cnf(2)*phi(385)
          work1(9) = -phi(517) - cnf(3)*phi(385)
          work1(10) = work1(10) - cnf(4)*phi(255)
          work1(12) = work1(12) - cnf(4)*phi(255)
          work1(15) = work1(15) - cnf(4)*phi(255)
          iii = 24
          ! BCTE
          work1(10) = work1(10) + cnf(1)*work1(7)
          work1(11) = work1(11) + cnf(2)*work1(7)
          work1(12) = work1(12) + cnf(2)*work1(8)
          work1(13) = work1(13) + cnf(3)*work1(7)
          work1(14) = work1(14) + cnf(3)*work1(8)
          work1(15) = work1(15) + cnf(3)*work1(9)
          work1(16) = work1(10) - cnf(1)*work1(2)
          work1(17) = work1(11) - cnf(1)*work1(3)
          work1(18) = work1(13) - cnf(1)*work1(4)
          work1(19) = work1(11) - cnf(2)*work1(2)
          work1(20) = work1(12) - cnf(2)*work1(3)
          work1(21) = work1(14) - cnf(2)*work1(4)
          work1(22) = work1(13) - cnf(3)*work1(2)
          work1(23) = work1(14) - cnf(3)*work1(3)
          work1(24) = work1(15) - cnf(3)*work1(4)
          j = 15
          do i = 1, 9
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

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

          work2(15, 2) = work2(15, 2) + cnf(1)*work2(9, 2)
          work2(15, 2) = work2(15, 2) + cnf(5)*work2(6, 2)
          work2(16, 2) = work2(16, 2) + cnf(1)*work2(10, 2)
          work2(16, 2) = work2(16, 2) + cnf(4)*work2(7, 2)
          work2(17, 2) = work2(17, 2) + cnf(1)*work2(11, 2)
          work2(18, 2) = work2(18, 2) + cnf(2)*work2(11, 2)
          work2(18, 2) = work2(18, 2) + cnf(5)*work2(7, 2)
          work2(19, 2) = work2(19, 2) + cnf(1)*work2(12, 2)
          work2(19, 2) = work2(19, 2) + cnf(4)*work2(8, 2)
          work2(20, 2) = work2(20, 2) + cnf(1)*work2(13, 2)
          work2(21, 2) = work2(21, 2) + cnf(2)*work2(13, 2)
          work2(21, 2) = work2(21, 2) + cnf(4)*work2(8, 2)
          work2(22, 2) = work2(22, 2) + cnf(1)*work2(14, 2)
          work2(23, 2) = work2(23, 2) + cnf(2)*work2(14, 2)
          work2(24, 2) = work2(24, 2) + cnf(3)*work2(14, 2)
          work2(24, 2) = work2(24, 2) + cnf(5)*work2(8, 2)
          work2(9, 2) = work2(9, 2) + cnf(1)*work2(3, 2)
          work2(9, 2) = work2(9, 2) + cnf(4)*work2(2, 2)
          work2(10, 2) = work2(10, 2) + cnf(1)*work2(4, 2)
          work2(11, 2) = work2(11, 2) + cnf(2)*work2(4, 2)
          work2(11, 2) = work2(11, 2) + cnf(4)*work2(2, 2)
          work2(12, 2) = work2(12, 2) + cnf(1)*work2(5, 2)
          work2(13, 2) = work2(13, 2) + cnf(2)*work2(5, 2)
          work2(14, 2) = work2(14, 2) + cnf(3)*work2(5, 2)
          work2(14, 2) = work2(14, 2) + cnf(4)*work2(2, 2)
          work2(3, 2) = work2(3, 2) + cnf(1)*work2(1, 2)
          work2(4, 2) = work2(4, 2) + cnf(2)*work2(1, 2)
          work2(5, 2) = work2(5, 2) + cnf(3)*work2(1, 2)
          work2(6, 2) = work2(6, 2) + cnf(1)*work2(2, 2)
          work2(7, 2) = work2(7, 2) + cnf(2)*work2(2, 2)
          work2(8, 2) = work2(8, 2) + cnf(3)*work2(2, 2)
          work2(15, 2) = work2(15, 2) + cnf(1)*work2(9, 2)
          work2(15, 2) = work2(15, 2) + cnf(4)*work2(6, 2)
          work2(16, 2) = work2(16, 2) + cnf(1)*work2(10, 2)
          work2(17, 2) = work2(17, 2) + cnf(2)*work2(10, 2)
          work2(17, 2) = work2(17, 2) + cnf(4)*work2(6, 2)
          work2(18, 2) = work2(18, 2) + cnf(2)*work2(11, 2)
          work2(18, 2) = work2(18, 2) + cnf(4)*work2(7, 2)
          work2(19, 2) = work2(19, 2) + cnf(1)*work2(12, 2)
          work2(20, 2) = work2(20, 2) + cnf(2)*work2(12, 2)
          work2(21, 2) = work2(21, 2) + cnf(2)*work2(13, 2)
          work2(22, 2) = work2(22, 2) + cnf(3)*work2(12, 2)
          work2(22, 2) = work2(22, 2) + cnf(4)*work2(6, 2)
          work2(23, 2) = work2(23, 2) + cnf(3)*work2(13, 2)
          work2(23, 2) = work2(23, 2) + cnf(4)*work2(7, 2)
          work2(24, 2) = work2(24, 2) + cnf(3)*work2(14, 2)
          work2(24, 2) = work2(24, 2) + cnf(4)*work2(8, 2)
          work2(9, 2) = work2(9, 2) + cnf(1)*work2(3, 2)
          work2(10, 2) = work2(10, 2) + cnf(2)*work2(3, 2)
          work2(11, 2) = work2(11, 2) + cnf(2)*work2(4, 2)
          work2(12, 2) = work2(12, 2) + cnf(3)*work2(3, 2)
          work2(13, 2) = work2(13, 2) + cnf(3)*work2(4, 2)
          work2(14, 2) = work2(14, 2) + cnf(3)*work2(5, 2)
          work2(15, 2) = work2(15, 2) + cnf(1)*work2(9, 2)
          work2(16, 2) = work2(16, 2) + cnf(2)*work2(9, 2)
          work2(17, 2) = work2(17, 2) + cnf(2)*work2(10, 2)
          work2(18, 2) = work2(18, 2) + cnf(2)*work2(11, 2)
          work2(19, 2) = work2(19, 2) + cnf(3)*work2(9, 2)
          work2(20, 2) = work2(20, 2) + cnf(3)*work2(10, 2)
          work2(21, 2) = work2(21, 2) + cnf(3)*work2(11, 2)
          work2(22, 2) = work2(22, 2) + cnf(3)*work2(12, 2)
          work2(23, 2) = work2(23, 2) + cnf(3)*work2(13, 2)
          work2(24, 2) = work2(24, 2) + cnf(3)*work2(14, 2)

          work2(15, 3) = work2(15, 3) + cnf(1)*work2(9, 3)
          work2(15, 3) = work2(15, 3) + cnf(5)*work2(6, 3)
          work2(16, 3) = work2(16, 3) + cnf(1)*work2(10, 3)
          work2(16, 3) = work2(16, 3) + cnf(4)*work2(7, 3)
          work2(17, 3) = work2(17, 3) + cnf(1)*work2(11, 3)
          work2(18, 3) = work2(18, 3) + cnf(2)*work2(11, 3)
          work2(18, 3) = work2(18, 3) + cnf(5)*work2(7, 3)
          work2(19, 3) = work2(19, 3) + cnf(1)*work2(12, 3)
          work2(19, 3) = work2(19, 3) + cnf(4)*work2(8, 3)
          work2(20, 3) = work2(20, 3) + cnf(1)*work2(13, 3)
          work2(21, 3) = work2(21, 3) + cnf(2)*work2(13, 3)
          work2(21, 3) = work2(21, 3) + cnf(4)*work2(8, 3)
          work2(22, 3) = work2(22, 3) + cnf(1)*work2(14, 3)
          work2(23, 3) = work2(23, 3) + cnf(2)*work2(14, 3)
          work2(24, 3) = work2(24, 3) + cnf(3)*work2(14, 3)
          work2(24, 3) = work2(24, 3) + cnf(5)*work2(8, 3)
          work2(9, 3) = work2(9, 3) + cnf(1)*work2(3, 3)
          work2(9, 3) = work2(9, 3) + cnf(4)*work2(2, 3)
          work2(10, 3) = work2(10, 3) + cnf(1)*work2(4, 3)
          work2(11, 3) = work2(11, 3) + cnf(2)*work2(4, 3)
          work2(11, 3) = work2(11, 3) + cnf(4)*work2(2, 3)
          work2(12, 3) = work2(12, 3) + cnf(1)*work2(5, 3)
          work2(13, 3) = work2(13, 3) + cnf(2)*work2(5, 3)
          work2(14, 3) = work2(14, 3) + cnf(3)*work2(5, 3)
          work2(14, 3) = work2(14, 3) + cnf(4)*work2(2, 3)
          work2(3, 3) = work2(3, 3) + cnf(1)*work2(1, 3)
          work2(4, 3) = work2(4, 3) + cnf(2)*work2(1, 3)
          work2(5, 3) = work2(5, 3) + cnf(3)*work2(1, 3)
          work2(6, 3) = work2(6, 3) + cnf(1)*work2(2, 3)
          work2(7, 3) = work2(7, 3) + cnf(2)*work2(2, 3)
          work2(8, 3) = work2(8, 3) + cnf(3)*work2(2, 3)
          work2(15, 3) = work2(15, 3) + cnf(1)*work2(9, 3)
          work2(15, 3) = work2(15, 3) + cnf(4)*work2(6, 3)
          work2(16, 3) = work2(16, 3) + cnf(1)*work2(10, 3)
          work2(17, 3) = work2(17, 3) + cnf(2)*work2(10, 3)
          work2(17, 3) = work2(17, 3) + cnf(4)*work2(6, 3)
          work2(18, 3) = work2(18, 3) + cnf(2)*work2(11, 3)
          work2(18, 3) = work2(18, 3) + cnf(4)*work2(7, 3)
          work2(19, 3) = work2(19, 3) + cnf(1)*work2(12, 3)
          work2(20, 3) = work2(20, 3) + cnf(2)*work2(12, 3)
          work2(21, 3) = work2(21, 3) + cnf(2)*work2(13, 3)
          work2(22, 3) = work2(22, 3) + cnf(3)*work2(12, 3)
          work2(22, 3) = work2(22, 3) + cnf(4)*work2(6, 3)
          work2(23, 3) = work2(23, 3) + cnf(3)*work2(13, 3)
          work2(23, 3) = work2(23, 3) + cnf(4)*work2(7, 3)
          work2(24, 3) = work2(24, 3) + cnf(3)*work2(14, 3)
          work2(24, 3) = work2(24, 3) + cnf(4)*work2(8, 3)
          work2(9, 3) = work2(9, 3) + cnf(1)*work2(3, 3)
          work2(10, 3) = work2(10, 3) + cnf(2)*work2(3, 3)
          work2(11, 3) = work2(11, 3) + cnf(2)*work2(4, 3)
          work2(12, 3) = work2(12, 3) + cnf(3)*work2(3, 3)
          work2(13, 3) = work2(13, 3) + cnf(3)*work2(4, 3)
          work2(14, 3) = work2(14, 3) + cnf(3)*work2(5, 3)
          work2(15, 3) = work2(15, 3) + cnf(1)*work2(9, 3)
          work2(16, 3) = work2(16, 3) + cnf(2)*work2(9, 3)
          work2(17, 3) = work2(17, 3) + cnf(2)*work2(10, 3)
          work2(18, 3) = work2(18, 3) + cnf(2)*work2(11, 3)
          work2(19, 3) = work2(19, 3) + cnf(3)*work2(9, 3)
          work2(20, 3) = work2(20, 3) + cnf(3)*work2(10, 3)
          work2(21, 3) = work2(21, 3) + cnf(3)*work2(11, 3)
          work2(22, 3) = work2(22, 3) + cnf(3)*work2(12, 3)
          work2(23, 3) = work2(23, 3) + cnf(3)*work2(13, 3)
          work2(24, 3) = work2(24, 3) + cnf(3)*work2(14, 3)

          work2(15, 4) = work2(15, 4) + cnf(1)*work2(9, 4)
          work2(15, 4) = work2(15, 4) + cnf(5)*work2(6, 4)
          work2(16, 4) = work2(16, 4) + cnf(1)*work2(10, 4)
          work2(16, 4) = work2(16, 4) + cnf(4)*work2(7, 4)
          work2(17, 4) = work2(17, 4) + cnf(1)*work2(11, 4)
          work2(18, 4) = work2(18, 4) + cnf(2)*work2(11, 4)
          work2(18, 4) = work2(18, 4) + cnf(5)*work2(7, 4)
          work2(19, 4) = work2(19, 4) + cnf(1)*work2(12, 4)
          work2(19, 4) = work2(19, 4) + cnf(4)*work2(8, 4)
          work2(20, 4) = work2(20, 4) + cnf(1)*work2(13, 4)
          work2(21, 4) = work2(21, 4) + cnf(2)*work2(13, 4)
          work2(21, 4) = work2(21, 4) + cnf(4)*work2(8, 4)
          work2(22, 4) = work2(22, 4) + cnf(1)*work2(14, 4)
          work2(23, 4) = work2(23, 4) + cnf(2)*work2(14, 4)
          work2(24, 4) = work2(24, 4) + cnf(3)*work2(14, 4)
          work2(24, 4) = work2(24, 4) + cnf(5)*work2(8, 4)
          work2(9, 4) = work2(9, 4) + cnf(1)*work2(3, 4)
          work2(9, 4) = work2(9, 4) + cnf(4)*work2(2, 4)
          work2(10, 4) = work2(10, 4) + cnf(1)*work2(4, 4)
          work2(11, 4) = work2(11, 4) + cnf(2)*work2(4, 4)
          work2(11, 4) = work2(11, 4) + cnf(4)*work2(2, 4)
          work2(12, 4) = work2(12, 4) + cnf(1)*work2(5, 4)
          work2(13, 4) = work2(13, 4) + cnf(2)*work2(5, 4)
          work2(14, 4) = work2(14, 4) + cnf(3)*work2(5, 4)
          work2(14, 4) = work2(14, 4) + cnf(4)*work2(2, 4)
          work2(3, 4) = work2(3, 4) + cnf(1)*work2(1, 4)
          work2(4, 4) = work2(4, 4) + cnf(2)*work2(1, 4)
          work2(5, 4) = work2(5, 4) + cnf(3)*work2(1, 4)
          work2(6, 4) = work2(6, 4) + cnf(1)*work2(2, 4)
          work2(7, 4) = work2(7, 4) + cnf(2)*work2(2, 4)
          work2(8, 4) = work2(8, 4) + cnf(3)*work2(2, 4)
          work2(15, 4) = work2(15, 4) + cnf(1)*work2(9, 4)
          work2(15, 4) = work2(15, 4) + cnf(4)*work2(6, 4)
          work2(16, 4) = work2(16, 4) + cnf(1)*work2(10, 4)
          work2(17, 4) = work2(17, 4) + cnf(2)*work2(10, 4)
          work2(17, 4) = work2(17, 4) + cnf(4)*work2(6, 4)
          work2(18, 4) = work2(18, 4) + cnf(2)*work2(11, 4)
          work2(18, 4) = work2(18, 4) + cnf(4)*work2(7, 4)
          work2(19, 4) = work2(19, 4) + cnf(1)*work2(12, 4)
          work2(20, 4) = work2(20, 4) + cnf(2)*work2(12, 4)
          work2(21, 4) = work2(21, 4) + cnf(2)*work2(13, 4)
          work2(22, 4) = work2(22, 4) + cnf(3)*work2(12, 4)
          work2(22, 4) = work2(22, 4) + cnf(4)*work2(6, 4)
          work2(23, 4) = work2(23, 4) + cnf(3)*work2(13, 4)
          work2(23, 4) = work2(23, 4) + cnf(4)*work2(7, 4)
          work2(24, 4) = work2(24, 4) + cnf(3)*work2(14, 4)
          work2(24, 4) = work2(24, 4) + cnf(4)*work2(8, 4)
          work2(9, 4) = work2(9, 4) + cnf(1)*work2(3, 4)
          work2(10, 4) = work2(10, 4) + cnf(2)*work2(3, 4)
          work2(11, 4) = work2(11, 4) + cnf(2)*work2(4, 4)
          work2(12, 4) = work2(12, 4) + cnf(3)*work2(3, 4)
          work2(13, 4) = work2(13, 4) + cnf(3)*work2(4, 4)
          work2(14, 4) = work2(14, 4) + cnf(3)*work2(5, 4)
          work2(15, 4) = work2(15, 4) + cnf(1)*work2(9, 4)
          work2(16, 4) = work2(16, 4) + cnf(2)*work2(9, 4)
          work2(17, 4) = work2(17, 4) + cnf(2)*work2(10, 4)
          work2(18, 4) = work2(18, 4) + cnf(2)*work2(11, 4)
          work2(19, 4) = work2(19, 4) + cnf(3)*work2(9, 4)
          work2(20, 4) = work2(20, 4) + cnf(3)*work2(10, 4)
          work2(21, 4) = work2(21, 4) + cnf(3)*work2(11, 4)
          work2(22, 4) = work2(22, 4) + cnf(3)*work2(12, 4)
          work2(23, 4) = work2(23, 4) + cnf(3)*work2(13, 4)
          work2(24, 4) = work2(24, 4) + cnf(3)*work2(14, 4)

          work2(15, 5) = work2(15, 5) + cnf(1)*work2(9, 5)
          work2(15, 5) = work2(15, 5) + cnf(5)*work2(6, 5)
          work2(16, 5) = work2(16, 5) + cnf(1)*work2(10, 5)
          work2(16, 5) = work2(16, 5) + cnf(4)*work2(7, 5)
          work2(17, 5) = work2(17, 5) + cnf(1)*work2(11, 5)
          work2(18, 5) = work2(18, 5) + cnf(2)*work2(11, 5)
          work2(18, 5) = work2(18, 5) + cnf(5)*work2(7, 5)
          work2(19, 5) = work2(19, 5) + cnf(1)*work2(12, 5)
          work2(19, 5) = work2(19, 5) + cnf(4)*work2(8, 5)
          work2(20, 5) = work2(20, 5) + cnf(1)*work2(13, 5)
          work2(21, 5) = work2(21, 5) + cnf(2)*work2(13, 5)
          work2(21, 5) = work2(21, 5) + cnf(4)*work2(8, 5)
          work2(22, 5) = work2(22, 5) + cnf(1)*work2(14, 5)
          work2(23, 5) = work2(23, 5) + cnf(2)*work2(14, 5)
          work2(24, 5) = work2(24, 5) + cnf(3)*work2(14, 5)
          work2(24, 5) = work2(24, 5) + cnf(5)*work2(8, 5)
          work2(9, 5) = work2(9, 5) + cnf(1)*work2(3, 5)
          work2(9, 5) = work2(9, 5) + cnf(4)*work2(2, 5)
          work2(10, 5) = work2(10, 5) + cnf(1)*work2(4, 5)
          work2(11, 5) = work2(11, 5) + cnf(2)*work2(4, 5)
          work2(11, 5) = work2(11, 5) + cnf(4)*work2(2, 5)
          work2(12, 5) = work2(12, 5) + cnf(1)*work2(5, 5)
          work2(13, 5) = work2(13, 5) + cnf(2)*work2(5, 5)
          work2(14, 5) = work2(14, 5) + cnf(3)*work2(5, 5)
          work2(14, 5) = work2(14, 5) + cnf(4)*work2(2, 5)
          work2(3, 5) = work2(3, 5) + cnf(1)*work2(1, 5)
          work2(4, 5) = work2(4, 5) + cnf(2)*work2(1, 5)
          work2(5, 5) = work2(5, 5) + cnf(3)*work2(1, 5)
          work2(6, 5) = work2(6, 5) + cnf(1)*work2(2, 5)
          work2(7, 5) = work2(7, 5) + cnf(2)*work2(2, 5)
          work2(8, 5) = work2(8, 5) + cnf(3)*work2(2, 5)
          work2(15, 5) = work2(15, 5) + cnf(1)*work2(9, 5)
          work2(15, 5) = work2(15, 5) + cnf(4)*work2(6, 5)
          work2(16, 5) = work2(16, 5) + cnf(1)*work2(10, 5)
          work2(17, 5) = work2(17, 5) + cnf(2)*work2(10, 5)
          work2(17, 5) = work2(17, 5) + cnf(4)*work2(6, 5)
          work2(18, 5) = work2(18, 5) + cnf(2)*work2(11, 5)
          work2(18, 5) = work2(18, 5) + cnf(4)*work2(7, 5)
          work2(19, 5) = work2(19, 5) + cnf(1)*work2(12, 5)
          work2(20, 5) = work2(20, 5) + cnf(2)*work2(12, 5)
          work2(21, 5) = work2(21, 5) + cnf(2)*work2(13, 5)
          work2(22, 5) = work2(22, 5) + cnf(3)*work2(12, 5)
          work2(22, 5) = work2(22, 5) + cnf(4)*work2(6, 5)
          work2(23, 5) = work2(23, 5) + cnf(3)*work2(13, 5)
          work2(23, 5) = work2(23, 5) + cnf(4)*work2(7, 5)
          work2(24, 5) = work2(24, 5) + cnf(3)*work2(14, 5)
          work2(24, 5) = work2(24, 5) + cnf(4)*work2(8, 5)
          work2(9, 5) = work2(9, 5) + cnf(1)*work2(3, 5)
          work2(10, 5) = work2(10, 5) + cnf(2)*work2(3, 5)
          work2(11, 5) = work2(11, 5) + cnf(2)*work2(4, 5)
          work2(12, 5) = work2(12, 5) + cnf(3)*work2(3, 5)
          work2(13, 5) = work2(13, 5) + cnf(3)*work2(4, 5)
          work2(14, 5) = work2(14, 5) + cnf(3)*work2(5, 5)
          work2(15, 5) = work2(15, 5) + cnf(1)*work2(9, 5)
          work2(16, 5) = work2(16, 5) + cnf(2)*work2(9, 5)
          work2(17, 5) = work2(17, 5) + cnf(2)*work2(10, 5)
          work2(18, 5) = work2(18, 5) + cnf(2)*work2(11, 5)
          work2(19, 5) = work2(19, 5) + cnf(3)*work2(9, 5)
          work2(20, 5) = work2(20, 5) + cnf(3)*work2(10, 5)
          work2(21, 5) = work2(21, 5) + cnf(3)*work2(11, 5)
          work2(22, 5) = work2(22, 5) + cnf(3)*work2(12, 5)
          work2(23, 5) = work2(23, 5) + cnf(3)*work2(13, 5)
          work2(24, 5) = work2(24, 5) + cnf(3)*work2(14, 5)

          work2(15, 6) = work2(15, 6) + cnf(1)*work2(9, 6)
          work2(15, 6) = work2(15, 6) + cnf(5)*work2(6, 6)
          work2(16, 6) = work2(16, 6) + cnf(1)*work2(10, 6)
          work2(16, 6) = work2(16, 6) + cnf(4)*work2(7, 6)
          work2(17, 6) = work2(17, 6) + cnf(1)*work2(11, 6)
          work2(18, 6) = work2(18, 6) + cnf(2)*work2(11, 6)
          work2(18, 6) = work2(18, 6) + cnf(5)*work2(7, 6)
          work2(19, 6) = work2(19, 6) + cnf(1)*work2(12, 6)
          work2(19, 6) = work2(19, 6) + cnf(4)*work2(8, 6)
          work2(20, 6) = work2(20, 6) + cnf(1)*work2(13, 6)
          work2(21, 6) = work2(21, 6) + cnf(2)*work2(13, 6)
          work2(21, 6) = work2(21, 6) + cnf(4)*work2(8, 6)
          work2(22, 6) = work2(22, 6) + cnf(1)*work2(14, 6)
          work2(23, 6) = work2(23, 6) + cnf(2)*work2(14, 6)
          work2(24, 6) = work2(24, 6) + cnf(3)*work2(14, 6)
          work2(24, 6) = work2(24, 6) + cnf(5)*work2(8, 6)
          work2(9, 6) = work2(9, 6) + cnf(1)*work2(3, 6)
          work2(9, 6) = work2(9, 6) + cnf(4)*work2(2, 6)
          work2(10, 6) = work2(10, 6) + cnf(1)*work2(4, 6)
          work2(11, 6) = work2(11, 6) + cnf(2)*work2(4, 6)
          work2(11, 6) = work2(11, 6) + cnf(4)*work2(2, 6)
          work2(12, 6) = work2(12, 6) + cnf(1)*work2(5, 6)
          work2(13, 6) = work2(13, 6) + cnf(2)*work2(5, 6)
          work2(14, 6) = work2(14, 6) + cnf(3)*work2(5, 6)
          work2(14, 6) = work2(14, 6) + cnf(4)*work2(2, 6)
          work2(3, 6) = work2(3, 6) + cnf(1)*work2(1, 6)
          work2(4, 6) = work2(4, 6) + cnf(2)*work2(1, 6)
          work2(5, 6) = work2(5, 6) + cnf(3)*work2(1, 6)
          work2(6, 6) = work2(6, 6) + cnf(1)*work2(2, 6)
          work2(7, 6) = work2(7, 6) + cnf(2)*work2(2, 6)
          work2(8, 6) = work2(8, 6) + cnf(3)*work2(2, 6)
          work2(15, 6) = work2(15, 6) + cnf(1)*work2(9, 6)
          work2(15, 6) = work2(15, 6) + cnf(4)*work2(6, 6)
          work2(16, 6) = work2(16, 6) + cnf(1)*work2(10, 6)
          work2(17, 6) = work2(17, 6) + cnf(2)*work2(10, 6)
          work2(17, 6) = work2(17, 6) + cnf(4)*work2(6, 6)
          work2(18, 6) = work2(18, 6) + cnf(2)*work2(11, 6)
          work2(18, 6) = work2(18, 6) + cnf(4)*work2(7, 6)
          work2(19, 6) = work2(19, 6) + cnf(1)*work2(12, 6)
          work2(20, 6) = work2(20, 6) + cnf(2)*work2(12, 6)
          work2(21, 6) = work2(21, 6) + cnf(2)*work2(13, 6)
          work2(22, 6) = work2(22, 6) + cnf(3)*work2(12, 6)
          work2(22, 6) = work2(22, 6) + cnf(4)*work2(6, 6)
          work2(23, 6) = work2(23, 6) + cnf(3)*work2(13, 6)
          work2(23, 6) = work2(23, 6) + cnf(4)*work2(7, 6)
          work2(24, 6) = work2(24, 6) + cnf(3)*work2(14, 6)
          work2(24, 6) = work2(24, 6) + cnf(4)*work2(8, 6)
          work2(9, 6) = work2(9, 6) + cnf(1)*work2(3, 6)
          work2(10, 6) = work2(10, 6) + cnf(2)*work2(3, 6)
          work2(11, 6) = work2(11, 6) + cnf(2)*work2(4, 6)
          work2(12, 6) = work2(12, 6) + cnf(3)*work2(3, 6)
          work2(13, 6) = work2(13, 6) + cnf(3)*work2(4, 6)
          work2(14, 6) = work2(14, 6) + cnf(3)*work2(5, 6)
          work2(15, 6) = work2(15, 6) + cnf(1)*work2(9, 6)
          work2(16, 6) = work2(16, 6) + cnf(2)*work2(9, 6)
          work2(17, 6) = work2(17, 6) + cnf(2)*work2(10, 6)
          work2(18, 6) = work2(18, 6) + cnf(2)*work2(11, 6)
          work2(19, 6) = work2(19, 6) + cnf(3)*work2(9, 6)
          work2(20, 6) = work2(20, 6) + cnf(3)*work2(10, 6)
          work2(21, 6) = work2(21, 6) + cnf(3)*work2(11, 6)
          work2(22, 6) = work2(22, 6) + cnf(3)*work2(12, 6)
          work2(23, 6) = work2(23, 6) + cnf(3)*work2(13, 6)
          work2(24, 6) = work2(24, 6) + cnf(3)*work2(14, 6)

          work2(15, 7) = work2(15, 7) + cnf(1)*work2(9, 7)
          work2(15, 7) = work2(15, 7) + cnf(5)*work2(6, 7)
          work2(16, 7) = work2(16, 7) + cnf(1)*work2(10, 7)
          work2(16, 7) = work2(16, 7) + cnf(4)*work2(7, 7)
          work2(17, 7) = work2(17, 7) + cnf(1)*work2(11, 7)
          work2(18, 7) = work2(18, 7) + cnf(2)*work2(11, 7)
          work2(18, 7) = work2(18, 7) + cnf(5)*work2(7, 7)
          work2(19, 7) = work2(19, 7) + cnf(1)*work2(12, 7)
          work2(19, 7) = work2(19, 7) + cnf(4)*work2(8, 7)
          work2(20, 7) = work2(20, 7) + cnf(1)*work2(13, 7)
          work2(21, 7) = work2(21, 7) + cnf(2)*work2(13, 7)
          work2(21, 7) = work2(21, 7) + cnf(4)*work2(8, 7)
          work2(22, 7) = work2(22, 7) + cnf(1)*work2(14, 7)
          work2(23, 7) = work2(23, 7) + cnf(2)*work2(14, 7)
          work2(24, 7) = work2(24, 7) + cnf(3)*work2(14, 7)
          work2(24, 7) = work2(24, 7) + cnf(5)*work2(8, 7)
          work2(9, 7) = work2(9, 7) + cnf(1)*work2(3, 7)
          work2(9, 7) = work2(9, 7) + cnf(4)*work2(2, 7)
          work2(10, 7) = work2(10, 7) + cnf(1)*work2(4, 7)
          work2(11, 7) = work2(11, 7) + cnf(2)*work2(4, 7)
          work2(11, 7) = work2(11, 7) + cnf(4)*work2(2, 7)
          work2(12, 7) = work2(12, 7) + cnf(1)*work2(5, 7)
          work2(13, 7) = work2(13, 7) + cnf(2)*work2(5, 7)
          work2(14, 7) = work2(14, 7) + cnf(3)*work2(5, 7)
          work2(14, 7) = work2(14, 7) + cnf(4)*work2(2, 7)
          work2(3, 7) = work2(3, 7) + cnf(1)*work2(1, 7)
          work2(4, 7) = work2(4, 7) + cnf(2)*work2(1, 7)
          work2(5, 7) = work2(5, 7) + cnf(3)*work2(1, 7)
          work2(6, 7) = work2(6, 7) + cnf(1)*work2(2, 7)
          work2(7, 7) = work2(7, 7) + cnf(2)*work2(2, 7)
          work2(8, 7) = work2(8, 7) + cnf(3)*work2(2, 7)
          work2(15, 7) = work2(15, 7) + cnf(1)*work2(9, 7)
          work2(15, 7) = work2(15, 7) + cnf(4)*work2(6, 7)
          work2(16, 7) = work2(16, 7) + cnf(1)*work2(10, 7)
          work2(17, 7) = work2(17, 7) + cnf(2)*work2(10, 7)
          work2(17, 7) = work2(17, 7) + cnf(4)*work2(6, 7)
          work2(18, 7) = work2(18, 7) + cnf(2)*work2(11, 7)
          work2(18, 7) = work2(18, 7) + cnf(4)*work2(7, 7)
          work2(19, 7) = work2(19, 7) + cnf(1)*work2(12, 7)
          work2(20, 7) = work2(20, 7) + cnf(2)*work2(12, 7)
          work2(21, 7) = work2(21, 7) + cnf(2)*work2(13, 7)
          work2(22, 7) = work2(22, 7) + cnf(3)*work2(12, 7)
          work2(22, 7) = work2(22, 7) + cnf(4)*work2(6, 7)
          work2(23, 7) = work2(23, 7) + cnf(3)*work2(13, 7)
          work2(23, 7) = work2(23, 7) + cnf(4)*work2(7, 7)
          work2(24, 7) = work2(24, 7) + cnf(3)*work2(14, 7)
          work2(24, 7) = work2(24, 7) + cnf(4)*work2(8, 7)
          work2(9, 7) = work2(9, 7) + cnf(1)*work2(3, 7)
          work2(10, 7) = work2(10, 7) + cnf(2)*work2(3, 7)
          work2(11, 7) = work2(11, 7) + cnf(2)*work2(4, 7)
          work2(12, 7) = work2(12, 7) + cnf(3)*work2(3, 7)
          work2(13, 7) = work2(13, 7) + cnf(3)*work2(4, 7)
          work2(14, 7) = work2(14, 7) + cnf(3)*work2(5, 7)
          work2(15, 7) = work2(15, 7) + cnf(1)*work2(9, 7)
          work2(16, 7) = work2(16, 7) + cnf(2)*work2(9, 7)
          work2(17, 7) = work2(17, 7) + cnf(2)*work2(10, 7)
          work2(18, 7) = work2(18, 7) + cnf(2)*work2(11, 7)
          work2(19, 7) = work2(19, 7) + cnf(3)*work2(9, 7)
          work2(20, 7) = work2(20, 7) + cnf(3)*work2(10, 7)
          work2(21, 7) = work2(21, 7) + cnf(3)*work2(11, 7)
          work2(22, 7) = work2(22, 7) + cnf(3)*work2(12, 7)
          work2(23, 7) = work2(23, 7) + cnf(3)*work2(13, 7)
          work2(24, 7) = work2(24, 7) + cnf(3)*work2(14, 7)

          work2(15, 8) = work2(15, 8) + cnf(1)*work2(9, 8)
          work2(15, 8) = work2(15, 8) + cnf(5)*work2(6, 8)
          work2(16, 8) = work2(16, 8) + cnf(1)*work2(10, 8)
          work2(16, 8) = work2(16, 8) + cnf(4)*work2(7, 8)
          work2(17, 8) = work2(17, 8) + cnf(1)*work2(11, 8)
          work2(18, 8) = work2(18, 8) + cnf(2)*work2(11, 8)
          work2(18, 8) = work2(18, 8) + cnf(5)*work2(7, 8)
          work2(19, 8) = work2(19, 8) + cnf(1)*work2(12, 8)
          work2(19, 8) = work2(19, 8) + cnf(4)*work2(8, 8)
          work2(20, 8) = work2(20, 8) + cnf(1)*work2(13, 8)
          work2(21, 8) = work2(21, 8) + cnf(2)*work2(13, 8)
          work2(21, 8) = work2(21, 8) + cnf(4)*work2(8, 8)
          work2(22, 8) = work2(22, 8) + cnf(1)*work2(14, 8)
          work2(23, 8) = work2(23, 8) + cnf(2)*work2(14, 8)
          work2(24, 8) = work2(24, 8) + cnf(3)*work2(14, 8)
          work2(24, 8) = work2(24, 8) + cnf(5)*work2(8, 8)
          work2(9, 8) = work2(9, 8) + cnf(1)*work2(3, 8)
          work2(9, 8) = work2(9, 8) + cnf(4)*work2(2, 8)
          work2(10, 8) = work2(10, 8) + cnf(1)*work2(4, 8)
          work2(11, 8) = work2(11, 8) + cnf(2)*work2(4, 8)
          work2(11, 8) = work2(11, 8) + cnf(4)*work2(2, 8)
          work2(12, 8) = work2(12, 8) + cnf(1)*work2(5, 8)
          work2(13, 8) = work2(13, 8) + cnf(2)*work2(5, 8)
          work2(14, 8) = work2(14, 8) + cnf(3)*work2(5, 8)
          work2(14, 8) = work2(14, 8) + cnf(4)*work2(2, 8)
          work2(3, 8) = work2(3, 8) + cnf(1)*work2(1, 8)
          work2(4, 8) = work2(4, 8) + cnf(2)*work2(1, 8)
          work2(5, 8) = work2(5, 8) + cnf(3)*work2(1, 8)
          work2(6, 8) = work2(6, 8) + cnf(1)*work2(2, 8)
          work2(7, 8) = work2(7, 8) + cnf(2)*work2(2, 8)
          work2(8, 8) = work2(8, 8) + cnf(3)*work2(2, 8)
          work2(15, 8) = work2(15, 8) + cnf(1)*work2(9, 8)
          work2(15, 8) = work2(15, 8) + cnf(4)*work2(6, 8)
          work2(16, 8) = work2(16, 8) + cnf(1)*work2(10, 8)
          work2(17, 8) = work2(17, 8) + cnf(2)*work2(10, 8)
          work2(17, 8) = work2(17, 8) + cnf(4)*work2(6, 8)
          work2(18, 8) = work2(18, 8) + cnf(2)*work2(11, 8)
          work2(18, 8) = work2(18, 8) + cnf(4)*work2(7, 8)
          work2(19, 8) = work2(19, 8) + cnf(1)*work2(12, 8)
          work2(20, 8) = work2(20, 8) + cnf(2)*work2(12, 8)
          work2(21, 8) = work2(21, 8) + cnf(2)*work2(13, 8)
          work2(22, 8) = work2(22, 8) + cnf(3)*work2(12, 8)
          work2(22, 8) = work2(22, 8) + cnf(4)*work2(6, 8)
          work2(23, 8) = work2(23, 8) + cnf(3)*work2(13, 8)
          work2(23, 8) = work2(23, 8) + cnf(4)*work2(7, 8)
          work2(24, 8) = work2(24, 8) + cnf(3)*work2(14, 8)
          work2(24, 8) = work2(24, 8) + cnf(4)*work2(8, 8)
          work2(9, 8) = work2(9, 8) + cnf(1)*work2(3, 8)
          work2(10, 8) = work2(10, 8) + cnf(2)*work2(3, 8)
          work2(11, 8) = work2(11, 8) + cnf(2)*work2(4, 8)
          work2(12, 8) = work2(12, 8) + cnf(3)*work2(3, 8)
          work2(13, 8) = work2(13, 8) + cnf(3)*work2(4, 8)
          work2(14, 8) = work2(14, 8) + cnf(3)*work2(5, 8)
          work2(15, 8) = work2(15, 8) + cnf(1)*work2(9, 8)
          work2(16, 8) = work2(16, 8) + cnf(2)*work2(9, 8)
          work2(17, 8) = work2(17, 8) + cnf(2)*work2(10, 8)
          work2(18, 8) = work2(18, 8) + cnf(2)*work2(11, 8)
          work2(19, 8) = work2(19, 8) + cnf(3)*work2(9, 8)
          work2(20, 8) = work2(20, 8) + cnf(3)*work2(10, 8)
          work2(21, 8) = work2(21, 8) + cnf(3)*work2(11, 8)
          work2(22, 8) = work2(22, 8) + cnf(3)*work2(12, 8)
          work2(23, 8) = work2(23, 8) + cnf(3)*work2(13, 8)
          work2(24, 8) = work2(24, 8) + cnf(3)*work2(14, 8)

          work2(15, 9) = work2(15, 9) + cnf(1)*work2(9, 9)
          work2(15, 9) = work2(15, 9) + cnf(5)*work2(6, 9)
          work2(16, 9) = work2(16, 9) + cnf(1)*work2(10, 9)
          work2(16, 9) = work2(16, 9) + cnf(4)*work2(7, 9)
          work2(17, 9) = work2(17, 9) + cnf(1)*work2(11, 9)
          work2(18, 9) = work2(18, 9) + cnf(2)*work2(11, 9)
          work2(18, 9) = work2(18, 9) + cnf(5)*work2(7, 9)
          work2(19, 9) = work2(19, 9) + cnf(1)*work2(12, 9)
          work2(19, 9) = work2(19, 9) + cnf(4)*work2(8, 9)
          work2(20, 9) = work2(20, 9) + cnf(1)*work2(13, 9)
          work2(21, 9) = work2(21, 9) + cnf(2)*work2(13, 9)
          work2(21, 9) = work2(21, 9) + cnf(4)*work2(8, 9)
          work2(22, 9) = work2(22, 9) + cnf(1)*work2(14, 9)
          work2(23, 9) = work2(23, 9) + cnf(2)*work2(14, 9)
          work2(24, 9) = work2(24, 9) + cnf(3)*work2(14, 9)
          work2(24, 9) = work2(24, 9) + cnf(5)*work2(8, 9)
          work2(9, 9) = work2(9, 9) + cnf(1)*work2(3, 9)
          work2(9, 9) = work2(9, 9) + cnf(4)*work2(2, 9)
          work2(10, 9) = work2(10, 9) + cnf(1)*work2(4, 9)
          work2(11, 9) = work2(11, 9) + cnf(2)*work2(4, 9)
          work2(11, 9) = work2(11, 9) + cnf(4)*work2(2, 9)
          work2(12, 9) = work2(12, 9) + cnf(1)*work2(5, 9)
          work2(13, 9) = work2(13, 9) + cnf(2)*work2(5, 9)
          work2(14, 9) = work2(14, 9) + cnf(3)*work2(5, 9)
          work2(14, 9) = work2(14, 9) + cnf(4)*work2(2, 9)
          work2(3, 9) = work2(3, 9) + cnf(1)*work2(1, 9)
          work2(4, 9) = work2(4, 9) + cnf(2)*work2(1, 9)
          work2(5, 9) = work2(5, 9) + cnf(3)*work2(1, 9)
          work2(6, 9) = work2(6, 9) + cnf(1)*work2(2, 9)
          work2(7, 9) = work2(7, 9) + cnf(2)*work2(2, 9)
          work2(8, 9) = work2(8, 9) + cnf(3)*work2(2, 9)
          work2(15, 9) = work2(15, 9) + cnf(1)*work2(9, 9)
          work2(15, 9) = work2(15, 9) + cnf(4)*work2(6, 9)
          work2(16, 9) = work2(16, 9) + cnf(1)*work2(10, 9)
          work2(17, 9) = work2(17, 9) + cnf(2)*work2(10, 9)
          work2(17, 9) = work2(17, 9) + cnf(4)*work2(6, 9)
          work2(18, 9) = work2(18, 9) + cnf(2)*work2(11, 9)
          work2(18, 9) = work2(18, 9) + cnf(4)*work2(7, 9)
          work2(19, 9) = work2(19, 9) + cnf(1)*work2(12, 9)
          work2(20, 9) = work2(20, 9) + cnf(2)*work2(12, 9)
          work2(21, 9) = work2(21, 9) + cnf(2)*work2(13, 9)
          work2(22, 9) = work2(22, 9) + cnf(3)*work2(12, 9)
          work2(22, 9) = work2(22, 9) + cnf(4)*work2(6, 9)
          work2(23, 9) = work2(23, 9) + cnf(3)*work2(13, 9)
          work2(23, 9) = work2(23, 9) + cnf(4)*work2(7, 9)
          work2(24, 9) = work2(24, 9) + cnf(3)*work2(14, 9)
          work2(24, 9) = work2(24, 9) + cnf(4)*work2(8, 9)
          work2(9, 9) = work2(9, 9) + cnf(1)*work2(3, 9)
          work2(10, 9) = work2(10, 9) + cnf(2)*work2(3, 9)
          work2(11, 9) = work2(11, 9) + cnf(2)*work2(4, 9)
          work2(12, 9) = work2(12, 9) + cnf(3)*work2(3, 9)
          work2(13, 9) = work2(13, 9) + cnf(3)*work2(4, 9)
          work2(14, 9) = work2(14, 9) + cnf(3)*work2(5, 9)
          work2(15, 9) = work2(15, 9) + cnf(1)*work2(9, 9)
          work2(16, 9) = work2(16, 9) + cnf(2)*work2(9, 9)
          work2(17, 9) = work2(17, 9) + cnf(2)*work2(10, 9)
          work2(18, 9) = work2(18, 9) + cnf(2)*work2(11, 9)
          work2(19, 9) = work2(19, 9) + cnf(3)*work2(9, 9)
          work2(20, 9) = work2(20, 9) + cnf(3)*work2(10, 9)
          work2(21, 9) = work2(21, 9) + cnf(3)*work2(11, 9)
          work2(22, 9) = work2(22, 9) + cnf(3)*work2(12, 9)
          work2(23, 9) = work2(23, 9) + cnf(3)*work2(13, 9)
          work2(24, 9) = work2(24, 9) + cnf(3)*work2(14, 9)

! ************************
! *                      *
! * Form final integrals *
! *                      *
! ************************

          eri_value(1) = work2(15, 1)
          eri_value(2) = work2(18, 1)
          eri_value(3) = work2(24, 1)
          eri_value(4) = work2(16, 1)*sqrt5
          eri_value(5) = work2(19, 1)*sqrt5
          eri_value(6) = work2(17, 1)*sqrt5
          eri_value(7) = work2(21, 1)*sqrt5
          eri_value(8) = work2(22, 1)*sqrt5
          eri_value(9) = work2(23, 1)*sqrt5
          eri_value(10) = work2(20, 1)*sqrt15
          eri_value(11) = work2(15, 4)
          eri_value(12) = work2(18, 4)
          eri_value(13) = work2(24, 4)
          eri_value(14) = work2(16, 4)*sqrt5
          eri_value(15) = work2(19, 4)*sqrt5
          eri_value(16) = work2(17, 4)*sqrt5
          eri_value(17) = work2(21, 4)*sqrt5
          eri_value(18) = work2(22, 4)*sqrt5
          eri_value(19) = work2(23, 4)*sqrt5
          eri_value(20) = work2(20, 4)*sqrt15
          eri_value(21) = work2(15, 7)
          eri_value(22) = work2(18, 7)
          eri_value(23) = work2(24, 7)
          eri_value(24) = work2(16, 7)*sqrt5
          eri_value(25) = work2(19, 7)*sqrt5
          eri_value(26) = work2(17, 7)*sqrt5
          eri_value(27) = work2(21, 7)*sqrt5
          eri_value(28) = work2(22, 7)*sqrt5
          eri_value(29) = work2(23, 7)*sqrt5
          eri_value(30) = work2(20, 7)*sqrt15
          eri_value(31) = work2(15, 2)
          eri_value(32) = work2(18, 2)
          eri_value(33) = work2(24, 2)
          eri_value(34) = work2(16, 2)*sqrt5
          eri_value(35) = work2(19, 2)*sqrt5
          eri_value(36) = work2(17, 2)*sqrt5
          eri_value(37) = work2(21, 2)*sqrt5
          eri_value(38) = work2(22, 2)*sqrt5
          eri_value(39) = work2(23, 2)*sqrt5
          eri_value(40) = work2(20, 2)*sqrt15
          eri_value(41) = work2(15, 5)
          eri_value(42) = work2(18, 5)
          eri_value(43) = work2(24, 5)
          eri_value(44) = work2(16, 5)*sqrt5
          eri_value(45) = work2(19, 5)*sqrt5
          eri_value(46) = work2(17, 5)*sqrt5
          eri_value(47) = work2(21, 5)*sqrt5
          eri_value(48) = work2(22, 5)*sqrt5
          eri_value(49) = work2(23, 5)*sqrt5
          eri_value(50) = work2(20, 5)*sqrt15
          eri_value(51) = work2(15, 8)
          eri_value(52) = work2(18, 8)
          eri_value(53) = work2(24, 8)
          eri_value(54) = work2(16, 8)*sqrt5
          eri_value(55) = work2(19, 8)*sqrt5
          eri_value(56) = work2(17, 8)*sqrt5
          eri_value(57) = work2(21, 8)*sqrt5
          eri_value(58) = work2(22, 8)*sqrt5
          eri_value(59) = work2(23, 8)*sqrt5
          eri_value(60) = work2(20, 8)*sqrt15
          eri_value(61) = work2(15, 3)
          eri_value(62) = work2(18, 3)
          eri_value(63) = work2(24, 3)
          eri_value(64) = work2(16, 3)*sqrt5
          eri_value(65) = work2(19, 3)*sqrt5
          eri_value(66) = work2(17, 3)*sqrt5
          eri_value(67) = work2(21, 3)*sqrt5
          eri_value(68) = work2(22, 3)*sqrt5
          eri_value(69) = work2(23, 3)*sqrt5
          eri_value(70) = work2(20, 3)*sqrt15
          eri_value(71) = work2(15, 6)
          eri_value(72) = work2(18, 6)
          eri_value(73) = work2(24, 6)
          eri_value(74) = work2(16, 6)*sqrt5
          eri_value(75) = work2(19, 6)*sqrt5
          eri_value(76) = work2(17, 6)*sqrt5
          eri_value(77) = work2(21, 6)*sqrt5
          eri_value(78) = work2(22, 6)*sqrt5
          eri_value(79) = work2(23, 6)*sqrt5
          eri_value(80) = work2(20, 6)*sqrt15
          eri_value(81) = work2(15, 9)
          eri_value(82) = work2(18, 9)
          eri_value(83) = work2(24, 9)
          eri_value(84) = work2(16, 9)*sqrt5
          eri_value(85) = work2(19, 9)*sqrt5
          eri_value(86) = work2(17, 9)*sqrt5
          eri_value(87) = work2(21, 9)*sqrt5
          eri_value(88) = work2(22, 9)*sqrt5
          eri_value(89) = work2(23, 9)*sqrt5
          eri_value(90) = work2(20, 9)*sqrt15

          maxi = 3
          maxj = 3
          maxk = 10

! ******************************
! *                            *
! * Digestion into Fock matrix *
! *                            *
! ******************************

          iandj = ish .eq. jsh

          loci = res%atom_loc(ish) - 1
          locj = res%atom_loc(jsh) - 1
          lock = res%atom_loc(ksh) - 1
          locl = res%atom_loc(lsh) - 1

          maxj2 = maxj
          nij = 0

          kstride = 1
          jstride = 10
          istride = 30

          ip = 1
          do i = 1, maxi
            if (iandj) maxj2 = i
            ii1 = i + loci
            ijp = ip
            ip = ip + istride

            do j = 1, maxj2
              nij = nij + 1
              jj1 = j + locj
              i2 = ii1
              j2 = jj1
              if (ii1 .lt. jj1) then ! sort <ij|
                i2 = jj1
                j2 = ii1
              end if

              ijkp = ijp
              ijp = ijp + jstride
              nkl = nij

              do k = 1, maxk
                kk1 = k + lock
                kk1 = k + lock

                ijklp = ijkp
                ijkp = ijkp + kstride

                l = 1
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
                !
              end do
            end do
          end do

        end if ! test.gt.cutoff_schwarz
      end do ! iquart
!$omp end target teams distribute parallel do

    end do ! itile

  end subroutine int1130
end submodule
