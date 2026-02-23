! The total angular momentum of this class is:           5
! The algorithm chosen is: PHR a.k.a. ERIC
! Writing an ERIC kernel
submodule(eric_kernels) int2030_impl
contains
  module subroutine int2030(sd_pair, sf_pair, density, fock, res)

    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: sd_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

! Variables for the class
    integer(kind=int64), allocatable :: n02bra(:), n03ket(:)
    real(dp), allocatable :: xint02bra(:), xint03ket(:)
    integer(kind=int64) :: nsdbra, nsfket
    real(dp) :: scutsdbra, scutsfket, test
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
    real(dp) :: ft(9), phi(456), c_factor(56), rxyz(3)
    real(dp) :: work2(24, 6), work1(11)
    real(dp) :: eri_value(60), angl(20)
    real(dp) :: ai, aij, aijk, aijkl
    integer(kind=int64) :: iord(20)
    integer(kind=int64) :: istride, kstride
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

    allocate (n02bra(res%n_s_shl*res%n_d_shl))
    allocate (xint02bra(res%n_s_shl*res%n_d_shl))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

! Start screening

    scutsdbra = cutoff_schwarz/maxval(sd_pair%xints)
    nsdbra = 0
    do ij = 1, res%n_s_shl*res%n_d_shl
      if (sd_pair%xints(ij) .ge. scutsdbra) then
        nsdbra = nsdbra + 1
        xint02bra(nsdbra) = sd_pair%xints(ij)
        n02bra(nsdbra) = ij
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

    if ((nsdbra*nsfket) .le. nchunksize_k10) nchunksize_k10 = nsdbra*nsfket
    ntile = int(nsdbra*nsfket/nchunksize_k10)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_k10 + 1
      iend = itile*nchunksize_k10
      if (itile .eq. ntile) iend = nsdbra*nsfket

! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

! Mappings to GPU

!$omp target teams distribute parallel do default(none) &
!$omp shared(res, density, fock, nquart_start, nquart_end, nsdbra, xint02bra, n02bra, xint03ket, n03ket, nsfket, sd_pair, sf_pair) &
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
!$omp private(loci,locj,lock,locl,nij,istride,kstride) &
!$omp private(ip,ii1,ijp,jj1,i2,j2,ijkp,nkl,kk1) &
!$omp private(ijklp,buff,ll1,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ii,jj,kk,ll)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsfket + 1
        kl_tmp = mod(iquart - 1, nsfket) + 1

        test = xint02bra(ij_tmp)*xint03ket(kl_tmp)

        if (test .gt. cutoff_schwarz) then

          ij = n02bra(ij_tmp)
          kl = n03ket(kl_tmp)

          ish_tmp = (ij - 1)/res%n_d_shl + 1
          jsh_tmp = mod(ij - 1, res%n_d_shl) + 1
          ksh_tmp = (kl - 1)/res%n_f_shl + 1
          lsh_tmp = mod(kl - 1, res%n_f_shl) + 1

          ish = res%i_d_shl(jsh_tmp)
          jsh = res%i_s_shl(ish_tmp)
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
          phi(90) = 0.0_dp
          phi(91) = 0.0_dp
          phi(92) = 0.0_dp
          phi(93) = 0.0_dp
          phi(94) = 0.0_dp
          phi(95) = 0.0_dp
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
          phi(121) = 0.0_dp
          phi(122) = 0.0_dp
          phi(123) = 0.0_dp
          phi(124) = 0.0_dp
          phi(125) = 0.0_dp
          phi(126) = 0.0_dp
          phi(127) = 0.0_dp
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
          phi(169) = 0.0_dp
          phi(177) = 0.0_dp
          phi(178) = 0.0_dp
          phi(179) = 0.0_dp
          phi(180) = 0.0_dp
          phi(181) = 0.0_dp
          phi(182) = 0.0_dp
          phi(183) = 0.0_dp
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
          phi(243) = 0.0_dp
          phi(244) = 0.0_dp
          phi(245) = 0.0_dp
          phi(246) = 0.0_dp
          phi(247) = 0.0_dp
          phi(248) = 0.0_dp
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
          phi(367) = 0.0_dp
          phi(368) = 0.0_dp
          phi(369) = 0.0_dp
          phi(370) = 0.0_dp
          phi(371) = 0.0_dp
          phi(372) = 0.0_dp
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
          phi(386) = 0.0_dp
          phi(387) = 0.0_dp
          phi(388) = 0.0_dp
          phi(423) = 0.0_dp
          phi(424) = 0.0_dp
          phi(425) = 0.0_dp
          phi(426) = 0.0_dp
          phi(427) = 0.0_dp
          phi(428) = 0.0_dp
          phi(429) = 0.0_dp
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
          phi(452) = 0.0_dp
          phi(453) = 0.0_dp
          phi(454) = 0.0_dp
          phi(455) = 0.0_dp
          phi(456) = 0.0_dp

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
            phi(88) = 0.0_dp
            phi(89) = 0.0_dp
            phi(96) = 0.0_dp
            phi(97) = 0.0_dp
            phi(98) = 0.0_dp
            phi(105) = 0.0_dp
            phi(106) = 0.0_dp
            phi(107) = 0.0_dp
            phi(114) = 0.0_dp
            phi(115) = 0.0_dp
            phi(116) = 0.0_dp
            phi(117) = 0.0_dp
            phi(118) = 0.0_dp
            phi(119) = 0.0_dp
            phi(120) = 0.0_dp
            phi(128) = 0.0_dp
            phi(129) = 0.0_dp
            phi(130) = 0.0_dp
            phi(131) = 0.0_dp
            phi(132) = 0.0_dp
            phi(133) = 0.0_dp
            phi(134) = 0.0_dp
            phi(149) = 0.0_dp
            phi(150) = 0.0_dp
            phi(151) = 0.0_dp
            phi(152) = 0.0_dp
            phi(153) = 0.0_dp
            phi(154) = 0.0_dp
            phi(155) = 0.0_dp
            phi(170) = 0.0_dp
            phi(171) = 0.0_dp
            phi(172) = 0.0_dp
            phi(173) = 0.0_dp
            phi(174) = 0.0_dp
            phi(175) = 0.0_dp
            phi(176) = 0.0_dp
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
            phi(221) = 0.0_dp
            phi(222) = 0.0_dp
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
            phi(287) = 0.0_dp
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
            phi(321) = 0.0_dp
            phi(322) = 0.0_dp
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
            phi(355) = 0.0_dp
            phi(356) = 0.0_dp
            phi(357) = 0.0_dp
            phi(358) = 0.0_dp
            phi(359) = 0.0_dp
            phi(360) = 0.0_dp
            phi(361) = 0.0_dp
            phi(362) = 0.0_dp
            phi(363) = 0.0_dp
            phi(364) = 0.0_dp
            phi(365) = 0.0_dp
            phi(366) = 0.0_dp
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

! Begin looping over ij primitives

            do i = 1, ijtop

              t_expon_ab = sd_pair%t_expon_ab(sd_pair%pair_loc(ij) + 1 + bra_loop) ! exp_c + exp_d
              t_expon_a = sd_pair%expon_b(sd_pair%pair_loc(ij) + 1 + bra_loop)
              t_expon_b = sd_pair%expon_a(sd_pair%pair_loc(ij) + 1 + bra_loop)

              t_inverse_expon_ab = 1.0_dp/t_expon_ab
              ccfbra = sd_pair%sq(sd_pair%pair_loc(ij) + 1 + bra_loop)*sqrt2_pi_5_4
              slbra = pi_1_4_div_sqrt2*sqrt(t_inverse_expon_ab)

              xij = ((t_expon_a*res%coord_sh(ish, 1)) + (t_expon_b*res%coord_sh(jsh, 1)))*t_inverse_expon_ab
              yij = ((t_expon_a*res%coord_sh(ish, 2)) + (t_expon_b*res%coord_sh(jsh, 2)))*t_inverse_expon_ab
              zij = ((t_expon_a*res%coord_sh(ish, 3)) + (t_expon_b*res%coord_sh(jsh, 3)))*t_inverse_expon_ab

              rxbra = sd_pair%t_inverse_expon_ab(sd_pair%pair_loc(ij) + 1 + bra_loop)! inverse_expon_ab*0.5_dp

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
              phi(84) = phi(84) + scale_factor2(3)*scale_factor1(3)*phi(1)

              phi(87) = phi(87) + scale_factor1(2)*phi(2)
              phi(88) = phi(88) + scale_factor1(2)*phi(3)
              phi(89) = phi(89) + scale_factor1(2)*phi(4)
              phi(96) = phi(96) + scale_factor2(2)*scale_factor1(3)*phi(2)
              phi(97) = phi(97) + scale_factor2(2)*scale_factor1(3)*phi(3)
              phi(98) = phi(98) + scale_factor2(2)*scale_factor1(3)*phi(4)
              phi(105) = phi(105) + scale_factor2(3)*scale_factor1(3)*phi(2)
              phi(106) = phi(106) + scale_factor2(3)*scale_factor1(3)*phi(3)
              phi(107) = phi(107) + scale_factor2(3)*scale_factor1(3)*phi(4)

              phi(114) = phi(114) + scale_factor1(2)*phi(5)
              phi(115) = phi(115) + scale_factor1(2)*phi(6)
              phi(116) = phi(116) + scale_factor1(2)*phi(7)
              phi(117) = phi(117) + scale_factor1(2)*phi(8)
              phi(118) = phi(118) + scale_factor1(2)*phi(9)
              phi(119) = phi(119) + scale_factor1(2)*phi(10)
              phi(120) = phi(120) + scale_factor1(2)*phi(11)
              phi(128) = phi(128) + scale_factor1(3)*phi(5)
              phi(129) = phi(129) + scale_factor1(3)*phi(6)
              phi(130) = phi(130) + scale_factor1(3)*phi(7)
              phi(131) = phi(131) + scale_factor1(3)*phi(8)
              phi(132) = phi(132) + scale_factor1(3)*phi(9)
              phi(133) = phi(133) + scale_factor1(3)*phi(10)
              phi(134) = phi(134) + scale_factor1(3)*phi(11)
              phi(149) = phi(149) + scale_factor2(2)*scale_factor1(3)*phi(5)
              phi(150) = phi(150) + scale_factor2(2)*scale_factor1(3)*phi(6)
              phi(151) = phi(151) + scale_factor2(2)*scale_factor1(3)*phi(7)
              phi(152) = phi(152) + scale_factor2(2)*scale_factor1(3)*phi(8)
              phi(153) = phi(153) + scale_factor2(2)*scale_factor1(3)*phi(9)
              phi(154) = phi(154) + scale_factor2(2)*scale_factor1(3)*phi(10)
              phi(155) = phi(155) + scale_factor2(2)*scale_factor1(3)*phi(11)
              phi(170) = phi(170) + scale_factor2(3)*scale_factor1(3)*phi(5)
              phi(171) = phi(171) + scale_factor2(3)*scale_factor1(3)*phi(6)
              phi(172) = phi(172) + scale_factor2(3)*scale_factor1(3)*phi(7)
              phi(173) = phi(173) + scale_factor2(3)*scale_factor1(3)*phi(8)
              phi(174) = phi(174) + scale_factor2(3)*scale_factor1(3)*phi(9)
              phi(175) = phi(175) + scale_factor2(3)*scale_factor1(3)*phi(10)
              phi(176) = phi(176) + scale_factor2(3)*scale_factor1(3)*phi(11)

              phi(184) = phi(184) + scale_factor1(2)*phi(12)
              phi(185) = phi(185) + scale_factor1(2)*phi(13)
              phi(186) = phi(186) + scale_factor1(2)*phi(14)
              phi(187) = phi(187) + scale_factor1(2)*phi(15)
              phi(188) = phi(188) + scale_factor1(2)*phi(16)
              phi(189) = phi(189) + scale_factor1(2)*phi(17)
              phi(190) = phi(190) + scale_factor1(2)*phi(18)
              phi(191) = phi(191) + scale_factor1(2)*phi(19)
              phi(192) = phi(192) + scale_factor1(2)*phi(20)
              phi(193) = phi(193) + scale_factor1(2)*phi(21)
              phi(194) = phi(194) + scale_factor1(2)*phi(22)
              phi(195) = phi(195) + scale_factor1(2)*phi(23)
              phi(196) = phi(196) + scale_factor1(2)*phi(24)
              phi(210) = phi(210) + scale_factor1(3)*phi(12)
              phi(211) = phi(211) + scale_factor1(3)*phi(13)
              phi(212) = phi(212) + scale_factor1(3)*phi(14)
              phi(213) = phi(213) + scale_factor1(3)*phi(15)
              phi(214) = phi(214) + scale_factor1(3)*phi(16)
              phi(215) = phi(215) + scale_factor1(3)*phi(17)
              phi(216) = phi(216) + scale_factor1(3)*phi(18)
              phi(217) = phi(217) + scale_factor1(3)*phi(19)
              phi(218) = phi(218) + scale_factor1(3)*phi(20)
              phi(219) = phi(219) + scale_factor1(3)*phi(21)
              phi(220) = phi(220) + scale_factor1(3)*phi(22)
              phi(221) = phi(221) + scale_factor1(3)*phi(23)
              phi(222) = phi(222) + scale_factor1(3)*phi(24)
              phi(249) = phi(249) + scale_factor2(2)*scale_factor1(3)*phi(12)
              phi(250) = phi(250) + scale_factor2(2)*scale_factor1(3)*phi(13)
              phi(251) = phi(251) + scale_factor2(2)*scale_factor1(3)*phi(14)
              phi(252) = phi(252) + scale_factor2(2)*scale_factor1(3)*phi(15)
              phi(253) = phi(253) + scale_factor2(2)*scale_factor1(3)*phi(16)
              phi(254) = phi(254) + scale_factor2(2)*scale_factor1(3)*phi(17)
              phi(255) = phi(255) + scale_factor2(2)*scale_factor1(3)*phi(18)
              phi(256) = phi(256) + scale_factor2(2)*scale_factor1(3)*phi(19)
              phi(257) = phi(257) + scale_factor2(2)*scale_factor1(3)*phi(20)
              phi(258) = phi(258) + scale_factor2(2)*scale_factor1(3)*phi(21)
              phi(259) = phi(259) + scale_factor2(2)*scale_factor1(3)*phi(22)
              phi(260) = phi(260) + scale_factor2(2)*scale_factor1(3)*phi(23)
              phi(261) = phi(261) + scale_factor2(2)*scale_factor1(3)*phi(24)
              phi(275) = phi(275) + scale_factor2(3)*scale_factor1(3)*phi(12)
              phi(276) = phi(276) + scale_factor2(3)*scale_factor1(3)*phi(13)
              phi(277) = phi(277) + scale_factor2(3)*scale_factor1(3)*phi(14)
              phi(278) = phi(278) + scale_factor2(3)*scale_factor1(3)*phi(15)
              phi(279) = phi(279) + scale_factor2(3)*scale_factor1(3)*phi(16)
              phi(280) = phi(280) + scale_factor2(3)*scale_factor1(3)*phi(17)
              phi(281) = phi(281) + scale_factor2(3)*scale_factor1(3)*phi(18)
              phi(282) = phi(282) + scale_factor2(3)*scale_factor1(3)*phi(19)
              phi(283) = phi(283) + scale_factor2(3)*scale_factor1(3)*phi(20)
              phi(284) = phi(284) + scale_factor2(3)*scale_factor1(3)*phi(21)
              phi(285) = phi(285) + scale_factor2(3)*scale_factor1(3)*phi(22)
              phi(286) = phi(286) + scale_factor2(3)*scale_factor1(3)*phi(23)
              phi(287) = phi(287) + scale_factor2(3)*scale_factor1(3)*phi(24)

              phi(301) = phi(301) + scale_factor1(3)*phi(25)
              phi(302) = phi(302) + scale_factor1(3)*phi(26)
              phi(303) = phi(303) + scale_factor1(3)*phi(27)
              phi(304) = phi(304) + scale_factor1(3)*phi(28)
              phi(305) = phi(305) + scale_factor1(3)*phi(29)
              phi(306) = phi(306) + scale_factor1(3)*phi(30)
              phi(307) = phi(307) + scale_factor1(3)*phi(31)
              phi(308) = phi(308) + scale_factor1(3)*phi(32)
              phi(309) = phi(309) + scale_factor1(3)*phi(33)
              phi(310) = phi(310) + scale_factor1(3)*phi(34)
              phi(311) = phi(311) + scale_factor1(3)*phi(35)
              phi(312) = phi(312) + scale_factor1(3)*phi(36)
              phi(313) = phi(313) + scale_factor1(3)*phi(37)
              phi(314) = phi(314) + scale_factor1(3)*phi(38)
              phi(315) = phi(315) + scale_factor1(3)*phi(39)
              phi(316) = phi(316) + scale_factor1(3)*phi(40)
              phi(317) = phi(317) + scale_factor1(3)*phi(41)
              phi(318) = phi(318) + scale_factor1(3)*phi(42)
              phi(319) = phi(319) + scale_factor1(3)*phi(43)
              phi(320) = phi(320) + scale_factor1(3)*phi(44)
              phi(321) = phi(321) + scale_factor1(3)*phi(45)
              phi(322) = phi(322) + scale_factor1(3)*phi(46)
              phi(345) = phi(345) + scale_factor2(2)*scale_factor1(3)*phi(25)
              phi(346) = phi(346) + scale_factor2(2)*scale_factor1(3)*phi(26)
              phi(347) = phi(347) + scale_factor2(2)*scale_factor1(3)*phi(27)
              phi(348) = phi(348) + scale_factor2(2)*scale_factor1(3)*phi(28)
              phi(349) = phi(349) + scale_factor2(2)*scale_factor1(3)*phi(29)
              phi(350) = phi(350) + scale_factor2(2)*scale_factor1(3)*phi(30)
              phi(351) = phi(351) + scale_factor2(2)*scale_factor1(3)*phi(31)
              phi(352) = phi(352) + scale_factor2(2)*scale_factor1(3)*phi(32)
              phi(353) = phi(353) + scale_factor2(2)*scale_factor1(3)*phi(33)
              phi(354) = phi(354) + scale_factor2(2)*scale_factor1(3)*phi(34)
              phi(355) = phi(355) + scale_factor2(2)*scale_factor1(3)*phi(35)
              phi(356) = phi(356) + scale_factor2(2)*scale_factor1(3)*phi(36)
              phi(357) = phi(357) + scale_factor2(2)*scale_factor1(3)*phi(37)
              phi(358) = phi(358) + scale_factor2(2)*scale_factor1(3)*phi(38)
              phi(359) = phi(359) + scale_factor2(2)*scale_factor1(3)*phi(39)
              phi(360) = phi(360) + scale_factor2(2)*scale_factor1(3)*phi(40)
              phi(361) = phi(361) + scale_factor2(2)*scale_factor1(3)*phi(41)
              phi(362) = phi(362) + scale_factor2(2)*scale_factor1(3)*phi(42)
              phi(363) = phi(363) + scale_factor2(2)*scale_factor1(3)*phi(43)
              phi(364) = phi(364) + scale_factor2(2)*scale_factor1(3)*phi(44)
              phi(365) = phi(365) + scale_factor2(2)*scale_factor1(3)*phi(45)
              phi(366) = phi(366) + scale_factor2(2)*scale_factor1(3)*phi(46)

              phi(389) = phi(389) + scale_factor1(3)*phi(47)
              phi(390) = phi(390) + scale_factor1(3)*phi(48)
              phi(391) = phi(391) + scale_factor1(3)*phi(49)
              phi(392) = phi(392) + scale_factor1(3)*phi(50)
              phi(393) = phi(393) + scale_factor1(3)*phi(51)
              phi(394) = phi(394) + scale_factor1(3)*phi(52)
              phi(395) = phi(395) + scale_factor1(3)*phi(53)
              phi(396) = phi(396) + scale_factor1(3)*phi(54)
              phi(397) = phi(397) + scale_factor1(3)*phi(55)
              phi(398) = phi(398) + scale_factor1(3)*phi(56)
              phi(399) = phi(399) + scale_factor1(3)*phi(57)
              phi(400) = phi(400) + scale_factor1(3)*phi(58)
              phi(401) = phi(401) + scale_factor1(3)*phi(59)
              phi(402) = phi(402) + scale_factor1(3)*phi(60)
              phi(403) = phi(403) + scale_factor1(3)*phi(61)
              phi(404) = phi(404) + scale_factor1(3)*phi(62)
              phi(405) = phi(405) + scale_factor1(3)*phi(63)
              phi(406) = phi(406) + scale_factor1(3)*phi(64)
              phi(407) = phi(407) + scale_factor1(3)*phi(65)
              phi(408) = phi(408) + scale_factor1(3)*phi(66)
              phi(409) = phi(409) + scale_factor1(3)*phi(67)
              phi(410) = phi(410) + scale_factor1(3)*phi(68)
              phi(411) = phi(411) + scale_factor1(3)*phi(69)
              phi(412) = phi(412) + scale_factor1(3)*phi(70)
              phi(413) = phi(413) + scale_factor1(3)*phi(71)
              phi(414) = phi(414) + scale_factor1(3)*phi(72)
              phi(415) = phi(415) + scale_factor1(3)*phi(73)
              phi(416) = phi(416) + scale_factor1(3)*phi(74)
              phi(417) = phi(417) + scale_factor1(3)*phi(75)
              phi(418) = phi(418) + scale_factor1(3)*phi(76)
              phi(419) = phi(419) + scale_factor1(3)*phi(77)
              phi(420) = phi(420) + scale_factor1(3)*phi(78)
              phi(421) = phi(421) + scale_factor1(3)*phi(79)
              phi(422) = phi(422) + scale_factor1(3)*phi(80)

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

            phi(90) = phi(90) + scale_factor1(3)*phi(87)
            phi(91) = phi(91) + scale_factor1(3)*phi(88)
            phi(92) = phi(92) + scale_factor1(3)*phi(89)
            phi(93) = phi(93) + scale_factor2(3)*scale_factor1(4)*phi(87)
            phi(94) = phi(94) + scale_factor2(3)*scale_factor1(4)*phi(88)
            phi(95) = phi(95) + scale_factor2(3)*scale_factor1(4)*phi(89)
            phi(99) = phi(99) + scale_factor2(2)*scale_factor1(3)*phi(96)
            phi(100) = phi(100) + scale_factor2(2)*scale_factor1(3)*phi(97)
            phi(101) = phi(101) + scale_factor2(2)*scale_factor1(3)*phi(98)
            phi(102) = phi(102) + scale_factor2(4)*scale_factor1(4)*phi(96)
            phi(103) = phi(103) + scale_factor2(4)*scale_factor1(4)*phi(97)
            phi(104) = phi(104) + scale_factor2(4)*scale_factor1(4)*phi(98)
            phi(108) = phi(108) + scale_factor1(3)*phi(105)
            phi(109) = phi(109) + scale_factor1(3)*phi(106)
            phi(110) = phi(110) + scale_factor1(3)*phi(107)
            phi(111) = phi(111) + scale_factor2(3)*scale_factor1(4)*phi(105)
            phi(112) = phi(112) + scale_factor2(3)*scale_factor1(4)*phi(106)
            phi(113) = phi(113) + scale_factor2(3)*scale_factor1(4)*phi(107)

            phi(121) = phi(121) + scale_factor2(2)*scale_factor1(4)*phi(114)
            phi(122) = phi(122) + scale_factor2(2)*scale_factor1(4)*phi(115)
            phi(123) = phi(123) + scale_factor2(2)*scale_factor1(4)*phi(116)
            phi(124) = phi(124) + scale_factor2(2)*scale_factor1(4)*phi(117)
            phi(125) = phi(125) + scale_factor2(2)*scale_factor1(4)*phi(118)
            phi(126) = phi(126) + scale_factor2(2)*scale_factor1(4)*phi(119)
            phi(127) = phi(127) + scale_factor2(2)*scale_factor1(4)*phi(120)
            phi(135) = phi(135) + scale_factor2(2)*scale_factor1(3)*phi(128)
            phi(136) = phi(136) + scale_factor2(2)*scale_factor1(3)*phi(129)
            phi(137) = phi(137) + scale_factor2(2)*scale_factor1(3)*phi(130)
            phi(138) = phi(138) + scale_factor2(2)*scale_factor1(3)*phi(131)
            phi(139) = phi(139) + scale_factor2(2)*scale_factor1(3)*phi(132)
            phi(140) = phi(140) + scale_factor2(2)*scale_factor1(3)*phi(133)
            phi(141) = phi(141) + scale_factor2(2)*scale_factor1(3)*phi(134)
            phi(142) = phi(142) + scale_factor2(4)*scale_factor1(4)*phi(128)
            phi(143) = phi(143) + scale_factor2(4)*scale_factor1(4)*phi(129)
            phi(144) = phi(144) + scale_factor2(4)*scale_factor1(4)*phi(130)
            phi(145) = phi(145) + scale_factor2(4)*scale_factor1(4)*phi(131)
            phi(146) = phi(146) + scale_factor2(4)*scale_factor1(4)*phi(132)
            phi(147) = phi(147) + scale_factor2(4)*scale_factor1(4)*phi(133)
            phi(148) = phi(148) + scale_factor2(4)*scale_factor1(4)*phi(134)
            phi(156) = phi(156) + scale_factor1(3)*phi(149)
            phi(157) = phi(157) + scale_factor1(3)*phi(150)
            phi(158) = phi(158) + scale_factor1(3)*phi(151)
            phi(159) = phi(159) + scale_factor1(3)*phi(152)
            phi(160) = phi(160) + scale_factor1(3)*phi(153)
            phi(161) = phi(161) + scale_factor1(3)*phi(154)
            phi(162) = phi(162) + scale_factor1(3)*phi(155)
            phi(163) = phi(163) + scale_factor2(3)*scale_factor1(4)*phi(149)
            phi(164) = phi(164) + scale_factor2(3)*scale_factor1(4)*phi(150)
            phi(165) = phi(165) + scale_factor2(3)*scale_factor1(4)*phi(151)
            phi(166) = phi(166) + scale_factor2(3)*scale_factor1(4)*phi(152)
            phi(167) = phi(167) + scale_factor2(3)*scale_factor1(4)*phi(153)
            phi(168) = phi(168) + scale_factor2(3)*scale_factor1(4)*phi(154)
            phi(169) = phi(169) + scale_factor2(3)*scale_factor1(4)*phi(155)
            phi(177) = phi(177) + scale_factor2(2)*scale_factor1(4)*phi(170)
            phi(178) = phi(178) + scale_factor2(2)*scale_factor1(4)*phi(171)
            phi(179) = phi(179) + scale_factor2(2)*scale_factor1(4)*phi(172)
            phi(180) = phi(180) + scale_factor2(2)*scale_factor1(4)*phi(173)
            phi(181) = phi(181) + scale_factor2(2)*scale_factor1(4)*phi(174)
            phi(182) = phi(182) + scale_factor2(2)*scale_factor1(4)*phi(175)
            phi(183) = phi(183) + scale_factor2(2)*scale_factor1(4)*phi(176)

            phi(197) = phi(197) + scale_factor1(4)*phi(184)
            phi(198) = phi(198) + scale_factor1(4)*phi(185)
            phi(199) = phi(199) + scale_factor1(4)*phi(186)
            phi(200) = phi(200) + scale_factor1(4)*phi(187)
            phi(201) = phi(201) + scale_factor1(4)*phi(188)
            phi(202) = phi(202) + scale_factor1(4)*phi(189)
            phi(203) = phi(203) + scale_factor1(4)*phi(190)
            phi(204) = phi(204) + scale_factor1(4)*phi(191)
            phi(205) = phi(205) + scale_factor1(4)*phi(192)
            phi(206) = phi(206) + scale_factor1(4)*phi(193)
            phi(207) = phi(207) + scale_factor1(4)*phi(194)
            phi(208) = phi(208) + scale_factor1(4)*phi(195)
            phi(209) = phi(209) + scale_factor1(4)*phi(196)
            phi(223) = phi(223) + scale_factor1(3)*phi(210)
            phi(224) = phi(224) + scale_factor1(3)*phi(211)
            phi(225) = phi(225) + scale_factor1(3)*phi(212)
            phi(226) = phi(226) + scale_factor1(3)*phi(213)
            phi(227) = phi(227) + scale_factor1(3)*phi(214)
            phi(228) = phi(228) + scale_factor1(3)*phi(215)
            phi(229) = phi(229) + scale_factor1(3)*phi(216)
            phi(230) = phi(230) + scale_factor1(3)*phi(217)
            phi(231) = phi(231) + scale_factor1(3)*phi(218)
            phi(232) = phi(232) + scale_factor1(3)*phi(219)
            phi(233) = phi(233) + scale_factor1(3)*phi(220)
            phi(234) = phi(234) + scale_factor1(3)*phi(221)
            phi(235) = phi(235) + scale_factor1(3)*phi(222)
            phi(236) = phi(236) + scale_factor2(3)*scale_factor1(4)*phi(210)
            phi(237) = phi(237) + scale_factor2(3)*scale_factor1(4)*phi(211)
            phi(238) = phi(238) + scale_factor2(3)*scale_factor1(4)*phi(212)
            phi(239) = phi(239) + scale_factor2(3)*scale_factor1(4)*phi(213)
            phi(240) = phi(240) + scale_factor2(3)*scale_factor1(4)*phi(214)
            phi(241) = phi(241) + scale_factor2(3)*scale_factor1(4)*phi(215)
            phi(242) = phi(242) + scale_factor2(3)*scale_factor1(4)*phi(216)
            phi(243) = phi(243) + scale_factor2(3)*scale_factor1(4)*phi(217)
            phi(244) = phi(244) + scale_factor2(3)*scale_factor1(4)*phi(218)
            phi(245) = phi(245) + scale_factor2(3)*scale_factor1(4)*phi(219)
            phi(246) = phi(246) + scale_factor2(3)*scale_factor1(4)*phi(220)
            phi(247) = phi(247) + scale_factor2(3)*scale_factor1(4)*phi(221)
            phi(248) = phi(248) + scale_factor2(3)*scale_factor1(4)*phi(222)
            phi(262) = phi(262) + scale_factor2(2)*scale_factor1(4)*phi(249)
            phi(263) = phi(263) + scale_factor2(2)*scale_factor1(4)*phi(250)
            phi(264) = phi(264) + scale_factor2(2)*scale_factor1(4)*phi(251)
            phi(265) = phi(265) + scale_factor2(2)*scale_factor1(4)*phi(252)
            phi(266) = phi(266) + scale_factor2(2)*scale_factor1(4)*phi(253)
            phi(267) = phi(267) + scale_factor2(2)*scale_factor1(4)*phi(254)
            phi(268) = phi(268) + scale_factor2(2)*scale_factor1(4)*phi(255)
            phi(269) = phi(269) + scale_factor2(2)*scale_factor1(4)*phi(256)
            phi(270) = phi(270) + scale_factor2(2)*scale_factor1(4)*phi(257)
            phi(271) = phi(271) + scale_factor2(2)*scale_factor1(4)*phi(258)
            phi(272) = phi(272) + scale_factor2(2)*scale_factor1(4)*phi(259)
            phi(273) = phi(273) + scale_factor2(2)*scale_factor1(4)*phi(260)
            phi(274) = phi(274) + scale_factor2(2)*scale_factor1(4)*phi(261)
            phi(288) = phi(288) + scale_factor1(4)*phi(275)
            phi(289) = phi(289) + scale_factor1(4)*phi(276)
            phi(290) = phi(290) + scale_factor1(4)*phi(277)
            phi(291) = phi(291) + scale_factor1(4)*phi(278)
            phi(292) = phi(292) + scale_factor1(4)*phi(279)
            phi(293) = phi(293) + scale_factor1(4)*phi(280)
            phi(294) = phi(294) + scale_factor1(4)*phi(281)
            phi(295) = phi(295) + scale_factor1(4)*phi(282)
            phi(296) = phi(296) + scale_factor1(4)*phi(283)
            phi(297) = phi(297) + scale_factor1(4)*phi(284)
            phi(298) = phi(298) + scale_factor1(4)*phi(285)
            phi(299) = phi(299) + scale_factor1(4)*phi(286)
            phi(300) = phi(300) + scale_factor1(4)*phi(287)

            phi(323) = phi(323) + scale_factor2(2)*scale_factor1(4)*phi(301)
            phi(324) = phi(324) + scale_factor2(2)*scale_factor1(4)*phi(302)
            phi(325) = phi(325) + scale_factor2(2)*scale_factor1(4)*phi(303)
            phi(326) = phi(326) + scale_factor2(2)*scale_factor1(4)*phi(304)
            phi(327) = phi(327) + scale_factor2(2)*scale_factor1(4)*phi(305)
            phi(328) = phi(328) + scale_factor2(2)*scale_factor1(4)*phi(306)
            phi(329) = phi(329) + scale_factor2(2)*scale_factor1(4)*phi(307)
            phi(330) = phi(330) + scale_factor2(2)*scale_factor1(4)*phi(308)
            phi(331) = phi(331) + scale_factor2(2)*scale_factor1(4)*phi(309)
            phi(332) = phi(332) + scale_factor2(2)*scale_factor1(4)*phi(310)
            phi(333) = phi(333) + scale_factor2(2)*scale_factor1(4)*phi(311)
            phi(334) = phi(334) + scale_factor2(2)*scale_factor1(4)*phi(312)
            phi(335) = phi(335) + scale_factor2(2)*scale_factor1(4)*phi(313)
            phi(336) = phi(336) + scale_factor2(2)*scale_factor1(4)*phi(314)
            phi(337) = phi(337) + scale_factor2(2)*scale_factor1(4)*phi(315)
            phi(338) = phi(338) + scale_factor2(2)*scale_factor1(4)*phi(316)
            phi(339) = phi(339) + scale_factor2(2)*scale_factor1(4)*phi(317)
            phi(340) = phi(340) + scale_factor2(2)*scale_factor1(4)*phi(318)
            phi(341) = phi(341) + scale_factor2(2)*scale_factor1(4)*phi(319)
            phi(342) = phi(342) + scale_factor2(2)*scale_factor1(4)*phi(320)
            phi(343) = phi(343) + scale_factor2(2)*scale_factor1(4)*phi(321)
            phi(344) = phi(344) + scale_factor2(2)*scale_factor1(4)*phi(322)
            phi(367) = phi(367) + scale_factor1(4)*phi(345)
            phi(368) = phi(368) + scale_factor1(4)*phi(346)
            phi(369) = phi(369) + scale_factor1(4)*phi(347)
            phi(370) = phi(370) + scale_factor1(4)*phi(348)
            phi(371) = phi(371) + scale_factor1(4)*phi(349)
            phi(372) = phi(372) + scale_factor1(4)*phi(350)
            phi(373) = phi(373) + scale_factor1(4)*phi(351)
            phi(374) = phi(374) + scale_factor1(4)*phi(352)
            phi(375) = phi(375) + scale_factor1(4)*phi(353)
            phi(376) = phi(376) + scale_factor1(4)*phi(354)
            phi(377) = phi(377) + scale_factor1(4)*phi(355)
            phi(378) = phi(378) + scale_factor1(4)*phi(356)
            phi(379) = phi(379) + scale_factor1(4)*phi(357)
            phi(380) = phi(380) + scale_factor1(4)*phi(358)
            phi(381) = phi(381) + scale_factor1(4)*phi(359)
            phi(382) = phi(382) + scale_factor1(4)*phi(360)
            phi(383) = phi(383) + scale_factor1(4)*phi(361)
            phi(384) = phi(384) + scale_factor1(4)*phi(362)
            phi(385) = phi(385) + scale_factor1(4)*phi(363)
            phi(386) = phi(386) + scale_factor1(4)*phi(364)
            phi(387) = phi(387) + scale_factor1(4)*phi(365)
            phi(388) = phi(388) + scale_factor1(4)*phi(366)

            phi(423) = phi(423) + scale_factor1(4)*phi(389)
            phi(424) = phi(424) + scale_factor1(4)*phi(390)
            phi(425) = phi(425) + scale_factor1(4)*phi(391)
            phi(426) = phi(426) + scale_factor1(4)*phi(392)
            phi(427) = phi(427) + scale_factor1(4)*phi(393)
            phi(428) = phi(428) + scale_factor1(4)*phi(394)
            phi(429) = phi(429) + scale_factor1(4)*phi(395)
            phi(430) = phi(430) + scale_factor1(4)*phi(396)
            phi(431) = phi(431) + scale_factor1(4)*phi(397)
            phi(432) = phi(432) + scale_factor1(4)*phi(398)
            phi(433) = phi(433) + scale_factor1(4)*phi(399)
            phi(434) = phi(434) + scale_factor1(4)*phi(400)
            phi(435) = phi(435) + scale_factor1(4)*phi(401)
            phi(436) = phi(436) + scale_factor1(4)*phi(402)
            phi(437) = phi(437) + scale_factor1(4)*phi(403)
            phi(438) = phi(438) + scale_factor1(4)*phi(404)
            phi(439) = phi(439) + scale_factor1(4)*phi(405)
            phi(440) = phi(440) + scale_factor1(4)*phi(406)
            phi(441) = phi(441) + scale_factor1(4)*phi(407)
            phi(442) = phi(442) + scale_factor1(4)*phi(408)
            phi(443) = phi(443) + scale_factor1(4)*phi(409)
            phi(444) = phi(444) + scale_factor1(4)*phi(410)
            phi(445) = phi(445) + scale_factor1(4)*phi(411)
            phi(446) = phi(446) + scale_factor1(4)*phi(412)
            phi(447) = phi(447) + scale_factor1(4)*phi(413)
            phi(448) = phi(448) + scale_factor1(4)*phi(414)
            phi(449) = phi(449) + scale_factor1(4)*phi(415)
            phi(450) = phi(450) + scale_factor1(4)*phi(416)
            phi(451) = phi(451) + scale_factor1(4)*phi(417)
            phi(452) = phi(452) + scale_factor1(4)*phi(418)
            phi(453) = phi(453) + scale_factor1(4)*phi(419)
            phi(454) = phi(454) + scale_factor1(4)*phi(420)
            phi(455) = phi(455) + scale_factor1(4)*phi(421)
            phi(456) = phi(456) + scale_factor1(4)*phi(422)

          end do ! k

!                   ******************************
!                   *                            *
!                   *   Post-contraction phase   *
!                   *                            *
!                   ******************************

          phi(122) = phi(122) - phi(121)
          phi(124) = phi(124) - phi(121)
          phi(127) = phi(127) - phi(121)

          phi(136) = phi(136) - phi(135)
          phi(138) = phi(138) - phi(135)
          phi(141) = phi(141) - phi(135)

          phi(143) = phi(143) - phi(142)
          phi(145) = phi(145) - phi(142)
          phi(148) = phi(148) - phi(142)

          phi(157) = phi(157) - phi(156)
          phi(159) = phi(159) - phi(156)
          phi(162) = phi(162) - phi(156)

          phi(164) = phi(164) - phi(163)
          phi(166) = phi(166) - phi(163)
          phi(169) = phi(169) - phi(163)

          phi(178) = phi(178) - phi(177)
          phi(180) = phi(180) - phi(177)
          phi(183) = phi(183) - phi(177)

          phi(200) = phi(200) - phi(197)
          phi(201) = phi(201) - phi(198)
          phi(203) = phi(203) - phi(198)
          phi(204) = phi(204) - phi(199)
          phi(206) = phi(206) - phi(199)
          phi(209) = phi(209) - phi(199)
          phi(200) = phi(200) - (2.0D+00)*phi(197)
          phi(202) = phi(202) - phi(197)
          phi(203) = phi(203) - (2.0D+00)*phi(198)
          phi(207) = phi(207) - phi(197)
          phi(208) = phi(208) - phi(198)
          phi(209) = phi(209) - (2.0D+00)*phi(199)

          phi(226) = phi(226) - phi(223)
          phi(227) = phi(227) - phi(224)
          phi(229) = phi(229) - phi(224)
          phi(230) = phi(230) - phi(225)
          phi(232) = phi(232) - phi(225)
          phi(235) = phi(235) - phi(225)
          phi(226) = phi(226) - (2.0D+00)*phi(223)
          phi(228) = phi(228) - phi(223)
          phi(229) = phi(229) - (2.0D+00)*phi(224)
          phi(233) = phi(233) - phi(223)
          phi(234) = phi(234) - phi(224)
          phi(235) = phi(235) - (2.0D+00)*phi(225)

          phi(239) = phi(239) - phi(236)
          phi(240) = phi(240) - phi(237)
          phi(242) = phi(242) - phi(237)
          phi(243) = phi(243) - phi(238)
          phi(245) = phi(245) - phi(238)
          phi(248) = phi(248) - phi(238)
          phi(239) = phi(239) - (2.0D+00)*phi(236)
          phi(241) = phi(241) - phi(236)
          phi(242) = phi(242) - (2.0D+00)*phi(237)
          phi(246) = phi(246) - phi(236)
          phi(247) = phi(247) - phi(237)
          phi(248) = phi(248) - (2.0D+00)*phi(238)

          phi(265) = phi(265) - phi(262)
          phi(266) = phi(266) - phi(263)
          phi(268) = phi(268) - phi(263)
          phi(269) = phi(269) - phi(264)
          phi(271) = phi(271) - phi(264)
          phi(274) = phi(274) - phi(264)
          phi(265) = phi(265) - (2.0D+00)*phi(262)
          phi(267) = phi(267) - phi(262)
          phi(268) = phi(268) - (2.0D+00)*phi(263)
          phi(272) = phi(272) - phi(262)
          phi(273) = phi(273) - phi(263)
          phi(274) = phi(274) - (2.0D+00)*phi(264)

          phi(291) = phi(291) - phi(288)
          phi(292) = phi(292) - phi(289)
          phi(294) = phi(294) - phi(289)
          phi(295) = phi(295) - phi(290)
          phi(297) = phi(297) - phi(290)
          phi(300) = phi(300) - phi(290)
          phi(291) = phi(291) - (2.0D+00)*phi(288)
          phi(293) = phi(293) - phi(288)
          phi(294) = phi(294) - (2.0D+00)*phi(289)
          phi(298) = phi(298) - phi(288)
          phi(299) = phi(299) - phi(289)
          phi(300) = phi(300) - (2.0D+00)*phi(290)

          phi(330) = phi(330) - phi(324)
          phi(331) = phi(331) - phi(325)
          phi(332) = phi(332) - phi(326)
          phi(334) = phi(334) - phi(326)
          phi(335) = phi(335) - phi(327)
          phi(336) = phi(336) - phi(328)
          phi(338) = phi(338) - phi(328)
          phi(339) = phi(339) - phi(329)
          phi(341) = phi(341) - phi(329)
          phi(344) = phi(344) - phi(329)
          phi(330) = phi(330) - (2.0D+00)*phi(324)
          phi(331) = phi(331) - (2.0D+00)*phi(325)
          phi(333) = phi(333) - phi(325)
          phi(334) = phi(334) - (2.0D+00)*phi(326)
          phi(335) = phi(335) - (2.0D+00)*phi(327)
          phi(337) = phi(337) - phi(327)
          phi(338) = phi(338) - (2.0D+00)*phi(328)
          phi(342) = phi(342) - phi(327)
          phi(343) = phi(343) - phi(328)
          phi(344) = phi(344) - (2.0D+00)*phi(329)
          phi(324) = phi(324) - phi(323)
          phi(326) = phi(326) - phi(323)
          phi(329) = phi(329) - phi(323)
          phi(330) = phi(330) - (3.0D+00)*phi(324)
          phi(332) = phi(332) - phi(324)
          phi(333) = phi(333) - (2.0D+00)*phi(325)
          phi(334) = phi(334) - (3.0D+00)*phi(326)
          phi(339) = phi(339) - phi(324)
          phi(340) = phi(340) - phi(325)
          phi(341) = phi(341) - phi(326)
          phi(342) = phi(342) - (2.0D+00)*phi(327)
          phi(343) = phi(343) - (2.0D+00)*phi(328)
          phi(344) = phi(344) - (3.0D+00)*phi(329)

          phi(374) = phi(374) - phi(368)
          phi(375) = phi(375) - phi(369)
          phi(376) = phi(376) - phi(370)
          phi(378) = phi(378) - phi(370)
          phi(379) = phi(379) - phi(371)
          phi(380) = phi(380) - phi(372)
          phi(382) = phi(382) - phi(372)
          phi(383) = phi(383) - phi(373)
          phi(385) = phi(385) - phi(373)
          phi(388) = phi(388) - phi(373)
          phi(374) = phi(374) - (2.0D+00)*phi(368)
          phi(375) = phi(375) - (2.0D+00)*phi(369)
          phi(377) = phi(377) - phi(369)
          phi(378) = phi(378) - (2.0D+00)*phi(370)
          phi(379) = phi(379) - (2.0D+00)*phi(371)
          phi(381) = phi(381) - phi(371)
          phi(382) = phi(382) - (2.0D+00)*phi(372)
          phi(386) = phi(386) - phi(371)
          phi(387) = phi(387) - phi(372)
          phi(388) = phi(388) - (2.0D+00)*phi(373)
          phi(368) = phi(368) - phi(367)
          phi(370) = phi(370) - phi(367)
          phi(373) = phi(373) - phi(367)
          phi(374) = phi(374) - (3.0D+00)*phi(368)
          phi(376) = phi(376) - phi(368)
          phi(377) = phi(377) - (2.0D+00)*phi(369)
          phi(378) = phi(378) - (3.0D+00)*phi(370)
          phi(383) = phi(383) - phi(368)
          phi(384) = phi(384) - phi(369)
          phi(385) = phi(385) - phi(370)
          phi(386) = phi(386) - (2.0D+00)*phi(371)
          phi(387) = phi(387) - (2.0D+00)*phi(372)
          phi(388) = phi(388) - (3.0D+00)*phi(373)

          phi(436) = phi(436) - phi(426)
          phi(437) = phi(437) - phi(427)
          phi(438) = phi(438) - phi(428)
          phi(439) = phi(439) - phi(429)
          phi(441) = phi(441) - phi(429)
          phi(442) = phi(442) - phi(430)
          phi(443) = phi(443) - phi(431)
          phi(444) = phi(444) - phi(432)
          phi(446) = phi(446) - phi(432)
          phi(447) = phi(447) - phi(433)
          phi(448) = phi(448) - phi(434)
          phi(450) = phi(450) - phi(434)
          phi(451) = phi(451) - phi(435)
          phi(453) = phi(453) - phi(435)
          phi(456) = phi(456) - phi(435)
          phi(436) = phi(436) - (2.0D+00)*phi(426)
          phi(437) = phi(437) - (2.0D+00)*phi(427)
          phi(438) = phi(438) - (2.0D+00)*phi(428)
          phi(440) = phi(440) - phi(428)
          phi(441) = phi(441) - (2.0D+00)*phi(429)
          phi(442) = phi(442) - (2.0D+00)*phi(430)
          phi(443) = phi(443) - (2.0D+00)*phi(431)
          phi(445) = phi(445) - phi(431)
          phi(446) = phi(446) - (2.0D+00)*phi(432)
          phi(447) = phi(447) - (2.0D+00)*phi(433)
          phi(449) = phi(449) - phi(433)
          phi(450) = phi(450) - (2.0D+00)*phi(434)
          phi(454) = phi(454) - phi(433)
          phi(455) = phi(455) - phi(434)
          phi(456) = phi(456) - (2.0D+00)*phi(435)
          phi(426) = phi(426) - phi(423)
          phi(427) = phi(427) - phi(424)
          phi(429) = phi(429) - phi(424)
          phi(430) = phi(430) - phi(425)
          phi(432) = phi(432) - phi(425)
          phi(435) = phi(435) - phi(425)
          phi(436) = phi(436) - (3.0D+00)*phi(426)
          phi(437) = phi(437) - (3.0D+00)*phi(427)
          phi(439) = phi(439) - phi(427)
          phi(440) = phi(440) - (2.0D+00)*phi(428)
          phi(441) = phi(441) - (3.0D+00)*phi(429)
          phi(442) = phi(442) - (3.0D+00)*phi(430)
          phi(444) = phi(444) - phi(430)
          phi(445) = phi(445) - (2.0D+00)*phi(431)
          phi(446) = phi(446) - (3.0D+00)*phi(432)
          phi(451) = phi(451) - phi(430)
          phi(452) = phi(452) - phi(431)
          phi(453) = phi(453) - phi(432)
          phi(454) = phi(454) - (2.0D+00)*phi(433)
          phi(455) = phi(455) - (2.0D+00)*phi(434)
          phi(456) = phi(456) - (3.0D+00)*phi(435)
          phi(426) = phi(426) - (2.0D+00)*phi(423)
          phi(428) = phi(428) - phi(423)
          phi(429) = phi(429) - (2.0D+00)*phi(424)
          phi(433) = phi(433) - phi(423)
          phi(434) = phi(434) - phi(424)
          phi(435) = phi(435) - (2.0D+00)*phi(425)
          phi(436) = phi(436) - (4.0D+00)*phi(426)
          phi(438) = phi(438) - phi(426)
          phi(439) = phi(439) - (2.0D+00)*phi(427)
          phi(440) = phi(440) - (3.0D+00)*phi(428)
          phi(441) = phi(441) - (4.0D+00)*phi(429)
          phi(447) = phi(447) - phi(426)
          phi(448) = phi(448) - phi(427)
          phi(449) = phi(449) - phi(428)
          phi(450) = phi(450) - phi(429)
          phi(451) = phi(451) - (2.0D+00)*phi(430)
          phi(452) = phi(452) - (2.0D+00)*phi(431)
          phi(453) = phi(453) - (2.0D+00)*phi(432)
          phi(454) = phi(454) - (3.0D+00)*phi(433)
          phi(455) = phi(455) - (3.0D+00)*phi(434)
          phi(456) = phi(456) - (4.0D+00)*phi(435)

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

          work1(6) = phi(143) + cnf(1)*phi(102)
          work1(7) = phi(144) + cnf(1)*phi(103)
          work1(8) = phi(145) + cnf(2)*phi(103)
          work1(9) = phi(146) + cnf(1)*phi(104)
          work1(10) = phi(147) + cnf(2)*phi(104)
          work1(11) = phi(148) + cnf(3)*phi(104)
          work1(3) = phi(102) + cnf(1)*phi(86)
          work1(4) = phi(103) + cnf(2)*phi(86)
          work1(5) = phi(104) + cnf(3)*phi(86)
          work1(6) = work1(6) + cnf(4)*phi(83)
          work1(8) = work1(8) + cnf(4)*phi(83)
          work1(11) = work1(11) + cnf(4)*phi(83)

          iii = 1
          ! Generating BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! bcte done

          work1(6) = phi(136) + cnf(1)*phi(99)
          work1(7) = phi(137) + cnf(1)*phi(100)
          work1(8) = phi(138) + cnf(2)*phi(100)
          work1(9) = phi(139) + cnf(1)*phi(101)
          work1(10) = phi(140) + cnf(2)*phi(101)
          work1(11) = phi(141) + cnf(3)*phi(101)
          work1(3) = phi(99) + cnf(1)*phi(85)
          work1(4) = phi(100) + cnf(2)*phi(85)
          work1(5) = phi(101) + cnf(3)*phi(85)
          work1(6) = work1(6) + cnf(4)*phi(82)
          work1(8) = work1(8) + cnf(4)*phi(82)
          work1(11) = work1(11) + cnf(4)*phi(82)

          iii = 2
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(239) - cnf(1)*phi(164)
          work1(7) = -phi(240) - cnf(1)*phi(165)
          work1(8) = -phi(241) - cnf(2)*phi(165)
          work1(9) = -phi(243) - cnf(1)*phi(167)
          work1(10) = -phi(244) - cnf(2)*phi(167)
          work1(11) = -phi(246) - cnf(3)*phi(167)
          work1(3) = -phi(164) - cnf(1)*phi(111)
          work1(4) = -phi(165) - cnf(2)*phi(111)
          work1(5) = -phi(167) - cnf(3)*phi(111)
          work1(6) = work1(6) - cnf(4)*phi(93)
          work1(8) = work1(8) - cnf(4)*phi(93)
          work1(11) = work1(11) - cnf(4)*phi(93)

          iii = 3
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(240) - cnf(1)*phi(165)
          work1(7) = -phi(241) - cnf(1)*phi(166)
          work1(8) = -phi(242) - cnf(2)*phi(166)
          work1(9) = -phi(244) - cnf(1)*phi(168)
          work1(10) = -phi(245) - cnf(2)*phi(168)
          work1(11) = -phi(247) - cnf(3)*phi(168)
          work1(3) = -phi(165) - cnf(1)*phi(112)
          work1(4) = -phi(166) - cnf(2)*phi(112)
          work1(5) = -phi(168) - cnf(3)*phi(112)
          work1(6) = work1(6) - cnf(4)*phi(94)
          work1(8) = work1(8) - cnf(4)*phi(94)
          work1(11) = work1(11) - cnf(4)*phi(94)

          iii = 4
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(243) - cnf(1)*phi(167)
          work1(7) = -phi(244) - cnf(1)*phi(168)
          work1(8) = -phi(245) - cnf(2)*phi(168)
          work1(9) = -phi(246) - cnf(1)*phi(169)
          work1(10) = -phi(247) - cnf(2)*phi(169)
          work1(11) = -phi(248) - cnf(3)*phi(169)
          work1(3) = -phi(167) - cnf(1)*phi(113)
          work1(4) = -phi(168) - cnf(2)*phi(113)
          work1(5) = -phi(169) - cnf(3)*phi(113)
          work1(6) = work1(6) - cnf(4)*phi(95)
          work1(8) = work1(8) - cnf(4)*phi(95)
          work1(11) = work1(11) - cnf(4)*phi(95)

          iii = 5
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(226) - cnf(1)*phi(157)
          work1(7) = -phi(227) - cnf(1)*phi(158)
          work1(8) = -phi(228) - cnf(2)*phi(158)
          work1(9) = -phi(230) - cnf(1)*phi(160)
          work1(10) = -phi(231) - cnf(2)*phi(160)
          work1(11) = -phi(233) - cnf(3)*phi(160)
          work1(3) = -phi(157) - cnf(1)*phi(108)
          work1(4) = -phi(158) - cnf(2)*phi(108)
          work1(5) = -phi(160) - cnf(3)*phi(108)
          work1(6) = work1(6) - cnf(4)*phi(90)
          work1(8) = work1(8) - cnf(4)*phi(90)
          work1(11) = work1(11) - cnf(4)*phi(90)

          iii = 6
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(227) - cnf(1)*phi(158)
          work1(7) = -phi(228) - cnf(1)*phi(159)
          work1(8) = -phi(229) - cnf(2)*phi(159)
          work1(9) = -phi(231) - cnf(1)*phi(161)
          work1(10) = -phi(232) - cnf(2)*phi(161)
          work1(11) = -phi(234) - cnf(3)*phi(161)
          work1(3) = -phi(158) - cnf(1)*phi(109)
          work1(4) = -phi(159) - cnf(2)*phi(109)
          work1(5) = -phi(161) - cnf(3)*phi(109)
          work1(6) = work1(6) - cnf(4)*phi(91)
          work1(8) = work1(8) - cnf(4)*phi(91)
          work1(11) = work1(11) - cnf(4)*phi(91)

          iii = 7
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(230) - cnf(1)*phi(160)
          work1(7) = -phi(231) - cnf(1)*phi(161)
          work1(8) = -phi(232) - cnf(2)*phi(161)
          work1(9) = -phi(233) - cnf(1)*phi(162)
          work1(10) = -phi(234) - cnf(2)*phi(162)
          work1(11) = -phi(235) - cnf(3)*phi(162)
          work1(3) = -phi(160) - cnf(1)*phi(110)
          work1(4) = -phi(161) - cnf(2)*phi(110)
          work1(5) = -phi(162) - cnf(3)*phi(110)
          work1(6) = work1(6) - cnf(4)*phi(92)
          work1(8) = work1(8) - cnf(4)*phi(92)
          work1(11) = work1(11) - cnf(4)*phi(92)

          iii = 8
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(330) + cnf(1)*phi(265)
          work1(7) = phi(331) + cnf(1)*phi(266)
          work1(8) = phi(332) + cnf(2)*phi(266)
          work1(9) = phi(335) + cnf(1)*phi(269)
          work1(10) = phi(336) + cnf(2)*phi(269)
          work1(11) = phi(339) + cnf(3)*phi(269)
          work1(3) = phi(265) + cnf(1)*phi(178)
          work1(4) = phi(266) + cnf(2)*phi(178)
          work1(5) = phi(269) + cnf(3)*phi(178)
          work1(6) = work1(6) + cnf(4)*phi(122)
          work1(8) = work1(8) + cnf(4)*phi(122)
          work1(11) = work1(11) + cnf(4)*phi(122)

          iii = 9
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(331) + cnf(1)*phi(266)
          work1(7) = phi(332) + cnf(1)*phi(267)
          work1(8) = phi(333) + cnf(2)*phi(267)
          work1(9) = phi(336) + cnf(1)*phi(270)
          work1(10) = phi(337) + cnf(2)*phi(270)
          work1(11) = phi(340) + cnf(3)*phi(270)
          work1(3) = phi(266) + cnf(1)*phi(179)
          work1(4) = phi(267) + cnf(2)*phi(179)
          work1(5) = phi(270) + cnf(3)*phi(179)
          work1(6) = work1(6) + cnf(4)*phi(123)
          work1(8) = work1(8) + cnf(4)*phi(123)
          work1(11) = work1(11) + cnf(4)*phi(123)

          iii = 10
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(332) + cnf(1)*phi(267)
          work1(7) = phi(333) + cnf(1)*phi(268)
          work1(8) = phi(334) + cnf(2)*phi(268)
          work1(9) = phi(337) + cnf(1)*phi(271)
          work1(10) = phi(338) + cnf(2)*phi(271)
          work1(11) = phi(341) + cnf(3)*phi(271)
          work1(3) = phi(267) + cnf(1)*phi(180)
          work1(4) = phi(268) + cnf(2)*phi(180)
          work1(5) = phi(271) + cnf(3)*phi(180)
          work1(6) = work1(6) + cnf(4)*phi(124)
          work1(8) = work1(8) + cnf(4)*phi(124)
          work1(11) = work1(11) + cnf(4)*phi(124)

          iii = 11
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(335) + cnf(1)*phi(269)
          work1(7) = phi(336) + cnf(1)*phi(270)
          work1(8) = phi(337) + cnf(2)*phi(270)
          work1(9) = phi(339) + cnf(1)*phi(272)
          work1(10) = phi(340) + cnf(2)*phi(272)
          work1(11) = phi(342) + cnf(3)*phi(272)
          work1(3) = phi(269) + cnf(1)*phi(181)
          work1(4) = phi(270) + cnf(2)*phi(181)
          work1(5) = phi(272) + cnf(3)*phi(181)
          work1(6) = work1(6) + cnf(4)*phi(125)
          work1(8) = work1(8) + cnf(4)*phi(125)
          work1(11) = work1(11) + cnf(4)*phi(125)

          iii = 12
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(336) + cnf(1)*phi(270)
          work1(7) = phi(337) + cnf(1)*phi(271)
          work1(8) = phi(338) + cnf(2)*phi(271)
          work1(9) = phi(340) + cnf(1)*phi(273)
          work1(10) = phi(341) + cnf(2)*phi(273)
          work1(11) = phi(343) + cnf(3)*phi(273)
          work1(3) = phi(270) + cnf(1)*phi(182)
          work1(4) = phi(271) + cnf(2)*phi(182)
          work1(5) = phi(273) + cnf(3)*phi(182)
          work1(6) = work1(6) + cnf(4)*phi(126)
          work1(8) = work1(8) + cnf(4)*phi(126)
          work1(11) = work1(11) + cnf(4)*phi(126)

          iii = 13
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = phi(339) + cnf(1)*phi(272)
          work1(7) = phi(340) + cnf(1)*phi(273)
          work1(8) = phi(341) + cnf(2)*phi(273)
          work1(9) = phi(342) + cnf(1)*phi(274)
          work1(10) = phi(343) + cnf(2)*phi(274)
          work1(11) = phi(344) + cnf(3)*phi(274)
          work1(3) = phi(272) + cnf(1)*phi(183)
          work1(4) = phi(273) + cnf(2)*phi(183)
          work1(5) = phi(274) + cnf(3)*phi(183)
          work1(6) = work1(6) + cnf(4)*phi(127)
          work1(8) = work1(8) + cnf(4)*phi(127)
          work1(11) = work1(11) + cnf(4)*phi(127)

          iii = 14
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(436) - cnf(1)*phi(374)
          work1(7) = -phi(437) - cnf(1)*phi(375)
          work1(8) = -phi(438) - cnf(2)*phi(375)
          work1(9) = -phi(442) - cnf(1)*phi(379)
          work1(10) = -phi(443) - cnf(2)*phi(379)
          work1(11) = -phi(447) - cnf(3)*phi(379)
          work1(3) = -phi(374) - cnf(1)*phi(291)
          work1(4) = -phi(375) - cnf(2)*phi(291)
          work1(5) = -phi(379) - cnf(3)*phi(291)
          work1(6) = work1(6) - cnf(4)*phi(200)
          work1(8) = work1(8) - cnf(4)*phi(200)
          work1(11) = work1(11) - cnf(4)*phi(200)

          iii = 15
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(437) - cnf(1)*phi(375)
          work1(7) = -phi(438) - cnf(1)*phi(376)
          work1(8) = -phi(439) - cnf(2)*phi(376)
          work1(9) = -phi(443) - cnf(1)*phi(380)
          work1(10) = -phi(444) - cnf(2)*phi(380)
          work1(11) = -phi(448) - cnf(3)*phi(380)
          work1(3) = -phi(375) - cnf(1)*phi(292)
          work1(4) = -phi(376) - cnf(2)*phi(292)
          work1(5) = -phi(380) - cnf(3)*phi(292)
          work1(6) = work1(6) - cnf(4)*phi(201)
          work1(8) = work1(8) - cnf(4)*phi(201)
          work1(11) = work1(11) - cnf(4)*phi(201)

          iii = 16
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(438) - cnf(1)*phi(376)
          work1(7) = -phi(439) - cnf(1)*phi(377)
          work1(8) = -phi(440) - cnf(2)*phi(377)
          work1(9) = -phi(444) - cnf(1)*phi(381)
          work1(10) = -phi(445) - cnf(2)*phi(381)
          work1(11) = -phi(449) - cnf(3)*phi(381)
          work1(3) = -phi(376) - cnf(1)*phi(293)
          work1(4) = -phi(377) - cnf(2)*phi(293)
          work1(5) = -phi(381) - cnf(3)*phi(293)
          work1(6) = work1(6) - cnf(4)*phi(202)
          work1(8) = work1(8) - cnf(4)*phi(202)
          work1(11) = work1(11) - cnf(4)*phi(202)

          iii = 17
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(439) - cnf(1)*phi(377)
          work1(7) = -phi(440) - cnf(1)*phi(378)
          work1(8) = -phi(441) - cnf(2)*phi(378)
          work1(9) = -phi(445) - cnf(1)*phi(382)
          work1(10) = -phi(446) - cnf(2)*phi(382)
          work1(11) = -phi(450) - cnf(3)*phi(382)
          work1(3) = -phi(377) - cnf(1)*phi(294)
          work1(4) = -phi(378) - cnf(2)*phi(294)
          work1(5) = -phi(382) - cnf(3)*phi(294)
          work1(6) = work1(6) - cnf(4)*phi(203)
          work1(8) = work1(8) - cnf(4)*phi(203)
          work1(11) = work1(11) - cnf(4)*phi(203)

          iii = 18
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(442) - cnf(1)*phi(379)
          work1(7) = -phi(443) - cnf(1)*phi(380)
          work1(8) = -phi(444) - cnf(2)*phi(380)
          work1(9) = -phi(447) - cnf(1)*phi(383)
          work1(10) = -phi(448) - cnf(2)*phi(383)
          work1(11) = -phi(451) - cnf(3)*phi(383)
          work1(3) = -phi(379) - cnf(1)*phi(295)
          work1(4) = -phi(380) - cnf(2)*phi(295)
          work1(5) = -phi(383) - cnf(3)*phi(295)
          work1(6) = work1(6) - cnf(4)*phi(204)
          work1(8) = work1(8) - cnf(4)*phi(204)
          work1(11) = work1(11) - cnf(4)*phi(204)

          iii = 19
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(443) - cnf(1)*phi(380)
          work1(7) = -phi(444) - cnf(1)*phi(381)
          work1(8) = -phi(445) - cnf(2)*phi(381)
          work1(9) = -phi(448) - cnf(1)*phi(384)
          work1(10) = -phi(449) - cnf(2)*phi(384)
          work1(11) = -phi(452) - cnf(3)*phi(384)
          work1(3) = -phi(380) - cnf(1)*phi(296)
          work1(4) = -phi(381) - cnf(2)*phi(296)
          work1(5) = -phi(384) - cnf(3)*phi(296)
          work1(6) = work1(6) - cnf(4)*phi(205)
          work1(8) = work1(8) - cnf(4)*phi(205)
          work1(11) = work1(11) - cnf(4)*phi(205)

          iii = 20
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(444) - cnf(1)*phi(381)
          work1(7) = -phi(445) - cnf(1)*phi(382)
          work1(8) = -phi(446) - cnf(2)*phi(382)
          work1(9) = -phi(449) - cnf(1)*phi(385)
          work1(10) = -phi(450) - cnf(2)*phi(385)
          work1(11) = -phi(453) - cnf(3)*phi(385)
          work1(3) = -phi(381) - cnf(1)*phi(297)
          work1(4) = -phi(382) - cnf(2)*phi(297)
          work1(5) = -phi(385) - cnf(3)*phi(297)
          work1(6) = work1(6) - cnf(4)*phi(206)
          work1(8) = work1(8) - cnf(4)*phi(206)
          work1(11) = work1(11) - cnf(4)*phi(206)

          iii = 21
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(447) - cnf(1)*phi(383)
          work1(7) = -phi(448) - cnf(1)*phi(384)
          work1(8) = -phi(449) - cnf(2)*phi(384)
          work1(9) = -phi(451) - cnf(1)*phi(386)
          work1(10) = -phi(452) - cnf(2)*phi(386)
          work1(11) = -phi(454) - cnf(3)*phi(386)
          work1(3) = -phi(383) - cnf(1)*phi(298)
          work1(4) = -phi(384) - cnf(2)*phi(298)
          work1(5) = -phi(386) - cnf(3)*phi(298)
          work1(6) = work1(6) - cnf(4)*phi(207)
          work1(8) = work1(8) - cnf(4)*phi(207)
          work1(11) = work1(11) - cnf(4)*phi(207)

          iii = 22
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(448) - cnf(1)*phi(384)
          work1(7) = -phi(449) - cnf(1)*phi(385)
          work1(8) = -phi(450) - cnf(2)*phi(385)
          work1(9) = -phi(452) - cnf(1)*phi(387)
          work1(10) = -phi(453) - cnf(2)*phi(387)
          work1(11) = -phi(455) - cnf(3)*phi(387)
          work1(3) = -phi(384) - cnf(1)*phi(299)
          work1(4) = -phi(385) - cnf(2)*phi(299)
          work1(5) = -phi(387) - cnf(3)*phi(299)
          work1(6) = work1(6) - cnf(4)*phi(208)
          work1(8) = work1(8) - cnf(4)*phi(208)
          work1(11) = work1(11) - cnf(4)*phi(208)

          iii = 23
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
            j = j + 1
            work2(iii, i) = work1(j)
          end do
          ! BCTE END

          work1(6) = -phi(451) - cnf(1)*phi(386)
          work1(7) = -phi(452) - cnf(1)*phi(387)
          work1(8) = -phi(453) - cnf(2)*phi(387)
          work1(9) = -phi(454) - cnf(1)*phi(388)
          work1(10) = -phi(455) - cnf(2)*phi(388)
          work1(11) = -phi(456) - cnf(3)*phi(388)
          work1(3) = -phi(386) - cnf(1)*phi(300)
          work1(4) = -phi(387) - cnf(2)*phi(300)
          work1(5) = -phi(388) - cnf(3)*phi(300)
          work1(6) = work1(6) - cnf(4)*phi(209)
          work1(8) = work1(8) - cnf(4)*phi(209)
          work1(11) = work1(11) - cnf(4)*phi(209)

          iii = 24
          ! BCTE
          work1(6) = work1(6) + cnf(1)*work1(3)
          work1(7) = work1(7) + cnf(2)*work1(3)
          work1(8) = work1(8) + cnf(2)*work1(4)
          work1(9) = work1(9) + cnf(3)*work1(3)
          work1(10) = work1(10) + cnf(3)*work1(4)
          work1(11) = work1(11) + cnf(3)*work1(5)
          j = 5
          do i = 1, 6
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
          eri_value(11) = work2(15, 3)
          eri_value(12) = work2(18, 3)
          eri_value(13) = work2(24, 3)
          eri_value(14) = work2(16, 3)*sqrt5
          eri_value(15) = work2(19, 3)*sqrt5
          eri_value(16) = work2(17, 3)*sqrt5
          eri_value(17) = work2(21, 3)*sqrt5
          eri_value(18) = work2(22, 3)*sqrt5
          eri_value(19) = work2(23, 3)*sqrt5
          eri_value(20) = work2(20, 3)*sqrt15
          eri_value(21) = work2(15, 6)
          eri_value(22) = work2(18, 6)
          eri_value(23) = work2(24, 6)
          eri_value(24) = work2(16, 6)*sqrt5
          eri_value(25) = work2(19, 6)*sqrt5
          eri_value(26) = work2(17, 6)*sqrt5
          eri_value(27) = work2(21, 6)*sqrt5
          eri_value(28) = work2(22, 6)*sqrt5
          eri_value(29) = work2(23, 6)*sqrt5
          eri_value(30) = work2(20, 6)*sqrt15
          eri_value(31) = work2(15, 2)*sqrt3
          eri_value(32) = work2(18, 2)*sqrt3
          eri_value(33) = work2(24, 2)*sqrt3
          eri_value(34) = work2(16, 2)*sqrt3*sqrt5
          eri_value(35) = work2(19, 2)*sqrt3*sqrt5
          eri_value(36) = work2(17, 2)*sqrt3*sqrt5
          eri_value(37) = work2(21, 2)*sqrt3*sqrt5
          eri_value(38) = work2(22, 2)*sqrt3*sqrt5
          eri_value(39) = work2(23, 2)*sqrt3*sqrt5
          eri_value(40) = work2(20, 2)*sqrt3*sqrt15
          eri_value(41) = work2(15, 4)*sqrt3
          eri_value(42) = work2(18, 4)*sqrt3
          eri_value(43) = work2(24, 4)*sqrt3
          eri_value(44) = work2(16, 4)*sqrt3*sqrt5
          eri_value(45) = work2(19, 4)*sqrt3*sqrt5
          eri_value(46) = work2(17, 4)*sqrt3*sqrt5
          eri_value(47) = work2(21, 4)*sqrt3*sqrt5
          eri_value(48) = work2(22, 4)*sqrt3*sqrt5
          eri_value(49) = work2(23, 4)*sqrt3*sqrt5
          eri_value(50) = work2(20, 4)*sqrt3*sqrt15
          eri_value(51) = work2(15, 5)*sqrt3
          eri_value(52) = work2(18, 5)*sqrt3
          eri_value(53) = work2(24, 5)*sqrt3
          eri_value(54) = work2(16, 5)*sqrt3*sqrt5
          eri_value(55) = work2(19, 5)*sqrt3*sqrt5
          eri_value(56) = work2(17, 5)*sqrt3*sqrt5
          eri_value(57) = work2(21, 5)*sqrt3*sqrt5
          eri_value(58) = work2(22, 5)*sqrt3*sqrt5
          eri_value(59) = work2(23, 5)*sqrt3*sqrt5
          eri_value(60) = work2(20, 5)*sqrt3*sqrt15

          maxi = 6
          maxk = 10

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

          kstride = 1
          istride = 10

          ip = 1
          do i = 1, maxi
            ii1 = i + loci
            ijp = ip
            ip = ip + istride

            j = 1
            nij = nij + 1
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
            !
          end do

        end if ! test.gt.cutoff_schwarz
      end do ! iquart
!$omp end target teams distribute parallel do

    end do ! itile

  end subroutine int2030
end submodule
