! The total angular momentum of this class is:           3
! The algorithm chosen is: PHR a.k.a. ERIC
! Writing an ERIC kernel
submodule(eric_kernels) int0030_impl
contains
  module subroutine int0030(ss_pair, sf_pair, density, fock, res)

    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: ss_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

! Variables for the class
    integer(kind=int64), allocatable :: n00bra(:), n03ket(:)
    real(dp), allocatable :: xint00bra(:), xint03ket(:)
    integer(kind=int64) :: nssbra, nsfket
    real(dp) :: scutssbra, scutsfket, test
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
    real(dp) :: ft(4), phi(76), c_factor(20), rxyz(3)
    real(dp) :: work2(24, 1)
    real(dp) :: eri_value(10), angl(20)
    real(dp) :: ai, aij, aijk, aijkl
    integer(kind=int64) :: iord(20)
    integer(kind=int64) :: kstride
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
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

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

    if ((nssbra*nsfket) .le. nchunksize_k10) nchunksize_k10 = nssbra*nsfket
    ntile = int(nssbra*nsfket/nchunksize_k10)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_k10 + 1
      iend = itile*nchunksize_k10
      if (itile .eq. ntile) iend = nssbra*nsfket

! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

! Mappings to GPU

!$omp target teams distribute parallel do default(none) &
!$omp shared(res, density, fock, nquart_start, nquart_end, nssbra, xint00bra, n00bra, xint03ket, n03ket, nsfket, sf_pair, ss_pair) &
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
!$omp private(loci,locj,lock,locl,nij,kstride) &
!$omp private(ip,ii1,ijp,jj1,i2,j2,ijkp,nkl,kk1) &
!$omp private(ijklp,buff,ll1,k2,l2,ii2,jj2,kk2,ik,il,jk,jl,ii,jj,kk,ll)
      do iquart = nquart_start, nquart_end

        ij_tmp = (iquart - 1)/nsfket + 1
        kl_tmp = mod(iquart - 1, nsfket) + 1

        test = xint00bra(ij_tmp)*xint03ket(kl_tmp)

        if (test .gt. cutoff_schwarz) then

          ij = n00bra(ij_tmp)
          kl = n03ket(kl_tmp)

          ish_tmp = (1 + sqrt(1.0 + 8.0*(ij - 1)))/2
          jsh_tmp = ij - ish_tmp*(ish_tmp - 1)/2
          ksh_tmp = (kl - 1)/res%n_f_shl + 1
          lsh_tmp = mod(kl - 1, res%n_f_shl) + 1

          ish = res%i_s_shl(ish_tmp)
          jsh = res%i_s_shl(jsh_tmp)
          ksh = res%i_f_shl(lsh_tmp)
          lsh = res%i_s_shl(ksh_tmp)

          ijtop = res%contr_num(ish)*res%contr_num(jsh)

          kltop = res%contr_num(ksh)*res%contr_num(lsh)

          eri_value = 0.0_dp
          ket_loop = 0

          ! 0.0_dp out phi elements to contract over kl
          ! PHI = "pre-Hermite integrals"

          phi(26) = 0.0_dp
          phi(27) = 0.0_dp
          phi(31) = 0.0_dp
          phi(32) = 0.0_dp
          phi(33) = 0.0_dp
          phi(34) = 0.0_dp
          phi(35) = 0.0_dp
          phi(36) = 0.0_dp
          phi(44) = 0.0_dp
          phi(45) = 0.0_dp
          phi(46) = 0.0_dp
          phi(47) = 0.0_dp
          phi(48) = 0.0_dp
          phi(49) = 0.0_dp
          phi(50) = 0.0_dp
          phi(64) = 0.0_dp
          phi(65) = 0.0_dp
          phi(66) = 0.0_dp
          phi(67) = 0.0_dp
          phi(68) = 0.0_dp
          phi(69) = 0.0_dp
          phi(70) = 0.0_dp
          phi(71) = 0.0_dp
          phi(72) = 0.0_dp
          phi(73) = 0.0_dp
          phi(74) = 0.0_dp
          phi(75) = 0.0_dp
          phi(76) = 0.0_dp

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

            phi(25) = 0.0_dp
            phi(28) = 0.0_dp
            phi(29) = 0.0_dp
            phi(30) = 0.0_dp
            phi(37) = 0.0_dp
            phi(38) = 0.0_dp
            phi(39) = 0.0_dp
            phi(40) = 0.0_dp
            phi(41) = 0.0_dp
            phi(42) = 0.0_dp
            phi(43) = 0.0_dp
            phi(51) = 0.0_dp
            phi(52) = 0.0_dp
            phi(53) = 0.0_dp
            phi(54) = 0.0_dp
            phi(55) = 0.0_dp
            phi(56) = 0.0_dp
            phi(57) = 0.0_dp
            phi(58) = 0.0_dp
            phi(59) = 0.0_dp
            phi(60) = 0.0_dp
            phi(61) = 0.0_dp
            phi(62) = 0.0_dp
            phi(63) = 0.0_dp

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

              n = 3
              if (tt .le. t_max) then

                ! Boys function evaluation with Chebyshev interpolation

                t_new = tt*0.1571158409626437D+02
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

                ft(4) = fmt
                t2 = tt + tt
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

              ! Begin contraction & scaling over ij
              ! i shell scaling factors
              fac1 = 1.0_dp
              scale_factor1(1) = 1.0_dp
              ! j shell scaling factors
              fac2 = 1.0_dp
              scale_factor2(1) = 1.0_dp

              ! ij contraction
              phi(25) = phi(25) + phi(1)

              phi(28) = phi(28) + phi(2)
              phi(29) = phi(29) + phi(3)
              phi(30) = phi(30) + phi(4)

              phi(37) = phi(37) + phi(5)
              phi(38) = phi(38) + phi(6)
              phi(39) = phi(39) + phi(7)
              phi(40) = phi(40) + phi(8)
              phi(41) = phi(41) + phi(9)
              phi(42) = phi(42) + phi(10)
              phi(43) = phi(43) + phi(11)

              phi(51) = phi(51) + phi(12)
              phi(52) = phi(52) + phi(13)
              phi(53) = phi(53) + phi(14)
              phi(54) = phi(54) + phi(15)
              phi(55) = phi(55) + phi(16)
              phi(56) = phi(56) + phi(17)
              phi(57) = phi(57) + phi(18)
              phi(58) = phi(58) + phi(19)
              phi(59) = phi(59) + phi(20)
              phi(60) = phi(60) + phi(21)
              phi(61) = phi(61) + phi(22)
              phi(62) = phi(62) + phi(23)
              phi(63) = phi(63) + phi(24)

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
            phi(26) = phi(26) + scale_factor2(2)*scale_factor1(3)*phi(25)
            phi(27) = phi(27) + scale_factor2(4)*scale_factor1(4)*phi(25)

            phi(31) = phi(31) + scale_factor1(3)*phi(28)
            phi(32) = phi(32) + scale_factor1(3)*phi(29)
            phi(33) = phi(33) + scale_factor1(3)*phi(30)
            phi(34) = phi(34) + scale_factor2(3)*scale_factor1(4)*phi(28)
            phi(35) = phi(35) + scale_factor2(3)*scale_factor1(4)*phi(29)
            phi(36) = phi(36) + scale_factor2(3)*scale_factor1(4)*phi(30)

            phi(44) = phi(44) + scale_factor2(2)*scale_factor1(4)*phi(37)
            phi(45) = phi(45) + scale_factor2(2)*scale_factor1(4)*phi(38)
            phi(46) = phi(46) + scale_factor2(2)*scale_factor1(4)*phi(39)
            phi(47) = phi(47) + scale_factor2(2)*scale_factor1(4)*phi(40)
            phi(48) = phi(48) + scale_factor2(2)*scale_factor1(4)*phi(41)
            phi(49) = phi(49) + scale_factor2(2)*scale_factor1(4)*phi(42)
            phi(50) = phi(50) + scale_factor2(2)*scale_factor1(4)*phi(43)

            phi(64) = phi(64) + scale_factor1(4)*phi(51)
            phi(65) = phi(65) + scale_factor1(4)*phi(52)
            phi(66) = phi(66) + scale_factor1(4)*phi(53)
            phi(67) = phi(67) + scale_factor1(4)*phi(54)
            phi(68) = phi(68) + scale_factor1(4)*phi(55)
            phi(69) = phi(69) + scale_factor1(4)*phi(56)
            phi(70) = phi(70) + scale_factor1(4)*phi(57)
            phi(71) = phi(71) + scale_factor1(4)*phi(58)
            phi(72) = phi(72) + scale_factor1(4)*phi(59)
            phi(73) = phi(73) + scale_factor1(4)*phi(60)
            phi(74) = phi(74) + scale_factor1(4)*phi(61)
            phi(75) = phi(75) + scale_factor1(4)*phi(62)
            phi(76) = phi(76) + scale_factor1(4)*phi(63)

          end do ! k

!                   ******************************
!                   *                            *
!                   *   Post-contraction phase   *
!                   *                            *
!                   ******************************

          phi(45) = phi(45) - phi(44)
          phi(47) = phi(47) - phi(44)
          phi(50) = phi(50) - phi(44)

          phi(67) = phi(67) - phi(64)
          phi(68) = phi(68) - phi(65)
          phi(70) = phi(70) - phi(65)
          phi(71) = phi(71) - phi(66)
          phi(73) = phi(73) - phi(66)
          phi(76) = phi(76) - phi(66)
          phi(67) = phi(67) - (2.0D+00)*phi(64)
          phi(69) = phi(69) - phi(64)
          phi(70) = phi(70) - (2.0D+00)*phi(65)
          phi(74) = phi(74) - phi(64)
          phi(75) = phi(75) - phi(65)
          phi(76) = phi(76) - (2.0D+00)*phi(66)

!           *************************************
!           *                                   *
!           * -- Angular momentum conversion -- *
!           *   1-center (r) to 2-center (p|q)  *
!           *  Bra contracted transfer equation *
!           *   Horizontal recurrence relation  *
!           *                                   *
!           *************************************

          work2(1, 1) = phi(27)
          work2(2, 1) = phi(26)
          work2(3, 1) = -phi(34)
          work2(4, 1) = -phi(35)
          work2(5, 1) = -phi(36)
          work2(6, 1) = -phi(31)
          work2(7, 1) = -phi(32)
          work2(8, 1) = -phi(33)
          work2(9, 1) = phi(45)
          work2(10, 1) = phi(46)
          work2(11, 1) = phi(47)
          work2(12, 1) = phi(48)
          work2(13, 1) = phi(49)
          work2(14, 1) = phi(50)
          work2(15, 1) = -phi(67)
          work2(16, 1) = -phi(68)
          work2(17, 1) = -phi(69)
          work2(18, 1) = -phi(70)
          work2(19, 1) = -phi(71)
          work2(20, 1) = -phi(72)
          work2(21, 1) = -phi(73)
          work2(22, 1) = -phi(74)
          work2(23, 1) = -phi(75)
          work2(24, 1) = -phi(76)

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

          ip = 1
          i = 1
          ii1 = i + loci
          ijp = ip

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
!

        end if ! test.gt.cutoff_schwarz
      end do ! iquart
!$omp end target teams distribute parallel do

    end do ! itile

  end subroutine int0030
end submodule
