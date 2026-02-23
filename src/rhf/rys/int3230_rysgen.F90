! The total angular momentum of this class is:           8
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3230_impl
contains
  module subroutine int3230(df_pair, sf_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, sf_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n03ket(:)
    real(dp), allocatable :: xint23bra(:), xint03ket(:)
    integer(kind=int64) :: ndfbra, nsfket
    real(dp) :: scutdfbra, scutsfket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2
    integer(kind=int64) :: n, i1, i3, i4, i5, k3, k4, nn, nm, km, nj, ni, nl, nk
    real(dp) :: cp10, c10, cp01, c01
    integer(kind=int64) :: nx, ny, nz, mx, my, mz
    integer(kind=int64) :: bra_loop, ket_loop, ijtop, kltop
    real(dp) :: t_expon_ab, t_expon_cd, t_inverse_expon_ab, t_inverse_expon_cd
    real(dp) :: t_expon_abcd_inverse, rho, expe, dum, rab, rcd
    real(dp) :: brrk, akxk, akyk, akzk, t_expon_c, t_expon_d, t_expon_a, t_expon_b
    real(dp) :: xa, ya, za, axak, ayak, azak, axai, ayai, azai, bbrrk, xb, yb, zb, bxbk
    real(dp) :: bybk, bzbk, bxbi, bybi, bzbi, xx, c1x, c2x, c3x, c4x, c1y, c2y, c3y, c4y
    real(dp) :: c1z, c2z, c3z, c4z, f00, u2, duminv, dm2inv, bp01, b00, b10, xcp00, xc00
    real(dp) :: ycp00, zcp00, zc00, yc00, dij, dxij, dyij, dzij
    real(dp) :: buff(9)
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(240), yin(240), zin(240)
    real(dp) :: eri_value(600)
    real(dp) :: d23bra(60), d03ket(10)
    integer(kind=int64) :: ix(10), jx(6), kx(10), lx(1)
    integer(kind=int64) :: iy(10), jy(6), ky(10), ly(1)
    integer(kind=int64) :: iz(10), jz(6), kz(10), lz(1)
    integer(kind=int64) :: in(6), in1(6), kn(4)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(10), kly(10), klz(10)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 13
    in1(3) = 25
    in1(4) = 37
    in1(5) = 41
    in1(6) = 45

    kn(1) = 0
    kn(2) = 1
    kn(3) = 2
    kn(4) = 3

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 0

    kx(1) = 3
    kx(2) = 0
    kx(3) = 0
    kx(4) = 2
    kx(5) = 2
    kx(6) = 1
    kx(7) = 0
    kx(8) = 1
    kx(9) = 0
    kx(10) = 1

    jx(1) = 8
    jx(2) = 0
    jx(3) = 0
    jx(4) = 4
    jx(5) = 4
    jx(6) = 0

    ix(1) = 37
    ix(2) = 1
    ix(3) = 1
    ix(4) = 25
    ix(5) = 25
    ix(6) = 13
    ix(7) = 1
    ix(8) = 13
    ix(9) = 1
    ix(10) = 13

    ! y-arrays

    ly(1) = 0

    ky(1) = 0
    ky(2) = 3
    ky(3) = 0
    ky(4) = 1
    ky(5) = 0
    ky(6) = 2
    ky(7) = 2
    ky(8) = 0
    ky(9) = 1
    ky(10) = 1

    jy(1) = 0
    jy(2) = 8
    jy(3) = 0
    jy(4) = 4
    jy(5) = 0
    jy(6) = 4

    iy(1) = 1
    iy(2) = 37
    iy(3) = 1
    iy(4) = 13
    iy(5) = 1
    iy(6) = 25
    iy(7) = 25
    iy(8) = 1
    iy(9) = 13
    iy(10) = 13

    ! z-arrays

    lz(1) = 0

    kz(1) = 0
    kz(2) = 0
    kz(3) = 3
    kz(4) = 0
    kz(5) = 1
    kz(6) = 0
    kz(7) = 1
    kz(8) = 2
    kz(9) = 2
    kz(10) = 1

    jz(1) = 0
    jz(2) = 0
    jz(3) = 8
    jz(4) = 0
    jz(5) = 4
    jz(6) = 4

    iz(1) = 1
    iz(2) = 1
    iz(3) = 37
    iz(4) = 1
    iz(5) = 13
    iz(6) = 1
    iz(7) = 13
    iz(8) = 25
    iz(9) = 25
    iz(10) = 13

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 45
    ijx(2) = 37
    ijx(3) = 37
    ijx(4) = 41
    ijx(5) = 41
    ijx(6) = 37
    ijx(7) = 9
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 5
    ijx(11) = 5
    ijx(12) = 1
    ijx(13) = 9
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 5
    ijx(17) = 5
    ijx(18) = 1
    ijx(19) = 33
    ijx(20) = 25
    ijx(21) = 25
    ijx(22) = 29
    ijx(23) = 29
    ijx(24) = 25
    ijx(25) = 33
    ijx(26) = 25
    ijx(27) = 25
    ijx(28) = 29
    ijx(29) = 29
    ijx(30) = 25
    ijx(31) = 21
    ijx(32) = 13
    ijx(33) = 13
    ijx(34) = 17
    ijx(35) = 17
    ijx(36) = 13
    ijx(37) = 9
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 5
    ijx(41) = 5
    ijx(42) = 1
    ijx(43) = 21
    ijx(44) = 13
    ijx(45) = 13
    ijx(46) = 17
    ijx(47) = 17
    ijx(48) = 13
    ijx(49) = 9
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 5
    ijx(53) = 5
    ijx(54) = 1
    ijx(55) = 21
    ijx(56) = 13
    ijx(57) = 13
    ijx(58) = 17
    ijx(59) = 17
    ijx(60) = 13

    ijy(1) = 1
    ijy(2) = 9
    ijy(3) = 1
    ijy(4) = 5
    ijy(5) = 1
    ijy(6) = 5
    ijy(7) = 37
    ijy(8) = 45
    ijy(9) = 37
    ijy(10) = 41
    ijy(11) = 37
    ijy(12) = 41
    ijy(13) = 1
    ijy(14) = 9
    ijy(15) = 1
    ijy(16) = 5
    ijy(17) = 1
    ijy(18) = 5
    ijy(19) = 13
    ijy(20) = 21
    ijy(21) = 13
    ijy(22) = 17
    ijy(23) = 13
    ijy(24) = 17
    ijy(25) = 1
    ijy(26) = 9
    ijy(27) = 1
    ijy(28) = 5
    ijy(29) = 1
    ijy(30) = 5
    ijy(31) = 25
    ijy(32) = 33
    ijy(33) = 25
    ijy(34) = 29
    ijy(35) = 25
    ijy(36) = 29
    ijy(37) = 25
    ijy(38) = 33
    ijy(39) = 25
    ijy(40) = 29
    ijy(41) = 25
    ijy(42) = 29
    ijy(43) = 1
    ijy(44) = 9
    ijy(45) = 1
    ijy(46) = 5
    ijy(47) = 1
    ijy(48) = 5
    ijy(49) = 13
    ijy(50) = 21
    ijy(51) = 13
    ijy(52) = 17
    ijy(53) = 13
    ijy(54) = 17
    ijy(55) = 13
    ijy(56) = 21
    ijy(57) = 13
    ijy(58) = 17
    ijy(59) = 13
    ijy(60) = 17

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 9
    ijz(4) = 1
    ijz(5) = 5
    ijz(6) = 5
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 9
    ijz(10) = 1
    ijz(11) = 5
    ijz(12) = 5
    ijz(13) = 37
    ijz(14) = 37
    ijz(15) = 45
    ijz(16) = 37
    ijz(17) = 41
    ijz(18) = 41
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 9
    ijz(22) = 1
    ijz(23) = 5
    ijz(24) = 5
    ijz(25) = 13
    ijz(26) = 13
    ijz(27) = 21
    ijz(28) = 13
    ijz(29) = 17
    ijz(30) = 17
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 9
    ijz(34) = 1
    ijz(35) = 5
    ijz(36) = 5
    ijz(37) = 13
    ijz(38) = 13
    ijz(39) = 21
    ijz(40) = 13
    ijz(41) = 17
    ijz(42) = 17
    ijz(43) = 25
    ijz(44) = 25
    ijz(45) = 33
    ijz(46) = 25
    ijz(47) = 29
    ijz(48) = 29
    ijz(49) = 25
    ijz(50) = 25
    ijz(51) = 33
    ijz(52) = 25
    ijz(53) = 29
    ijz(54) = 29
    ijz(55) = 13
    ijz(56) = 13
    ijz(57) = 21
    ijz(58) = 13
    ijz(59) = 17
    ijz(60) = 17

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 3
    klx(2) = 0
    klx(3) = 0
    klx(4) = 2
    klx(5) = 2
    klx(6) = 1
    klx(7) = 0
    klx(8) = 1
    klx(9) = 0
    klx(10) = 1

    kly(1) = 0
    kly(2) = 3
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 2
    kly(7) = 2
    kly(8) = 0
    kly(9) = 1
    kly(10) = 1

    klz(1) = 0
    klz(2) = 0
    klz(3) = 3
    klz(4) = 0
    klz(5) = 1
    klz(6) = 0
    klz(7) = 1
    klz(8) = 2
    klz(9) = 2
    klz(10) = 1

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n03ket(res%n_s_shl*res%n_f_shl))
    allocate (xint03ket(res%n_s_shl*res%n_f_shl))

    ! Start screening

    scutdfbra = cutoff_schwarz/maxval(df_pair%xints)
    ndfbra = 0
    do ij = 1, res%n_d_shl*res%n_f_shl
      if (df_pair%xints(ij) .ge. scutdfbra) then
        ndfbra = ndfbra + 1
        xint23bra(ndfbra) = df_pair%xints(ij)
        n23bra(ndfbra) = ij
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

    nchunksize_int64 = 375000000

    if ((ndfbra*nsfket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nsfket
    ntile = int(ndfbra*nsfket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nsfket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint03ket, n03ket, df_pair, sf_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d03ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint03ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n03ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                lsh_tmp = (kl - 1)/res%n_f_shl + 1

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
                ksh = res%i_f_shl(ksh_tmp)
                lsh = res%i_s_shl(lsh_tmp)

                ijtop = res%contr_num(ish)*res%contr_num(jsh)

                kltop = res%contr_num(ksh)*res%contr_num(lsh)

                ! --- Distances between shells ----

                rab = ((res%coord_sh(ish, 1) - res%coord_sh(jsh, 1))*(res%coord_sh(ish, 1) - res%coord_sh(jsh, 1)) + &
                       (res%coord_sh(ish, 2) - res%coord_sh(jsh, 2))*(res%coord_sh(ish, 2) - res%coord_sh(jsh, 2)) + &
                       (res%coord_sh(ish, 3) - res%coord_sh(jsh, 3))*(res%coord_sh(ish, 3) - res%coord_sh(jsh, 3)))

                rcd = ((res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1))*(res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1)) + &
                       (res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2))*(res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2)) + &
                       (res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3))*(res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3)))

                ! Distances used in the bra transfer equation

                dxij = res%coord_sh(ish, 1) - res%coord_sh(jsh, 1)
                dyij = res%coord_sh(ish, 2) - res%coord_sh(jsh, 2)
                dzij = res%coord_sh(ish, 3) - res%coord_sh(jsh, 3)

                ! Prepare to begin looping

                eri_value = 0.0_dp

                ket_loop = 0

                ! --- Start looping over primitives ---

                do k = 1, kltop

                  ket_loop = ket_loop + 1

                  t_expon_cd = sf_pair%t_expon_ab(sf_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = sf_pair%expon_b(sf_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = sf_pair%expon_a(sf_pair%pair_loc(kl) + ket_loop)
                  t_inverse_expon_cd = 1.0_dp/t_expon_cd

                  brrk = t_expon_c*rcd                ! Same variable names
                  akxk = t_expon_c*res%coord_sh(ksh, 1)    ! as in the original
                  akyk = t_expon_c*res%coord_sh(ksh, 2)    ! GAMESS subroutine
                  akzk = t_expon_c*res%coord_sh(ksh, 3)

                  bbrrk = t_expon_d*brrk*t_inverse_expon_cd
                  xb = (akxk + t_expon_d*res%coord_sh(lsh, 1))*t_inverse_expon_cd
                  yb = (akyk + t_expon_d*res%coord_sh(lsh, 2))*t_inverse_expon_cd
                  zb = (akzk + t_expon_d*res%coord_sh(lsh, 3))*t_inverse_expon_cd

                  bxbk = t_expon_cd*(xb - res%coord_sh(ksh, 1))
                  bybk = t_expon_cd*(yb - res%coord_sh(ksh, 2))
                  bzbk = t_expon_cd*(zb - res%coord_sh(ksh, 3))

                  bxbi = t_expon_cd*(xb - res%coord_sh(ish, 1))
                  bybi = t_expon_cd*(yb - res%coord_sh(ish, 2))
                  bzbi = t_expon_cd*(zb - res%coord_sh(ish, 3))

                  ! Contraction coefficient factor + angular momentum normalization
                  ! This could be optimized by moving twopi_5_2 to the bra
                  ! Whether this would have much impact or be benfictial at all
                  ! is debatable

                  d03ket(1) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(2) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(3) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d03ket(4) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(5) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(6) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(7) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(8) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(9) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                  d03ket(10) = sf_pair%d_coeff_alt(sf_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3

                  bra_loop = 0

                  do i = 1, ijtop

                    bra_loop = bra_loop + 1

                    t_expon_ab = df_pair%t_expon_ab(df_pair%pair_loc(ij) + bra_loop)
                    t_expon_a = df_pair%expon_b(df_pair%pair_loc(ij) + bra_loop)
                    t_expon_b = df_pair%expon_a(df_pair%pair_loc(ij) + bra_loop)
                    t_inverse_expon_ab = 1.0_dp/t_expon_ab

                    dum = bbrrk + t_expon_a*t_expon_b*rab*t_inverse_expon_ab
                    t_expon_abcd_inverse = 1.0_dp/(t_expon_ab + t_expon_cd)
                    expe = exp(-dum)*sqrt(t_expon_abcd_inverse)
                    rho = t_expon_ab*t_expon_cd*t_expon_abcd_inverse

                    xa = (t_expon_a*res%coord_sh(ish, 1) + t_expon_b*res%coord_sh(jsh, 1))
                    ya = (t_expon_a*res%coord_sh(ish, 2) + t_expon_b*res%coord_sh(jsh, 2))
                    za = (t_expon_a*res%coord_sh(ish, 3) + t_expon_b*res%coord_sh(jsh, 3))

                    xa = xa*t_inverse_expon_ab
                    ya = ya*t_inverse_expon_ab
                    za = za*t_inverse_expon_ab

                    axak = t_expon_ab*(xa - res%coord_sh(ksh, 1))
                    ayak = t_expon_ab*(ya - res%coord_sh(ksh, 2))
                    azak = t_expon_ab*(za - res%coord_sh(ksh, 3))

                    axai = t_expon_ab*(xa - res%coord_sh(ish, 1))
                    ayai = t_expon_ab*(ya - res%coord_sh(ish, 2))
                    azai = t_expon_ab*(za - res%coord_sh(ish, 3))

                    c1x = bxbk + axak
                    c2x = t_expon_ab*bxbk
                    c3x = bxbi + axai
                    c4x = t_expon_cd*axai

                    c1y = bybk + ayak
                    c2y = t_expon_ab*bybk
                    c3y = bybi + ayai
                    c4y = t_expon_cd*ayai

                    c1z = bzbk + azak
                    c2z = t_expon_ab*bzbk
                    c3z = bzbi + azai
                    c4z = t_expon_cd*azai

                    ! Contraction coefficient factor + angular momentum normalization
                    ! Multiply only when the array element is not 1.0_dp to begin with

                    d23bra(1) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(2) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(3) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(4) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(5) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(6) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(7) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(8) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(9) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(10) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(11) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(12) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(13) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(14) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(15) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)
                    d23bra(16) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(17) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(18) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt3
                    d23bra(19) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(20) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(21) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(22) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(23) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(24) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(25) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(26) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(27) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(28) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(29) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(30) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(31) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(32) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(33) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(34) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(35) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(36) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(37) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(38) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(39) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(40) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(41) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(42) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(43) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(44) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(45) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(46) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(47) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(48) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(49) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(50) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(51) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5
                    d23bra(52) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(53) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(54) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(55) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(56) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(57) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3
                    d23bra(58) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3
                    d23bra(59) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3
                    d23bra(60) = df_pair%d_coeff_alt(df_pair%pair_loc(ij) + bra_loop)*sqrt5*sqrt3*sqrt3

                    ! Compute xx and evaluate for the respective number of roots and weights

                    xx = rho*((xa - xb)*(xa - xb) + (ya - yb)*(ya - yb) + (za - zb)*(za - zb))

                    ! Evaluate roots and weights

                    if (xx .ge. 55.0D+00) then ! Asymptotic form

                      factr = 1.0_dp/xx
                      factw = sqrt(factr)

                      rts(1) = factr*0.1175813202117781D+00
                      rts(2) = factr*0.1074562012436902D+01
                      rts(3) = factr*0.3085937443717543D+01
                      rts(4) = factr*0.6414729733662035D+01
                      rts(5) = factr*0.1180718948997174D+02

                      wts(1) = factw*0.6108626337353259D+00
                      wts(2) = factw*0.2401386110823149D+00
                      wts(3) = factw*0.3387439445548109D-01
                      wts(4) = factw*0.1343645746781227D-02
                      wts(5) = factw*0.7640432855232614D-05

                    else ! "regular" evaluation

                      rgrid(1) = 0.1314956323727580D-05
                      rgrid(2) = 0.3638644488856184D-04
                      rgrid(3) = 0.2184836363610052D-03
                      rgrid(4) = 0.7467882061060832D-03
                      rgrid(5) = 0.1898594940807536D-02
                      rgrid(6) = 0.4018347564842567D-02
                      rgrid(7) = 0.7503899366138534D-02
                      rgrid(8) = 0.1279045049262753D-01
                      rgrid(9) = 0.2033269213059855D-01
                      rgrid(10) = 0.3058574876834499D-01
                      rgrid(11) = 0.4398555195277359D-01
                      rgrid(12) = 0.6092930105499368D-01
                      rgrid(13) = 0.8175666785532240D-01
                      rgrid(14) = 0.1067323821680571D+00
                      rgrid(15) = 0.1360307958356441D+00
                      rgrid(16) = 0.1697229634514230D+00
                      rgrid(17) = 0.2077667019402563D+00
                      rgrid(18) = 0.2499999999999998D+00
                      rgrid(19) = 0.2961380452159158D+00
                      rgrid(20) = 0.3457740246174122D+00
                      rgrid(21) = 0.3983837370449402D+00
                      rgrid(22) = 0.4533339365988709D+00
                      rgrid(23) = 0.5098942093731368D+00
                      rgrid(24) = 0.5672520742964827D+00
                      rgrid(25) = 0.6245308967025374D+00
                      rgrid(26) = 0.6808101134342349D+00
                      rgrid(27) = 0.7351471936872270D+00
                      rgrid(28) = 0.7866007027795398D+00
                      rgrid(29) = 0.8342537984583634D+00
                      rgrid(30) = 0.8772374725900651D+00
                      rgrid(31) = 0.9147528563001249D+00
                      rgrid(32) = 0.9460919364139329D+00
                      rgrid(33) = 0.9706560996755904D+00
                      rgrid(34) = 0.9879721508887392D+00
                      rgrid(35) = 0.9977078840559234D+00

                      wgrid(1) = 0.2941716710221880D-02*exp(-xx*0.1314956323727580D-05)
                      wgrid(2) = 0.6825414174180503D-02*exp(-xx*0.3638644488856184D-04)
                      wgrid(3) = 0.1066148995574242D-01*exp(-xx*0.2184836363610052D-03)
                      wgrid(4) = 0.1441463005444643D-01*exp(-xx*0.7467882061060832D-03)
                      wgrid(5) = 0.1805505793173160D-01*exp(-xx*0.1898594940807536D-02)
                      wgrid(6) = 0.2155421116308512D-01*exp(-xx*0.4018347564842567D-02)
                      wgrid(7) = 0.2488468520067677D-01*exp(-xx*0.7503899366138534D-02)
                      wgrid(8) = 0.2802040810618503D-01*exp(-xx*0.1279045049262753D-01)
                      wgrid(9) = 0.3093683598304008D-01*exp(-xx*0.2033269213059855D-01)
                      wgrid(10) = 0.3361114263454359D-01*exp(-xx*0.3058574876834499D-01)
                      wgrid(11) = 0.3602239738628011D-01*exp(-xx*0.4398555195277359D-01)
                      wgrid(12) = 0.3815172857772094D-01*exp(-xx*0.6092930105499368D-01)
                      wgrid(13) = 0.3998247112116213D-01*exp(-xx*0.8175666785532240D-01)
                      wgrid(14) = 0.4150029686442836D-01*exp(-xx*0.1067323821680571D+00)
                      wgrid(15) = 0.4269332669604951D-01*exp(-xx*0.1360307958356441D+00)
                      wgrid(16) = 0.4355222349859141D-01*exp(-xx*0.1697229634514230D+00)
                      wgrid(17) = 0.4407026521513794D-01*exp(-xx*0.2077667019402563D+00)
                      wgrid(18) = 0.4424339745355207D-01*exp(-xx*0.2499999999999998D+00)
                      wgrid(19) = 0.4407026521513754D-01*exp(-xx*0.2961380452159158D+00)
                      wgrid(20) = 0.4355222349859152D-01*exp(-xx*0.3457740246174122D+00)
                      wgrid(21) = 0.4269332669604975D-01*exp(-xx*0.3983837370449402D+00)
                      wgrid(22) = 0.4150029686442849D-01*exp(-xx*0.4533339365988709D+00)
                      wgrid(23) = 0.3998247112116278D-01*exp(-xx*0.5098942093731368D+00)
                      wgrid(24) = 0.3815172857772100D-01*exp(-xx*0.5672520742964827D+00)
                      wgrid(25) = 0.3602239738628051D-01*exp(-xx*0.6245308967025374D+00)
                      wgrid(26) = 0.3361114263454346D-01*exp(-xx*0.6808101134342349D+00)
                      wgrid(27) = 0.3093683598304035D-01*exp(-xx*0.7351471936872270D+00)
                      wgrid(28) = 0.2802040810618515D-01*exp(-xx*0.7866007027795398D+00)
                      wgrid(29) = 0.2488468520067681D-01*exp(-xx*0.8342537984583634D+00)
                      wgrid(30) = 0.2155421116308484D-01*exp(-xx*0.8772374725900651D+00)
                      wgrid(31) = 0.1805505793173228D-01*exp(-xx*0.9147528563001249D+00)
                      wgrid(32) = 0.1441463005444692D-01*exp(-xx*0.9460919364139329D+00)
                      wgrid(33) = 0.1066148995574185D-01*exp(-xx*0.9706560996755904D+00)
                      wgrid(34) = 0.6825414174181070D-02*exp(-xx*0.9879721508887392D+00)
                      wgrid(35) = 0.2941716710221409D-02*exp(-xx*0.9977078840559234D+00)

                      ! Call to RYSDS

                      sum0 = 0.0D+00
                      sum1 = 0.0D+00

                      do m = 1, 35
                        sum0 = sum0 + wgrid(m)
                        sum1 = sum1 + wgrid(m)*rgrid(m)
                      end do

                      alpha(1) = sum1/sum0
                      beta(1) = sum0

                      do m = 1, 35
                        p1(m) = 0.0D+00
                        p2(m) = 1.0D+00
                      end do

                      do kk = 1, 4

                        sum1 = 0.0D+00
                        sum2 = 0.0D+00

                        do 30 m = 1, 35

                          if (wgrid(m) .eq. 0.0d+00) goto 30

                          p0(m) = p1(m)
                          p1(m) = p2(m)
                          p2(m) = (rgrid(m) - alpha(kk))*p1(m) - beta(kk)*p0(m)
                          t = wgrid(m)*p2(m)*p2(m)
                          sum1 = sum1 + t
                          sum2 = sum2 + t*rgrid(m)

30                        continue

                          alpha(kk + 1) = sum2/sum1
                          beta(kk + 1) = sum1/sum0
                          sum0 = sum1

                        end do

                        ! End of RYSDS

                        ! Call to RYSGW

                        rts(1) = alpha(1)
                        wts(1) = 1.0D+00
                        wrk(5) = 0.0D+00
                        do 100 kk = 2, 5

                          rts(kk) = alpha(kk)
                          wrk(kk - 1) = sqrt(beta(kk))
                          wts(kk) = 0.0D+00

100                       continue

                          do 240 l = 1, 5

                            jj = 0

105                         do 110 m = l, 5
                              if (m .eq. 5) go to 120
                              if (abs(wrk(m)) .le. (1.0D-14)*(abs(rts(m)) + abs(rts(m + 1)))) go to 120
110                           continue

120                           dpp = rts(l)
                              if (m .eq. l) go to 240
                              jj = jj + 1

                              dg = (rts(l + 1) - dpp)/(2.0D+00*wrk(l))
                              dr = sqrt(dg*dg + 1.0D+00)
                              dg = rts(m) - dpp + wrk(l)/(dg + sign(dr, dg))
                              ds = 1.0D+00
                              dc = 1.0D+00
                              dpp = 0.0D+00
                              mml = m - l

                              do 200 ii = 1, mml

                                mmii = m - ii
                                df = ds*wrk(mmii)
                                db = dc*wrk(mmii)
                                if (abs(df) .lt. abs(dg)) go to 150
                                dc = dg/df
                                dr = sqrt(dc*dc + 1.0D+00)
                                wrk(mmii + 1) = df*dr
                                ds = 1.0D+00/dr
                                dc = dc*ds
                                go to 160
150                             ds = df/dg
                                dr = sqrt(ds*ds + 1.0D+00)
                                wrk(mmii + 1) = dg*dr
                                dc = 1.0D+00/dr
                                ds = ds*dc
160                             dg = rts(mmii + 1) - dpp
                                dr = (rts(mmii) - dg)*ds + 2.0D+00*dc*db
                                dpp = ds*dr
                                rts(mmii + 1) = dg + dpp
                                dg = dc*dr - db
                                df = wts(mmii + 1)
                                wts(mmii + 1) = ds*wts(mmii) + dc*df
                                wts(mmii) = dc*wts(mmii) - ds*df

200                             continue

                                rts(l) = rts(l) - dpp
                                wrk(l) = dg
                                wrk(m) = 0.0D+00
                                go to 105

240                             continue

                                do 300 ii = 2, 5

                                  iim1 = ii - 1
                                  kk = iim1
                                  dpp = rts(iim1)

                                  do 260 jj = ii, 5
                                    if (rts(jj) .ge. dpp) go to 260
                                    kk = jj
                                    dpp = rts(jj)
260                                 continue

                                    if (kk .eq. iim1) go to 300

                                    rts(kk) = rts(iim1)
                                    rts(iim1) = dpp
                                    dpp = wts(iim1)
                                    wts(iim1) = wts(kk)
                                    wts(kk) = dpp

300                                 continue

                                    do 310 kk = 1, 5
                                      wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                   continue

                                      end if

                                      do kk = 1, 5
                                        roots(kk) = rts(kk)/(1.0_dp - rts(kk))
                                        wghts(kk) = wts(kk)
                                      end do

                                      ! Start computation of 2.0_dp-electron integrals for each root

                                      ! mm = 0
                                      ! do m = 1, nroots

                                      u2 = roots(1)*rho
                                      f00 = expe*wghts(1)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =    1

                                      xin(1) = 1.0_dp
                                      yin(1) = 1.0_dp
                                      zin(1) = f00

                                      ! i2 = in(2) =   13
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(13) = xc00
                                      yin(13) = yc00
                                      zin(13) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    2

                                      xin(2) = xcp00
                                      yin(2) = ycp00
                                      zin(2) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   14
                                      ! i2 =   13

                                      xin(14) = xcp00*xin(13) + cp10
                                      yin(14) = ycp00*yin(13) + cp10
                                      zin(14) = zcp00*zin(13) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   25
                                      ! i3 =    1
                                      ! i4 =   13

                                      xin(25) = c10*xin(1) + xc00*xin(13)
                                      yin(25) = c10*yin(1) + yc00*yin(13)
                                      zin(25) = c10*zin(1) + zc00*zin(13)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   26
                                      ! i5 =   25
                                      ! i4 =   13

                                      xin(26) = xcp00*xin(25) + cp10*xin(13)
                                      yin(26) = ycp00*yin(25) + cp10*yin(13)
                                      zin(26) = zcp00*zin(25) + cp10*zin(13)

                                      ! ------------------

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   37
                                      ! i3 =   13
                                      ! i4 =   25

                                      xin(37) = c10*xin(13) + xc00*xin(25)
                                      yin(37) = c10*yin(13) + yc00*yin(25)
                                      zin(37) = c10*zin(13) + zc00*zin(25)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   38
                                      ! i5 =   37
                                      ! i4 =   25

                                      xin(38) = xcp00*xin(37) + cp10*xin(25)
                                      yin(38) = ycp00*yin(37) + cp10*yin(25)
                                      zin(38) = zcp00*zin(37) + cp10*zin(25)

                                      ! ------------------

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   41
                                      ! i3 =   25
                                      ! i4 =   37

                                      xin(41) = c10*xin(25) + xc00*xin(37)
                                      yin(41) = c10*yin(25) + yc00*yin(37)
                                      zin(41) = c10*zin(25) + zc00*zin(37)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   42
                                      ! i5 =   41
                                      ! i4 =   37

                                      xin(42) = xcp00*xin(41) + cp10*xin(37)
                                      yin(42) = ycp00*yin(41) + cp10*yin(37)
                                      zin(42) = zcp00*zin(41) + cp10*zin(37)

                                      ! ------------------

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   45
                                      ! i3 =   37
                                      ! i4 =   41

                                      xin(45) = c10*xin(37) + xc00*xin(41)
                                      yin(45) = c10*yin(37) + yc00*yin(41)
                                      zin(45) = c10*zin(37) + zc00*zin(41)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   46
                                      ! i5 =   45
                                      ! i4 =   41

                                      xin(46) = xcp00*xin(45) + cp10*xin(41)
                                      yin(46) = ycp00*yin(45) + cp10*yin(41)
                                      zin(46) = zcp00*zin(45) + cp10*zin(41)

                                      ! ------------------

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    2

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    3
                                      ! i3 =    1
                                      ! i4 =    2

                                      xin(3) = cp01*xin(1) + xcp00*xin(2)
                                      yin(3) = cp01*yin(1) + ycp00*yin(2)
                                      zin(3) = cp01*zin(1) + zcp00*zin(2)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   15

                                      xin(15) = xc00*xin(3) + c01*xin(2)
                                      yin(15) = yc00*yin(3) + c01*yin(2)
                                      zin(15) = zc00*zin(3) + c01*zin(2)

                                      ! ------------------

                                      ! i3 = i4 =    2
                                      ! i4 = i5 =    3

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    4
                                      ! i3 =    2
                                      ! i4 =    3

                                      xin(4) = cp01*xin(2) + xcp00*xin(3)
                                      yin(4) = cp01*yin(2) + ycp00*yin(3)
                                      zin(4) = cp01*zin(2) + zcp00*zin(3)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   16

                                      xin(16) = xc00*xin(4) + c01*xin(3)
                                      yin(16) = yc00*yin(4) + c01*yin(3)
                                      zin(16) = zc00*zin(4) + c01*zin(3)

                                      ! ------------------

                                      ! i3 = i4 =    3
                                      ! i4 = i5 =    4

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   25

                                      xin(27) = c10*xin(3) + xc00*xin(15) + c01*xin(14)
                                      yin(27) = c10*yin(3) + yc00*yin(15) + c01*yin(14)
                                      zin(27) = c10*zin(3) + zc00*zin(15) + c01*zin(14)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(39) = c10*xin(15) + xc00*xin(27) + c01*xin(26)
                                      yin(39) = c10*yin(15) + yc00*yin(27) + c01*yin(26)
                                      zin(39) = c10*zin(15) + zc00*zin(27) + c01*zin(26)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   41

                                      xin(43) = c10*xin(27) + xc00*xin(39) + c01*xin(38)
                                      yin(43) = c10*yin(27) + yc00*yin(39) + c01*yin(38)
                                      zin(43) = c10*zin(27) + zc00*zin(39) + c01*zin(38)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   45

                                      xin(47) = c10*xin(39) + xc00*xin(43) + c01*xin(42)
                                      yin(47) = c10*yin(39) + yc00*yin(43) + c01*yin(42)
                                      zin(47) = c10*zin(39) + zc00*zin(43) + c01*zin(42)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   13

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   25

                                      xin(28) = c10*xin(4) + xc00*xin(16) + c01*xin(15)
                                      yin(28) = c10*yin(4) + yc00*yin(16) + c01*yin(15)
                                      zin(28) = c10*zin(4) + zc00*zin(16) + c01*zin(15)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   13
                                      ! i4 = i5 =   25

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   37

                                      xin(40) = c10*xin(16) + xc00*xin(28) + c01*xin(27)
                                      yin(40) = c10*yin(16) + yc00*yin(28) + c01*yin(27)
                                      zin(40) = c10*zin(16) + zc00*zin(28) + c01*zin(27)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   25
                                      ! i4 = i5 =   37

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   41

                                      xin(44) = c10*xin(28) + xc00*xin(40) + c01*xin(39)
                                      yin(44) = c10*yin(28) + yc00*yin(40) + c01*yin(39)
                                      zin(44) = c10*zin(28) + zc00*zin(40) + c01*zin(39)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   37
                                      ! i4 = i5 =   41

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   45

                                      xin(48) = c10*xin(40) + xc00*xin(44) + c01*xin(43)
                                      yin(48) = c10*yin(40) + yc00*yin(44) + c01*yin(43)
                                      zin(48) = c10*zin(40) + zc00*zin(44) + c01*zin(43)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   41
                                      ! i4 = i5 =   45

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   45

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(45) = xin(45) + dxij*xin(41)
                                      yin(45) = yin(45) + dyij*yin(41)
                                      zin(45) = zin(45) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   37

                                      xin(41) = xin(41) + dxij*xin(37)
                                      yin(41) = yin(41) + dyij*yin(37)
                                      zin(41) = zin(41) + dzij*zin(37)

                                      ! i3 = i4 =   37
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   45

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   41

                                      xin(45) = xin(45) + dxij*xin(41)
                                      yin(45) = yin(45) + dyij*yin(41)
                                      zin(45) = zin(45) + dzij*zin(41)

                                      ! i3 = i4 =   41
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    5

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    5

                                      ! do ni = 1,    3

                                      xin(5) = xin(13) + dxij*xin(1)
                                      yin(5) = yin(13) + dyij*yin(1)
                                      zin(5) = zin(13) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   17

                                      ! ni =    2

                                      xin(17) = xin(25) + dxij*xin(13)
                                      yin(17) = yin(25) + dyij*yin(13)
                                      zin(17) = zin(25) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   29

                                      ! ni =    3

                                      xin(29) = xin(37) + dxij*xin(25)
                                      yin(29) = yin(37) + dyij*yin(25)
                                      zin(29) = zin(37) + dzij*zin(25)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   41

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =    9

                                      ! nj =    2

                                      ! i4 = i3 =    9

                                      ! do ni = 1,    3

                                      xin(9) = xin(17) + dxij*xin(5)
                                      yin(9) = yin(17) + dyij*yin(5)
                                      zin(9) = zin(17) + dzij*zin(5)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   21

                                      ! ni =    2

                                      xin(21) = xin(29) + dxij*xin(17)
                                      yin(21) = yin(29) + dyij*yin(17)
                                      zin(21) = zin(29) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   33

                                      ! ni =    3

                                      xin(33) = xin(41) + dxij*xin(29)
                                      yin(33) = yin(41) + dyij*yin(29)
                                      zin(33) = zin(41) + dzij*zin(29)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   13

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   42

                                      xin(46) = xin(46) + dxij*xin(42)
                                      yin(46) = yin(46) + dyij*yin(42)
                                      zin(46) = zin(46) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   38

                                      xin(42) = xin(42) + dxij*xin(38)
                                      yin(42) = yin(42) + dyij*yin(38)
                                      zin(42) = zin(42) + dzij*zin(38)

                                      ! i3 = i4 =   38
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   46

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   42

                                      xin(46) = xin(46) + dxij*xin(42)
                                      yin(46) = yin(46) + dyij*yin(42)
                                      zin(46) = zin(46) + dzij*zin(42)

                                      ! i3 = i4 =   42
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    6

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    6

                                      ! do ni = 1,    3

                                      xin(6) = xin(14) + dxij*xin(2)
                                      yin(6) = yin(14) + dyij*yin(2)
                                      zin(6) = zin(14) + dzij*zin(2)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   18

                                      ! ni =    2

                                      xin(18) = xin(26) + dxij*xin(14)
                                      yin(18) = yin(26) + dyij*yin(14)
                                      zin(18) = zin(26) + dzij*zin(14)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   30

                                      ! ni =    3

                                      xin(30) = xin(38) + dxij*xin(26)
                                      yin(30) = yin(38) + dyij*yin(26)
                                      zin(30) = zin(38) + dzij*zin(26)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   42

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   10

                                      ! nj =    2

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(18) + dxij*xin(6)
                                      yin(10) = yin(18) + dyij*yin(6)
                                      zin(10) = zin(18) + dzij*zin(6)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   22

                                      ! ni =    2

                                      xin(22) = xin(30) + dxij*xin(18)
                                      yin(22) = yin(30) + dyij*yin(18)
                                      zin(22) = zin(30) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   34

                                      ! ni =    3

                                      xin(34) = xin(42) + dxij*xin(30)
                                      yin(34) = yin(42) + dyij*yin(30)
                                      zin(34) = zin(42) + dzij*zin(30)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   14

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(47) = xin(47) + dxij*xin(43)
                                      yin(47) = yin(47) + dyij*yin(43)
                                      zin(47) = zin(47) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   39

                                      xin(43) = xin(43) + dxij*xin(39)
                                      yin(43) = yin(43) + dyij*yin(39)
                                      zin(43) = zin(43) + dzij*zin(39)

                                      ! i3 = i4 =   39
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   47

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   43

                                      xin(47) = xin(47) + dxij*xin(43)
                                      yin(47) = yin(47) + dyij*yin(43)
                                      zin(47) = zin(47) + dzij*zin(43)

                                      ! i3 = i4 =   43
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    7

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    7

                                      ! do ni = 1,    3

                                      xin(7) = xin(15) + dxij*xin(3)
                                      yin(7) = yin(15) + dyij*yin(3)
                                      zin(7) = zin(15) + dzij*zin(3)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   19

                                      ! ni =    2

                                      xin(19) = xin(27) + dxij*xin(15)
                                      yin(19) = yin(27) + dyij*yin(15)
                                      zin(19) = zin(27) + dzij*zin(15)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   31

                                      ! ni =    3

                                      xin(31) = xin(39) + dxij*xin(27)
                                      yin(31) = yin(39) + dyij*yin(27)
                                      zin(31) = zin(39) + dzij*zin(27)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   11

                                      ! nj =    2

                                      ! i4 = i3 =   11

                                      ! do ni = 1,    3

                                      xin(11) = xin(19) + dxij*xin(7)
                                      yin(11) = yin(19) + dyij*yin(7)
                                      zin(11) = zin(19) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   23

                                      ! ni =    2

                                      xin(23) = xin(31) + dxij*xin(19)
                                      yin(23) = yin(31) + dyij*yin(19)
                                      zin(23) = zin(31) + dzij*zin(19)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   35

                                      ! ni =    3

                                      xin(35) = xin(43) + dxij*xin(31)
                                      yin(35) = yin(43) + dyij*yin(31)
                                      zin(35) = zin(43) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   47

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   15

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(48) = xin(48) + dxij*xin(44)
                                      yin(48) = yin(48) + dyij*yin(44)
                                      zin(48) = zin(48) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   40

                                      xin(44) = xin(44) + dxij*xin(40)
                                      yin(44) = yin(44) + dyij*yin(40)
                                      zin(44) = zin(44) + dzij*zin(40)

                                      ! i3 = i4 =   40
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   48

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   44

                                      xin(48) = xin(48) + dxij*xin(44)
                                      yin(48) = yin(48) + dyij*yin(44)
                                      zin(48) = zin(48) + dzij*zin(44)

                                      ! i3 = i4 =   44
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =    8

                                      ! do nj = 1,    2

                                      ! i4 = i3 =    8

                                      ! do ni = 1,    3

                                      xin(8) = xin(16) + dxij*xin(4)
                                      yin(8) = yin(16) + dyij*yin(4)
                                      zin(8) = zin(16) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   20

                                      ! ni =    2

                                      xin(20) = xin(28) + dxij*xin(16)
                                      yin(20) = yin(28) + dyij*yin(16)
                                      zin(20) = zin(28) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   32

                                      ! ni =    3

                                      xin(32) = xin(40) + dxij*xin(28)
                                      yin(32) = yin(40) + dyij*yin(28)
                                      zin(32) = zin(40) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   12

                                      ! nj =    2

                                      ! i4 = i3 =   12

                                      ! do ni = 1,    3

                                      xin(12) = xin(20) + dxij*xin(8)
                                      yin(12) = yin(20) + dyij*yin(8)
                                      zin(12) = zin(20) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   24

                                      ! ni =    2

                                      xin(24) = xin(32) + dxij*xin(20)
                                      yin(24) = yin(32) + dyij*yin(20)
                                      zin(24) = zin(32) + dzij*zin(20)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   36

                                      ! ni =    3

                                      xin(36) = xin(44) + dxij*xin(32)
                                      yin(36) = yin(44) + dyij*yin(32)
                                      zin(36) = zin(44) + dzij*zin(32)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   48

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   16

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   48

                                      u2 = roots(2)*rho
                                      f00 = expe*wghts(2)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =   49

                                      xin(49) = 1.0_dp
                                      yin(49) = 1.0_dp
                                      zin(49) = f00

                                      ! i2 = in(2) =   61
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(61) = xc00
                                      yin(61) = yc00
                                      zin(61) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   50

                                      xin(50) = xcp00
                                      yin(50) = ycp00
                                      zin(50) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   62
                                      ! i2 =   61

                                      xin(62) = xcp00*xin(61) + cp10
                                      yin(62) = ycp00*yin(61) + cp10
                                      zin(62) = zcp00*zin(61) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   73
                                      ! i3 =   49
                                      ! i4 =   61

                                      xin(73) = c10*xin(49) + xc00*xin(61)
                                      yin(73) = c10*yin(49) + yc00*yin(61)
                                      zin(73) = c10*zin(49) + zc00*zin(61)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   74
                                      ! i5 =   73
                                      ! i4 =   61

                                      xin(74) = xcp00*xin(73) + cp10*xin(61)
                                      yin(74) = ycp00*yin(73) + cp10*yin(61)
                                      zin(74) = zcp00*zin(73) + cp10*zin(61)

                                      ! ------------------

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   85
                                      ! i3 =   61
                                      ! i4 =   73

                                      xin(85) = c10*xin(61) + xc00*xin(73)
                                      yin(85) = c10*yin(61) + yc00*yin(73)
                                      zin(85) = c10*zin(61) + zc00*zin(73)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   86
                                      ! i5 =   85
                                      ! i4 =   73

                                      xin(86) = xcp00*xin(85) + cp10*xin(73)
                                      yin(86) = ycp00*yin(85) + cp10*yin(73)
                                      zin(86) = zcp00*zin(85) + cp10*zin(73)

                                      ! ------------------

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   89
                                      ! i3 =   73
                                      ! i4 =   85

                                      xin(89) = c10*xin(73) + xc00*xin(85)
                                      yin(89) = c10*yin(73) + yc00*yin(85)
                                      zin(89) = c10*zin(73) + zc00*zin(85)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   90
                                      ! i5 =   89
                                      ! i4 =   85

                                      xin(90) = xcp00*xin(89) + cp10*xin(85)
                                      yin(90) = ycp00*yin(89) + cp10*yin(85)
                                      zin(90) = zcp00*zin(89) + cp10*zin(85)

                                      ! ------------------

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   89

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   93
                                      ! i3 =   85
                                      ! i4 =   89

                                      xin(93) = c10*xin(85) + xc00*xin(89)
                                      yin(93) = c10*yin(85) + yc00*yin(89)
                                      zin(93) = c10*zin(85) + zc00*zin(89)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   94
                                      ! i5 =   93
                                      ! i4 =   89

                                      xin(94) = xcp00*xin(93) + cp10*xin(89)
                                      yin(94) = ycp00*yin(93) + cp10*yin(89)
                                      zin(94) = zcp00*zin(93) + cp10*zin(89)

                                      ! ------------------

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   49
                                      ! i4 = i1+k2 =   50

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   51
                                      ! i3 =   49
                                      ! i4 =   50

                                      xin(51) = cp01*xin(49) + xcp00*xin(50)
                                      yin(51) = cp01*yin(49) + ycp00*yin(50)
                                      zin(51) = cp01*zin(49) + zcp00*zin(50)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   63

                                      xin(63) = xc00*xin(51) + c01*xin(50)
                                      yin(63) = yc00*yin(51) + c01*yin(50)
                                      zin(63) = zc00*zin(51) + c01*zin(50)

                                      ! ------------------

                                      ! i3 = i4 =   50
                                      ! i4 = i5 =   51

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   52
                                      ! i3 =   50
                                      ! i4 =   51

                                      xin(52) = cp01*xin(50) + xcp00*xin(51)
                                      yin(52) = cp01*yin(50) + ycp00*yin(51)
                                      zin(52) = cp01*zin(50) + zcp00*zin(51)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   64

                                      xin(64) = xc00*xin(52) + c01*xin(51)
                                      yin(64) = yc00*yin(52) + c01*yin(51)
                                      zin(64) = zc00*zin(52) + c01*zin(51)

                                      ! ------------------

                                      ! i3 = i4 =   51
                                      ! i4 = i5 =   52

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   73

                                      xin(75) = c10*xin(51) + xc00*xin(63) + c01*xin(62)
                                      yin(75) = c10*yin(51) + yc00*yin(63) + c01*yin(62)
                                      zin(75) = c10*zin(51) + zc00*zin(63) + c01*zin(62)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(87) = c10*xin(63) + xc00*xin(75) + c01*xin(74)
                                      yin(87) = c10*yin(63) + yc00*yin(75) + c01*yin(74)
                                      zin(87) = c10*zin(63) + zc00*zin(75) + c01*zin(74)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   89

                                      xin(91) = c10*xin(75) + xc00*xin(87) + c01*xin(86)
                                      yin(91) = c10*yin(75) + yc00*yin(87) + c01*yin(86)
                                      zin(91) = c10*zin(75) + zc00*zin(87) + c01*zin(86)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   89

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   93

                                      xin(95) = c10*xin(87) + xc00*xin(91) + c01*xin(90)
                                      yin(95) = c10*yin(87) + yc00*yin(91) + c01*yin(90)
                                      zin(95) = c10*zin(87) + zc00*zin(91) + c01*zin(90)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   49
                                      ! i4 = i2 =   61

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   73

                                      xin(76) = c10*xin(52) + xc00*xin(64) + c01*xin(63)
                                      yin(76) = c10*yin(52) + yc00*yin(64) + c01*yin(63)
                                      zin(76) = c10*zin(52) + zc00*zin(64) + c01*zin(63)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   61
                                      ! i4 = i5 =   73

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   85

                                      xin(88) = c10*xin(64) + xc00*xin(76) + c01*xin(75)
                                      yin(88) = c10*yin(64) + yc00*yin(76) + c01*yin(75)
                                      zin(88) = c10*zin(64) + zc00*zin(76) + c01*zin(75)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   73
                                      ! i4 = i5 =   85

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   89

                                      xin(92) = c10*xin(76) + xc00*xin(88) + c01*xin(87)
                                      yin(92) = c10*yin(76) + yc00*yin(88) + c01*yin(87)
                                      zin(92) = c10*zin(76) + zc00*zin(88) + c01*zin(87)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   85
                                      ! i4 = i5 =   89

                                      ! nn =    5

                                      ! i5 = in(nn+1) =   93

                                      xin(96) = c10*xin(88) + xc00*xin(92) + c01*xin(91)
                                      yin(96) = c10*yin(88) + yc00*yin(92) + c01*yin(91)
                                      zin(96) = c10*zin(88) + zc00*zin(92) + c01*zin(91)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   89
                                      ! i4 = i5 =   93

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =   93

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   85

                                      xin(89) = xin(89) + dxij*xin(85)
                                      yin(89) = yin(89) + dyij*yin(85)
                                      zin(89) = zin(89) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   93

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   89

                                      xin(93) = xin(93) + dxij*xin(89)
                                      yin(93) = yin(93) + dyij*yin(89)
                                      zin(93) = zin(93) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   53

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   53

                                      ! do ni = 1,    3

                                      xin(53) = xin(61) + dxij*xin(49)
                                      yin(53) = yin(61) + dyij*yin(49)
                                      zin(53) = zin(61) + dzij*zin(49)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   65

                                      ! ni =    2

                                      xin(65) = xin(73) + dxij*xin(61)
                                      yin(65) = yin(73) + dyij*yin(61)
                                      zin(65) = zin(73) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   77

                                      ! ni =    3

                                      xin(77) = xin(85) + dxij*xin(73)
                                      yin(77) = yin(85) + dyij*yin(73)
                                      zin(77) = zin(85) + dzij*zin(73)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   89

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   57

                                      ! nj =    2

                                      ! i4 = i3 =   57

                                      ! do ni = 1,    3

                                      xin(57) = xin(65) + dxij*xin(53)
                                      yin(57) = yin(65) + dyij*yin(53)
                                      zin(57) = zin(65) + dzij*zin(53)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   69

                                      ! ni =    2

                                      xin(69) = xin(77) + dxij*xin(65)
                                      yin(69) = yin(77) + dyij*yin(65)
                                      zin(69) = zin(77) + dzij*zin(65)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(89) + dxij*xin(77)
                                      yin(81) = yin(89) + dyij*yin(77)
                                      zin(81) = zin(89) + dzij*zin(77)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   93

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   61

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(94) = xin(94) + dxij*xin(90)
                                      yin(94) = yin(94) + dyij*yin(90)
                                      zin(94) = zin(94) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   86

                                      xin(90) = xin(90) + dxij*xin(86)
                                      yin(90) = yin(90) + dyij*yin(86)
                                      zin(90) = zin(90) + dzij*zin(86)

                                      ! i3 = i4 =   86
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   94

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   90

                                      xin(94) = xin(94) + dxij*xin(90)
                                      yin(94) = yin(94) + dyij*yin(90)
                                      zin(94) = zin(94) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   54

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   54

                                      ! do ni = 1,    3

                                      xin(54) = xin(62) + dxij*xin(50)
                                      yin(54) = yin(62) + dyij*yin(50)
                                      zin(54) = zin(62) + dzij*zin(50)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   66

                                      ! ni =    2

                                      xin(66) = xin(74) + dxij*xin(62)
                                      yin(66) = yin(74) + dyij*yin(62)
                                      zin(66) = zin(74) + dzij*zin(62)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   78

                                      ! ni =    3

                                      xin(78) = xin(86) + dxij*xin(74)
                                      yin(78) = yin(86) + dyij*yin(74)
                                      zin(78) = zin(86) + dzij*zin(74)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   90

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   58

                                      ! nj =    2

                                      ! i4 = i3 =   58

                                      ! do ni = 1,    3

                                      xin(58) = xin(66) + dxij*xin(54)
                                      yin(58) = yin(66) + dyij*yin(54)
                                      zin(58) = zin(66) + dzij*zin(54)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    2

                                      xin(70) = xin(78) + dxij*xin(66)
                                      yin(70) = yin(78) + dyij*yin(66)
                                      zin(70) = zin(78) + dzij*zin(66)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   82

                                      ! ni =    3

                                      xin(82) = xin(90) + dxij*xin(78)
                                      yin(82) = yin(90) + dyij*yin(78)
                                      zin(82) = zin(90) + dzij*zin(78)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   62

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   87

                                      xin(91) = xin(91) + dxij*xin(87)
                                      yin(91) = yin(91) + dyij*yin(87)
                                      zin(91) = zin(91) + dzij*zin(87)

                                      ! i3 = i4 =   87
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   95

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(95) = xin(95) + dxij*xin(91)
                                      yin(95) = yin(95) + dyij*yin(91)
                                      zin(95) = zin(95) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   55

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   55

                                      ! do ni = 1,    3

                                      xin(55) = xin(63) + dxij*xin(51)
                                      yin(55) = yin(63) + dyij*yin(51)
                                      zin(55) = zin(63) + dzij*zin(51)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    2

                                      xin(67) = xin(75) + dxij*xin(63)
                                      yin(67) = yin(75) + dyij*yin(63)
                                      zin(67) = zin(75) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(87) + dxij*xin(75)
                                      yin(79) = yin(87) + dyij*yin(75)
                                      zin(79) = zin(87) + dzij*zin(75)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   59

                                      ! nj =    2

                                      ! i4 = i3 =   59

                                      ! do ni = 1,    3

                                      xin(59) = xin(67) + dxij*xin(55)
                                      yin(59) = yin(67) + dyij*yin(55)
                                      zin(59) = zin(67) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    2

                                      xin(71) = xin(79) + dxij*xin(67)
                                      yin(71) = yin(79) + dyij*yin(67)
                                      zin(71) = zin(79) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   83

                                      ! ni =    3

                                      xin(83) = xin(91) + dxij*xin(79)
                                      yin(83) = yin(91) + dyij*yin(79)
                                      zin(83) = zin(91) + dzij*zin(79)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   63

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   88

                                      xin(92) = xin(92) + dxij*xin(88)
                                      yin(92) = yin(92) + dyij*yin(88)
                                      zin(92) = zin(92) + dzij*zin(88)

                                      ! i3 = i4 =   88
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =   96

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   92

                                      xin(96) = xin(96) + dxij*xin(92)
                                      yin(96) = yin(96) + dyij*yin(92)
                                      zin(96) = zin(96) + dzij*zin(92)

                                      ! i3 = i4 =   92
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   56

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   56

                                      ! do ni = 1,    3

                                      xin(56) = xin(64) + dxij*xin(52)
                                      yin(56) = yin(64) + dyij*yin(52)
                                      zin(56) = zin(64) + dzij*zin(52)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   68

                                      ! ni =    2

                                      xin(68) = xin(76) + dxij*xin(64)
                                      yin(68) = yin(76) + dyij*yin(64)
                                      zin(68) = zin(76) + dzij*zin(64)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    3

                                      xin(80) = xin(88) + dxij*xin(76)
                                      yin(80) = yin(88) + dyij*yin(76)
                                      zin(80) = zin(88) + dzij*zin(76)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   92

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   60

                                      ! nj =    2

                                      ! i4 = i3 =   60

                                      ! do ni = 1,    3

                                      xin(60) = xin(68) + dxij*xin(56)
                                      yin(60) = yin(68) + dyij*yin(56)
                                      zin(60) = zin(68) + dzij*zin(56)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    2

                                      xin(72) = xin(80) + dxij*xin(68)
                                      yin(72) = yin(80) + dyij*yin(68)
                                      zin(72) = zin(80) + dzij*zin(68)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   84

                                      ! ni =    3

                                      xin(84) = xin(92) + dxij*xin(80)
                                      yin(84) = yin(92) + dyij*yin(80)
                                      zin(84) = zin(92) + dzij*zin(80)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   64

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =   96

                                      u2 = roots(3)*rho
                                      f00 = expe*wghts(3)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =   97

                                      xin(97) = 1.0_dp
                                      yin(97) = 1.0_dp
                                      zin(97) = f00

                                      ! i2 = in(2) =  109
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(109) = xc00
                                      yin(109) = yc00
                                      zin(109) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =   98

                                      xin(98) = xcp00
                                      yin(98) = ycp00
                                      zin(98) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  110
                                      ! i2 =  109

                                      xin(110) = xcp00*xin(109) + cp10
                                      yin(110) = ycp00*yin(109) + cp10
                                      zin(110) = zcp00*zin(109) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  121
                                      ! i3 =   97
                                      ! i4 =  109

                                      xin(121) = c10*xin(97) + xc00*xin(109)
                                      yin(121) = c10*yin(97) + yc00*yin(109)
                                      zin(121) = c10*zin(97) + zc00*zin(109)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  122
                                      ! i5 =  121
                                      ! i4 =  109

                                      xin(122) = xcp00*xin(121) + cp10*xin(109)
                                      yin(122) = ycp00*yin(121) + cp10*yin(109)
                                      zin(122) = zcp00*zin(121) + cp10*zin(109)

                                      ! ------------------

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  133
                                      ! i3 =  109
                                      ! i4 =  121

                                      xin(133) = c10*xin(109) + xc00*xin(121)
                                      yin(133) = c10*yin(109) + yc00*yin(121)
                                      zin(133) = c10*zin(109) + zc00*zin(121)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  134
                                      ! i5 =  133
                                      ! i4 =  121

                                      xin(134) = xcp00*xin(133) + cp10*xin(121)
                                      yin(134) = ycp00*yin(133) + cp10*yin(121)
                                      zin(134) = zcp00*zin(133) + cp10*zin(121)

                                      ! ------------------

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  137
                                      ! i3 =  121
                                      ! i4 =  133

                                      xin(137) = c10*xin(121) + xc00*xin(133)
                                      yin(137) = c10*yin(121) + yc00*yin(133)
                                      zin(137) = c10*zin(121) + zc00*zin(133)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  138
                                      ! i5 =  137
                                      ! i4 =  133

                                      xin(138) = xcp00*xin(137) + cp10*xin(133)
                                      yin(138) = ycp00*yin(137) + cp10*yin(133)
                                      zin(138) = zcp00*zin(137) + cp10*zin(133)

                                      ! ------------------

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  141
                                      ! i3 =  133
                                      ! i4 =  137

                                      xin(141) = c10*xin(133) + xc00*xin(137)
                                      yin(141) = c10*yin(133) + yc00*yin(137)
                                      zin(141) = c10*zin(133) + zc00*zin(137)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  142
                                      ! i5 =  141
                                      ! i4 =  137

                                      xin(142) = xcp00*xin(141) + cp10*xin(137)
                                      yin(142) = ycp00*yin(141) + cp10*yin(137)
                                      zin(142) = zcp00*zin(141) + cp10*zin(137)

                                      ! ------------------

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =   97
                                      ! i4 = i1+k2 =   98

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =   99
                                      ! i3 =   97
                                      ! i4 =   98

                                      xin(99) = cp01*xin(97) + xcp00*xin(98)
                                      yin(99) = cp01*yin(97) + ycp00*yin(98)
                                      zin(99) = cp01*zin(97) + zcp00*zin(98)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  111

                                      xin(111) = xc00*xin(99) + c01*xin(98)
                                      yin(111) = yc00*yin(99) + c01*yin(98)
                                      zin(111) = zc00*zin(99) + c01*zin(98)

                                      ! ------------------

                                      ! i3 = i4 =   98
                                      ! i4 = i5 =   99

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  100
                                      ! i3 =   98
                                      ! i4 =   99

                                      xin(100) = cp01*xin(98) + xcp00*xin(99)
                                      yin(100) = cp01*yin(98) + ycp00*yin(99)
                                      zin(100) = cp01*zin(98) + zcp00*zin(99)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  112

                                      xin(112) = xc00*xin(100) + c01*xin(99)
                                      yin(112) = yc00*yin(100) + c01*yin(99)
                                      zin(112) = zc00*zin(100) + c01*zin(99)

                                      ! ------------------

                                      ! i3 = i4 =   99
                                      ! i4 = i5 =  100

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  121

                                      xin(123) = c10*xin(99) + xc00*xin(111) + c01*xin(110)
                                      yin(123) = c10*yin(99) + yc00*yin(111) + c01*yin(110)
                                      zin(123) = c10*zin(99) + zc00*zin(111) + c01*zin(110)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(135) = c10*xin(111) + xc00*xin(123) + c01*xin(122)
                                      yin(135) = c10*yin(111) + yc00*yin(123) + c01*yin(122)
                                      zin(135) = c10*zin(111) + zc00*zin(123) + c01*zin(122)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(139) = c10*xin(123) + xc00*xin(135) + c01*xin(134)
                                      yin(139) = c10*yin(123) + yc00*yin(135) + c01*yin(134)
                                      zin(139) = c10*zin(123) + zc00*zin(135) + c01*zin(134)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  141

                                      xin(143) = c10*xin(135) + xc00*xin(139) + c01*xin(138)
                                      yin(143) = c10*yin(135) + yc00*yin(139) + c01*yin(138)
                                      zin(143) = c10*zin(135) + zc00*zin(139) + c01*zin(138)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =   97
                                      ! i4 = i2 =  109

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  121

                                      xin(124) = c10*xin(100) + xc00*xin(112) + c01*xin(111)
                                      yin(124) = c10*yin(100) + yc00*yin(112) + c01*yin(111)
                                      zin(124) = c10*zin(100) + zc00*zin(112) + c01*zin(111)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  109
                                      ! i4 = i5 =  121

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  133

                                      xin(136) = c10*xin(112) + xc00*xin(124) + c01*xin(123)
                                      yin(136) = c10*yin(112) + yc00*yin(124) + c01*yin(123)
                                      zin(136) = c10*zin(112) + zc00*zin(124) + c01*zin(123)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  121
                                      ! i4 = i5 =  133

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  137

                                      xin(140) = c10*xin(124) + xc00*xin(136) + c01*xin(135)
                                      yin(140) = c10*yin(124) + yc00*yin(136) + c01*yin(135)
                                      zin(140) = c10*zin(124) + zc00*zin(136) + c01*zin(135)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  133
                                      ! i4 = i5 =  137

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  141

                                      xin(144) = c10*xin(136) + xc00*xin(140) + c01*xin(139)
                                      yin(144) = c10*yin(136) + yc00*yin(140) + c01*yin(139)
                                      zin(144) = c10*zin(136) + zc00*zin(140) + c01*zin(139)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  137
                                      ! i4 = i5 =  141

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  141

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  133

                                      xin(137) = xin(137) + dxij*xin(133)
                                      yin(137) = yin(137) + dyij*yin(133)
                                      zin(137) = zin(137) + dzij*zin(133)

                                      ! i3 = i4 =  133
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  141

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  137

                                      xin(141) = xin(141) + dxij*xin(137)
                                      yin(141) = yin(141) + dyij*yin(137)
                                      zin(141) = zin(141) + dzij*zin(137)

                                      ! i3 = i4 =  137
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  101

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  101

                                      ! do ni = 1,    3

                                      xin(101) = xin(109) + dxij*xin(97)
                                      yin(101) = yin(109) + dyij*yin(97)
                                      zin(101) = zin(109) + dzij*zin(97)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  113

                                      ! ni =    2

                                      xin(113) = xin(121) + dxij*xin(109)
                                      yin(113) = yin(121) + dyij*yin(109)
                                      zin(113) = zin(121) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  125

                                      ! ni =    3

                                      xin(125) = xin(133) + dxij*xin(121)
                                      yin(125) = yin(133) + dyij*yin(121)
                                      zin(125) = zin(133) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  137

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  105

                                      ! nj =    2

                                      ! i4 = i3 =  105

                                      ! do ni = 1,    3

                                      xin(105) = xin(113) + dxij*xin(101)
                                      yin(105) = yin(113) + dyij*yin(101)
                                      zin(105) = zin(113) + dzij*zin(101)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  117

                                      ! ni =    2

                                      xin(117) = xin(125) + dxij*xin(113)
                                      yin(117) = yin(125) + dyij*yin(113)
                                      zin(117) = zin(125) + dzij*zin(113)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  129

                                      ! ni =    3

                                      xin(129) = xin(137) + dxij*xin(125)
                                      yin(129) = yin(137) + dyij*yin(125)
                                      zin(129) = zin(137) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  141

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  109

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(142) = xin(142) + dxij*xin(138)
                                      yin(142) = yin(142) + dyij*yin(138)
                                      zin(142) = zin(142) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  134

                                      xin(138) = xin(138) + dxij*xin(134)
                                      yin(138) = yin(138) + dyij*yin(134)
                                      zin(138) = zin(138) + dzij*zin(134)

                                      ! i3 = i4 =  134
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  142

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  138

                                      xin(142) = xin(142) + dxij*xin(138)
                                      yin(142) = yin(142) + dyij*yin(138)
                                      zin(142) = zin(142) + dzij*zin(138)

                                      ! i3 = i4 =  138
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  102

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  102

                                      ! do ni = 1,    3

                                      xin(102) = xin(110) + dxij*xin(98)
                                      yin(102) = yin(110) + dyij*yin(98)
                                      zin(102) = zin(110) + dzij*zin(98)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  114

                                      ! ni =    2

                                      xin(114) = xin(122) + dxij*xin(110)
                                      yin(114) = yin(122) + dyij*yin(110)
                                      zin(114) = zin(122) + dzij*zin(110)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  126

                                      ! ni =    3

                                      xin(126) = xin(134) + dxij*xin(122)
                                      yin(126) = yin(134) + dyij*yin(122)
                                      zin(126) = zin(134) + dzij*zin(122)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  138

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  106

                                      ! nj =    2

                                      ! i4 = i3 =  106

                                      ! do ni = 1,    3

                                      xin(106) = xin(114) + dxij*xin(102)
                                      yin(106) = yin(114) + dyij*yin(102)
                                      zin(106) = zin(114) + dzij*zin(102)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  118

                                      ! ni =    2

                                      xin(118) = xin(126) + dxij*xin(114)
                                      yin(118) = yin(126) + dyij*yin(114)
                                      zin(118) = zin(126) + dzij*zin(114)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                      ! ni =    3

                                      xin(130) = xin(138) + dxij*xin(126)
                                      yin(130) = yin(138) + dyij*yin(126)
                                      zin(130) = zin(138) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  110

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  135

                                      xin(139) = xin(139) + dxij*xin(135)
                                      yin(139) = yin(139) + dyij*yin(135)
                                      zin(139) = zin(139) + dzij*zin(135)

                                      ! i3 = i4 =  135
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  143

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  139

                                      xin(143) = xin(143) + dxij*xin(139)
                                      yin(143) = yin(143) + dyij*yin(139)
                                      zin(143) = zin(143) + dzij*zin(139)

                                      ! i3 = i4 =  139
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  103

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  103

                                      ! do ni = 1,    3

                                      xin(103) = xin(111) + dxij*xin(99)
                                      yin(103) = yin(111) + dyij*yin(99)
                                      zin(103) = zin(111) + dzij*zin(99)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  115

                                      ! ni =    2

                                      xin(115) = xin(123) + dxij*xin(111)
                                      yin(115) = yin(123) + dyij*yin(111)
                                      zin(115) = zin(123) + dzij*zin(111)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                      ! ni =    3

                                      xin(127) = xin(135) + dxij*xin(123)
                                      yin(127) = yin(135) + dyij*yin(123)
                                      zin(127) = zin(135) + dzij*zin(123)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  107

                                      ! nj =    2

                                      ! i4 = i3 =  107

                                      ! do ni = 1,    3

                                      xin(107) = xin(115) + dxij*xin(103)
                                      yin(107) = yin(115) + dyij*yin(103)
                                      zin(107) = zin(115) + dzij*zin(103)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  119

                                      ! ni =    2

                                      xin(119) = xin(127) + dxij*xin(115)
                                      yin(119) = yin(127) + dyij*yin(115)
                                      zin(119) = zin(127) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                      ! ni =    3

                                      xin(131) = xin(139) + dxij*xin(127)
                                      yin(131) = yin(139) + dyij*yin(127)
                                      zin(131) = zin(139) + dzij*zin(127)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  111

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  136

                                      xin(140) = xin(140) + dxij*xin(136)
                                      yin(140) = yin(140) + dyij*yin(136)
                                      zin(140) = zin(140) + dzij*zin(136)

                                      ! i3 = i4 =  136
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  144

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  140

                                      xin(144) = xin(144) + dxij*xin(140)
                                      yin(144) = yin(144) + dyij*yin(140)
                                      zin(144) = zin(144) + dzij*zin(140)

                                      ! i3 = i4 =  140
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  104

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  104

                                      ! do ni = 1,    3

                                      xin(104) = xin(112) + dxij*xin(100)
                                      yin(104) = yin(112) + dyij*yin(100)
                                      zin(104) = zin(112) + dzij*zin(100)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  116

                                      ! ni =    2

                                      xin(116) = xin(124) + dxij*xin(112)
                                      yin(116) = yin(124) + dyij*yin(112)
                                      zin(116) = zin(124) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  128

                                      ! ni =    3

                                      xin(128) = xin(136) + dxij*xin(124)
                                      yin(128) = yin(136) + dyij*yin(124)
                                      zin(128) = zin(136) + dzij*zin(124)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  140

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  108

                                      ! nj =    2

                                      ! i4 = i3 =  108

                                      ! do ni = 1,    3

                                      xin(108) = xin(116) + dxij*xin(104)
                                      yin(108) = yin(116) + dyij*yin(104)
                                      zin(108) = zin(116) + dzij*zin(104)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  120

                                      ! ni =    2

                                      xin(120) = xin(128) + dxij*xin(116)
                                      yin(120) = yin(128) + dyij*yin(116)
                                      zin(120) = zin(128) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                      ! ni =    3

                                      xin(132) = xin(140) + dxij*xin(128)
                                      yin(132) = yin(140) + dyij*yin(128)
                                      zin(132) = zin(140) + dzij*zin(128)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  112

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

                                      u2 = roots(4)*rho
                                      f00 = expe*wghts(4)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  145

                                      xin(145) = 1.0_dp
                                      yin(145) = 1.0_dp
                                      zin(145) = f00

                                      ! i2 = in(2) =  157
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(157) = xc00
                                      yin(157) = yc00
                                      zin(157) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  146

                                      xin(146) = xcp00
                                      yin(146) = ycp00
                                      zin(146) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  158
                                      ! i2 =  157

                                      xin(158) = xcp00*xin(157) + cp10
                                      yin(158) = ycp00*yin(157) + cp10
                                      zin(158) = zcp00*zin(157) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  169
                                      ! i3 =  145
                                      ! i4 =  157

                                      xin(169) = c10*xin(145) + xc00*xin(157)
                                      yin(169) = c10*yin(145) + yc00*yin(157)
                                      zin(169) = c10*zin(145) + zc00*zin(157)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  170
                                      ! i5 =  169
                                      ! i4 =  157

                                      xin(170) = xcp00*xin(169) + cp10*xin(157)
                                      yin(170) = ycp00*yin(169) + cp10*yin(157)
                                      zin(170) = zcp00*zin(169) + cp10*zin(157)

                                      ! ------------------

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  181
                                      ! i3 =  157
                                      ! i4 =  169

                                      xin(181) = c10*xin(157) + xc00*xin(169)
                                      yin(181) = c10*yin(157) + yc00*yin(169)
                                      zin(181) = c10*zin(157) + zc00*zin(169)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  182
                                      ! i5 =  181
                                      ! i4 =  169

                                      xin(182) = xcp00*xin(181) + cp10*xin(169)
                                      yin(182) = ycp00*yin(181) + cp10*yin(169)
                                      zin(182) = zcp00*zin(181) + cp10*zin(169)

                                      ! ------------------

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  185
                                      ! i3 =  169
                                      ! i4 =  181

                                      xin(185) = c10*xin(169) + xc00*xin(181)
                                      yin(185) = c10*yin(169) + yc00*yin(181)
                                      zin(185) = c10*zin(169) + zc00*zin(181)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  186
                                      ! i5 =  185
                                      ! i4 =  181

                                      xin(186) = xcp00*xin(185) + cp10*xin(181)
                                      yin(186) = ycp00*yin(185) + cp10*yin(181)
                                      zin(186) = zcp00*zin(185) + cp10*zin(181)

                                      ! ------------------

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  189
                                      ! i3 =  181
                                      ! i4 =  185

                                      xin(189) = c10*xin(181) + xc00*xin(185)
                                      yin(189) = c10*yin(181) + yc00*yin(185)
                                      zin(189) = c10*zin(181) + zc00*zin(185)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  190
                                      ! i5 =  189
                                      ! i4 =  185

                                      xin(190) = xcp00*xin(189) + cp10*xin(185)
                                      yin(190) = ycp00*yin(189) + cp10*yin(185)
                                      zin(190) = zcp00*zin(189) + cp10*zin(185)

                                      ! ------------------

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  145
                                      ! i4 = i1+k2 =  146

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  147
                                      ! i3 =  145
                                      ! i4 =  146

                                      xin(147) = cp01*xin(145) + xcp00*xin(146)
                                      yin(147) = cp01*yin(145) + ycp00*yin(146)
                                      zin(147) = cp01*zin(145) + zcp00*zin(146)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  159

                                      xin(159) = xc00*xin(147) + c01*xin(146)
                                      yin(159) = yc00*yin(147) + c01*yin(146)
                                      zin(159) = zc00*zin(147) + c01*zin(146)

                                      ! ------------------

                                      ! i3 = i4 =  146
                                      ! i4 = i5 =  147

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  148
                                      ! i3 =  146
                                      ! i4 =  147

                                      xin(148) = cp01*xin(146) + xcp00*xin(147)
                                      yin(148) = cp01*yin(146) + ycp00*yin(147)
                                      zin(148) = cp01*zin(146) + zcp00*zin(147)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  160

                                      xin(160) = xc00*xin(148) + c01*xin(147)
                                      yin(160) = yc00*yin(148) + c01*yin(147)
                                      zin(160) = zc00*zin(148) + c01*zin(147)

                                      ! ------------------

                                      ! i3 = i4 =  147
                                      ! i4 = i5 =  148

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  169

                                      xin(171) = c10*xin(147) + xc00*xin(159) + c01*xin(158)
                                      yin(171) = c10*yin(147) + yc00*yin(159) + c01*yin(158)
                                      zin(171) = c10*zin(147) + zc00*zin(159) + c01*zin(158)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(183) = c10*xin(159) + xc00*xin(171) + c01*xin(170)
                                      yin(183) = c10*yin(159) + yc00*yin(171) + c01*yin(170)
                                      zin(183) = c10*zin(159) + zc00*zin(171) + c01*zin(170)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  185

                                      xin(187) = c10*xin(171) + xc00*xin(183) + c01*xin(182)
                                      yin(187) = c10*yin(171) + yc00*yin(183) + c01*yin(182)
                                      zin(187) = c10*zin(171) + zc00*zin(183) + c01*zin(182)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  189

                                      xin(191) = c10*xin(183) + xc00*xin(187) + c01*xin(186)
                                      yin(191) = c10*yin(183) + yc00*yin(187) + c01*yin(186)
                                      zin(191) = c10*zin(183) + zc00*zin(187) + c01*zin(186)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  145
                                      ! i4 = i2 =  157

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  169

                                      xin(172) = c10*xin(148) + xc00*xin(160) + c01*xin(159)
                                      yin(172) = c10*yin(148) + yc00*yin(160) + c01*yin(159)
                                      zin(172) = c10*zin(148) + zc00*zin(160) + c01*zin(159)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  157
                                      ! i4 = i5 =  169

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  181

                                      xin(184) = c10*xin(160) + xc00*xin(172) + c01*xin(171)
                                      yin(184) = c10*yin(160) + yc00*yin(172) + c01*yin(171)
                                      zin(184) = c10*zin(160) + zc00*zin(172) + c01*zin(171)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  169
                                      ! i4 = i5 =  181

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  185

                                      xin(188) = c10*xin(172) + xc00*xin(184) + c01*xin(183)
                                      yin(188) = c10*yin(172) + yc00*yin(184) + c01*yin(183)
                                      zin(188) = c10*zin(172) + zc00*zin(184) + c01*zin(183)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  181
                                      ! i4 = i5 =  185

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  189

                                      xin(192) = c10*xin(184) + xc00*xin(188) + c01*xin(187)
                                      yin(192) = c10*yin(184) + yc00*yin(188) + c01*yin(187)
                                      zin(192) = c10*zin(184) + zc00*zin(188) + c01*zin(187)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  185
                                      ! i4 = i5 =  189

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  189

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  181

                                      xin(185) = xin(185) + dxij*xin(181)
                                      yin(185) = yin(185) + dyij*yin(181)
                                      zin(185) = zin(185) + dzij*zin(181)

                                      ! i3 = i4 =  181
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  189

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  185

                                      xin(189) = xin(189) + dxij*xin(185)
                                      yin(189) = yin(189) + dyij*yin(185)
                                      zin(189) = zin(189) + dzij*zin(185)

                                      ! i3 = i4 =  185
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  149

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  149

                                      ! do ni = 1,    3

                                      xin(149) = xin(157) + dxij*xin(145)
                                      yin(149) = yin(157) + dyij*yin(145)
                                      zin(149) = zin(157) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    2

                                      xin(161) = xin(169) + dxij*xin(157)
                                      yin(161) = yin(169) + dyij*yin(157)
                                      zin(161) = zin(169) + dzij*zin(157)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  173

                                      ! ni =    3

                                      xin(173) = xin(181) + dxij*xin(169)
                                      yin(173) = yin(181) + dyij*yin(169)
                                      zin(173) = zin(181) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  185

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  153

                                      ! nj =    2

                                      ! i4 = i3 =  153

                                      ! do ni = 1,    3

                                      xin(153) = xin(161) + dxij*xin(149)
                                      yin(153) = yin(161) + dyij*yin(149)
                                      zin(153) = zin(161) + dzij*zin(149)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  165

                                      ! ni =    2

                                      xin(165) = xin(173) + dxij*xin(161)
                                      yin(165) = yin(173) + dyij*yin(161)
                                      zin(165) = zin(173) + dzij*zin(161)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  177

                                      ! ni =    3

                                      xin(177) = xin(185) + dxij*xin(173)
                                      yin(177) = yin(185) + dyij*yin(173)
                                      zin(177) = zin(185) + dzij*zin(173)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  157

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(190) = xin(190) + dxij*xin(186)
                                      yin(190) = yin(190) + dyij*yin(186)
                                      zin(190) = zin(190) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  182

                                      xin(186) = xin(186) + dxij*xin(182)
                                      yin(186) = yin(186) + dyij*yin(182)
                                      zin(186) = zin(186) + dzij*zin(182)

                                      ! i3 = i4 =  182
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  190

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  186

                                      xin(190) = xin(190) + dxij*xin(186)
                                      yin(190) = yin(190) + dyij*yin(186)
                                      zin(190) = zin(190) + dzij*zin(186)

                                      ! i3 = i4 =  186
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  150

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  150

                                      ! do ni = 1,    3

                                      xin(150) = xin(158) + dxij*xin(146)
                                      yin(150) = yin(158) + dyij*yin(146)
                                      zin(150) = zin(158) + dzij*zin(146)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  162

                                      ! ni =    2

                                      xin(162) = xin(170) + dxij*xin(158)
                                      yin(162) = yin(170) + dyij*yin(158)
                                      zin(162) = zin(170) + dzij*zin(158)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  174

                                      ! ni =    3

                                      xin(174) = xin(182) + dxij*xin(170)
                                      yin(174) = yin(182) + dyij*yin(170)
                                      zin(174) = zin(182) + dzij*zin(170)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  186

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  154

                                      ! nj =    2

                                      ! i4 = i3 =  154

                                      ! do ni = 1,    3

                                      xin(154) = xin(162) + dxij*xin(150)
                                      yin(154) = yin(162) + dyij*yin(150)
                                      zin(154) = zin(162) + dzij*zin(150)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  166

                                      ! ni =    2

                                      xin(166) = xin(174) + dxij*xin(162)
                                      yin(166) = yin(174) + dyij*yin(162)
                                      zin(166) = zin(174) + dzij*zin(162)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  178

                                      ! ni =    3

                                      xin(178) = xin(186) + dxij*xin(174)
                                      yin(178) = yin(186) + dyij*yin(174)
                                      zin(178) = zin(186) + dzij*zin(174)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  190

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  158

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  183

                                      xin(187) = xin(187) + dxij*xin(183)
                                      yin(187) = yin(187) + dyij*yin(183)
                                      zin(187) = zin(187) + dzij*zin(183)

                                      ! i3 = i4 =  183
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  191

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  187

                                      xin(191) = xin(191) + dxij*xin(187)
                                      yin(191) = yin(191) + dyij*yin(187)
                                      zin(191) = zin(191) + dzij*zin(187)

                                      ! i3 = i4 =  187
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  151

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  151

                                      ! do ni = 1,    3

                                      xin(151) = xin(159) + dxij*xin(147)
                                      yin(151) = yin(159) + dyij*yin(147)
                                      zin(151) = zin(159) + dzij*zin(147)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! ni =    2

                                      xin(163) = xin(171) + dxij*xin(159)
                                      yin(163) = yin(171) + dyij*yin(159)
                                      zin(163) = zin(171) + dzij*zin(159)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(183) + dxij*xin(171)
                                      yin(175) = yin(183) + dyij*yin(171)
                                      zin(175) = zin(183) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  155

                                      ! nj =    2

                                      ! i4 = i3 =  155

                                      ! do ni = 1,    3

                                      xin(155) = xin(163) + dxij*xin(151)
                                      yin(155) = yin(163) + dyij*yin(151)
                                      zin(155) = zin(163) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  167

                                      ! ni =    2

                                      xin(167) = xin(175) + dxij*xin(163)
                                      yin(167) = yin(175) + dyij*yin(163)
                                      zin(167) = zin(175) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    3

                                      xin(179) = xin(187) + dxij*xin(175)
                                      yin(179) = yin(187) + dyij*yin(175)
                                      zin(179) = zin(187) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  191

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  159

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  184

                                      xin(188) = xin(188) + dxij*xin(184)
                                      yin(188) = yin(188) + dyij*yin(184)
                                      zin(188) = zin(188) + dzij*zin(184)

                                      ! i3 = i4 =  184
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  192

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  188

                                      xin(192) = xin(192) + dxij*xin(188)
                                      yin(192) = yin(192) + dyij*yin(188)
                                      zin(192) = zin(192) + dzij*zin(188)

                                      ! i3 = i4 =  188
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  152

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  152

                                      ! do ni = 1,    3

                                      xin(152) = xin(160) + dxij*xin(148)
                                      yin(152) = yin(160) + dyij*yin(148)
                                      zin(152) = zin(160) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  164

                                      ! ni =    2

                                      xin(164) = xin(172) + dxij*xin(160)
                                      yin(164) = yin(172) + dyij*yin(160)
                                      zin(164) = zin(172) + dzij*zin(160)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  176

                                      ! ni =    3

                                      xin(176) = xin(184) + dxij*xin(172)
                                      yin(176) = yin(184) + dyij*yin(172)
                                      zin(176) = zin(184) + dzij*zin(172)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  188

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  156

                                      ! nj =    2

                                      ! i4 = i3 =  156

                                      ! do ni = 1,    3

                                      xin(156) = xin(164) + dxij*xin(152)
                                      yin(156) = yin(164) + dyij*yin(152)
                                      zin(156) = zin(164) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  168

                                      ! ni =    2

                                      xin(168) = xin(176) + dxij*xin(164)
                                      yin(168) = yin(176) + dyij*yin(164)
                                      zin(168) = zin(176) + dzij*zin(164)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    3

                                      xin(180) = xin(188) + dxij*xin(176)
                                      yin(180) = yin(188) + dyij*yin(176)
                                      zin(180) = zin(188) + dzij*zin(176)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  192

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  160

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  192

                                      u2 = roots(5)*rho
                                      f00 = expe*wghts(5)

                                      ! do iii = 1, nroot
                                      !     in(iii) = in1(iii) + mm ! Indices for current root
                                      ! end do

                                      duminv = 1.0_dp/(t_expon_ab*t_expon_cd + u2*(t_expon_ab + t_expon_cd))
                                      dm2inv = 0.5_dp*duminv
                                      bp01 = (t_expon_ab + u2)*dm2inv
                                      b00 = u2*dm2inv
                                      b10 = (t_expon_cd + u2)*dm2inv
                                      xcp00 = (u2*c1x + c2x)*duminv
                                      xc00 = (u2*c3x + c4x)*duminv
                                      ycp00 = (u2*c1y + c2y)*duminv
                                      yc00 = (u2*c3y + c4y)*duminv
                                      zcp00 = (u2*c1z + c2z)*duminv
                                      zc00 = (u2*c3z + c4z)*duminv

                                      !                       --- XYZINT ---
                                      ! Calculate intermediate 2D integrals with the coefficients from above
                                      !
                                      ! Ix,Iy,Iz are the intermediate 2D integrals in the x,y,z directions
                                      ! they are evaluated via recurrence relations and transfer equations

                                      ! ----- I(0,0) -----

                                      ! i1 = in(1) =  193

                                      xin(193) = 1.0_dp
                                      yin(193) = 1.0_dp
                                      zin(193) = f00

                                      ! i2 = in(2) =  205
                                      ! k2 = kn(2) =    1
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(205) = xc00
                                      yin(205) = yc00
                                      zin(205) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  194

                                      xin(194) = xcp00
                                      yin(194) = ycp00
                                      zin(194) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  206
                                      ! i2 =  205

                                      xin(206) = xcp00*xin(205) + cp10
                                      yin(206) = ycp00*yin(205) + cp10
                                      zin(206) = zcp00*zin(205) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  205

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  217
                                      ! i3 =  193
                                      ! i4 =  205

                                      xin(217) = c10*xin(193) + xc00*xin(205)
                                      yin(217) = c10*yin(193) + yc00*yin(205)
                                      zin(217) = c10*zin(193) + zc00*zin(205)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  218
                                      ! i5 =  217
                                      ! i4 =  205

                                      xin(218) = xcp00*xin(217) + cp10*xin(205)
                                      yin(218) = ycp00*yin(217) + cp10*yin(205)
                                      zin(218) = zcp00*zin(217) + cp10*zin(205)

                                      ! ------------------

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  217

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  229
                                      ! i3 =  205
                                      ! i4 =  217

                                      xin(229) = c10*xin(205) + xc00*xin(217)
                                      yin(229) = c10*yin(205) + yc00*yin(217)
                                      zin(229) = c10*zin(205) + zc00*zin(217)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  230
                                      ! i5 =  229
                                      ! i4 =  217

                                      xin(230) = xcp00*xin(229) + cp10*xin(217)
                                      yin(230) = ycp00*yin(229) + cp10*yin(217)
                                      zin(230) = zcp00*zin(229) + cp10*zin(217)

                                      ! ------------------

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  229

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  233
                                      ! i3 =  217
                                      ! i4 =  229

                                      xin(233) = c10*xin(217) + xc00*xin(229)
                                      yin(233) = c10*yin(217) + yc00*yin(229)
                                      zin(233) = c10*zin(217) + zc00*zin(229)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  234
                                      ! i5 =  233
                                      ! i4 =  229

                                      xin(234) = xcp00*xin(233) + cp10*xin(229)
                                      yin(234) = ycp00*yin(233) + cp10*yin(229)
                                      zin(234) = zcp00*zin(233) + cp10*zin(229)

                                      ! ------------------

                                      ! i3 = i4 =  229
                                      ! i4 = i5 =  233

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  237
                                      ! i3 =  229
                                      ! i4 =  233

                                      xin(237) = c10*xin(229) + xc00*xin(233)
                                      yin(237) = c10*yin(229) + yc00*yin(233)
                                      zin(237) = c10*zin(229) + zc00*zin(233)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  238
                                      ! i5 =  237
                                      ! i4 =  233

                                      xin(238) = xcp00*xin(237) + cp10*xin(233)
                                      yin(238) = ycp00*yin(237) + cp10*yin(233)
                                      zin(238) = zcp00*zin(237) + cp10*zin(233)

                                      ! ------------------

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  237

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  193
                                      ! i4 = i1+k2 =  194

                                      ! do n = 2,    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  195
                                      ! i3 =  193
                                      ! i4 =  194

                                      xin(195) = cp01*xin(193) + xcp00*xin(194)
                                      yin(195) = cp01*yin(193) + ycp00*yin(194)
                                      zin(195) = cp01*zin(193) + zcp00*zin(194)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  207

                                      xin(207) = xc00*xin(195) + c01*xin(194)
                                      yin(207) = yc00*yin(195) + c01*yin(194)
                                      zin(207) = zc00*zin(195) + c01*zin(194)

                                      ! ------------------

                                      ! i3 = i4 =  194
                                      ! i4 = i5 =  195

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  196
                                      ! i3 =  194
                                      ! i4 =  195

                                      xin(196) = cp01*xin(194) + xcp00*xin(195)
                                      yin(196) = cp01*yin(194) + ycp00*yin(195)
                                      zin(196) = cp01*zin(194) + zcp00*zin(195)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  208

                                      xin(208) = xc00*xin(196) + c01*xin(195)
                                      yin(208) = yc00*yin(196) + c01*yin(195)
                                      zin(208) = zc00*zin(196) + c01*zin(195)

                                      ! ------------------

                                      ! i3 = i4 =  195
                                      ! i4 = i5 =  196

                                      ! n =    4

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    1

                                      ! do n = 2,    3

                                      ! k4 = kn(n+1) =    2
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  205

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  217

                                      xin(219) = c10*xin(195) + xc00*xin(207) + c01*xin(206)
                                      yin(219) = c10*yin(195) + yc00*yin(207) + c01*yin(206)
                                      zin(219) = c10*zin(195) + zc00*zin(207) + c01*zin(206)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  229

                                      xin(231) = c10*xin(207) + xc00*xin(219) + c01*xin(218)
                                      yin(231) = c10*yin(207) + yc00*yin(219) + c01*yin(218)
                                      zin(231) = c10*zin(207) + zc00*zin(219) + c01*zin(218)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  229

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  233

                                      xin(235) = c10*xin(219) + xc00*xin(231) + c01*xin(230)
                                      yin(235) = c10*yin(219) + yc00*yin(231) + c01*yin(230)
                                      zin(235) = c10*zin(219) + zc00*zin(231) + c01*zin(230)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  229
                                      ! i4 = i5 =  233

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  237

                                      xin(239) = c10*xin(231) + xc00*xin(235) + c01*xin(234)
                                      yin(239) = c10*yin(231) + yc00*yin(235) + c01*yin(234)
                                      zin(239) = c10*zin(231) + zc00*zin(235) + c01*zin(234)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  237

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   2

                                      ! n =    3

                                      ! k4 = kn(n+1) =    3
                                      ! i3 = i1 =  193
                                      ! i4 = i2 =  205

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  217

                                      xin(220) = c10*xin(196) + xc00*xin(208) + c01*xin(207)
                                      yin(220) = c10*yin(196) + yc00*yin(208) + c01*yin(207)
                                      zin(220) = c10*zin(196) + zc00*zin(208) + c01*zin(207)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  205
                                      ! i4 = i5 =  217

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  229

                                      xin(232) = c10*xin(208) + xc00*xin(220) + c01*xin(219)
                                      yin(232) = c10*yin(208) + yc00*yin(220) + c01*yin(219)
                                      zin(232) = c10*zin(208) + zc00*zin(220) + c01*zin(219)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  217
                                      ! i4 = i5 =  229

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  233

                                      xin(236) = c10*xin(220) + xc00*xin(232) + c01*xin(231)
                                      yin(236) = c10*yin(220) + yc00*yin(232) + c01*yin(231)
                                      zin(236) = c10*zin(220) + zc00*zin(232) + c01*zin(231)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  229
                                      ! i4 = i5 =  233

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  237

                                      xin(240) = c10*xin(232) + xc00*xin(236) + c01*xin(235)
                                      yin(240) = c10*yin(232) + yc00*yin(236) + c01*yin(235)
                                      zin(240) = c10*zin(232) + zc00*zin(236) + c01*zin(235)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  233
                                      ! i4 = i5 =  237

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   3

                                      ! n =    4

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  237

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  237

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  233

                                      xin(237) = xin(237) + dxij*xin(233)
                                      yin(237) = yin(237) + dyij*yin(233)
                                      zin(237) = zin(237) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  229

                                      xin(233) = xin(233) + dxij*xin(229)
                                      yin(233) = yin(233) + dyij*yin(229)
                                      zin(233) = zin(233) + dzij*zin(229)

                                      ! i3 = i4 =  229
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  237

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  233

                                      xin(237) = xin(237) + dxij*xin(233)
                                      yin(237) = yin(237) + dyij*yin(233)
                                      zin(237) = zin(237) + dzij*zin(233)

                                      ! i3 = i4 =  233
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  197

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  197

                                      ! do ni = 1,    3

                                      xin(197) = xin(205) + dxij*xin(193)
                                      yin(197) = yin(205) + dyij*yin(193)
                                      zin(197) = zin(205) + dzij*zin(193)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  209

                                      ! ni =    2

                                      xin(209) = xin(217) + dxij*xin(205)
                                      yin(209) = yin(217) + dyij*yin(205)
                                      zin(209) = zin(217) + dzij*zin(205)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  221

                                      ! ni =    3

                                      xin(221) = xin(229) + dxij*xin(217)
                                      yin(221) = yin(229) + dyij*yin(217)
                                      zin(221) = zin(229) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  233

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  201

                                      ! nj =    2

                                      ! i4 = i3 =  201

                                      ! do ni = 1,    3

                                      xin(201) = xin(209) + dxij*xin(197)
                                      yin(201) = yin(209) + dyij*yin(197)
                                      zin(201) = zin(209) + dzij*zin(197)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  213

                                      ! ni =    2

                                      xin(213) = xin(221) + dxij*xin(209)
                                      yin(213) = yin(221) + dyij*yin(209)
                                      zin(213) = zin(221) + dzij*zin(209)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  225

                                      ! ni =    3

                                      xin(225) = xin(233) + dxij*xin(221)
                                      yin(225) = yin(233) + dyij*yin(221)
                                      zin(225) = zin(233) + dzij*zin(221)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  237

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  205

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    1

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  238

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  234

                                      xin(238) = xin(238) + dxij*xin(234)
                                      yin(238) = yin(238) + dyij*yin(234)
                                      zin(238) = zin(238) + dzij*zin(234)

                                      ! i3 = i4 =  234
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  230

                                      xin(234) = xin(234) + dxij*xin(230)
                                      yin(234) = yin(234) + dyij*yin(230)
                                      zin(234) = zin(234) + dzij*zin(230)

                                      ! i3 = i4 =  230
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  238

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  234

                                      xin(238) = xin(238) + dxij*xin(234)
                                      yin(238) = yin(238) + dyij*yin(234)
                                      zin(238) = zin(238) + dzij*zin(234)

                                      ! i3 = i4 =  234
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  198

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  198

                                      ! do ni = 1,    3

                                      xin(198) = xin(206) + dxij*xin(194)
                                      yin(198) = yin(206) + dyij*yin(194)
                                      zin(198) = zin(206) + dzij*zin(194)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  210

                                      ! ni =    2

                                      xin(210) = xin(218) + dxij*xin(206)
                                      yin(210) = yin(218) + dyij*yin(206)
                                      zin(210) = zin(218) + dzij*zin(206)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  222

                                      ! ni =    3

                                      xin(222) = xin(230) + dxij*xin(218)
                                      yin(222) = yin(230) + dyij*yin(218)
                                      zin(222) = zin(230) + dzij*zin(218)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  234

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  202

                                      ! nj =    2

                                      ! i4 = i3 =  202

                                      ! do ni = 1,    3

                                      xin(202) = xin(210) + dxij*xin(198)
                                      yin(202) = yin(210) + dyij*yin(198)
                                      zin(202) = zin(210) + dzij*zin(198)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                      ! ni =    2

                                      xin(214) = xin(222) + dxij*xin(210)
                                      yin(214) = yin(222) + dyij*yin(210)
                                      zin(214) = zin(222) + dzij*zin(210)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  226

                                      ! ni =    3

                                      xin(226) = xin(234) + dxij*xin(222)
                                      yin(226) = yin(234) + dyij*yin(222)
                                      zin(226) = zin(234) + dzij*zin(222)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  238

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  206

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    2

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  239

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  235

                                      xin(239) = xin(239) + dxij*xin(235)
                                      yin(239) = yin(239) + dyij*yin(235)
                                      zin(239) = zin(239) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  231

                                      xin(235) = xin(235) + dxij*xin(231)
                                      yin(235) = yin(235) + dyij*yin(231)
                                      zin(235) = zin(235) + dzij*zin(231)

                                      ! i3 = i4 =  231
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  239

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  235

                                      xin(239) = xin(239) + dxij*xin(235)
                                      yin(239) = yin(239) + dyij*yin(235)
                                      zin(239) = zin(239) + dzij*zin(235)

                                      ! i3 = i4 =  235
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  199

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  199

                                      ! do ni = 1,    3

                                      xin(199) = xin(207) + dxij*xin(195)
                                      yin(199) = yin(207) + dyij*yin(195)
                                      zin(199) = zin(207) + dzij*zin(195)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    2

                                      xin(211) = xin(219) + dxij*xin(207)
                                      yin(211) = yin(219) + dyij*yin(207)
                                      zin(211) = zin(219) + dzij*zin(207)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  223

                                      ! ni =    3

                                      xin(223) = xin(231) + dxij*xin(219)
                                      yin(223) = yin(231) + dyij*yin(219)
                                      zin(223) = zin(231) + dzij*zin(219)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  203

                                      ! nj =    2

                                      ! i4 = i3 =  203

                                      ! do ni = 1,    3

                                      xin(203) = xin(211) + dxij*xin(199)
                                      yin(203) = yin(211) + dyij*yin(199)
                                      zin(203) = zin(211) + dzij*zin(199)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    2

                                      xin(215) = xin(223) + dxij*xin(211)
                                      yin(215) = yin(223) + dyij*yin(211)
                                      zin(215) = zin(223) + dzij*zin(211)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  227

                                      ! ni =    3

                                      xin(227) = xin(235) + dxij*xin(223)
                                      yin(227) = yin(235) + dyij*yin(223)
                                      zin(227) = zin(235) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  207

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  240

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  236

                                      xin(240) = xin(240) + dxij*xin(236)
                                      yin(240) = yin(240) + dyij*yin(236)
                                      zin(240) = zin(240) + dzij*zin(236)

                                      ! i3 = i4 =  236
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  232

                                      xin(236) = xin(236) + dxij*xin(232)
                                      yin(236) = yin(236) + dyij*yin(232)
                                      zin(236) = zin(236) + dzij*zin(232)

                                      ! i3 = i4 =  232
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  240

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  236

                                      xin(240) = xin(240) + dxij*xin(236)
                                      yin(240) = yin(240) + dyij*yin(236)
                                      zin(240) = zin(240) + dzij*zin(236)

                                      ! i3 = i4 =  236
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  200

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  200

                                      ! do ni = 1,    3

                                      xin(200) = xin(208) + dxij*xin(196)
                                      yin(200) = yin(208) + dyij*yin(196)
                                      zin(200) = zin(208) + dzij*zin(196)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  212

                                      ! ni =    2

                                      xin(212) = xin(220) + dxij*xin(208)
                                      yin(212) = yin(220) + dyij*yin(208)
                                      zin(212) = zin(220) + dzij*zin(208)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  224

                                      ! ni =    3

                                      xin(224) = xin(232) + dxij*xin(220)
                                      yin(224) = yin(232) + dyij*yin(220)
                                      zin(224) = zin(232) + dzij*zin(220)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  236

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  204

                                      ! nj =    2

                                      ! i4 = i3 =  204

                                      ! do ni = 1,    3

                                      xin(204) = xin(212) + dxij*xin(200)
                                      yin(204) = yin(212) + dyij*yin(200)
                                      zin(204) = zin(212) + dzij*zin(200)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    2

                                      xin(216) = xin(224) + dxij*xin(212)
                                      yin(216) = yin(224) + dyij*yin(212)
                                      zin(216) = zin(224) + dzij*zin(212)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  228

                                      ! ni =    3

                                      xin(228) = xin(236) + dxij*xin(224)
                                      yin(228) = yin(236) + dyij*yin(224)
                                      zin(228) = zin(236) + dzij*zin(224)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  208

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  240

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 600! loop over all integrals

                                        l = n - 10*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d23bra(j)*d03ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 48)*yin(my + 48)*zin(mz + 48) & ! root  2
                                                        + xin(mx + 96)*yin(my + 96)*zin(mz + 96) & ! root  3
                                                        + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  4
                                                        + xin(mx + 192)*yin(my + 192)*zin(mz + 192)) ! root  5

                                        j = int(n/10) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*60 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*10 + ip ! Add stride between functions in j

                                      do k = 1, 10 ! # of cartesians in k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*1 + ijp ! Add stride between functions in k

                                        do l = 1, 1! # of cartesians in l

                                          ijklp = ijkp + l ! No stride between functions in l

                                          buff(1) = eri_value(ijklp)

                                          if (abs(buff(1)) .lt. 5.0D-11) cycle ! Skip small integrals

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

                                          ! Account for identical permutations

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

                                        end do ! l
                                      end do ! k
                                    end do ! j
                                  end do ! i

                                  !                  --- END DIRFCK_RHF ---

                                  end if ! Screening if

                                end do ! iquart

                                !$omp end target teams distribute parallel do


                              end do ! itile

                              deallocate (n23bra)
                              deallocate (xint23bra)
                              deallocate (n03ket)
                              deallocate (xint03ket)

                              end subroutine int3230
                              end submodule
