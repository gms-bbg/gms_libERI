! The total angular momentum of this class is:           9
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3222_impl
contains
  module subroutine int3222(df_pair, dd_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair, dd_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:), n22ket(:)
    real(dp), allocatable :: xint23bra(:), xint22ket(:)
    integer(kind=int64) :: ndfbra, nddket
    real(dp) :: scutdfbra, scutddket, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iquart
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2
    integer(kind=int64) :: n, i1, i3, i4, i5, k3, k4, nn, nm, km, iaa, ib, nj, ni, nl, nk
    real(dp) :: cp10, c10, cp01, c01
    integer(kind=int64) :: nx, ny, nz, mx, my, mz
    integer(kind=int64) :: bra_loop, ket_loop, ijtop, kltop
    real(dp) :: t_expon_ab, t_expon_cd, t_inverse_expon_ab, t_inverse_expon_cd
    real(dp) :: t_expon_abcd_inverse, rho, expe, dum, rab, rcd
    real(dp) :: brrk, akxk, akyk, akzk, t_expon_c, t_expon_d, t_expon_a, t_expon_b
    real(dp) :: xa, ya, za, axak, ayak, azak, axai, ayai, azai, bbrrk, xb, yb, zb, bxbk
    real(dp) :: bybk, bzbk, bxbi, bybi, bzbi, xx, c1x, c2x, c3x, c4x, c1y, c2y, c3y, c4y
    real(dp) :: c1z, c2z, c3z, c4z, f00, u2, duminv, dm2inv, bp01, b00, b10, xcp00, xc00
    real(dp) :: ycp00, zcp00, zc00, yc00, dij, dxij, dyij, dzij, dxkl, dykl, dzkl
    real(dp) :: buff(9)
    real(dp) :: roots(5), wghts(5)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(35), wgrid(35), p0(35), p1(35), p2(35)
    real(dp) :: rts(5), wts(5), alpha(5), beta(5), wrk(5)
    real(dp) :: xin(540), yin(540), zin(540)
    real(dp) :: eri_value(2160)
    real(dp) :: d23bra(60), d22ket(36)
    integer(kind=int64) :: ix(10), jx(6), kx(6), lx(6)
    integer(kind=int64) :: iy(10), jy(6), ky(6), ly(6)
    integer(kind=int64) :: iz(10), jz(6), kz(6), lz(6)
    integer(kind=int64) :: in(6), in1(6), kn(5)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(36), kly(36), klz(36)
    integer(kind=int64) :: nchunksize_int64, istart, iend, itile, ntile
    logical :: kandl

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 28
    in1(3) = 55
    in1(4) = 82
    in1(5) = 91
    in1(6) = 100

    kn(1) = 0
    kn(2) = 3
    kn(3) = 6
    kn(4) = 7
    kn(5) = 8

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 2
    lx(2) = 0
    lx(3) = 0
    lx(4) = 1
    lx(5) = 1
    lx(6) = 0

    kx(1) = 6
    kx(2) = 0
    kx(3) = 0
    kx(4) = 3
    kx(5) = 3
    kx(6) = 0

    jx(1) = 18
    jx(2) = 0
    jx(3) = 0
    jx(4) = 9
    jx(5) = 9
    jx(6) = 0

    ix(1) = 82
    ix(2) = 1
    ix(3) = 1
    ix(4) = 55
    ix(5) = 55
    ix(6) = 28
    ix(7) = 1
    ix(8) = 28
    ix(9) = 1
    ix(10) = 28

    ! y-arrays

    ly(1) = 0
    ly(2) = 2
    ly(3) = 0
    ly(4) = 1
    ly(5) = 0
    ly(6) = 1

    ky(1) = 0
    ky(2) = 6
    ky(3) = 0
    ky(4) = 3
    ky(5) = 0
    ky(6) = 3

    jy(1) = 0
    jy(2) = 18
    jy(3) = 0
    jy(4) = 9
    jy(5) = 0
    jy(6) = 9

    iy(1) = 1
    iy(2) = 82
    iy(3) = 1
    iy(4) = 28
    iy(5) = 1
    iy(6) = 55
    iy(7) = 55
    iy(8) = 1
    iy(9) = 28
    iy(10) = 28

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 2
    lz(4) = 0
    lz(5) = 1
    lz(6) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 6
    kz(4) = 0
    kz(5) = 3
    kz(6) = 3

    jz(1) = 0
    jz(2) = 0
    jz(3) = 18
    jz(4) = 0
    jz(5) = 9
    jz(6) = 9

    iz(1) = 1
    iz(2) = 1
    iz(3) = 82
    iz(4) = 1
    iz(5) = 28
    iz(6) = 1
    iz(7) = 28
    iz(8) = 55
    iz(9) = 55
    iz(10) = 28

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 100
    ijx(2) = 82
    ijx(3) = 82
    ijx(4) = 91
    ijx(5) = 91
    ijx(6) = 82
    ijx(7) = 19
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 10
    ijx(11) = 10
    ijx(12) = 1
    ijx(13) = 19
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 10
    ijx(17) = 10
    ijx(18) = 1
    ijx(19) = 73
    ijx(20) = 55
    ijx(21) = 55
    ijx(22) = 64
    ijx(23) = 64
    ijx(24) = 55
    ijx(25) = 73
    ijx(26) = 55
    ijx(27) = 55
    ijx(28) = 64
    ijx(29) = 64
    ijx(30) = 55
    ijx(31) = 46
    ijx(32) = 28
    ijx(33) = 28
    ijx(34) = 37
    ijx(35) = 37
    ijx(36) = 28
    ijx(37) = 19
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 10
    ijx(41) = 10
    ijx(42) = 1
    ijx(43) = 46
    ijx(44) = 28
    ijx(45) = 28
    ijx(46) = 37
    ijx(47) = 37
    ijx(48) = 28
    ijx(49) = 19
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 10
    ijx(53) = 10
    ijx(54) = 1
    ijx(55) = 46
    ijx(56) = 28
    ijx(57) = 28
    ijx(58) = 37
    ijx(59) = 37
    ijx(60) = 28

    ijy(1) = 1
    ijy(2) = 19
    ijy(3) = 1
    ijy(4) = 10
    ijy(5) = 1
    ijy(6) = 10
    ijy(7) = 82
    ijy(8) = 100
    ijy(9) = 82
    ijy(10) = 91
    ijy(11) = 82
    ijy(12) = 91
    ijy(13) = 1
    ijy(14) = 19
    ijy(15) = 1
    ijy(16) = 10
    ijy(17) = 1
    ijy(18) = 10
    ijy(19) = 28
    ijy(20) = 46
    ijy(21) = 28
    ijy(22) = 37
    ijy(23) = 28
    ijy(24) = 37
    ijy(25) = 1
    ijy(26) = 19
    ijy(27) = 1
    ijy(28) = 10
    ijy(29) = 1
    ijy(30) = 10
    ijy(31) = 55
    ijy(32) = 73
    ijy(33) = 55
    ijy(34) = 64
    ijy(35) = 55
    ijy(36) = 64
    ijy(37) = 55
    ijy(38) = 73
    ijy(39) = 55
    ijy(40) = 64
    ijy(41) = 55
    ijy(42) = 64
    ijy(43) = 1
    ijy(44) = 19
    ijy(45) = 1
    ijy(46) = 10
    ijy(47) = 1
    ijy(48) = 10
    ijy(49) = 28
    ijy(50) = 46
    ijy(51) = 28
    ijy(52) = 37
    ijy(53) = 28
    ijy(54) = 37
    ijy(55) = 28
    ijy(56) = 46
    ijy(57) = 28
    ijy(58) = 37
    ijy(59) = 28
    ijy(60) = 37

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 19
    ijz(4) = 1
    ijz(5) = 10
    ijz(6) = 10
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 19
    ijz(10) = 1
    ijz(11) = 10
    ijz(12) = 10
    ijz(13) = 82
    ijz(14) = 82
    ijz(15) = 100
    ijz(16) = 82
    ijz(17) = 91
    ijz(18) = 91
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 19
    ijz(22) = 1
    ijz(23) = 10
    ijz(24) = 10
    ijz(25) = 28
    ijz(26) = 28
    ijz(27) = 46
    ijz(28) = 28
    ijz(29) = 37
    ijz(30) = 37
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 19
    ijz(34) = 1
    ijz(35) = 10
    ijz(36) = 10
    ijz(37) = 28
    ijz(38) = 28
    ijz(39) = 46
    ijz(40) = 28
    ijz(41) = 37
    ijz(42) = 37
    ijz(43) = 55
    ijz(44) = 55
    ijz(45) = 73
    ijz(46) = 55
    ijz(47) = 64
    ijz(48) = 64
    ijz(49) = 55
    ijz(50) = 55
    ijz(51) = 73
    ijz(52) = 55
    ijz(53) = 64
    ijz(54) = 64
    ijz(55) = 28
    ijz(56) = 28
    ijz(57) = 46
    ijz(58) = 28
    ijz(59) = 37
    ijz(60) = 37

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 8
    klx(2) = 6
    klx(3) = 6
    klx(4) = 7
    klx(5) = 7
    klx(6) = 6
    klx(7) = 2
    klx(8) = 0
    klx(9) = 0
    klx(10) = 1
    klx(11) = 1
    klx(12) = 0
    klx(13) = 2
    klx(14) = 0
    klx(15) = 0
    klx(16) = 1
    klx(17) = 1
    klx(18) = 0
    klx(19) = 5
    klx(20) = 3
    klx(21) = 3
    klx(22) = 4
    klx(23) = 4
    klx(24) = 3
    klx(25) = 5
    klx(26) = 3
    klx(27) = 3
    klx(28) = 4
    klx(29) = 4
    klx(30) = 3
    klx(31) = 2
    klx(32) = 0
    klx(33) = 0
    klx(34) = 1
    klx(35) = 1
    klx(36) = 0

    kly(1) = 0
    kly(2) = 2
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 1
    kly(7) = 6
    kly(8) = 8
    kly(9) = 6
    kly(10) = 7
    kly(11) = 6
    kly(12) = 7
    kly(13) = 0
    kly(14) = 2
    kly(15) = 0
    kly(16) = 1
    kly(17) = 0
    kly(18) = 1
    kly(19) = 3
    kly(20) = 5
    kly(21) = 3
    kly(22) = 4
    kly(23) = 3
    kly(24) = 4
    kly(25) = 0
    kly(26) = 2
    kly(27) = 0
    kly(28) = 1
    kly(29) = 0
    kly(30) = 1
    kly(31) = 3
    kly(32) = 5
    kly(33) = 3
    kly(34) = 4
    kly(35) = 3
    kly(36) = 4

    klz(1) = 0
    klz(2) = 0
    klz(3) = 2
    klz(4) = 0
    klz(5) = 1
    klz(6) = 1
    klz(7) = 0
    klz(8) = 0
    klz(9) = 2
    klz(10) = 0
    klz(11) = 1
    klz(12) = 1
    klz(13) = 6
    klz(14) = 6
    klz(15) = 8
    klz(16) = 6
    klz(17) = 7
    klz(18) = 7
    klz(19) = 0
    klz(20) = 0
    klz(21) = 2
    klz(22) = 0
    klz(23) = 1
    klz(24) = 1
    klz(25) = 3
    klz(26) = 3
    klz(27) = 5
    klz(28) = 3
    klz(29) = 4
    klz(30) = 4
    klz(31) = 3
    klz(32) = 3
    klz(33) = 5
    klz(34) = 3
    klz(35) = 4
    klz(36) = 4

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))
    allocate (n22ket(res%n_d_shl*(res%n_d_shl + 1)/2))
    allocate (xint22ket(res%n_d_shl*(res%n_d_shl + 1)/2))

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

    scutddket = cutoff_schwarz/maxval(dd_pair%xints)
    nddket = 0
    do ij = 1, res%n_d_shl*(res%n_d_shl + 1)/2
      if (dd_pair%xints(ij) .ge. scutddket) then
        nddket = nddket + 1
        xint22ket(nddket) = dd_pair%xints(ij)
        n22ket(nddket) = ij
      end if
    end do

    nchunksize_int64 = 375000000

    if ((ndfbra*nddket) .le. nchunksize_int64) nchunksize_int64 = ndfbra*nddket
    ntile = int(ndfbra*nddket/nchunksize_int64)
    do itile = 1, ntile
      istart = (itile - 1)*nchunksize_int64 + 1
      iend = itile*nchunksize_int64
      if (itile .eq. ntile) iend = ndfbra*nddket

      ! --multi-gpu--work
      nchunk = (iend - istart)/res%n_size
      nquart_start = nchunk*res%n_rank + istart
      nquart_end = nquart_start + nchunk - 1
      if (res%n_rank .eq. res%n_size - 1) nquart_end = iend

      ! Mappings to GPU

 !$omp target teams distribute parallel do default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, xint22ket, n22ket, dd_pair, df_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(ij_tmp,kl_tmp,test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d22ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,kandl)
            do iquart = nquart_start, nquart_end

              ij_tmp = mod(iquart - 1, ndfbra) + 1
              kl_tmp = (iquart - 1)/ndfbra + 1

              test = xint23bra(ij_tmp)*xint22ket(kl_tmp)

              if (test .gt. cutoff_schwarz) then

                ij = n23bra(ij_tmp)
                kl = n22ket(kl_tmp)

                ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                jsh_tmp = (ij - 1)/res%n_f_shl + 1
                ksh_tmp = (1 + sqrt(1.0 + 8.0*(kl - 1)))/2
                lsh_tmp = kl - ksh_tmp*(ksh_tmp - 1)/2

                ish = res%i_f_shl(ish_tmp)
                jsh = res%i_d_shl(jsh_tmp)
                ksh = res%i_d_shl(ksh_tmp)
                lsh = res%i_d_shl(lsh_tmp)

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

                ! Distances used in the ket transfer equation

                dxkl = res%coord_sh(ksh, 1) - res%coord_sh(lsh, 1)
                dykl = res%coord_sh(ksh, 2) - res%coord_sh(lsh, 2)
                dzkl = res%coord_sh(ksh, 3) - res%coord_sh(lsh, 3)

                ! Prepare to begin looping

                eri_value = 0.0_dp

                ket_loop = 0

                ! --- Start looping over primitives ---

                do k = 1, kltop

                  ket_loop = ket_loop + 1

                  t_expon_cd = dd_pair%t_expon_ab(dd_pair%pair_loc(kl) + ket_loop)
                  t_expon_c = dd_pair%expon_a(dd_pair%pair_loc(kl) + ket_loop)
                  t_expon_d = dd_pair%expon_b(dd_pair%pair_loc(kl) + ket_loop)
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

                  d22ket(1) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(2) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(3) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(4) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(5) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(6) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(7) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(8) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(9) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(10) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(11) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(12) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(13) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(14) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(15) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                  d22ket(16) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(17) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(18) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(19) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(20) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(21) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(22) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(23) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(24) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(25) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(26) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(27) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(28) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(29) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(30) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(31) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(32) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(33) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                  d22ket(34) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(35) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3
                  d22ket(36) = dd_pair%d_coeff_alt(dd_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3*sqrt3

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

                                      ! i2 = in(2) =   28
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(28) = xc00
                                      yin(28) = yc00
                                      zin(28) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =    4

                                      xin(4) = xcp00
                                      yin(4) = ycp00
                                      zin(4) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =   31
                                      ! i2 =   28

                                      xin(31) = xcp00*xin(28) + cp10
                                      yin(31) = ycp00*yin(28) + cp10
                                      zin(31) = zcp00*zin(28) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   28

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   55
                                      ! i3 =    1
                                      ! i4 =   28

                                      xin(55) = c10*xin(1) + xc00*xin(28)
                                      yin(55) = c10*yin(1) + yc00*yin(28)
                                      zin(55) = c10*zin(1) + zc00*zin(28)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   58
                                      ! i5 =   55
                                      ! i4 =   28

                                      xin(58) = xcp00*xin(55) + cp10*xin(28)
                                      yin(58) = ycp00*yin(55) + cp10*yin(28)
                                      zin(58) = zcp00*zin(55) + cp10*zin(28)

                                      ! ------------------

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   55

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   82
                                      ! i3 =   28
                                      ! i4 =   55

                                      xin(82) = c10*xin(28) + xc00*xin(55)
                                      yin(82) = c10*yin(28) + yc00*yin(55)
                                      zin(82) = c10*zin(28) + zc00*zin(55)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   85
                                      ! i5 =   82
                                      ! i4 =   55

                                      xin(85) = xcp00*xin(82) + cp10*xin(55)
                                      yin(85) = ycp00*yin(82) + cp10*yin(55)
                                      zin(85) = zcp00*zin(82) + cp10*zin(55)

                                      ! ------------------

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   82

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =   91
                                      ! i3 =   55
                                      ! i4 =   82

                                      xin(91) = c10*xin(55) + xc00*xin(82)
                                      yin(91) = c10*yin(55) + yc00*yin(82)
                                      zin(91) = c10*zin(55) + zc00*zin(82)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =   94
                                      ! i5 =   91
                                      ! i4 =   82

                                      xin(94) = xcp00*xin(91) + cp10*xin(82)
                                      yin(94) = ycp00*yin(91) + cp10*yin(82)
                                      zin(94) = zcp00*zin(91) + cp10*zin(82)

                                      ! ------------------

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  100
                                      ! i3 =   82
                                      ! i4 =   91

                                      xin(100) = c10*xin(82) + xc00*xin(91)
                                      yin(100) = c10*yin(82) + yc00*yin(91)
                                      zin(100) = c10*zin(82) + zc00*zin(91)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  103
                                      ! i5 =  100
                                      ! i4 =   91

                                      xin(103) = xcp00*xin(100) + cp10*xin(91)
                                      yin(103) = ycp00*yin(100) + cp10*yin(91)
                                      zin(103) = zcp00*zin(100) + cp10*zin(91)

                                      ! ------------------

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =    1
                                      ! i4 = i1+k2 =    4

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    7
                                      ! i3 =    1
                                      ! i4 =    4

                                      xin(7) = cp01*xin(1) + xcp00*xin(4)
                                      yin(7) = cp01*yin(1) + ycp00*yin(4)
                                      zin(7) = cp01*zin(1) + zcp00*zin(4)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   34

                                      xin(34) = xc00*xin(7) + c01*xin(4)
                                      yin(34) = yc00*yin(7) + c01*yin(4)
                                      zin(34) = zc00*zin(7) + c01*zin(4)

                                      ! ------------------

                                      ! i3 = i4 =    4
                                      ! i4 = i5 =    7

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    8
                                      ! i3 =    4
                                      ! i4 =    7

                                      xin(8) = cp01*xin(4) + xcp00*xin(7)
                                      yin(8) = cp01*yin(4) + ycp00*yin(7)
                                      zin(8) = cp01*zin(4) + zcp00*zin(7)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   35

                                      xin(35) = xc00*xin(8) + c01*xin(7)
                                      yin(35) = yc00*yin(8) + c01*yin(7)
                                      zin(35) = zc00*zin(8) + c01*zin(7)

                                      ! ------------------

                                      ! i3 = i4 =    7
                                      ! i4 = i5 =    8

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =    9
                                      ! i3 =    7
                                      ! i4 =    8

                                      xin(9) = cp01*xin(7) + xcp00*xin(8)
                                      yin(9) = cp01*yin(7) + ycp00*yin(8)
                                      zin(9) = cp01*zin(7) + zcp00*zin(8)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =   36

                                      xin(36) = xc00*xin(9) + c01*xin(8)
                                      yin(36) = yc00*yin(9) + c01*yin(8)
                                      zin(36) = zc00*zin(9) + c01*zin(8)

                                      ! ------------------

                                      ! i3 = i4 =    8
                                      ! i4 = i5 =    9

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   28

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   55

                                      xin(61) = c10*xin(7) + xc00*xin(34) + c01*xin(31)
                                      yin(61) = c10*yin(7) + yc00*yin(34) + c01*yin(31)
                                      zin(61) = c10*zin(7) + zc00*zin(34) + c01*zin(31)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   55

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   82

                                      xin(88) = c10*xin(34) + xc00*xin(61) + c01*xin(58)
                                      yin(88) = c10*yin(34) + yc00*yin(61) + c01*yin(58)
                                      zin(88) = c10*zin(34) + zc00*zin(61) + c01*zin(58)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   82

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   91

                                      xin(97) = c10*xin(61) + xc00*xin(88) + c01*xin(85)
                                      yin(97) = c10*yin(61) + yc00*yin(88) + c01*yin(85)
                                      zin(97) = c10*zin(61) + zc00*zin(88) + c01*zin(85)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  100

                                      xin(106) = c10*xin(88) + xc00*xin(97) + c01*xin(94)
                                      yin(106) = c10*yin(88) + yc00*yin(97) + c01*yin(94)
                                      zin(106) = c10*zin(88) + zc00*zin(97) + c01*zin(94)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   28

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   55

                                      xin(62) = c10*xin(8) + xc00*xin(35) + c01*xin(34)
                                      yin(62) = c10*yin(8) + yc00*yin(35) + c01*yin(34)
                                      zin(62) = c10*zin(8) + zc00*zin(35) + c01*zin(34)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   55

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   82

                                      xin(89) = c10*xin(35) + xc00*xin(62) + c01*xin(61)
                                      yin(89) = c10*yin(35) + yc00*yin(62) + c01*yin(61)
                                      zin(89) = c10*zin(35) + zc00*zin(62) + c01*zin(61)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   82

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   91

                                      xin(98) = c10*xin(62) + xc00*xin(89) + c01*xin(88)
                                      yin(98) = c10*yin(62) + yc00*yin(89) + c01*yin(88)
                                      zin(98) = c10*zin(62) + zc00*zin(89) + c01*zin(88)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  100

                                      xin(107) = c10*xin(89) + xc00*xin(98) + c01*xin(97)
                                      yin(107) = c10*yin(89) + yc00*yin(98) + c01*yin(97)
                                      zin(107) = c10*zin(89) + zc00*zin(98) + c01*zin(97)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =    1
                                      ! i4 = i2 =   28

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =   55

                                      xin(63) = c10*xin(9) + xc00*xin(36) + c01*xin(35)
                                      yin(63) = c10*yin(9) + yc00*yin(36) + c01*yin(35)
                                      zin(63) = c10*zin(9) + zc00*zin(36) + c01*zin(35)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   28
                                      ! i4 = i5 =   55

                                      ! nn =    3

                                      ! i5 = in(nn+1) =   82

                                      xin(90) = c10*xin(36) + xc00*xin(63) + c01*xin(62)
                                      yin(90) = c10*yin(36) + yc00*yin(63) + c01*yin(62)
                                      zin(90) = c10*zin(36) + zc00*zin(63) + c01*zin(62)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   55
                                      ! i4 = i5 =   82

                                      ! nn =    4

                                      ! i5 = in(nn+1) =   91

                                      xin(99) = c10*xin(63) + xc00*xin(90) + c01*xin(89)
                                      yin(99) = c10*yin(63) + yc00*yin(90) + c01*yin(89)
                                      zin(99) = c10*zin(63) + zc00*zin(90) + c01*zin(89)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   82
                                      ! i4 = i5 =   91

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  100

                                      xin(108) = c10*xin(90) + xc00*xin(99) + c01*xin(98)
                                      yin(108) = c10*yin(90) + yc00*yin(99) + c01*yin(98)
                                      zin(108) = c10*zin(90) + zc00*zin(99) + c01*zin(98)

                                      c10 = c10 + b10

                                      ! i3 = i4 =   91
                                      ! i4 = i5 =  100

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  100

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  100

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(100) = xin(100) + dxij*xin(91)
                                      yin(100) = yin(100) + dyij*yin(91)
                                      zin(100) = zin(100) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   82

                                      xin(91) = xin(91) + dxij*xin(82)
                                      yin(91) = yin(91) + dyij*yin(82)
                                      zin(91) = zin(91) + dzij*zin(82)

                                      ! i3 = i4 =   82
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  100

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   91

                                      xin(100) = xin(100) + dxij*xin(91)
                                      yin(100) = yin(100) + dyij*yin(91)
                                      zin(100) = zin(100) + dzij*zin(91)

                                      ! i3 = i4 =   91
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   10

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   10

                                      ! do ni = 1,    3

                                      xin(10) = xin(28) + dxij*xin(1)
                                      yin(10) = yin(28) + dyij*yin(1)
                                      zin(10) = zin(28) + dzij*zin(1)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   37

                                      ! ni =    2

                                      xin(37) = xin(55) + dxij*xin(28)
                                      yin(37) = yin(55) + dyij*yin(28)
                                      zin(37) = zin(55) + dzij*zin(28)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                      ! ni =    3

                                      xin(64) = xin(82) + dxij*xin(55)
                                      yin(64) = yin(82) + dyij*yin(55)
                                      zin(64) = zin(82) + dzij*zin(55)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   19

                                      ! nj =    2

                                      ! i4 = i3 =   19

                                      ! do ni = 1,    3

                                      xin(19) = xin(37) + dxij*xin(10)
                                      yin(19) = yin(37) + dyij*yin(10)
                                      zin(19) = zin(37) + dzij*zin(10)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   46

                                      ! ni =    2

                                      xin(46) = xin(64) + dxij*xin(37)
                                      yin(46) = yin(64) + dyij*yin(37)
                                      zin(46) = zin(64) + dzij*zin(37)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   73

                                      ! ni =    3

                                      xin(73) = xin(91) + dxij*xin(64)
                                      yin(73) = yin(91) + dyij*yin(64)
                                      zin(73) = zin(91) + dzij*zin(64)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  100

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   28

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  103

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   94

                                      xin(103) = xin(103) + dxij*xin(94)
                                      yin(103) = yin(103) + dyij*yin(94)
                                      zin(103) = zin(103) + dzij*zin(94)

                                      ! i3 = i4 =   94
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   85

                                      xin(94) = xin(94) + dxij*xin(85)
                                      yin(94) = yin(94) + dyij*yin(85)
                                      zin(94) = zin(94) + dzij*zin(85)

                                      ! i3 = i4 =   85
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  103

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   94

                                      xin(103) = xin(103) + dxij*xin(94)
                                      yin(103) = yin(103) + dyij*yin(94)
                                      zin(103) = zin(103) + dzij*zin(94)

                                      ! i3 = i4 =   94
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   13

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   13

                                      ! do ni = 1,    3

                                      xin(13) = xin(31) + dxij*xin(4)
                                      yin(13) = yin(31) + dyij*yin(4)
                                      zin(13) = zin(31) + dzij*zin(4)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   40

                                      ! ni =    2

                                      xin(40) = xin(58) + dxij*xin(31)
                                      yin(40) = yin(58) + dyij*yin(31)
                                      zin(40) = zin(58) + dzij*zin(31)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                      ! ni =    3

                                      xin(67) = xin(85) + dxij*xin(58)
                                      yin(67) = yin(85) + dyij*yin(58)
                                      zin(67) = zin(85) + dzij*zin(58)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   22

                                      ! nj =    2

                                      ! i4 = i3 =   22

                                      ! do ni = 1,    3

                                      xin(22) = xin(40) + dxij*xin(13)
                                      yin(22) = yin(40) + dyij*yin(13)
                                      zin(22) = zin(40) + dzij*zin(13)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                      ! ni =    2

                                      xin(49) = xin(67) + dxij*xin(40)
                                      yin(49) = yin(67) + dyij*yin(40)
                                      zin(49) = zin(67) + dzij*zin(40)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   76

                                      ! ni =    3

                                      xin(76) = xin(94) + dxij*xin(67)
                                      yin(76) = yin(94) + dyij*yin(67)
                                      zin(76) = zin(94) + dzij*zin(67)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   31

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   97

                                      xin(106) = xin(106) + dxij*xin(97)
                                      yin(106) = yin(106) + dyij*yin(97)
                                      zin(106) = zin(106) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   88

                                      xin(97) = xin(97) + dxij*xin(88)
                                      yin(97) = yin(97) + dyij*yin(88)
                                      zin(97) = zin(97) + dzij*zin(88)

                                      ! i3 = i4 =   88
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  106

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   97

                                      xin(106) = xin(106) + dxij*xin(97)
                                      yin(106) = yin(106) + dyij*yin(97)
                                      zin(106) = zin(106) + dzij*zin(97)

                                      ! i3 = i4 =   97
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   16

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   16

                                      ! do ni = 1,    3

                                      xin(16) = xin(34) + dxij*xin(7)
                                      yin(16) = yin(34) + dyij*yin(7)
                                      zin(16) = zin(34) + dzij*zin(7)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   43

                                      ! ni =    2

                                      xin(43) = xin(61) + dxij*xin(34)
                                      yin(43) = yin(61) + dyij*yin(34)
                                      zin(43) = zin(61) + dzij*zin(34)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                      ! ni =    3

                                      xin(70) = xin(88) + dxij*xin(61)
                                      yin(70) = yin(88) + dyij*yin(61)
                                      zin(70) = zin(88) + dzij*zin(61)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   25

                                      ! nj =    2

                                      ! i4 = i3 =   25

                                      ! do ni = 1,    3

                                      xin(25) = xin(43) + dxij*xin(16)
                                      yin(25) = yin(43) + dyij*yin(16)
                                      zin(25) = zin(43) + dzij*zin(16)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   52

                                      ! ni =    2

                                      xin(52) = xin(70) + dxij*xin(43)
                                      yin(52) = yin(70) + dyij*yin(43)
                                      zin(52) = zin(70) + dzij*zin(43)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   79

                                      ! ni =    3

                                      xin(79) = xin(97) + dxij*xin(70)
                                      yin(79) = yin(97) + dyij*yin(70)
                                      zin(79) = zin(97) + dzij*zin(70)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   34

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   98

                                      xin(107) = xin(107) + dxij*xin(98)
                                      yin(107) = yin(107) + dyij*yin(98)
                                      zin(107) = zin(107) + dzij*zin(98)

                                      ! i3 = i4 =   98
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   89

                                      xin(98) = xin(98) + dxij*xin(89)
                                      yin(98) = yin(98) + dyij*yin(89)
                                      zin(98) = zin(98) + dzij*zin(89)

                                      ! i3 = i4 =   89
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  107

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   98

                                      xin(107) = xin(107) + dxij*xin(98)
                                      yin(107) = yin(107) + dyij*yin(98)
                                      zin(107) = zin(107) + dzij*zin(98)

                                      ! i3 = i4 =   98
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   17

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   17

                                      ! do ni = 1,    3

                                      xin(17) = xin(35) + dxij*xin(8)
                                      yin(17) = yin(35) + dyij*yin(8)
                                      zin(17) = zin(35) + dzij*zin(8)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   44

                                      ! ni =    2

                                      xin(44) = xin(62) + dxij*xin(35)
                                      yin(44) = yin(62) + dyij*yin(35)
                                      zin(44) = zin(62) + dzij*zin(35)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                      ! ni =    3

                                      xin(71) = xin(89) + dxij*xin(62)
                                      yin(71) = yin(89) + dyij*yin(62)
                                      zin(71) = zin(89) + dzij*zin(62)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   98

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   26

                                      ! nj =    2

                                      ! i4 = i3 =   26

                                      ! do ni = 1,    3

                                      xin(26) = xin(44) + dxij*xin(17)
                                      yin(26) = yin(44) + dyij*yin(17)
                                      zin(26) = zin(44) + dzij*zin(17)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   53

                                      ! ni =    2

                                      xin(53) = xin(71) + dxij*xin(44)
                                      yin(53) = yin(71) + dyij*yin(44)
                                      zin(53) = zin(71) + dzij*zin(44)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   80

                                      ! ni =    3

                                      xin(80) = xin(98) + dxij*xin(71)
                                      yin(80) = yin(98) + dyij*yin(71)
                                      zin(80) = zin(98) + dzij*zin(71)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   35

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   99

                                      xin(108) = xin(108) + dxij*xin(99)
                                      yin(108) = yin(108) + dyij*yin(99)
                                      zin(108) = zin(108) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =   90

                                      xin(99) = xin(99) + dxij*xin(90)
                                      yin(99) = yin(99) + dyij*yin(90)
                                      zin(99) = zin(99) + dzij*zin(90)

                                      ! i3 = i4 =   90
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  108

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =   99

                                      xin(108) = xin(108) + dxij*xin(99)
                                      yin(108) = yin(108) + dyij*yin(99)
                                      zin(108) = zin(108) + dzij*zin(99)

                                      ! i3 = i4 =   99
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =   18

                                      ! do nj = 1,    2

                                      ! i4 = i3 =   18

                                      ! do ni = 1,    3

                                      xin(18) = xin(36) + dxij*xin(9)
                                      yin(18) = yin(36) + dyij*yin(9)
                                      zin(18) = zin(36) + dzij*zin(9)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   45

                                      ! ni =    2

                                      xin(45) = xin(63) + dxij*xin(36)
                                      yin(45) = yin(63) + dyij*yin(36)
                                      zin(45) = zin(63) + dzij*zin(36)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                      ! ni =    3

                                      xin(72) = xin(90) + dxij*xin(63)
                                      yin(72) = yin(90) + dyij*yin(63)
                                      zin(72) = zin(90) + dzij*zin(63)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   99

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   27

                                      ! nj =    2

                                      ! i4 = i3 =   27

                                      ! do ni = 1,    3

                                      xin(27) = xin(45) + dxij*xin(18)
                                      yin(27) = yin(45) + dyij*yin(18)
                                      zin(27) = zin(45) + dzij*zin(18)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   54

                                      ! ni =    2

                                      xin(54) = xin(72) + dxij*xin(45)
                                      yin(54) = yin(72) + dyij*yin(45)
                                      zin(54) = zin(72) + dzij*zin(45)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   81

                                      ! ni =    3

                                      xin(81) = xin(99) + dxij*xin(72)
                                      yin(81) = yin(99) + dyij*yin(72)
                                      zin(81) = zin(99) + dzij*zin(72)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =   36

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =    1

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    9

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    8

                                      xin(9) = xin(9) + dxkl*xin(8)
                                      yin(9) = yin(9) + dykl*yin(8)
                                      zin(9) = zin(9) + dzkl*zin(8)

                                      ! i3 = i4 =    8
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =    7

                                      xin(8) = xin(8) + dxkl*xin(7)
                                      yin(8) = yin(8) + dykl*yin(7)
                                      zin(8) = zin(8) + dzkl*zin(7)

                                      ! i3 = i4 =    7
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =    9

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =    8

                                      xin(9) = xin(9) + dxkl*xin(8)
                                      yin(9) = yin(9) + dykl*yin(8)
                                      zin(9) = zin(9) + dzkl*zin(8)

                                      ! i3 = i4 =    8
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =    2

                                      ! do nl = 1,    2

                                      ! i4 = i3 =    2

                                      ! do nk = 1,    2

                                      xin(2) = xin(4) + dxkl*xin(1)
                                      yin(2) = yin(4) + dykl*yin(1)
                                      zin(2) = zin(4) + dzkl*zin(1)
                                      ! i4 = i4 + lang+1 =    5

                                      ! nk =    2

                                      xin(5) = xin(7) + dxkl*xin(4)
                                      yin(5) = yin(7) + dykl*yin(4)
                                      zin(5) = zin(7) + dzkl*zin(4)
                                      ! i4 = i4 + lang+1 =    8

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =    3

                                      ! nl =    2

                                      ! i4 = i3 =    3

                                      ! do nk = 1,    2

                                      xin(3) = xin(5) + dxkl*xin(2)
                                      yin(3) = yin(5) + dykl*yin(2)
                                      zin(3) = zin(5) + dzkl*zin(2)
                                      ! i4 = i4 + lang+1 =    6

                                      ! nk =    2

                                      xin(6) = xin(8) + dxkl*xin(5)
                                      yin(6) = yin(8) + dykl*yin(5)
                                      zin(6) = zin(8) + dzkl*zin(5)
                                      ! i4 = i4 + lang+1 =    9

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =    4

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   10

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   18

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   17

                                      xin(18) = xin(18) + dxkl*xin(17)
                                      yin(18) = yin(18) + dykl*yin(17)
                                      zin(18) = zin(18) + dzkl*zin(17)

                                      ! i3 = i4 =   17
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   16

                                      xin(17) = xin(17) + dxkl*xin(16)
                                      yin(17) = yin(17) + dykl*yin(16)
                                      zin(17) = zin(17) + dzkl*zin(16)

                                      ! i3 = i4 =   16
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   18

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   17

                                      xin(18) = xin(18) + dxkl*xin(17)
                                      yin(18) = yin(18) + dykl*yin(17)
                                      zin(18) = zin(18) + dzkl*zin(17)

                                      ! i3 = i4 =   17
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   11

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   11

                                      ! do nk = 1,    2

                                      xin(11) = xin(13) + dxkl*xin(10)
                                      yin(11) = yin(13) + dykl*yin(10)
                                      zin(11) = zin(13) + dzkl*zin(10)
                                      ! i4 = i4 + lang+1 =   14

                                      ! nk =    2

                                      xin(14) = xin(16) + dxkl*xin(13)
                                      yin(14) = yin(16) + dykl*yin(13)
                                      zin(14) = zin(16) + dzkl*zin(13)
                                      ! i4 = i4 + lang+1 =   17

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   12

                                      ! nl =    2

                                      ! i4 = i3 =   12

                                      ! do nk = 1,    2

                                      xin(12) = xin(14) + dxkl*xin(11)
                                      yin(12) = yin(14) + dykl*yin(11)
                                      zin(12) = zin(14) + dzkl*zin(11)
                                      ! i4 = i4 + lang+1 =   15

                                      ! nk =    2

                                      xin(15) = xin(17) + dxkl*xin(14)
                                      yin(15) = yin(17) + dykl*yin(14)
                                      zin(15) = zin(17) + dzkl*zin(14)
                                      ! i4 = i4 + lang+1 =   18

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   13

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   19

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   27

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   26

                                      xin(27) = xin(27) + dxkl*xin(26)
                                      yin(27) = yin(27) + dykl*yin(26)
                                      zin(27) = zin(27) + dzkl*zin(26)

                                      ! i3 = i4 =   26
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   25

                                      xin(26) = xin(26) + dxkl*xin(25)
                                      yin(26) = yin(26) + dykl*yin(25)
                                      zin(26) = zin(26) + dzkl*zin(25)

                                      ! i3 = i4 =   25
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   27

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   26

                                      xin(27) = xin(27) + dxkl*xin(26)
                                      yin(27) = yin(27) + dykl*yin(26)
                                      zin(27) = zin(27) + dzkl*zin(26)

                                      ! i3 = i4 =   26
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   20

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   20

                                      ! do nk = 1,    2

                                      xin(20) = xin(22) + dxkl*xin(19)
                                      yin(20) = yin(22) + dykl*yin(19)
                                      zin(20) = zin(22) + dzkl*zin(19)
                                      ! i4 = i4 + lang+1 =   23

                                      ! nk =    2

                                      xin(23) = xin(25) + dxkl*xin(22)
                                      yin(23) = yin(25) + dykl*yin(22)
                                      zin(23) = zin(25) + dzkl*zin(22)
                                      ! i4 = i4 + lang+1 =   26

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   21

                                      ! nl =    2

                                      ! i4 = i3 =   21

                                      ! do nk = 1,    2

                                      xin(21) = xin(23) + dxkl*xin(20)
                                      yin(21) = yin(23) + dykl*yin(20)
                                      zin(21) = zin(23) + dzkl*zin(20)
                                      ! i4 = i4 + lang+1 =   24

                                      ! nk =    2

                                      xin(24) = xin(26) + dxkl*xin(23)
                                      yin(24) = yin(26) + dykl*yin(23)
                                      zin(24) = zin(26) + dzkl*zin(23)
                                      ! i4 = i4 + lang+1 =   27

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   22

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   28

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   28

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   34

                                      xin(35) = xin(35) + dxkl*xin(34)
                                      yin(35) = yin(35) + dykl*yin(34)
                                      zin(35) = zin(35) + dzkl*zin(34)

                                      ! i3 = i4 =   34
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   36

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   35

                                      xin(36) = xin(36) + dxkl*xin(35)
                                      yin(36) = yin(36) + dykl*yin(35)
                                      zin(36) = zin(36) + dzkl*zin(35)

                                      ! i3 = i4 =   35
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   29

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   29

                                      ! do nk = 1,    2

                                      xin(29) = xin(31) + dxkl*xin(28)
                                      yin(29) = yin(31) + dykl*yin(28)
                                      zin(29) = zin(31) + dzkl*zin(28)
                                      ! i4 = i4 + lang+1 =   32

                                      ! nk =    2

                                      xin(32) = xin(34) + dxkl*xin(31)
                                      yin(32) = yin(34) + dykl*yin(31)
                                      zin(32) = zin(34) + dzkl*zin(31)
                                      ! i4 = i4 + lang+1 =   35

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   30

                                      ! nl =    2

                                      ! i4 = i3 =   30

                                      ! do nk = 1,    2

                                      xin(30) = xin(32) + dxkl*xin(29)
                                      yin(30) = yin(32) + dykl*yin(29)
                                      zin(30) = zin(32) + dzkl*zin(29)
                                      ! i4 = i4 + lang+1 =   33

                                      ! nk =    2

                                      xin(33) = xin(35) + dxkl*xin(32)
                                      yin(33) = yin(35) + dykl*yin(32)
                                      zin(33) = zin(35) + dzkl*zin(32)
                                      ! i4 = i4 + lang+1 =   36

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   31

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   37

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   45

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   44

                                      xin(45) = xin(45) + dxkl*xin(44)
                                      yin(45) = yin(45) + dykl*yin(44)
                                      zin(45) = zin(45) + dzkl*zin(44)

                                      ! i3 = i4 =   44
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   43

                                      xin(44) = xin(44) + dxkl*xin(43)
                                      yin(44) = yin(44) + dykl*yin(43)
                                      zin(44) = zin(44) + dzkl*zin(43)

                                      ! i3 = i4 =   43
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   45

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   44

                                      xin(45) = xin(45) + dxkl*xin(44)
                                      yin(45) = yin(45) + dykl*yin(44)
                                      zin(45) = zin(45) + dzkl*zin(44)

                                      ! i3 = i4 =   44
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   38

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   38

                                      ! do nk = 1,    2

                                      xin(38) = xin(40) + dxkl*xin(37)
                                      yin(38) = yin(40) + dykl*yin(37)
                                      zin(38) = zin(40) + dzkl*zin(37)
                                      ! i4 = i4 + lang+1 =   41

                                      ! nk =    2

                                      xin(41) = xin(43) + dxkl*xin(40)
                                      yin(41) = yin(43) + dykl*yin(40)
                                      zin(41) = zin(43) + dzkl*zin(40)
                                      ! i4 = i4 + lang+1 =   44

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   39

                                      ! nl =    2

                                      ! i4 = i3 =   39

                                      ! do nk = 1,    2

                                      xin(39) = xin(41) + dxkl*xin(38)
                                      yin(39) = yin(41) + dykl*yin(38)
                                      zin(39) = zin(41) + dzkl*zin(38)
                                      ! i4 = i4 + lang+1 =   42

                                      ! nk =    2

                                      xin(42) = xin(44) + dxkl*xin(41)
                                      yin(42) = yin(44) + dykl*yin(41)
                                      zin(42) = zin(44) + dzkl*zin(41)
                                      ! i4 = i4 + lang+1 =   45

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   40

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   46

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   54

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   53

                                      xin(54) = xin(54) + dxkl*xin(53)
                                      yin(54) = yin(54) + dykl*yin(53)
                                      zin(54) = zin(54) + dzkl*zin(53)

                                      ! i3 = i4 =   53
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   52

                                      xin(53) = xin(53) + dxkl*xin(52)
                                      yin(53) = yin(53) + dykl*yin(52)
                                      zin(53) = zin(53) + dzkl*zin(52)

                                      ! i3 = i4 =   52
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   54

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   53

                                      xin(54) = xin(54) + dxkl*xin(53)
                                      yin(54) = yin(54) + dykl*yin(53)
                                      zin(54) = zin(54) + dzkl*zin(53)

                                      ! i3 = i4 =   53
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   47

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   47

                                      ! do nk = 1,    2

                                      xin(47) = xin(49) + dxkl*xin(46)
                                      yin(47) = yin(49) + dykl*yin(46)
                                      zin(47) = zin(49) + dzkl*zin(46)
                                      ! i4 = i4 + lang+1 =   50

                                      ! nk =    2

                                      xin(50) = xin(52) + dxkl*xin(49)
                                      yin(50) = yin(52) + dykl*yin(49)
                                      zin(50) = zin(52) + dzkl*zin(49)
                                      ! i4 = i4 + lang+1 =   53

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   48

                                      ! nl =    2

                                      ! i4 = i3 =   48

                                      ! do nk = 1,    2

                                      xin(48) = xin(50) + dxkl*xin(47)
                                      yin(48) = yin(50) + dykl*yin(47)
                                      zin(48) = zin(50) + dzkl*zin(47)
                                      ! i4 = i4 + lang+1 =   51

                                      ! nk =    2

                                      xin(51) = xin(53) + dxkl*xin(50)
                                      yin(51) = yin(53) + dykl*yin(50)
                                      zin(51) = zin(53) + dzkl*zin(50)
                                      ! i4 = i4 + lang+1 =   54

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   49

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   55

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   55

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   63

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   62

                                      xin(63) = xin(63) + dxkl*xin(62)
                                      yin(63) = yin(63) + dykl*yin(62)
                                      zin(63) = zin(63) + dzkl*zin(62)

                                      ! i3 = i4 =   62
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   61

                                      xin(62) = xin(62) + dxkl*xin(61)
                                      yin(62) = yin(62) + dykl*yin(61)
                                      zin(62) = zin(62) + dzkl*zin(61)

                                      ! i3 = i4 =   61
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   63

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   62

                                      xin(63) = xin(63) + dxkl*xin(62)
                                      yin(63) = yin(63) + dykl*yin(62)
                                      zin(63) = zin(63) + dzkl*zin(62)

                                      ! i3 = i4 =   62
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   56

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   56

                                      ! do nk = 1,    2

                                      xin(56) = xin(58) + dxkl*xin(55)
                                      yin(56) = yin(58) + dykl*yin(55)
                                      zin(56) = zin(58) + dzkl*zin(55)
                                      ! i4 = i4 + lang+1 =   59

                                      ! nk =    2

                                      xin(59) = xin(61) + dxkl*xin(58)
                                      yin(59) = yin(61) + dykl*yin(58)
                                      zin(59) = zin(61) + dzkl*zin(58)
                                      ! i4 = i4 + lang+1 =   62

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   57

                                      ! nl =    2

                                      ! i4 = i3 =   57

                                      ! do nk = 1,    2

                                      xin(57) = xin(59) + dxkl*xin(56)
                                      yin(57) = yin(59) + dykl*yin(56)
                                      zin(57) = zin(59) + dzkl*zin(56)
                                      ! i4 = i4 + lang+1 =   60

                                      ! nk =    2

                                      xin(60) = xin(62) + dxkl*xin(59)
                                      yin(60) = yin(62) + dykl*yin(59)
                                      zin(60) = zin(62) + dzkl*zin(59)
                                      ! i4 = i4 + lang+1 =   63

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   58

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   64

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   70

                                      xin(71) = xin(71) + dxkl*xin(70)
                                      yin(71) = yin(71) + dykl*yin(70)
                                      zin(71) = zin(71) + dzkl*zin(70)

                                      ! i3 = i4 =   70
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   72

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   71

                                      xin(72) = xin(72) + dxkl*xin(71)
                                      yin(72) = yin(72) + dykl*yin(71)
                                      zin(72) = zin(72) + dzkl*zin(71)

                                      ! i3 = i4 =   71
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   65

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   65

                                      ! do nk = 1,    2

                                      xin(65) = xin(67) + dxkl*xin(64)
                                      yin(65) = yin(67) + dykl*yin(64)
                                      zin(65) = zin(67) + dzkl*zin(64)
                                      ! i4 = i4 + lang+1 =   68

                                      ! nk =    2

                                      xin(68) = xin(70) + dxkl*xin(67)
                                      yin(68) = yin(70) + dykl*yin(67)
                                      zin(68) = zin(70) + dzkl*zin(67)
                                      ! i4 = i4 + lang+1 =   71

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   66

                                      ! nl =    2

                                      ! i4 = i3 =   66

                                      ! do nk = 1,    2

                                      xin(66) = xin(68) + dxkl*xin(65)
                                      yin(66) = yin(68) + dykl*yin(65)
                                      zin(66) = zin(68) + dzkl*zin(65)
                                      ! i4 = i4 + lang+1 =   69

                                      ! nk =    2

                                      xin(69) = xin(71) + dxkl*xin(68)
                                      yin(69) = yin(71) + dykl*yin(68)
                                      zin(69) = zin(71) + dzkl*zin(68)
                                      ! i4 = i4 + lang+1 =   72

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   67

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   73

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   81

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   80

                                      xin(81) = xin(81) + dxkl*xin(80)
                                      yin(81) = yin(81) + dykl*yin(80)
                                      zin(81) = zin(81) + dzkl*zin(80)

                                      ! i3 = i4 =   80
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   79

                                      xin(80) = xin(80) + dxkl*xin(79)
                                      yin(80) = yin(80) + dykl*yin(79)
                                      zin(80) = zin(80) + dzkl*zin(79)

                                      ! i3 = i4 =   79
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   81

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   80

                                      xin(81) = xin(81) + dxkl*xin(80)
                                      yin(81) = yin(81) + dykl*yin(80)
                                      zin(81) = zin(81) + dzkl*zin(80)

                                      ! i3 = i4 =   80
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   74

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   74

                                      ! do nk = 1,    2

                                      xin(74) = xin(76) + dxkl*xin(73)
                                      yin(74) = yin(76) + dykl*yin(73)
                                      zin(74) = zin(76) + dzkl*zin(73)
                                      ! i4 = i4 + lang+1 =   77

                                      ! nk =    2

                                      xin(77) = xin(79) + dxkl*xin(76)
                                      yin(77) = yin(79) + dykl*yin(76)
                                      zin(77) = zin(79) + dzkl*zin(76)
                                      ! i4 = i4 + lang+1 =   80

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   75

                                      ! nl =    2

                                      ! i4 = i3 =   75

                                      ! do nk = 1,    2

                                      xin(75) = xin(77) + dxkl*xin(74)
                                      yin(75) = yin(77) + dykl*yin(74)
                                      zin(75) = zin(77) + dzkl*zin(74)
                                      ! i4 = i4 + lang+1 =   78

                                      ! nk =    2

                                      xin(78) = xin(80) + dxkl*xin(77)
                                      yin(78) = yin(80) + dykl*yin(77)
                                      zin(78) = zin(80) + dzkl*zin(77)
                                      ! i4 = i4 + lang+1 =   81

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   76

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   82

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   82

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   90

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   89

                                      xin(90) = xin(90) + dxkl*xin(89)
                                      yin(90) = yin(90) + dykl*yin(89)
                                      zin(90) = zin(90) + dzkl*zin(89)

                                      ! i3 = i4 =   89
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   88

                                      xin(89) = xin(89) + dxkl*xin(88)
                                      yin(89) = yin(89) + dykl*yin(88)
                                      zin(89) = zin(89) + dzkl*zin(88)

                                      ! i3 = i4 =   88
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   90

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   89

                                      xin(90) = xin(90) + dxkl*xin(89)
                                      yin(90) = yin(90) + dykl*yin(89)
                                      zin(90) = zin(90) + dzkl*zin(89)

                                      ! i3 = i4 =   89
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   83

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   83

                                      ! do nk = 1,    2

                                      xin(83) = xin(85) + dxkl*xin(82)
                                      yin(83) = yin(85) + dykl*yin(82)
                                      zin(83) = zin(85) + dzkl*zin(82)
                                      ! i4 = i4 + lang+1 =   86

                                      ! nk =    2

                                      xin(86) = xin(88) + dxkl*xin(85)
                                      yin(86) = yin(88) + dykl*yin(85)
                                      zin(86) = zin(88) + dzkl*zin(85)
                                      ! i4 = i4 + lang+1 =   89

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   84

                                      ! nl =    2

                                      ! i4 = i3 =   84

                                      ! do nk = 1,    2

                                      xin(84) = xin(86) + dxkl*xin(83)
                                      yin(84) = yin(86) + dykl*yin(83)
                                      zin(84) = zin(86) + dzkl*zin(83)
                                      ! i4 = i4 + lang+1 =   87

                                      ! nk =    2

                                      xin(87) = xin(89) + dxkl*xin(86)
                                      yin(87) = yin(89) + dykl*yin(86)
                                      zin(87) = zin(89) + dzkl*zin(86)
                                      ! i4 = i4 + lang+1 =   90

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   85

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =   91

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   99

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   98

                                      xin(99) = xin(99) + dxkl*xin(98)
                                      yin(99) = yin(99) + dykl*yin(98)
                                      zin(99) = zin(99) + dzkl*zin(98)

                                      ! i3 = i4 =   98
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =   97

                                      xin(98) = xin(98) + dxkl*xin(97)
                                      yin(98) = yin(98) + dykl*yin(97)
                                      zin(98) = zin(98) + dzkl*zin(97)

                                      ! i3 = i4 =   97
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =   99

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =   98

                                      xin(99) = xin(99) + dxkl*xin(98)
                                      yin(99) = yin(99) + dykl*yin(98)
                                      zin(99) = zin(99) + dzkl*zin(98)

                                      ! i3 = i4 =   98
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =   92

                                      ! do nl = 1,    2

                                      ! i4 = i3 =   92

                                      ! do nk = 1,    2

                                      xin(92) = xin(94) + dxkl*xin(91)
                                      yin(92) = yin(94) + dykl*yin(91)
                                      zin(92) = zin(94) + dzkl*zin(91)
                                      ! i4 = i4 + lang+1 =   95

                                      ! nk =    2

                                      xin(95) = xin(97) + dxkl*xin(94)
                                      yin(95) = yin(97) + dykl*yin(94)
                                      zin(95) = zin(97) + dzkl*zin(94)
                                      ! i4 = i4 + lang+1 =   98

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   93

                                      ! nl =    2

                                      ! i4 = i3 =   93

                                      ! do nk = 1,    2

                                      xin(93) = xin(95) + dxkl*xin(92)
                                      yin(93) = yin(95) + dykl*yin(92)
                                      zin(93) = zin(95) + dzkl*zin(92)
                                      ! i4 = i4 + lang+1 =   96

                                      ! nk =    2

                                      xin(96) = xin(98) + dxkl*xin(95)
                                      yin(96) = yin(98) + dykl*yin(95)
                                      zin(96) = zin(98) + dzkl*zin(95)
                                      ! i4 = i4 + lang+1 =   99

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =   94

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  100

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  106

                                      xin(107) = xin(107) + dxkl*xin(106)
                                      yin(107) = yin(107) + dykl*yin(106)
                                      zin(107) = zin(107) + dzkl*zin(106)

                                      ! i3 = i4 =  106
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  108

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  107

                                      xin(108) = xin(108) + dxkl*xin(107)
                                      yin(108) = yin(108) + dykl*yin(107)
                                      zin(108) = zin(108) + dzkl*zin(107)

                                      ! i3 = i4 =  107
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  101

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  101

                                      ! do nk = 1,    2

                                      xin(101) = xin(103) + dxkl*xin(100)
                                      yin(101) = yin(103) + dykl*yin(100)
                                      zin(101) = zin(103) + dzkl*zin(100)
                                      ! i4 = i4 + lang+1 =  104

                                      ! nk =    2

                                      xin(104) = xin(106) + dxkl*xin(103)
                                      yin(104) = yin(106) + dykl*yin(103)
                                      zin(104) = zin(106) + dzkl*zin(103)
                                      ! i4 = i4 + lang+1 =  107

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  102

                                      ! nl =    2

                                      ! i4 = i3 =  102

                                      ! do nk = 1,    2

                                      xin(102) = xin(104) + dxkl*xin(101)
                                      yin(102) = yin(104) + dykl*yin(101)
                                      zin(102) = zin(104) + dzkl*zin(101)
                                      ! i4 = i4 + lang+1 =  105

                                      ! nk =    2

                                      xin(105) = xin(107) + dxkl*xin(104)
                                      yin(105) = yin(107) + dykl*yin(104)
                                      zin(105) = zin(107) + dzkl*zin(104)
                                      ! i4 = i4 + lang+1 =  108

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  103

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  109

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                      ! end do

                                      ! *** Now root =    2

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  108

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

                                      ! i1 = in(1) =  109

                                      xin(109) = 1.0_dp
                                      yin(109) = 1.0_dp
                                      zin(109) = f00

                                      ! i2 = in(2) =  136
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(136) = xc00
                                      yin(136) = yc00
                                      zin(136) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  112

                                      xin(112) = xcp00
                                      yin(112) = ycp00
                                      zin(112) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  139
                                      ! i2 =  136

                                      xin(139) = xcp00*xin(136) + cp10
                                      yin(139) = ycp00*yin(136) + cp10
                                      zin(139) = zcp00*zin(136) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  136

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  163
                                      ! i3 =  109
                                      ! i4 =  136

                                      xin(163) = c10*xin(109) + xc00*xin(136)
                                      yin(163) = c10*yin(109) + yc00*yin(136)
                                      zin(163) = c10*zin(109) + zc00*zin(136)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  166
                                      ! i5 =  163
                                      ! i4 =  136

                                      xin(166) = xcp00*xin(163) + cp10*xin(136)
                                      yin(166) = ycp00*yin(163) + cp10*yin(136)
                                      zin(166) = zcp00*zin(163) + cp10*zin(136)

                                      ! ------------------

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  163

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  190
                                      ! i3 =  136
                                      ! i4 =  163

                                      xin(190) = c10*xin(136) + xc00*xin(163)
                                      yin(190) = c10*yin(136) + yc00*yin(163)
                                      zin(190) = c10*zin(136) + zc00*zin(163)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  193
                                      ! i5 =  190
                                      ! i4 =  163

                                      xin(193) = xcp00*xin(190) + cp10*xin(163)
                                      yin(193) = ycp00*yin(190) + cp10*yin(163)
                                      zin(193) = zcp00*zin(190) + cp10*zin(163)

                                      ! ------------------

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  190

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  199
                                      ! i3 =  163
                                      ! i4 =  190

                                      xin(199) = c10*xin(163) + xc00*xin(190)
                                      yin(199) = c10*yin(163) + yc00*yin(190)
                                      zin(199) = c10*zin(163) + zc00*zin(190)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  202
                                      ! i5 =  199
                                      ! i4 =  190

                                      xin(202) = xcp00*xin(199) + cp10*xin(190)
                                      yin(202) = ycp00*yin(199) + cp10*yin(190)
                                      zin(202) = zcp00*zin(199) + cp10*zin(190)

                                      ! ------------------

                                      ! i3 = i4 =  190
                                      ! i4 = i5 =  199

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  208
                                      ! i3 =  190
                                      ! i4 =  199

                                      xin(208) = c10*xin(190) + xc00*xin(199)
                                      yin(208) = c10*yin(190) + yc00*yin(199)
                                      zin(208) = c10*zin(190) + zc00*zin(199)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  211
                                      ! i5 =  208
                                      ! i4 =  199

                                      xin(211) = xcp00*xin(208) + cp10*xin(199)
                                      yin(211) = ycp00*yin(208) + cp10*yin(199)
                                      zin(211) = zcp00*zin(208) + cp10*zin(199)

                                      ! ------------------

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  208

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  109
                                      ! i4 = i1+k2 =  112

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  115
                                      ! i3 =  109
                                      ! i4 =  112

                                      xin(115) = cp01*xin(109) + xcp00*xin(112)
                                      yin(115) = cp01*yin(109) + ycp00*yin(112)
                                      zin(115) = cp01*zin(109) + zcp00*zin(112)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  142

                                      xin(142) = xc00*xin(115) + c01*xin(112)
                                      yin(142) = yc00*yin(115) + c01*yin(112)
                                      zin(142) = zc00*zin(115) + c01*zin(112)

                                      ! ------------------

                                      ! i3 = i4 =  112
                                      ! i4 = i5 =  115

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  116
                                      ! i3 =  112
                                      ! i4 =  115

                                      xin(116) = cp01*xin(112) + xcp00*xin(115)
                                      yin(116) = cp01*yin(112) + ycp00*yin(115)
                                      zin(116) = cp01*zin(112) + zcp00*zin(115)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  143

                                      xin(143) = xc00*xin(116) + c01*xin(115)
                                      yin(143) = yc00*yin(116) + c01*yin(115)
                                      zin(143) = zc00*zin(116) + c01*zin(115)

                                      ! ------------------

                                      ! i3 = i4 =  115
                                      ! i4 = i5 =  116

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  117
                                      ! i3 =  115
                                      ! i4 =  116

                                      xin(117) = cp01*xin(115) + xcp00*xin(116)
                                      yin(117) = cp01*yin(115) + ycp00*yin(116)
                                      zin(117) = cp01*zin(115) + zcp00*zin(116)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  144

                                      xin(144) = xc00*xin(117) + c01*xin(116)
                                      yin(144) = yc00*yin(117) + c01*yin(116)
                                      zin(144) = zc00*zin(117) + c01*zin(116)

                                      ! ------------------

                                      ! i3 = i4 =  116
                                      ! i4 = i5 =  117

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  136

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  163

                                      xin(169) = c10*xin(115) + xc00*xin(142) + c01*xin(139)
                                      yin(169) = c10*yin(115) + yc00*yin(142) + c01*yin(139)
                                      zin(169) = c10*zin(115) + zc00*zin(142) + c01*zin(139)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  163

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  190

                                      xin(196) = c10*xin(142) + xc00*xin(169) + c01*xin(166)
                                      yin(196) = c10*yin(142) + yc00*yin(169) + c01*yin(166)
                                      zin(196) = c10*zin(142) + zc00*zin(169) + c01*zin(166)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  190

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  199

                                      xin(205) = c10*xin(169) + xc00*xin(196) + c01*xin(193)
                                      yin(205) = c10*yin(169) + yc00*yin(196) + c01*yin(193)
                                      zin(205) = c10*zin(169) + zc00*zin(196) + c01*zin(193)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  190
                                      ! i4 = i5 =  199

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  208

                                      xin(214) = c10*xin(196) + xc00*xin(205) + c01*xin(202)
                                      yin(214) = c10*yin(196) + yc00*yin(205) + c01*yin(202)
                                      zin(214) = c10*zin(196) + zc00*zin(205) + c01*zin(202)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  208

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  136

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  163

                                      xin(170) = c10*xin(116) + xc00*xin(143) + c01*xin(142)
                                      yin(170) = c10*yin(116) + yc00*yin(143) + c01*yin(142)
                                      zin(170) = c10*zin(116) + zc00*zin(143) + c01*zin(142)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  163

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  190

                                      xin(197) = c10*xin(143) + xc00*xin(170) + c01*xin(169)
                                      yin(197) = c10*yin(143) + yc00*yin(170) + c01*yin(169)
                                      zin(197) = c10*zin(143) + zc00*zin(170) + c01*zin(169)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  190

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  199

                                      xin(206) = c10*xin(170) + xc00*xin(197) + c01*xin(196)
                                      yin(206) = c10*yin(170) + yc00*yin(197) + c01*yin(196)
                                      zin(206) = c10*zin(170) + zc00*zin(197) + c01*zin(196)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  190
                                      ! i4 = i5 =  199

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  208

                                      xin(215) = c10*xin(197) + xc00*xin(206) + c01*xin(205)
                                      yin(215) = c10*yin(197) + yc00*yin(206) + c01*yin(205)
                                      zin(215) = c10*zin(197) + zc00*zin(206) + c01*zin(205)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  208

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  109
                                      ! i4 = i2 =  136

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  163

                                      xin(171) = c10*xin(117) + xc00*xin(144) + c01*xin(143)
                                      yin(171) = c10*yin(117) + yc00*yin(144) + c01*yin(143)
                                      zin(171) = c10*zin(117) + zc00*zin(144) + c01*zin(143)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  136
                                      ! i4 = i5 =  163

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  190

                                      xin(198) = c10*xin(144) + xc00*xin(171) + c01*xin(170)
                                      yin(198) = c10*yin(144) + yc00*yin(171) + c01*yin(170)
                                      zin(198) = c10*zin(144) + zc00*zin(171) + c01*zin(170)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  163
                                      ! i4 = i5 =  190

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  199

                                      xin(207) = c10*xin(171) + xc00*xin(198) + c01*xin(197)
                                      yin(207) = c10*yin(171) + yc00*yin(198) + c01*yin(197)
                                      zin(207) = c10*zin(171) + zc00*zin(198) + c01*zin(197)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  190
                                      ! i4 = i5 =  199

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  208

                                      xin(216) = c10*xin(198) + xc00*xin(207) + c01*xin(206)
                                      yin(216) = c10*yin(198) + yc00*yin(207) + c01*yin(206)
                                      zin(216) = c10*zin(198) + zc00*zin(207) + c01*zin(206)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  199
                                      ! i4 = i5 =  208

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  208

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  208

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  199

                                      xin(208) = xin(208) + dxij*xin(199)
                                      yin(208) = yin(208) + dyij*yin(199)
                                      zin(208) = zin(208) + dzij*zin(199)

                                      ! i3 = i4 =  199
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  190

                                      xin(199) = xin(199) + dxij*xin(190)
                                      yin(199) = yin(199) + dyij*yin(190)
                                      zin(199) = zin(199) + dzij*zin(190)

                                      ! i3 = i4 =  190
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  208

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  199

                                      xin(208) = xin(208) + dxij*xin(199)
                                      yin(208) = yin(208) + dyij*yin(199)
                                      zin(208) = zin(208) + dzij*zin(199)

                                      ! i3 = i4 =  199
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  118

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  118

                                      ! do ni = 1,    3

                                      xin(118) = xin(136) + dxij*xin(109)
                                      yin(118) = yin(136) + dyij*yin(109)
                                      zin(118) = zin(136) + dzij*zin(109)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  145

                                      ! ni =    2

                                      xin(145) = xin(163) + dxij*xin(136)
                                      yin(145) = yin(163) + dyij*yin(136)
                                      zin(145) = zin(163) + dzij*zin(136)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  172

                                      ! ni =    3

                                      xin(172) = xin(190) + dxij*xin(163)
                                      yin(172) = yin(190) + dyij*yin(163)
                                      zin(172) = zin(190) + dzij*zin(163)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  199

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  127

                                      ! nj =    2

                                      ! i4 = i3 =  127

                                      ! do ni = 1,    3

                                      xin(127) = xin(145) + dxij*xin(118)
                                      yin(127) = yin(145) + dyij*yin(118)
                                      zin(127) = zin(145) + dzij*zin(118)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  154

                                      ! ni =    2

                                      xin(154) = xin(172) + dxij*xin(145)
                                      yin(154) = yin(172) + dyij*yin(145)
                                      zin(154) = zin(172) + dzij*zin(145)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  181

                                      ! ni =    3

                                      xin(181) = xin(199) + dxij*xin(172)
                                      yin(181) = yin(199) + dyij*yin(172)
                                      zin(181) = zin(199) + dzij*zin(172)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  136

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  202

                                      xin(211) = xin(211) + dxij*xin(202)
                                      yin(211) = yin(211) + dyij*yin(202)
                                      zin(211) = zin(211) + dzij*zin(202)

                                      ! i3 = i4 =  202
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  193

                                      xin(202) = xin(202) + dxij*xin(193)
                                      yin(202) = yin(202) + dyij*yin(193)
                                      zin(202) = zin(202) + dzij*zin(193)

                                      ! i3 = i4 =  193
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  211

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  202

                                      xin(211) = xin(211) + dxij*xin(202)
                                      yin(211) = yin(211) + dyij*yin(202)
                                      zin(211) = zin(211) + dzij*zin(202)

                                      ! i3 = i4 =  202
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  121

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  121

                                      ! do ni = 1,    3

                                      xin(121) = xin(139) + dxij*xin(112)
                                      yin(121) = yin(139) + dyij*yin(112)
                                      zin(121) = zin(139) + dzij*zin(112)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  148

                                      ! ni =    2

                                      xin(148) = xin(166) + dxij*xin(139)
                                      yin(148) = yin(166) + dyij*yin(139)
                                      zin(148) = zin(166) + dzij*zin(139)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  175

                                      ! ni =    3

                                      xin(175) = xin(193) + dxij*xin(166)
                                      yin(175) = yin(193) + dyij*yin(166)
                                      zin(175) = zin(193) + dzij*zin(166)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  202

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  130

                                      ! nj =    2

                                      ! i4 = i3 =  130

                                      ! do ni = 1,    3

                                      xin(130) = xin(148) + dxij*xin(121)
                                      yin(130) = yin(148) + dyij*yin(121)
                                      zin(130) = zin(148) + dzij*zin(121)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  157

                                      ! ni =    2

                                      xin(157) = xin(175) + dxij*xin(148)
                                      yin(157) = yin(175) + dyij*yin(148)
                                      zin(157) = zin(175) + dzij*zin(148)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  184

                                      ! ni =    3

                                      xin(184) = xin(202) + dxij*xin(175)
                                      yin(184) = yin(202) + dyij*yin(175)
                                      zin(184) = zin(202) + dzij*zin(175)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  139

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  214

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(214) = xin(214) + dxij*xin(205)
                                      yin(214) = yin(214) + dyij*yin(205)
                                      zin(214) = zin(214) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  196

                                      xin(205) = xin(205) + dxij*xin(196)
                                      yin(205) = yin(205) + dyij*yin(196)
                                      zin(205) = zin(205) + dzij*zin(196)

                                      ! i3 = i4 =  196
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  214

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  205

                                      xin(214) = xin(214) + dxij*xin(205)
                                      yin(214) = yin(214) + dyij*yin(205)
                                      zin(214) = zin(214) + dzij*zin(205)

                                      ! i3 = i4 =  205
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  124

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  124

                                      ! do ni = 1,    3

                                      xin(124) = xin(142) + dxij*xin(115)
                                      yin(124) = yin(142) + dyij*yin(115)
                                      zin(124) = zin(142) + dzij*zin(115)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  151

                                      ! ni =    2

                                      xin(151) = xin(169) + dxij*xin(142)
                                      yin(151) = yin(169) + dyij*yin(142)
                                      zin(151) = zin(169) + dzij*zin(142)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  178

                                      ! ni =    3

                                      xin(178) = xin(196) + dxij*xin(169)
                                      yin(178) = yin(196) + dyij*yin(169)
                                      zin(178) = zin(196) + dzij*zin(169)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  133

                                      ! nj =    2

                                      ! i4 = i3 =  133

                                      ! do ni = 1,    3

                                      xin(133) = xin(151) + dxij*xin(124)
                                      yin(133) = yin(151) + dyij*yin(124)
                                      zin(133) = zin(151) + dzij*zin(124)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  160

                                      ! ni =    2

                                      xin(160) = xin(178) + dxij*xin(151)
                                      yin(160) = yin(178) + dyij*yin(151)
                                      zin(160) = zin(178) + dzij*zin(151)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  187

                                      ! ni =    3

                                      xin(187) = xin(205) + dxij*xin(178)
                                      yin(187) = yin(205) + dyij*yin(178)
                                      zin(187) = zin(205) + dzij*zin(178)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  142

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  206

                                      xin(215) = xin(215) + dxij*xin(206)
                                      yin(215) = yin(215) + dyij*yin(206)
                                      zin(215) = zin(215) + dzij*zin(206)

                                      ! i3 = i4 =  206
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  197

                                      xin(206) = xin(206) + dxij*xin(197)
                                      yin(206) = yin(206) + dyij*yin(197)
                                      zin(206) = zin(206) + dzij*zin(197)

                                      ! i3 = i4 =  197
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  215

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  206

                                      xin(215) = xin(215) + dxij*xin(206)
                                      yin(215) = yin(215) + dyij*yin(206)
                                      zin(215) = zin(215) + dzij*zin(206)

                                      ! i3 = i4 =  206
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  125

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  125

                                      ! do ni = 1,    3

                                      xin(125) = xin(143) + dxij*xin(116)
                                      yin(125) = yin(143) + dyij*yin(116)
                                      zin(125) = zin(143) + dzij*zin(116)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  152

                                      ! ni =    2

                                      xin(152) = xin(170) + dxij*xin(143)
                                      yin(152) = yin(170) + dyij*yin(143)
                                      zin(152) = zin(170) + dzij*zin(143)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  179

                                      ! ni =    3

                                      xin(179) = xin(197) + dxij*xin(170)
                                      yin(179) = yin(197) + dyij*yin(170)
                                      zin(179) = zin(197) + dzij*zin(170)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  206

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  134

                                      ! nj =    2

                                      ! i4 = i3 =  134

                                      ! do ni = 1,    3

                                      xin(134) = xin(152) + dxij*xin(125)
                                      yin(134) = yin(152) + dyij*yin(125)
                                      zin(134) = zin(152) + dzij*zin(125)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  161

                                      ! ni =    2

                                      xin(161) = xin(179) + dxij*xin(152)
                                      yin(161) = yin(179) + dyij*yin(152)
                                      zin(161) = zin(179) + dzij*zin(152)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  188

                                      ! ni =    3

                                      xin(188) = xin(206) + dxij*xin(179)
                                      yin(188) = yin(206) + dyij*yin(179)
                                      zin(188) = zin(206) + dzij*zin(179)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  143

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(216) = xin(216) + dxij*xin(207)
                                      yin(216) = yin(216) + dyij*yin(207)
                                      zin(216) = zin(216) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  198

                                      xin(207) = xin(207) + dxij*xin(198)
                                      yin(207) = yin(207) + dyij*yin(198)
                                      zin(207) = zin(207) + dzij*zin(198)

                                      ! i3 = i4 =  198
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  216

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  207

                                      xin(216) = xin(216) + dxij*xin(207)
                                      yin(216) = yin(216) + dyij*yin(207)
                                      zin(216) = zin(216) + dzij*zin(207)

                                      ! i3 = i4 =  207
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  126

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  126

                                      ! do ni = 1,    3

                                      xin(126) = xin(144) + dxij*xin(117)
                                      yin(126) = yin(144) + dyij*yin(117)
                                      zin(126) = zin(144) + dzij*zin(117)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  153

                                      ! ni =    2

                                      xin(153) = xin(171) + dxij*xin(144)
                                      yin(153) = yin(171) + dyij*yin(144)
                                      zin(153) = zin(171) + dzij*zin(144)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  180

                                      ! ni =    3

                                      xin(180) = xin(198) + dxij*xin(171)
                                      yin(180) = yin(198) + dyij*yin(171)
                                      zin(180) = zin(198) + dzij*zin(171)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  207

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  135

                                      ! nj =    2

                                      ! i4 = i3 =  135

                                      ! do ni = 1,    3

                                      xin(135) = xin(153) + dxij*xin(126)
                                      yin(135) = yin(153) + dyij*yin(126)
                                      zin(135) = zin(153) + dzij*zin(126)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  162

                                      ! ni =    2

                                      xin(162) = xin(180) + dxij*xin(153)
                                      yin(162) = yin(180) + dyij*yin(153)
                                      zin(162) = zin(180) + dzij*zin(153)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  189

                                      ! ni =    3

                                      xin(189) = xin(207) + dxij*xin(180)
                                      yin(189) = yin(207) + dyij*yin(180)
                                      zin(189) = zin(207) + dzij*zin(180)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  144

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  109

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  117

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  116

                                      xin(117) = xin(117) + dxkl*xin(116)
                                      yin(117) = yin(117) + dykl*yin(116)
                                      zin(117) = zin(117) + dzkl*zin(116)

                                      ! i3 = i4 =  116
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  115

                                      xin(116) = xin(116) + dxkl*xin(115)
                                      yin(116) = yin(116) + dykl*yin(115)
                                      zin(116) = zin(116) + dzkl*zin(115)

                                      ! i3 = i4 =  115
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  117

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  116

                                      xin(117) = xin(117) + dxkl*xin(116)
                                      yin(117) = yin(117) + dykl*yin(116)
                                      zin(117) = zin(117) + dzkl*zin(116)

                                      ! i3 = i4 =  116
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  110

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  110

                                      ! do nk = 1,    2

                                      xin(110) = xin(112) + dxkl*xin(109)
                                      yin(110) = yin(112) + dykl*yin(109)
                                      zin(110) = zin(112) + dzkl*zin(109)
                                      ! i4 = i4 + lang+1 =  113

                                      ! nk =    2

                                      xin(113) = xin(115) + dxkl*xin(112)
                                      yin(113) = yin(115) + dykl*yin(112)
                                      zin(113) = zin(115) + dzkl*zin(112)
                                      ! i4 = i4 + lang+1 =  116

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  111

                                      ! nl =    2

                                      ! i4 = i3 =  111

                                      ! do nk = 1,    2

                                      xin(111) = xin(113) + dxkl*xin(110)
                                      yin(111) = yin(113) + dykl*yin(110)
                                      zin(111) = zin(113) + dzkl*zin(110)
                                      ! i4 = i4 + lang+1 =  114

                                      ! nk =    2

                                      xin(114) = xin(116) + dxkl*xin(113)
                                      yin(114) = yin(116) + dykl*yin(113)
                                      zin(114) = zin(116) + dzkl*zin(113)
                                      ! i4 = i4 + lang+1 =  117

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  112

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  118

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  126

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  125

                                      xin(126) = xin(126) + dxkl*xin(125)
                                      yin(126) = yin(126) + dykl*yin(125)
                                      zin(126) = zin(126) + dzkl*zin(125)

                                      ! i3 = i4 =  125
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  124

                                      xin(125) = xin(125) + dxkl*xin(124)
                                      yin(125) = yin(125) + dykl*yin(124)
                                      zin(125) = zin(125) + dzkl*zin(124)

                                      ! i3 = i4 =  124
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  126

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  125

                                      xin(126) = xin(126) + dxkl*xin(125)
                                      yin(126) = yin(126) + dykl*yin(125)
                                      zin(126) = zin(126) + dzkl*zin(125)

                                      ! i3 = i4 =  125
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  119

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  119

                                      ! do nk = 1,    2

                                      xin(119) = xin(121) + dxkl*xin(118)
                                      yin(119) = yin(121) + dykl*yin(118)
                                      zin(119) = zin(121) + dzkl*zin(118)
                                      ! i4 = i4 + lang+1 =  122

                                      ! nk =    2

                                      xin(122) = xin(124) + dxkl*xin(121)
                                      yin(122) = yin(124) + dykl*yin(121)
                                      zin(122) = zin(124) + dzkl*zin(121)
                                      ! i4 = i4 + lang+1 =  125

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  120

                                      ! nl =    2

                                      ! i4 = i3 =  120

                                      ! do nk = 1,    2

                                      xin(120) = xin(122) + dxkl*xin(119)
                                      yin(120) = yin(122) + dykl*yin(119)
                                      zin(120) = zin(122) + dzkl*zin(119)
                                      ! i4 = i4 + lang+1 =  123

                                      ! nk =    2

                                      xin(123) = xin(125) + dxkl*xin(122)
                                      yin(123) = yin(125) + dykl*yin(122)
                                      zin(123) = zin(125) + dzkl*zin(122)
                                      ! i4 = i4 + lang+1 =  126

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  121

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  127

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  135

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  134

                                      xin(135) = xin(135) + dxkl*xin(134)
                                      yin(135) = yin(135) + dykl*yin(134)
                                      zin(135) = zin(135) + dzkl*zin(134)

                                      ! i3 = i4 =  134
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  133

                                      xin(134) = xin(134) + dxkl*xin(133)
                                      yin(134) = yin(134) + dykl*yin(133)
                                      zin(134) = zin(134) + dzkl*zin(133)

                                      ! i3 = i4 =  133
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  135

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  134

                                      xin(135) = xin(135) + dxkl*xin(134)
                                      yin(135) = yin(135) + dykl*yin(134)
                                      zin(135) = zin(135) + dzkl*zin(134)

                                      ! i3 = i4 =  134
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  128

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  128

                                      ! do nk = 1,    2

                                      xin(128) = xin(130) + dxkl*xin(127)
                                      yin(128) = yin(130) + dykl*yin(127)
                                      zin(128) = zin(130) + dzkl*zin(127)
                                      ! i4 = i4 + lang+1 =  131

                                      ! nk =    2

                                      xin(131) = xin(133) + dxkl*xin(130)
                                      yin(131) = yin(133) + dykl*yin(130)
                                      zin(131) = zin(133) + dzkl*zin(130)
                                      ! i4 = i4 + lang+1 =  134

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  129

                                      ! nl =    2

                                      ! i4 = i3 =  129

                                      ! do nk = 1,    2

                                      xin(129) = xin(131) + dxkl*xin(128)
                                      yin(129) = yin(131) + dykl*yin(128)
                                      zin(129) = zin(131) + dzkl*zin(128)
                                      ! i4 = i4 + lang+1 =  132

                                      ! nk =    2

                                      xin(132) = xin(134) + dxkl*xin(131)
                                      yin(132) = yin(134) + dykl*yin(131)
                                      zin(132) = zin(134) + dzkl*zin(131)
                                      ! i4 = i4 + lang+1 =  135

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  130

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  136

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  136

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  142

                                      xin(143) = xin(143) + dxkl*xin(142)
                                      yin(143) = yin(143) + dykl*yin(142)
                                      zin(143) = zin(143) + dzkl*zin(142)

                                      ! i3 = i4 =  142
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  144

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  143

                                      xin(144) = xin(144) + dxkl*xin(143)
                                      yin(144) = yin(144) + dykl*yin(143)
                                      zin(144) = zin(144) + dzkl*zin(143)

                                      ! i3 = i4 =  143
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  137

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  137

                                      ! do nk = 1,    2

                                      xin(137) = xin(139) + dxkl*xin(136)
                                      yin(137) = yin(139) + dykl*yin(136)
                                      zin(137) = zin(139) + dzkl*zin(136)
                                      ! i4 = i4 + lang+1 =  140

                                      ! nk =    2

                                      xin(140) = xin(142) + dxkl*xin(139)
                                      yin(140) = yin(142) + dykl*yin(139)
                                      zin(140) = zin(142) + dzkl*zin(139)
                                      ! i4 = i4 + lang+1 =  143

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  138

                                      ! nl =    2

                                      ! i4 = i3 =  138

                                      ! do nk = 1,    2

                                      xin(138) = xin(140) + dxkl*xin(137)
                                      yin(138) = yin(140) + dykl*yin(137)
                                      zin(138) = zin(140) + dzkl*zin(137)
                                      ! i4 = i4 + lang+1 =  141

                                      ! nk =    2

                                      xin(141) = xin(143) + dxkl*xin(140)
                                      yin(141) = yin(143) + dykl*yin(140)
                                      zin(141) = zin(143) + dzkl*zin(140)
                                      ! i4 = i4 + lang+1 =  144

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  139

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  145

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  153

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  152

                                      xin(153) = xin(153) + dxkl*xin(152)
                                      yin(153) = yin(153) + dykl*yin(152)
                                      zin(153) = zin(153) + dzkl*zin(152)

                                      ! i3 = i4 =  152
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  151

                                      xin(152) = xin(152) + dxkl*xin(151)
                                      yin(152) = yin(152) + dykl*yin(151)
                                      zin(152) = zin(152) + dzkl*zin(151)

                                      ! i3 = i4 =  151
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  153

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  152

                                      xin(153) = xin(153) + dxkl*xin(152)
                                      yin(153) = yin(153) + dykl*yin(152)
                                      zin(153) = zin(153) + dzkl*zin(152)

                                      ! i3 = i4 =  152
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  146

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  146

                                      ! do nk = 1,    2

                                      xin(146) = xin(148) + dxkl*xin(145)
                                      yin(146) = yin(148) + dykl*yin(145)
                                      zin(146) = zin(148) + dzkl*zin(145)
                                      ! i4 = i4 + lang+1 =  149

                                      ! nk =    2

                                      xin(149) = xin(151) + dxkl*xin(148)
                                      yin(149) = yin(151) + dykl*yin(148)
                                      zin(149) = zin(151) + dzkl*zin(148)
                                      ! i4 = i4 + lang+1 =  152

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  147

                                      ! nl =    2

                                      ! i4 = i3 =  147

                                      ! do nk = 1,    2

                                      xin(147) = xin(149) + dxkl*xin(146)
                                      yin(147) = yin(149) + dykl*yin(146)
                                      zin(147) = zin(149) + dzkl*zin(146)
                                      ! i4 = i4 + lang+1 =  150

                                      ! nk =    2

                                      xin(150) = xin(152) + dxkl*xin(149)
                                      yin(150) = yin(152) + dykl*yin(149)
                                      zin(150) = zin(152) + dzkl*zin(149)
                                      ! i4 = i4 + lang+1 =  153

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  148

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  154

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  162

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  161

                                      xin(162) = xin(162) + dxkl*xin(161)
                                      yin(162) = yin(162) + dykl*yin(161)
                                      zin(162) = zin(162) + dzkl*zin(161)

                                      ! i3 = i4 =  161
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  160

                                      xin(161) = xin(161) + dxkl*xin(160)
                                      yin(161) = yin(161) + dykl*yin(160)
                                      zin(161) = zin(161) + dzkl*zin(160)

                                      ! i3 = i4 =  160
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  162

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  161

                                      xin(162) = xin(162) + dxkl*xin(161)
                                      yin(162) = yin(162) + dykl*yin(161)
                                      zin(162) = zin(162) + dzkl*zin(161)

                                      ! i3 = i4 =  161
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  155

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  155

                                      ! do nk = 1,    2

                                      xin(155) = xin(157) + dxkl*xin(154)
                                      yin(155) = yin(157) + dykl*yin(154)
                                      zin(155) = zin(157) + dzkl*zin(154)
                                      ! i4 = i4 + lang+1 =  158

                                      ! nk =    2

                                      xin(158) = xin(160) + dxkl*xin(157)
                                      yin(158) = yin(160) + dykl*yin(157)
                                      zin(158) = zin(160) + dzkl*zin(157)
                                      ! i4 = i4 + lang+1 =  161

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  156

                                      ! nl =    2

                                      ! i4 = i3 =  156

                                      ! do nk = 1,    2

                                      xin(156) = xin(158) + dxkl*xin(155)
                                      yin(156) = yin(158) + dykl*yin(155)
                                      zin(156) = zin(158) + dzkl*zin(155)
                                      ! i4 = i4 + lang+1 =  159

                                      ! nk =    2

                                      xin(159) = xin(161) + dxkl*xin(158)
                                      yin(159) = yin(161) + dykl*yin(158)
                                      zin(159) = zin(161) + dzkl*zin(158)
                                      ! i4 = i4 + lang+1 =  162

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  157

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  163

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  163

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  171

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  170

                                      xin(171) = xin(171) + dxkl*xin(170)
                                      yin(171) = yin(171) + dykl*yin(170)
                                      zin(171) = zin(171) + dzkl*zin(170)

                                      ! i3 = i4 =  170
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  169

                                      xin(170) = xin(170) + dxkl*xin(169)
                                      yin(170) = yin(170) + dykl*yin(169)
                                      zin(170) = zin(170) + dzkl*zin(169)

                                      ! i3 = i4 =  169
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  171

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  170

                                      xin(171) = xin(171) + dxkl*xin(170)
                                      yin(171) = yin(171) + dykl*yin(170)
                                      zin(171) = zin(171) + dzkl*zin(170)

                                      ! i3 = i4 =  170
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  164

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  164

                                      ! do nk = 1,    2

                                      xin(164) = xin(166) + dxkl*xin(163)
                                      yin(164) = yin(166) + dykl*yin(163)
                                      zin(164) = zin(166) + dzkl*zin(163)
                                      ! i4 = i4 + lang+1 =  167

                                      ! nk =    2

                                      xin(167) = xin(169) + dxkl*xin(166)
                                      yin(167) = yin(169) + dykl*yin(166)
                                      zin(167) = zin(169) + dzkl*zin(166)
                                      ! i4 = i4 + lang+1 =  170

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  165

                                      ! nl =    2

                                      ! i4 = i3 =  165

                                      ! do nk = 1,    2

                                      xin(165) = xin(167) + dxkl*xin(164)
                                      yin(165) = yin(167) + dykl*yin(164)
                                      zin(165) = zin(167) + dzkl*zin(164)
                                      ! i4 = i4 + lang+1 =  168

                                      ! nk =    2

                                      xin(168) = xin(170) + dxkl*xin(167)
                                      yin(168) = yin(170) + dykl*yin(167)
                                      zin(168) = zin(170) + dzkl*zin(167)
                                      ! i4 = i4 + lang+1 =  171

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  166

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  172

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  178

                                      xin(179) = xin(179) + dxkl*xin(178)
                                      yin(179) = yin(179) + dykl*yin(178)
                                      zin(179) = zin(179) + dzkl*zin(178)

                                      ! i3 = i4 =  178
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  180

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  179

                                      xin(180) = xin(180) + dxkl*xin(179)
                                      yin(180) = yin(180) + dykl*yin(179)
                                      zin(180) = zin(180) + dzkl*zin(179)

                                      ! i3 = i4 =  179
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  173

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  173

                                      ! do nk = 1,    2

                                      xin(173) = xin(175) + dxkl*xin(172)
                                      yin(173) = yin(175) + dykl*yin(172)
                                      zin(173) = zin(175) + dzkl*zin(172)
                                      ! i4 = i4 + lang+1 =  176

                                      ! nk =    2

                                      xin(176) = xin(178) + dxkl*xin(175)
                                      yin(176) = yin(178) + dykl*yin(175)
                                      zin(176) = zin(178) + dzkl*zin(175)
                                      ! i4 = i4 + lang+1 =  179

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  174

                                      ! nl =    2

                                      ! i4 = i3 =  174

                                      ! do nk = 1,    2

                                      xin(174) = xin(176) + dxkl*xin(173)
                                      yin(174) = yin(176) + dykl*yin(173)
                                      zin(174) = zin(176) + dzkl*zin(173)
                                      ! i4 = i4 + lang+1 =  177

                                      ! nk =    2

                                      xin(177) = xin(179) + dxkl*xin(176)
                                      yin(177) = yin(179) + dykl*yin(176)
                                      zin(177) = zin(179) + dzkl*zin(176)
                                      ! i4 = i4 + lang+1 =  180

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  175

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  181

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  189

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  188

                                      xin(189) = xin(189) + dxkl*xin(188)
                                      yin(189) = yin(189) + dykl*yin(188)
                                      zin(189) = zin(189) + dzkl*zin(188)

                                      ! i3 = i4 =  188
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  187

                                      xin(188) = xin(188) + dxkl*xin(187)
                                      yin(188) = yin(188) + dykl*yin(187)
                                      zin(188) = zin(188) + dzkl*zin(187)

                                      ! i3 = i4 =  187
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  189

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  188

                                      xin(189) = xin(189) + dxkl*xin(188)
                                      yin(189) = yin(189) + dykl*yin(188)
                                      zin(189) = zin(189) + dzkl*zin(188)

                                      ! i3 = i4 =  188
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  182

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  182

                                      ! do nk = 1,    2

                                      xin(182) = xin(184) + dxkl*xin(181)
                                      yin(182) = yin(184) + dykl*yin(181)
                                      zin(182) = zin(184) + dzkl*zin(181)
                                      ! i4 = i4 + lang+1 =  185

                                      ! nk =    2

                                      xin(185) = xin(187) + dxkl*xin(184)
                                      yin(185) = yin(187) + dykl*yin(184)
                                      zin(185) = zin(187) + dzkl*zin(184)
                                      ! i4 = i4 + lang+1 =  188

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  183

                                      ! nl =    2

                                      ! i4 = i3 =  183

                                      ! do nk = 1,    2

                                      xin(183) = xin(185) + dxkl*xin(182)
                                      yin(183) = yin(185) + dykl*yin(182)
                                      zin(183) = zin(185) + dzkl*zin(182)
                                      ! i4 = i4 + lang+1 =  186

                                      ! nk =    2

                                      xin(186) = xin(188) + dxkl*xin(185)
                                      yin(186) = yin(188) + dykl*yin(185)
                                      zin(186) = zin(188) + dzkl*zin(185)
                                      ! i4 = i4 + lang+1 =  189

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  184

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  190

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  190

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  198

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  197

                                      xin(198) = xin(198) + dxkl*xin(197)
                                      yin(198) = yin(198) + dykl*yin(197)
                                      zin(198) = zin(198) + dzkl*zin(197)

                                      ! i3 = i4 =  197
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  196

                                      xin(197) = xin(197) + dxkl*xin(196)
                                      yin(197) = yin(197) + dykl*yin(196)
                                      zin(197) = zin(197) + dzkl*zin(196)

                                      ! i3 = i4 =  196
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  198

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  197

                                      xin(198) = xin(198) + dxkl*xin(197)
                                      yin(198) = yin(198) + dykl*yin(197)
                                      zin(198) = zin(198) + dzkl*zin(197)

                                      ! i3 = i4 =  197
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  191

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  191

                                      ! do nk = 1,    2

                                      xin(191) = xin(193) + dxkl*xin(190)
                                      yin(191) = yin(193) + dykl*yin(190)
                                      zin(191) = zin(193) + dzkl*zin(190)
                                      ! i4 = i4 + lang+1 =  194

                                      ! nk =    2

                                      xin(194) = xin(196) + dxkl*xin(193)
                                      yin(194) = yin(196) + dykl*yin(193)
                                      zin(194) = zin(196) + dzkl*zin(193)
                                      ! i4 = i4 + lang+1 =  197

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  192

                                      ! nl =    2

                                      ! i4 = i3 =  192

                                      ! do nk = 1,    2

                                      xin(192) = xin(194) + dxkl*xin(191)
                                      yin(192) = yin(194) + dykl*yin(191)
                                      zin(192) = zin(194) + dzkl*zin(191)
                                      ! i4 = i4 + lang+1 =  195

                                      ! nk =    2

                                      xin(195) = xin(197) + dxkl*xin(194)
                                      yin(195) = yin(197) + dykl*yin(194)
                                      zin(195) = zin(197) + dzkl*zin(194)
                                      ! i4 = i4 + lang+1 =  198

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  193

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  199

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  207

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  206

                                      xin(207) = xin(207) + dxkl*xin(206)
                                      yin(207) = yin(207) + dykl*yin(206)
                                      zin(207) = zin(207) + dzkl*zin(206)

                                      ! i3 = i4 =  206
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  205

                                      xin(206) = xin(206) + dxkl*xin(205)
                                      yin(206) = yin(206) + dykl*yin(205)
                                      zin(206) = zin(206) + dzkl*zin(205)

                                      ! i3 = i4 =  205
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  207

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  206

                                      xin(207) = xin(207) + dxkl*xin(206)
                                      yin(207) = yin(207) + dykl*yin(206)
                                      zin(207) = zin(207) + dzkl*zin(206)

                                      ! i3 = i4 =  206
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  200

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  200

                                      ! do nk = 1,    2

                                      xin(200) = xin(202) + dxkl*xin(199)
                                      yin(200) = yin(202) + dykl*yin(199)
                                      zin(200) = zin(202) + dzkl*zin(199)
                                      ! i4 = i4 + lang+1 =  203

                                      ! nk =    2

                                      xin(203) = xin(205) + dxkl*xin(202)
                                      yin(203) = yin(205) + dykl*yin(202)
                                      zin(203) = zin(205) + dzkl*zin(202)
                                      ! i4 = i4 + lang+1 =  206

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  201

                                      ! nl =    2

                                      ! i4 = i3 =  201

                                      ! do nk = 1,    2

                                      xin(201) = xin(203) + dxkl*xin(200)
                                      yin(201) = yin(203) + dykl*yin(200)
                                      zin(201) = zin(203) + dzkl*zin(200)
                                      ! i4 = i4 + lang+1 =  204

                                      ! nk =    2

                                      xin(204) = xin(206) + dxkl*xin(203)
                                      yin(204) = yin(206) + dykl*yin(203)
                                      zin(204) = zin(206) + dzkl*zin(203)
                                      ! i4 = i4 + lang+1 =  207

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  202

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  208

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  214

                                      xin(215) = xin(215) + dxkl*xin(214)
                                      yin(215) = yin(215) + dykl*yin(214)
                                      zin(215) = zin(215) + dzkl*zin(214)

                                      ! i3 = i4 =  214
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  216

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  215

                                      xin(216) = xin(216) + dxkl*xin(215)
                                      yin(216) = yin(216) + dykl*yin(215)
                                      zin(216) = zin(216) + dzkl*zin(215)

                                      ! i3 = i4 =  215
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  209

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  209

                                      ! do nk = 1,    2

                                      xin(209) = xin(211) + dxkl*xin(208)
                                      yin(209) = yin(211) + dykl*yin(208)
                                      zin(209) = zin(211) + dzkl*zin(208)
                                      ! i4 = i4 + lang+1 =  212

                                      ! nk =    2

                                      xin(212) = xin(214) + dxkl*xin(211)
                                      yin(212) = yin(214) + dykl*yin(211)
                                      zin(212) = zin(214) + dzkl*zin(211)
                                      ! i4 = i4 + lang+1 =  215

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  210

                                      ! nl =    2

                                      ! i4 = i3 =  210

                                      ! do nk = 1,    2

                                      xin(210) = xin(212) + dxkl*xin(209)
                                      yin(210) = yin(212) + dykl*yin(209)
                                      zin(210) = zin(212) + dzkl*zin(209)
                                      ! i4 = i4 + lang+1 =  213

                                      ! nk =    2

                                      xin(213) = xin(215) + dxkl*xin(212)
                                      yin(213) = yin(215) + dykl*yin(212)
                                      zin(213) = zin(215) + dzkl*zin(212)
                                      ! i4 = i4 + lang+1 =  216

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  211

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  217

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                      ! end do

                                      ! *** Now root =    3

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  216

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

                                      ! i1 = in(1) =  217

                                      xin(217) = 1.0_dp
                                      yin(217) = 1.0_dp
                                      zin(217) = f00

                                      ! i2 = in(2) =  244
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(244) = xc00
                                      yin(244) = yc00
                                      zin(244) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  220

                                      xin(220) = xcp00
                                      yin(220) = ycp00
                                      zin(220) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  247
                                      ! i2 =  244

                                      xin(247) = xcp00*xin(244) + cp10
                                      yin(247) = ycp00*yin(244) + cp10
                                      zin(247) = zcp00*zin(244) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  244

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  271
                                      ! i3 =  217
                                      ! i4 =  244

                                      xin(271) = c10*xin(217) + xc00*xin(244)
                                      yin(271) = c10*yin(217) + yc00*yin(244)
                                      zin(271) = c10*zin(217) + zc00*zin(244)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  274
                                      ! i5 =  271
                                      ! i4 =  244

                                      xin(274) = xcp00*xin(271) + cp10*xin(244)
                                      yin(274) = ycp00*yin(271) + cp10*yin(244)
                                      zin(274) = zcp00*zin(271) + cp10*zin(244)

                                      ! ------------------

                                      ! i3 = i4 =  244
                                      ! i4 = i5 =  271

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  298
                                      ! i3 =  244
                                      ! i4 =  271

                                      xin(298) = c10*xin(244) + xc00*xin(271)
                                      yin(298) = c10*yin(244) + yc00*yin(271)
                                      zin(298) = c10*zin(244) + zc00*zin(271)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  301
                                      ! i5 =  298
                                      ! i4 =  271

                                      xin(301) = xcp00*xin(298) + cp10*xin(271)
                                      yin(301) = ycp00*yin(298) + cp10*yin(271)
                                      zin(301) = zcp00*zin(298) + cp10*zin(271)

                                      ! ------------------

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  298

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  307
                                      ! i3 =  271
                                      ! i4 =  298

                                      xin(307) = c10*xin(271) + xc00*xin(298)
                                      yin(307) = c10*yin(271) + yc00*yin(298)
                                      zin(307) = c10*zin(271) + zc00*zin(298)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  310
                                      ! i5 =  307
                                      ! i4 =  298

                                      xin(310) = xcp00*xin(307) + cp10*xin(298)
                                      yin(310) = ycp00*yin(307) + cp10*yin(298)
                                      zin(310) = zcp00*zin(307) + cp10*zin(298)

                                      ! ------------------

                                      ! i3 = i4 =  298
                                      ! i4 = i5 =  307

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  316
                                      ! i3 =  298
                                      ! i4 =  307

                                      xin(316) = c10*xin(298) + xc00*xin(307)
                                      yin(316) = c10*yin(298) + yc00*yin(307)
                                      zin(316) = c10*zin(298) + zc00*zin(307)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  319
                                      ! i5 =  316
                                      ! i4 =  307

                                      xin(319) = xcp00*xin(316) + cp10*xin(307)
                                      yin(319) = ycp00*yin(316) + cp10*yin(307)
                                      zin(319) = zcp00*zin(316) + cp10*zin(307)

                                      ! ------------------

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  316

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  217
                                      ! i4 = i1+k2 =  220

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  223
                                      ! i3 =  217
                                      ! i4 =  220

                                      xin(223) = cp01*xin(217) + xcp00*xin(220)
                                      yin(223) = cp01*yin(217) + ycp00*yin(220)
                                      zin(223) = cp01*zin(217) + zcp00*zin(220)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  250

                                      xin(250) = xc00*xin(223) + c01*xin(220)
                                      yin(250) = yc00*yin(223) + c01*yin(220)
                                      zin(250) = zc00*zin(223) + c01*zin(220)

                                      ! ------------------

                                      ! i3 = i4 =  220
                                      ! i4 = i5 =  223

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  224
                                      ! i3 =  220
                                      ! i4 =  223

                                      xin(224) = cp01*xin(220) + xcp00*xin(223)
                                      yin(224) = cp01*yin(220) + ycp00*yin(223)
                                      zin(224) = cp01*zin(220) + zcp00*zin(223)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  251

                                      xin(251) = xc00*xin(224) + c01*xin(223)
                                      yin(251) = yc00*yin(224) + c01*yin(223)
                                      zin(251) = zc00*zin(224) + c01*zin(223)

                                      ! ------------------

                                      ! i3 = i4 =  223
                                      ! i4 = i5 =  224

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  225
                                      ! i3 =  223
                                      ! i4 =  224

                                      xin(225) = cp01*xin(223) + xcp00*xin(224)
                                      yin(225) = cp01*yin(223) + ycp00*yin(224)
                                      zin(225) = cp01*zin(223) + zcp00*zin(224)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  252

                                      xin(252) = xc00*xin(225) + c01*xin(224)
                                      yin(252) = yc00*yin(225) + c01*yin(224)
                                      zin(252) = zc00*zin(225) + c01*zin(224)

                                      ! ------------------

                                      ! i3 = i4 =  224
                                      ! i4 = i5 =  225

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  244

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  271

                                      xin(277) = c10*xin(223) + xc00*xin(250) + c01*xin(247)
                                      yin(277) = c10*yin(223) + yc00*yin(250) + c01*yin(247)
                                      zin(277) = c10*zin(223) + zc00*zin(250) + c01*zin(247)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  244
                                      ! i4 = i5 =  271

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  298

                                      xin(304) = c10*xin(250) + xc00*xin(277) + c01*xin(274)
                                      yin(304) = c10*yin(250) + yc00*yin(277) + c01*yin(274)
                                      zin(304) = c10*zin(250) + zc00*zin(277) + c01*zin(274)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  298

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  307

                                      xin(313) = c10*xin(277) + xc00*xin(304) + c01*xin(301)
                                      yin(313) = c10*yin(277) + yc00*yin(304) + c01*yin(301)
                                      zin(313) = c10*zin(277) + zc00*zin(304) + c01*zin(301)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  298
                                      ! i4 = i5 =  307

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  316

                                      xin(322) = c10*xin(304) + xc00*xin(313) + c01*xin(310)
                                      yin(322) = c10*yin(304) + yc00*yin(313) + c01*yin(310)
                                      zin(322) = c10*zin(304) + zc00*zin(313) + c01*zin(310)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  316

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  244

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  271

                                      xin(278) = c10*xin(224) + xc00*xin(251) + c01*xin(250)
                                      yin(278) = c10*yin(224) + yc00*yin(251) + c01*yin(250)
                                      zin(278) = c10*zin(224) + zc00*zin(251) + c01*zin(250)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  244
                                      ! i4 = i5 =  271

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  298

                                      xin(305) = c10*xin(251) + xc00*xin(278) + c01*xin(277)
                                      yin(305) = c10*yin(251) + yc00*yin(278) + c01*yin(277)
                                      zin(305) = c10*zin(251) + zc00*zin(278) + c01*zin(277)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  298

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  307

                                      xin(314) = c10*xin(278) + xc00*xin(305) + c01*xin(304)
                                      yin(314) = c10*yin(278) + yc00*yin(305) + c01*yin(304)
                                      zin(314) = c10*zin(278) + zc00*zin(305) + c01*zin(304)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  298
                                      ! i4 = i5 =  307

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  316

                                      xin(323) = c10*xin(305) + xc00*xin(314) + c01*xin(313)
                                      yin(323) = c10*yin(305) + yc00*yin(314) + c01*yin(313)
                                      zin(323) = c10*zin(305) + zc00*zin(314) + c01*zin(313)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  316

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  217
                                      ! i4 = i2 =  244

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  271

                                      xin(279) = c10*xin(225) + xc00*xin(252) + c01*xin(251)
                                      yin(279) = c10*yin(225) + yc00*yin(252) + c01*yin(251)
                                      zin(279) = c10*zin(225) + zc00*zin(252) + c01*zin(251)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  244
                                      ! i4 = i5 =  271

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  298

                                      xin(306) = c10*xin(252) + xc00*xin(279) + c01*xin(278)
                                      yin(306) = c10*yin(252) + yc00*yin(279) + c01*yin(278)
                                      zin(306) = c10*zin(252) + zc00*zin(279) + c01*zin(278)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  271
                                      ! i4 = i5 =  298

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  307

                                      xin(315) = c10*xin(279) + xc00*xin(306) + c01*xin(305)
                                      yin(315) = c10*yin(279) + yc00*yin(306) + c01*yin(305)
                                      zin(315) = c10*zin(279) + zc00*zin(306) + c01*zin(305)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  298
                                      ! i4 = i5 =  307

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  316

                                      xin(324) = c10*xin(306) + xc00*xin(315) + c01*xin(314)
                                      yin(324) = c10*yin(306) + yc00*yin(315) + c01*yin(314)
                                      zin(324) = c10*zin(306) + zc00*zin(315) + c01*zin(314)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  307
                                      ! i4 = i5 =  316

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  316

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  316

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  307

                                      xin(316) = xin(316) + dxij*xin(307)
                                      yin(316) = yin(316) + dyij*yin(307)
                                      zin(316) = zin(316) + dzij*zin(307)

                                      ! i3 = i4 =  307
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  298

                                      xin(307) = xin(307) + dxij*xin(298)
                                      yin(307) = yin(307) + dyij*yin(298)
                                      zin(307) = zin(307) + dzij*zin(298)

                                      ! i3 = i4 =  298
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  316

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  307

                                      xin(316) = xin(316) + dxij*xin(307)
                                      yin(316) = yin(316) + dyij*yin(307)
                                      zin(316) = zin(316) + dzij*zin(307)

                                      ! i3 = i4 =  307
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  226

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  226

                                      ! do ni = 1,    3

                                      xin(226) = xin(244) + dxij*xin(217)
                                      yin(226) = yin(244) + dyij*yin(217)
                                      zin(226) = zin(244) + dzij*zin(217)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  253

                                      ! ni =    2

                                      xin(253) = xin(271) + dxij*xin(244)
                                      yin(253) = yin(271) + dyij*yin(244)
                                      zin(253) = zin(271) + dzij*zin(244)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  280

                                      ! ni =    3

                                      xin(280) = xin(298) + dxij*xin(271)
                                      yin(280) = yin(298) + dyij*yin(271)
                                      zin(280) = zin(298) + dzij*zin(271)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  307

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  235

                                      ! nj =    2

                                      ! i4 = i3 =  235

                                      ! do ni = 1,    3

                                      xin(235) = xin(253) + dxij*xin(226)
                                      yin(235) = yin(253) + dyij*yin(226)
                                      zin(235) = zin(253) + dzij*zin(226)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  262

                                      ! ni =    2

                                      xin(262) = xin(280) + dxij*xin(253)
                                      yin(262) = yin(280) + dyij*yin(253)
                                      zin(262) = zin(280) + dzij*zin(253)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  289

                                      ! ni =    3

                                      xin(289) = xin(307) + dxij*xin(280)
                                      yin(289) = yin(307) + dyij*yin(280)
                                      zin(289) = zin(307) + dzij*zin(280)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  316

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  244

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  319

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  310

                                      xin(319) = xin(319) + dxij*xin(310)
                                      yin(319) = yin(319) + dyij*yin(310)
                                      zin(319) = zin(319) + dzij*zin(310)

                                      ! i3 = i4 =  310
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  301

                                      xin(310) = xin(310) + dxij*xin(301)
                                      yin(310) = yin(310) + dyij*yin(301)
                                      zin(310) = zin(310) + dzij*zin(301)

                                      ! i3 = i4 =  301
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  319

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  310

                                      xin(319) = xin(319) + dxij*xin(310)
                                      yin(319) = yin(319) + dyij*yin(310)
                                      zin(319) = zin(319) + dzij*zin(310)

                                      ! i3 = i4 =  310
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  229

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  229

                                      ! do ni = 1,    3

                                      xin(229) = xin(247) + dxij*xin(220)
                                      yin(229) = yin(247) + dyij*yin(220)
                                      zin(229) = zin(247) + dzij*zin(220)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  256

                                      ! ni =    2

                                      xin(256) = xin(274) + dxij*xin(247)
                                      yin(256) = yin(274) + dyij*yin(247)
                                      zin(256) = zin(274) + dzij*zin(247)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                      ! ni =    3

                                      xin(283) = xin(301) + dxij*xin(274)
                                      yin(283) = yin(301) + dyij*yin(274)
                                      zin(283) = zin(301) + dzij*zin(274)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  310

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  238

                                      ! nj =    2

                                      ! i4 = i3 =  238

                                      ! do ni = 1,    3

                                      xin(238) = xin(256) + dxij*xin(229)
                                      yin(238) = yin(256) + dyij*yin(229)
                                      zin(238) = zin(256) + dzij*zin(229)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  265

                                      ! ni =    2

                                      xin(265) = xin(283) + dxij*xin(256)
                                      yin(265) = yin(283) + dyij*yin(256)
                                      zin(265) = zin(283) + dzij*zin(256)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  292

                                      ! ni =    3

                                      xin(292) = xin(310) + dxij*xin(283)
                                      yin(292) = yin(310) + dyij*yin(283)
                                      zin(292) = zin(310) + dzij*zin(283)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  319

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  247

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  322

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  313

                                      xin(322) = xin(322) + dxij*xin(313)
                                      yin(322) = yin(322) + dyij*yin(313)
                                      zin(322) = zin(322) + dzij*zin(313)

                                      ! i3 = i4 =  313
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  304

                                      xin(313) = xin(313) + dxij*xin(304)
                                      yin(313) = yin(313) + dyij*yin(304)
                                      zin(313) = zin(313) + dzij*zin(304)

                                      ! i3 = i4 =  304
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  322

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  313

                                      xin(322) = xin(322) + dxij*xin(313)
                                      yin(322) = yin(322) + dyij*yin(313)
                                      zin(322) = zin(322) + dzij*zin(313)

                                      ! i3 = i4 =  313
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  232

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  232

                                      ! do ni = 1,    3

                                      xin(232) = xin(250) + dxij*xin(223)
                                      yin(232) = yin(250) + dyij*yin(223)
                                      zin(232) = zin(250) + dzij*zin(223)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  259

                                      ! ni =    2

                                      xin(259) = xin(277) + dxij*xin(250)
                                      yin(259) = yin(277) + dyij*yin(250)
                                      zin(259) = zin(277) + dzij*zin(250)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  286

                                      ! ni =    3

                                      xin(286) = xin(304) + dxij*xin(277)
                                      yin(286) = yin(304) + dyij*yin(277)
                                      zin(286) = zin(304) + dzij*zin(277)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  313

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  241

                                      ! nj =    2

                                      ! i4 = i3 =  241

                                      ! do ni = 1,    3

                                      xin(241) = xin(259) + dxij*xin(232)
                                      yin(241) = yin(259) + dyij*yin(232)
                                      zin(241) = zin(259) + dzij*zin(232)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  268

                                      ! ni =    2

                                      xin(268) = xin(286) + dxij*xin(259)
                                      yin(268) = yin(286) + dyij*yin(259)
                                      zin(268) = zin(286) + dzij*zin(259)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  295

                                      ! ni =    3

                                      xin(295) = xin(313) + dxij*xin(286)
                                      yin(295) = yin(313) + dyij*yin(286)
                                      zin(295) = zin(313) + dzij*zin(286)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  322

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  250

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  323

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  314

                                      xin(323) = xin(323) + dxij*xin(314)
                                      yin(323) = yin(323) + dyij*yin(314)
                                      zin(323) = zin(323) + dzij*zin(314)

                                      ! i3 = i4 =  314
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  305

                                      xin(314) = xin(314) + dxij*xin(305)
                                      yin(314) = yin(314) + dyij*yin(305)
                                      zin(314) = zin(314) + dzij*zin(305)

                                      ! i3 = i4 =  305
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  323

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  314

                                      xin(323) = xin(323) + dxij*xin(314)
                                      yin(323) = yin(323) + dyij*yin(314)
                                      zin(323) = zin(323) + dzij*zin(314)

                                      ! i3 = i4 =  314
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  233

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  233

                                      ! do ni = 1,    3

                                      xin(233) = xin(251) + dxij*xin(224)
                                      yin(233) = yin(251) + dyij*yin(224)
                                      zin(233) = zin(251) + dzij*zin(224)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  260

                                      ! ni =    2

                                      xin(260) = xin(278) + dxij*xin(251)
                                      yin(260) = yin(278) + dyij*yin(251)
                                      zin(260) = zin(278) + dzij*zin(251)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                      ! ni =    3

                                      xin(287) = xin(305) + dxij*xin(278)
                                      yin(287) = yin(305) + dyij*yin(278)
                                      zin(287) = zin(305) + dzij*zin(278)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  314

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  242

                                      ! nj =    2

                                      ! i4 = i3 =  242

                                      ! do ni = 1,    3

                                      xin(242) = xin(260) + dxij*xin(233)
                                      yin(242) = yin(260) + dyij*yin(233)
                                      zin(242) = zin(260) + dzij*zin(233)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  269

                                      ! ni =    2

                                      xin(269) = xin(287) + dxij*xin(260)
                                      yin(269) = yin(287) + dyij*yin(260)
                                      zin(269) = zin(287) + dzij*zin(260)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  296

                                      ! ni =    3

                                      xin(296) = xin(314) + dxij*xin(287)
                                      yin(296) = yin(314) + dyij*yin(287)
                                      zin(296) = zin(314) + dzij*zin(287)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  323

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  251

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  324

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  315

                                      xin(324) = xin(324) + dxij*xin(315)
                                      yin(324) = yin(324) + dyij*yin(315)
                                      zin(324) = zin(324) + dzij*zin(315)

                                      ! i3 = i4 =  315
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  306

                                      xin(315) = xin(315) + dxij*xin(306)
                                      yin(315) = yin(315) + dyij*yin(306)
                                      zin(315) = zin(315) + dzij*zin(306)

                                      ! i3 = i4 =  306
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  324

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  315

                                      xin(324) = xin(324) + dxij*xin(315)
                                      yin(324) = yin(324) + dyij*yin(315)
                                      zin(324) = zin(324) + dzij*zin(315)

                                      ! i3 = i4 =  315
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  234

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  234

                                      ! do ni = 1,    3

                                      xin(234) = xin(252) + dxij*xin(225)
                                      yin(234) = yin(252) + dyij*yin(225)
                                      zin(234) = zin(252) + dzij*zin(225)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  261

                                      ! ni =    2

                                      xin(261) = xin(279) + dxij*xin(252)
                                      yin(261) = yin(279) + dyij*yin(252)
                                      zin(261) = zin(279) + dzij*zin(252)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                      ! ni =    3

                                      xin(288) = xin(306) + dxij*xin(279)
                                      yin(288) = yin(306) + dyij*yin(279)
                                      zin(288) = zin(306) + dzij*zin(279)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  315

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  243

                                      ! nj =    2

                                      ! i4 = i3 =  243

                                      ! do ni = 1,    3

                                      xin(243) = xin(261) + dxij*xin(234)
                                      yin(243) = yin(261) + dyij*yin(234)
                                      zin(243) = zin(261) + dzij*zin(234)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  270

                                      ! ni =    2

                                      xin(270) = xin(288) + dxij*xin(261)
                                      yin(270) = yin(288) + dyij*yin(261)
                                      zin(270) = zin(288) + dzij*zin(261)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  297

                                      ! ni =    3

                                      xin(297) = xin(315) + dxij*xin(288)
                                      yin(297) = yin(315) + dyij*yin(288)
                                      zin(297) = zin(315) + dzij*zin(288)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  324

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  252

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  217

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  225

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  224

                                      xin(225) = xin(225) + dxkl*xin(224)
                                      yin(225) = yin(225) + dykl*yin(224)
                                      zin(225) = zin(225) + dzkl*zin(224)

                                      ! i3 = i4 =  224
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  223

                                      xin(224) = xin(224) + dxkl*xin(223)
                                      yin(224) = yin(224) + dykl*yin(223)
                                      zin(224) = zin(224) + dzkl*zin(223)

                                      ! i3 = i4 =  223
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  225

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  224

                                      xin(225) = xin(225) + dxkl*xin(224)
                                      yin(225) = yin(225) + dykl*yin(224)
                                      zin(225) = zin(225) + dzkl*zin(224)

                                      ! i3 = i4 =  224
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  218

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  218

                                      ! do nk = 1,    2

                                      xin(218) = xin(220) + dxkl*xin(217)
                                      yin(218) = yin(220) + dykl*yin(217)
                                      zin(218) = zin(220) + dzkl*zin(217)
                                      ! i4 = i4 + lang+1 =  221

                                      ! nk =    2

                                      xin(221) = xin(223) + dxkl*xin(220)
                                      yin(221) = yin(223) + dykl*yin(220)
                                      zin(221) = zin(223) + dzkl*zin(220)
                                      ! i4 = i4 + lang+1 =  224

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  219

                                      ! nl =    2

                                      ! i4 = i3 =  219

                                      ! do nk = 1,    2

                                      xin(219) = xin(221) + dxkl*xin(218)
                                      yin(219) = yin(221) + dykl*yin(218)
                                      zin(219) = zin(221) + dzkl*zin(218)
                                      ! i4 = i4 + lang+1 =  222

                                      ! nk =    2

                                      xin(222) = xin(224) + dxkl*xin(221)
                                      yin(222) = yin(224) + dykl*yin(221)
                                      zin(222) = zin(224) + dzkl*zin(221)
                                      ! i4 = i4 + lang+1 =  225

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  220

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  226

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  234

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  233

                                      xin(234) = xin(234) + dxkl*xin(233)
                                      yin(234) = yin(234) + dykl*yin(233)
                                      zin(234) = zin(234) + dzkl*zin(233)

                                      ! i3 = i4 =  233
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  232

                                      xin(233) = xin(233) + dxkl*xin(232)
                                      yin(233) = yin(233) + dykl*yin(232)
                                      zin(233) = zin(233) + dzkl*zin(232)

                                      ! i3 = i4 =  232
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  234

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  233

                                      xin(234) = xin(234) + dxkl*xin(233)
                                      yin(234) = yin(234) + dykl*yin(233)
                                      zin(234) = zin(234) + dzkl*zin(233)

                                      ! i3 = i4 =  233
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  227

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  227

                                      ! do nk = 1,    2

                                      xin(227) = xin(229) + dxkl*xin(226)
                                      yin(227) = yin(229) + dykl*yin(226)
                                      zin(227) = zin(229) + dzkl*zin(226)
                                      ! i4 = i4 + lang+1 =  230

                                      ! nk =    2

                                      xin(230) = xin(232) + dxkl*xin(229)
                                      yin(230) = yin(232) + dykl*yin(229)
                                      zin(230) = zin(232) + dzkl*zin(229)
                                      ! i4 = i4 + lang+1 =  233

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  228

                                      ! nl =    2

                                      ! i4 = i3 =  228

                                      ! do nk = 1,    2

                                      xin(228) = xin(230) + dxkl*xin(227)
                                      yin(228) = yin(230) + dykl*yin(227)
                                      zin(228) = zin(230) + dzkl*zin(227)
                                      ! i4 = i4 + lang+1 =  231

                                      ! nk =    2

                                      xin(231) = xin(233) + dxkl*xin(230)
                                      yin(231) = yin(233) + dykl*yin(230)
                                      zin(231) = zin(233) + dzkl*zin(230)
                                      ! i4 = i4 + lang+1 =  234

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  229

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  235

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  243

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  242

                                      xin(243) = xin(243) + dxkl*xin(242)
                                      yin(243) = yin(243) + dykl*yin(242)
                                      zin(243) = zin(243) + dzkl*zin(242)

                                      ! i3 = i4 =  242
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  241

                                      xin(242) = xin(242) + dxkl*xin(241)
                                      yin(242) = yin(242) + dykl*yin(241)
                                      zin(242) = zin(242) + dzkl*zin(241)

                                      ! i3 = i4 =  241
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  243

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  242

                                      xin(243) = xin(243) + dxkl*xin(242)
                                      yin(243) = yin(243) + dykl*yin(242)
                                      zin(243) = zin(243) + dzkl*zin(242)

                                      ! i3 = i4 =  242
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  236

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  236

                                      ! do nk = 1,    2

                                      xin(236) = xin(238) + dxkl*xin(235)
                                      yin(236) = yin(238) + dykl*yin(235)
                                      zin(236) = zin(238) + dzkl*zin(235)
                                      ! i4 = i4 + lang+1 =  239

                                      ! nk =    2

                                      xin(239) = xin(241) + dxkl*xin(238)
                                      yin(239) = yin(241) + dykl*yin(238)
                                      zin(239) = zin(241) + dzkl*zin(238)
                                      ! i4 = i4 + lang+1 =  242

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  237

                                      ! nl =    2

                                      ! i4 = i3 =  237

                                      ! do nk = 1,    2

                                      xin(237) = xin(239) + dxkl*xin(236)
                                      yin(237) = yin(239) + dykl*yin(236)
                                      zin(237) = zin(239) + dzkl*zin(236)
                                      ! i4 = i4 + lang+1 =  240

                                      ! nk =    2

                                      xin(240) = xin(242) + dxkl*xin(239)
                                      yin(240) = yin(242) + dykl*yin(239)
                                      zin(240) = zin(242) + dzkl*zin(239)
                                      ! i4 = i4 + lang+1 =  243

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  238

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  244

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  244

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  250

                                      xin(251) = xin(251) + dxkl*xin(250)
                                      yin(251) = yin(251) + dykl*yin(250)
                                      zin(251) = zin(251) + dzkl*zin(250)

                                      ! i3 = i4 =  250
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  252

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  251

                                      xin(252) = xin(252) + dxkl*xin(251)
                                      yin(252) = yin(252) + dykl*yin(251)
                                      zin(252) = zin(252) + dzkl*zin(251)

                                      ! i3 = i4 =  251
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  245

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  245

                                      ! do nk = 1,    2

                                      xin(245) = xin(247) + dxkl*xin(244)
                                      yin(245) = yin(247) + dykl*yin(244)
                                      zin(245) = zin(247) + dzkl*zin(244)
                                      ! i4 = i4 + lang+1 =  248

                                      ! nk =    2

                                      xin(248) = xin(250) + dxkl*xin(247)
                                      yin(248) = yin(250) + dykl*yin(247)
                                      zin(248) = zin(250) + dzkl*zin(247)
                                      ! i4 = i4 + lang+1 =  251

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  246

                                      ! nl =    2

                                      ! i4 = i3 =  246

                                      ! do nk = 1,    2

                                      xin(246) = xin(248) + dxkl*xin(245)
                                      yin(246) = yin(248) + dykl*yin(245)
                                      zin(246) = zin(248) + dzkl*zin(245)
                                      ! i4 = i4 + lang+1 =  249

                                      ! nk =    2

                                      xin(249) = xin(251) + dxkl*xin(248)
                                      yin(249) = yin(251) + dykl*yin(248)
                                      zin(249) = zin(251) + dzkl*zin(248)
                                      ! i4 = i4 + lang+1 =  252

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  247

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  253

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  261

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  260

                                      xin(261) = xin(261) + dxkl*xin(260)
                                      yin(261) = yin(261) + dykl*yin(260)
                                      zin(261) = zin(261) + dzkl*zin(260)

                                      ! i3 = i4 =  260
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  259

                                      xin(260) = xin(260) + dxkl*xin(259)
                                      yin(260) = yin(260) + dykl*yin(259)
                                      zin(260) = zin(260) + dzkl*zin(259)

                                      ! i3 = i4 =  259
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  261

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  260

                                      xin(261) = xin(261) + dxkl*xin(260)
                                      yin(261) = yin(261) + dykl*yin(260)
                                      zin(261) = zin(261) + dzkl*zin(260)

                                      ! i3 = i4 =  260
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  254

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  254

                                      ! do nk = 1,    2

                                      xin(254) = xin(256) + dxkl*xin(253)
                                      yin(254) = yin(256) + dykl*yin(253)
                                      zin(254) = zin(256) + dzkl*zin(253)
                                      ! i4 = i4 + lang+1 =  257

                                      ! nk =    2

                                      xin(257) = xin(259) + dxkl*xin(256)
                                      yin(257) = yin(259) + dykl*yin(256)
                                      zin(257) = zin(259) + dzkl*zin(256)
                                      ! i4 = i4 + lang+1 =  260

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  255

                                      ! nl =    2

                                      ! i4 = i3 =  255

                                      ! do nk = 1,    2

                                      xin(255) = xin(257) + dxkl*xin(254)
                                      yin(255) = yin(257) + dykl*yin(254)
                                      zin(255) = zin(257) + dzkl*zin(254)
                                      ! i4 = i4 + lang+1 =  258

                                      ! nk =    2

                                      xin(258) = xin(260) + dxkl*xin(257)
                                      yin(258) = yin(260) + dykl*yin(257)
                                      zin(258) = zin(260) + dzkl*zin(257)
                                      ! i4 = i4 + lang+1 =  261

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  256

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  262

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  270

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  269

                                      xin(270) = xin(270) + dxkl*xin(269)
                                      yin(270) = yin(270) + dykl*yin(269)
                                      zin(270) = zin(270) + dzkl*zin(269)

                                      ! i3 = i4 =  269
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  268

                                      xin(269) = xin(269) + dxkl*xin(268)
                                      yin(269) = yin(269) + dykl*yin(268)
                                      zin(269) = zin(269) + dzkl*zin(268)

                                      ! i3 = i4 =  268
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  270

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  269

                                      xin(270) = xin(270) + dxkl*xin(269)
                                      yin(270) = yin(270) + dykl*yin(269)
                                      zin(270) = zin(270) + dzkl*zin(269)

                                      ! i3 = i4 =  269
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  263

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  263

                                      ! do nk = 1,    2

                                      xin(263) = xin(265) + dxkl*xin(262)
                                      yin(263) = yin(265) + dykl*yin(262)
                                      zin(263) = zin(265) + dzkl*zin(262)
                                      ! i4 = i4 + lang+1 =  266

                                      ! nk =    2

                                      xin(266) = xin(268) + dxkl*xin(265)
                                      yin(266) = yin(268) + dykl*yin(265)
                                      zin(266) = zin(268) + dzkl*zin(265)
                                      ! i4 = i4 + lang+1 =  269

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  264

                                      ! nl =    2

                                      ! i4 = i3 =  264

                                      ! do nk = 1,    2

                                      xin(264) = xin(266) + dxkl*xin(263)
                                      yin(264) = yin(266) + dykl*yin(263)
                                      zin(264) = zin(266) + dzkl*zin(263)
                                      ! i4 = i4 + lang+1 =  267

                                      ! nk =    2

                                      xin(267) = xin(269) + dxkl*xin(266)
                                      yin(267) = yin(269) + dykl*yin(266)
                                      zin(267) = zin(269) + dzkl*zin(266)
                                      ! i4 = i4 + lang+1 =  270

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  265

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  271

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  271

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  279

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  278

                                      xin(279) = xin(279) + dxkl*xin(278)
                                      yin(279) = yin(279) + dykl*yin(278)
                                      zin(279) = zin(279) + dzkl*zin(278)

                                      ! i3 = i4 =  278
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  277

                                      xin(278) = xin(278) + dxkl*xin(277)
                                      yin(278) = yin(278) + dykl*yin(277)
                                      zin(278) = zin(278) + dzkl*zin(277)

                                      ! i3 = i4 =  277
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  279

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  278

                                      xin(279) = xin(279) + dxkl*xin(278)
                                      yin(279) = yin(279) + dykl*yin(278)
                                      zin(279) = zin(279) + dzkl*zin(278)

                                      ! i3 = i4 =  278
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  272

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  272

                                      ! do nk = 1,    2

                                      xin(272) = xin(274) + dxkl*xin(271)
                                      yin(272) = yin(274) + dykl*yin(271)
                                      zin(272) = zin(274) + dzkl*zin(271)
                                      ! i4 = i4 + lang+1 =  275

                                      ! nk =    2

                                      xin(275) = xin(277) + dxkl*xin(274)
                                      yin(275) = yin(277) + dykl*yin(274)
                                      zin(275) = zin(277) + dzkl*zin(274)
                                      ! i4 = i4 + lang+1 =  278

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  273

                                      ! nl =    2

                                      ! i4 = i3 =  273

                                      ! do nk = 1,    2

                                      xin(273) = xin(275) + dxkl*xin(272)
                                      yin(273) = yin(275) + dykl*yin(272)
                                      zin(273) = zin(275) + dzkl*zin(272)
                                      ! i4 = i4 + lang+1 =  276

                                      ! nk =    2

                                      xin(276) = xin(278) + dxkl*xin(275)
                                      yin(276) = yin(278) + dykl*yin(275)
                                      zin(276) = zin(278) + dzkl*zin(275)
                                      ! i4 = i4 + lang+1 =  279

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  274

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  280

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  286

                                      xin(287) = xin(287) + dxkl*xin(286)
                                      yin(287) = yin(287) + dykl*yin(286)
                                      zin(287) = zin(287) + dzkl*zin(286)

                                      ! i3 = i4 =  286
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  288

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  287

                                      xin(288) = xin(288) + dxkl*xin(287)
                                      yin(288) = yin(288) + dykl*yin(287)
                                      zin(288) = zin(288) + dzkl*zin(287)

                                      ! i3 = i4 =  287
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  281

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  281

                                      ! do nk = 1,    2

                                      xin(281) = xin(283) + dxkl*xin(280)
                                      yin(281) = yin(283) + dykl*yin(280)
                                      zin(281) = zin(283) + dzkl*zin(280)
                                      ! i4 = i4 + lang+1 =  284

                                      ! nk =    2

                                      xin(284) = xin(286) + dxkl*xin(283)
                                      yin(284) = yin(286) + dykl*yin(283)
                                      zin(284) = zin(286) + dzkl*zin(283)
                                      ! i4 = i4 + lang+1 =  287

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  282

                                      ! nl =    2

                                      ! i4 = i3 =  282

                                      ! do nk = 1,    2

                                      xin(282) = xin(284) + dxkl*xin(281)
                                      yin(282) = yin(284) + dykl*yin(281)
                                      zin(282) = zin(284) + dzkl*zin(281)
                                      ! i4 = i4 + lang+1 =  285

                                      ! nk =    2

                                      xin(285) = xin(287) + dxkl*xin(284)
                                      yin(285) = yin(287) + dykl*yin(284)
                                      zin(285) = zin(287) + dzkl*zin(284)
                                      ! i4 = i4 + lang+1 =  288

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  283

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  289

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  297

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  296

                                      xin(297) = xin(297) + dxkl*xin(296)
                                      yin(297) = yin(297) + dykl*yin(296)
                                      zin(297) = zin(297) + dzkl*zin(296)

                                      ! i3 = i4 =  296
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  295

                                      xin(296) = xin(296) + dxkl*xin(295)
                                      yin(296) = yin(296) + dykl*yin(295)
                                      zin(296) = zin(296) + dzkl*zin(295)

                                      ! i3 = i4 =  295
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  297

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  296

                                      xin(297) = xin(297) + dxkl*xin(296)
                                      yin(297) = yin(297) + dykl*yin(296)
                                      zin(297) = zin(297) + dzkl*zin(296)

                                      ! i3 = i4 =  296
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  290

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  290

                                      ! do nk = 1,    2

                                      xin(290) = xin(292) + dxkl*xin(289)
                                      yin(290) = yin(292) + dykl*yin(289)
                                      zin(290) = zin(292) + dzkl*zin(289)
                                      ! i4 = i4 + lang+1 =  293

                                      ! nk =    2

                                      xin(293) = xin(295) + dxkl*xin(292)
                                      yin(293) = yin(295) + dykl*yin(292)
                                      zin(293) = zin(295) + dzkl*zin(292)
                                      ! i4 = i4 + lang+1 =  296

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  291

                                      ! nl =    2

                                      ! i4 = i3 =  291

                                      ! do nk = 1,    2

                                      xin(291) = xin(293) + dxkl*xin(290)
                                      yin(291) = yin(293) + dykl*yin(290)
                                      zin(291) = zin(293) + dzkl*zin(290)
                                      ! i4 = i4 + lang+1 =  294

                                      ! nk =    2

                                      xin(294) = xin(296) + dxkl*xin(293)
                                      yin(294) = yin(296) + dykl*yin(293)
                                      zin(294) = zin(296) + dzkl*zin(293)
                                      ! i4 = i4 + lang+1 =  297

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  292

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  298

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  298

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  306

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  305

                                      xin(306) = xin(306) + dxkl*xin(305)
                                      yin(306) = yin(306) + dykl*yin(305)
                                      zin(306) = zin(306) + dzkl*zin(305)

                                      ! i3 = i4 =  305
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  304

                                      xin(305) = xin(305) + dxkl*xin(304)
                                      yin(305) = yin(305) + dykl*yin(304)
                                      zin(305) = zin(305) + dzkl*zin(304)

                                      ! i3 = i4 =  304
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  306

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  305

                                      xin(306) = xin(306) + dxkl*xin(305)
                                      yin(306) = yin(306) + dykl*yin(305)
                                      zin(306) = zin(306) + dzkl*zin(305)

                                      ! i3 = i4 =  305
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  299

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  299

                                      ! do nk = 1,    2

                                      xin(299) = xin(301) + dxkl*xin(298)
                                      yin(299) = yin(301) + dykl*yin(298)
                                      zin(299) = zin(301) + dzkl*zin(298)
                                      ! i4 = i4 + lang+1 =  302

                                      ! nk =    2

                                      xin(302) = xin(304) + dxkl*xin(301)
                                      yin(302) = yin(304) + dykl*yin(301)
                                      zin(302) = zin(304) + dzkl*zin(301)
                                      ! i4 = i4 + lang+1 =  305

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  300

                                      ! nl =    2

                                      ! i4 = i3 =  300

                                      ! do nk = 1,    2

                                      xin(300) = xin(302) + dxkl*xin(299)
                                      yin(300) = yin(302) + dykl*yin(299)
                                      zin(300) = zin(302) + dzkl*zin(299)
                                      ! i4 = i4 + lang+1 =  303

                                      ! nk =    2

                                      xin(303) = xin(305) + dxkl*xin(302)
                                      yin(303) = yin(305) + dykl*yin(302)
                                      zin(303) = zin(305) + dzkl*zin(302)
                                      ! i4 = i4 + lang+1 =  306

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  301

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  307

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  315

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  314

                                      xin(315) = xin(315) + dxkl*xin(314)
                                      yin(315) = yin(315) + dykl*yin(314)
                                      zin(315) = zin(315) + dzkl*zin(314)

                                      ! i3 = i4 =  314
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  313

                                      xin(314) = xin(314) + dxkl*xin(313)
                                      yin(314) = yin(314) + dykl*yin(313)
                                      zin(314) = zin(314) + dzkl*zin(313)

                                      ! i3 = i4 =  313
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  315

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  314

                                      xin(315) = xin(315) + dxkl*xin(314)
                                      yin(315) = yin(315) + dykl*yin(314)
                                      zin(315) = zin(315) + dzkl*zin(314)

                                      ! i3 = i4 =  314
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  308

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  308

                                      ! do nk = 1,    2

                                      xin(308) = xin(310) + dxkl*xin(307)
                                      yin(308) = yin(310) + dykl*yin(307)
                                      zin(308) = zin(310) + dzkl*zin(307)
                                      ! i4 = i4 + lang+1 =  311

                                      ! nk =    2

                                      xin(311) = xin(313) + dxkl*xin(310)
                                      yin(311) = yin(313) + dykl*yin(310)
                                      zin(311) = zin(313) + dzkl*zin(310)
                                      ! i4 = i4 + lang+1 =  314

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  309

                                      ! nl =    2

                                      ! i4 = i3 =  309

                                      ! do nk = 1,    2

                                      xin(309) = xin(311) + dxkl*xin(308)
                                      yin(309) = yin(311) + dykl*yin(308)
                                      zin(309) = zin(311) + dzkl*zin(308)
                                      ! i4 = i4 + lang+1 =  312

                                      ! nk =    2

                                      xin(312) = xin(314) + dxkl*xin(311)
                                      yin(312) = yin(314) + dykl*yin(311)
                                      zin(312) = zin(314) + dzkl*zin(311)
                                      ! i4 = i4 + lang+1 =  315

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  310

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  316

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  324

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  323

                                      xin(324) = xin(324) + dxkl*xin(323)
                                      yin(324) = yin(324) + dykl*yin(323)
                                      zin(324) = zin(324) + dzkl*zin(323)

                                      ! i3 = i4 =  323
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  322

                                      xin(323) = xin(323) + dxkl*xin(322)
                                      yin(323) = yin(323) + dykl*yin(322)
                                      zin(323) = zin(323) + dzkl*zin(322)

                                      ! i3 = i4 =  322
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  324

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  323

                                      xin(324) = xin(324) + dxkl*xin(323)
                                      yin(324) = yin(324) + dykl*yin(323)
                                      zin(324) = zin(324) + dzkl*zin(323)

                                      ! i3 = i4 =  323
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  317

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  317

                                      ! do nk = 1,    2

                                      xin(317) = xin(319) + dxkl*xin(316)
                                      yin(317) = yin(319) + dykl*yin(316)
                                      zin(317) = zin(319) + dzkl*zin(316)
                                      ! i4 = i4 + lang+1 =  320

                                      ! nk =    2

                                      xin(320) = xin(322) + dxkl*xin(319)
                                      yin(320) = yin(322) + dykl*yin(319)
                                      zin(320) = zin(322) + dzkl*zin(319)
                                      ! i4 = i4 + lang+1 =  323

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  318

                                      ! nl =    2

                                      ! i4 = i3 =  318

                                      ! do nk = 1,    2

                                      xin(318) = xin(320) + dxkl*xin(317)
                                      yin(318) = yin(320) + dykl*yin(317)
                                      zin(318) = zin(320) + dzkl*zin(317)
                                      ! i4 = i4 + lang+1 =  321

                                      ! nk =    2

                                      xin(321) = xin(323) + dxkl*xin(320)
                                      yin(321) = yin(323) + dykl*yin(320)
                                      zin(321) = zin(323) + dzkl*zin(320)
                                      ! i4 = i4 + lang+1 =  324

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  319

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  325

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  325

                                      ! end do

                                      ! *** Now root =    4

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  324

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

                                      ! i1 = in(1) =  325

                                      xin(325) = 1.0_dp
                                      yin(325) = 1.0_dp
                                      zin(325) = f00

                                      ! i2 = in(2) =  352
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(352) = xc00
                                      yin(352) = yc00
                                      zin(352) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  328

                                      xin(328) = xcp00
                                      yin(328) = ycp00
                                      zin(328) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  355
                                      ! i2 =  352

                                      xin(355) = xcp00*xin(352) + cp10
                                      yin(355) = ycp00*yin(352) + cp10
                                      zin(355) = zcp00*zin(352) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  325
                                      ! i4 = i2 =  352

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  379
                                      ! i3 =  325
                                      ! i4 =  352

                                      xin(379) = c10*xin(325) + xc00*xin(352)
                                      yin(379) = c10*yin(325) + yc00*yin(352)
                                      zin(379) = c10*zin(325) + zc00*zin(352)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  382
                                      ! i5 =  379
                                      ! i4 =  352

                                      xin(382) = xcp00*xin(379) + cp10*xin(352)
                                      yin(382) = ycp00*yin(379) + cp10*yin(352)
                                      zin(382) = zcp00*zin(379) + cp10*zin(352)

                                      ! ------------------

                                      ! i3 = i4 =  352
                                      ! i4 = i5 =  379

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  406
                                      ! i3 =  352
                                      ! i4 =  379

                                      xin(406) = c10*xin(352) + xc00*xin(379)
                                      yin(406) = c10*yin(352) + yc00*yin(379)
                                      zin(406) = c10*zin(352) + zc00*zin(379)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  409
                                      ! i5 =  406
                                      ! i4 =  379

                                      xin(409) = xcp00*xin(406) + cp10*xin(379)
                                      yin(409) = ycp00*yin(406) + cp10*yin(379)
                                      zin(409) = zcp00*zin(406) + cp10*zin(379)

                                      ! ------------------

                                      ! i3 = i4 =  379
                                      ! i4 = i5 =  406

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  415
                                      ! i3 =  379
                                      ! i4 =  406

                                      xin(415) = c10*xin(379) + xc00*xin(406)
                                      yin(415) = c10*yin(379) + yc00*yin(406)
                                      zin(415) = c10*zin(379) + zc00*zin(406)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  418
                                      ! i5 =  415
                                      ! i4 =  406

                                      xin(418) = xcp00*xin(415) + cp10*xin(406)
                                      yin(418) = ycp00*yin(415) + cp10*yin(406)
                                      zin(418) = zcp00*zin(415) + cp10*zin(406)

                                      ! ------------------

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  424
                                      ! i3 =  406
                                      ! i4 =  415

                                      xin(424) = c10*xin(406) + xc00*xin(415)
                                      yin(424) = c10*yin(406) + yc00*yin(415)
                                      zin(424) = c10*zin(406) + zc00*zin(415)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  427
                                      ! i5 =  424
                                      ! i4 =  415

                                      xin(427) = xcp00*xin(424) + cp10*xin(415)
                                      yin(427) = ycp00*yin(424) + cp10*yin(415)
                                      zin(427) = zcp00*zin(424) + cp10*zin(415)

                                      ! ------------------

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  325
                                      ! i4 = i1+k2 =  328

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  331
                                      ! i3 =  325
                                      ! i4 =  328

                                      xin(331) = cp01*xin(325) + xcp00*xin(328)
                                      yin(331) = cp01*yin(325) + ycp00*yin(328)
                                      zin(331) = cp01*zin(325) + zcp00*zin(328)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  358

                                      xin(358) = xc00*xin(331) + c01*xin(328)
                                      yin(358) = yc00*yin(331) + c01*yin(328)
                                      zin(358) = zc00*zin(331) + c01*zin(328)

                                      ! ------------------

                                      ! i3 = i4 =  328
                                      ! i4 = i5 =  331

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  332
                                      ! i3 =  328
                                      ! i4 =  331

                                      xin(332) = cp01*xin(328) + xcp00*xin(331)
                                      yin(332) = cp01*yin(328) + ycp00*yin(331)
                                      zin(332) = cp01*zin(328) + zcp00*zin(331)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  359

                                      xin(359) = xc00*xin(332) + c01*xin(331)
                                      yin(359) = yc00*yin(332) + c01*yin(331)
                                      zin(359) = zc00*zin(332) + c01*zin(331)

                                      ! ------------------

                                      ! i3 = i4 =  331
                                      ! i4 = i5 =  332

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  333
                                      ! i3 =  331
                                      ! i4 =  332

                                      xin(333) = cp01*xin(331) + xcp00*xin(332)
                                      yin(333) = cp01*yin(331) + ycp00*yin(332)
                                      zin(333) = cp01*zin(331) + zcp00*zin(332)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  360

                                      xin(360) = xc00*xin(333) + c01*xin(332)
                                      yin(360) = yc00*yin(333) + c01*yin(332)
                                      zin(360) = zc00*zin(333) + c01*zin(332)

                                      ! ------------------

                                      ! i3 = i4 =  332
                                      ! i4 = i5 =  333

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  325
                                      ! i4 = i2 =  352

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  379

                                      xin(385) = c10*xin(331) + xc00*xin(358) + c01*xin(355)
                                      yin(385) = c10*yin(331) + yc00*yin(358) + c01*yin(355)
                                      zin(385) = c10*zin(331) + zc00*zin(358) + c01*zin(355)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  352
                                      ! i4 = i5 =  379

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  406

                                      xin(412) = c10*xin(358) + xc00*xin(385) + c01*xin(382)
                                      yin(412) = c10*yin(358) + yc00*yin(385) + c01*yin(382)
                                      zin(412) = c10*zin(358) + zc00*zin(385) + c01*zin(382)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  379
                                      ! i4 = i5 =  406

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  415

                                      xin(421) = c10*xin(385) + xc00*xin(412) + c01*xin(409)
                                      yin(421) = c10*yin(385) + yc00*yin(412) + c01*yin(409)
                                      zin(421) = c10*zin(385) + zc00*zin(412) + c01*zin(409)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  424

                                      xin(430) = c10*xin(412) + xc00*xin(421) + c01*xin(418)
                                      yin(430) = c10*yin(412) + yc00*yin(421) + c01*yin(418)
                                      zin(430) = c10*zin(412) + zc00*zin(421) + c01*zin(418)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  325
                                      ! i4 = i2 =  352

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  379

                                      xin(386) = c10*xin(332) + xc00*xin(359) + c01*xin(358)
                                      yin(386) = c10*yin(332) + yc00*yin(359) + c01*yin(358)
                                      zin(386) = c10*zin(332) + zc00*zin(359) + c01*zin(358)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  352
                                      ! i4 = i5 =  379

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  406

                                      xin(413) = c10*xin(359) + xc00*xin(386) + c01*xin(385)
                                      yin(413) = c10*yin(359) + yc00*yin(386) + c01*yin(385)
                                      zin(413) = c10*zin(359) + zc00*zin(386) + c01*zin(385)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  379
                                      ! i4 = i5 =  406

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  415

                                      xin(422) = c10*xin(386) + xc00*xin(413) + c01*xin(412)
                                      yin(422) = c10*yin(386) + yc00*yin(413) + c01*yin(412)
                                      zin(422) = c10*zin(386) + zc00*zin(413) + c01*zin(412)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  424

                                      xin(431) = c10*xin(413) + xc00*xin(422) + c01*xin(421)
                                      yin(431) = c10*yin(413) + yc00*yin(422) + c01*yin(421)
                                      zin(431) = c10*zin(413) + zc00*zin(422) + c01*zin(421)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  325
                                      ! i4 = i2 =  352

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  379

                                      xin(387) = c10*xin(333) + xc00*xin(360) + c01*xin(359)
                                      yin(387) = c10*yin(333) + yc00*yin(360) + c01*yin(359)
                                      zin(387) = c10*zin(333) + zc00*zin(360) + c01*zin(359)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  352
                                      ! i4 = i5 =  379

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  406

                                      xin(414) = c10*xin(360) + xc00*xin(387) + c01*xin(386)
                                      yin(414) = c10*yin(360) + yc00*yin(387) + c01*yin(386)
                                      zin(414) = c10*zin(360) + zc00*zin(387) + c01*zin(386)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  379
                                      ! i4 = i5 =  406

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  415

                                      xin(423) = c10*xin(387) + xc00*xin(414) + c01*xin(413)
                                      yin(423) = c10*yin(387) + yc00*yin(414) + c01*yin(413)
                                      zin(423) = c10*zin(387) + zc00*zin(414) + c01*zin(413)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  406
                                      ! i4 = i5 =  415

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  424

                                      xin(432) = c10*xin(414) + xc00*xin(423) + c01*xin(422)
                                      yin(432) = c10*yin(414) + yc00*yin(423) + c01*yin(422)
                                      zin(432) = c10*zin(414) + zc00*zin(423) + c01*zin(422)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  415
                                      ! i4 = i5 =  424

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  424

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  424

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  415

                                      xin(424) = xin(424) + dxij*xin(415)
                                      yin(424) = yin(424) + dyij*yin(415)
                                      zin(424) = zin(424) + dzij*zin(415)

                                      ! i3 = i4 =  415
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  406

                                      xin(415) = xin(415) + dxij*xin(406)
                                      yin(415) = yin(415) + dyij*yin(406)
                                      zin(415) = zin(415) + dzij*zin(406)

                                      ! i3 = i4 =  406
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  424

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  415

                                      xin(424) = xin(424) + dxij*xin(415)
                                      yin(424) = yin(424) + dyij*yin(415)
                                      zin(424) = zin(424) + dzij*zin(415)

                                      ! i3 = i4 =  415
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  334

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  334

                                      ! do ni = 1,    3

                                      xin(334) = xin(352) + dxij*xin(325)
                                      yin(334) = yin(352) + dyij*yin(325)
                                      zin(334) = zin(352) + dzij*zin(325)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  361

                                      ! ni =    2

                                      xin(361) = xin(379) + dxij*xin(352)
                                      yin(361) = yin(379) + dyij*yin(352)
                                      zin(361) = zin(379) + dzij*zin(352)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  388

                                      ! ni =    3

                                      xin(388) = xin(406) + dxij*xin(379)
                                      yin(388) = yin(406) + dyij*yin(379)
                                      zin(388) = zin(406) + dzij*zin(379)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  415

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  343

                                      ! nj =    2

                                      ! i4 = i3 =  343

                                      ! do ni = 1,    3

                                      xin(343) = xin(361) + dxij*xin(334)
                                      yin(343) = yin(361) + dyij*yin(334)
                                      zin(343) = zin(361) + dzij*zin(334)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  370

                                      ! ni =    2

                                      xin(370) = xin(388) + dxij*xin(361)
                                      yin(370) = yin(388) + dyij*yin(361)
                                      zin(370) = zin(388) + dzij*zin(361)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  397

                                      ! ni =    3

                                      xin(397) = xin(415) + dxij*xin(388)
                                      yin(397) = yin(415) + dyij*yin(388)
                                      zin(397) = zin(415) + dzij*zin(388)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  424

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  352

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  427

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  418

                                      xin(427) = xin(427) + dxij*xin(418)
                                      yin(427) = yin(427) + dyij*yin(418)
                                      zin(427) = zin(427) + dzij*zin(418)

                                      ! i3 = i4 =  418
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  409

                                      xin(418) = xin(418) + dxij*xin(409)
                                      yin(418) = yin(418) + dyij*yin(409)
                                      zin(418) = zin(418) + dzij*zin(409)

                                      ! i3 = i4 =  409
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  427

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  418

                                      xin(427) = xin(427) + dxij*xin(418)
                                      yin(427) = yin(427) + dyij*yin(418)
                                      zin(427) = zin(427) + dzij*zin(418)

                                      ! i3 = i4 =  418
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  337

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  337

                                      ! do ni = 1,    3

                                      xin(337) = xin(355) + dxij*xin(328)
                                      yin(337) = yin(355) + dyij*yin(328)
                                      zin(337) = zin(355) + dzij*zin(328)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  364

                                      ! ni =    2

                                      xin(364) = xin(382) + dxij*xin(355)
                                      yin(364) = yin(382) + dyij*yin(355)
                                      zin(364) = zin(382) + dzij*zin(355)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  391

                                      ! ni =    3

                                      xin(391) = xin(409) + dxij*xin(382)
                                      yin(391) = yin(409) + dyij*yin(382)
                                      zin(391) = zin(409) + dzij*zin(382)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  418

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  346

                                      ! nj =    2

                                      ! i4 = i3 =  346

                                      ! do ni = 1,    3

                                      xin(346) = xin(364) + dxij*xin(337)
                                      yin(346) = yin(364) + dyij*yin(337)
                                      zin(346) = zin(364) + dzij*zin(337)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  373

                                      ! ni =    2

                                      xin(373) = xin(391) + dxij*xin(364)
                                      yin(373) = yin(391) + dyij*yin(364)
                                      zin(373) = zin(391) + dzij*zin(364)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  400

                                      ! ni =    3

                                      xin(400) = xin(418) + dxij*xin(391)
                                      yin(400) = yin(418) + dyij*yin(391)
                                      zin(400) = zin(418) + dzij*zin(391)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  427

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  355

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  430

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  421

                                      xin(430) = xin(430) + dxij*xin(421)
                                      yin(430) = yin(430) + dyij*yin(421)
                                      zin(430) = zin(430) + dzij*zin(421)

                                      ! i3 = i4 =  421
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  412

                                      xin(421) = xin(421) + dxij*xin(412)
                                      yin(421) = yin(421) + dyij*yin(412)
                                      zin(421) = zin(421) + dzij*zin(412)

                                      ! i3 = i4 =  412
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  430

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  421

                                      xin(430) = xin(430) + dxij*xin(421)
                                      yin(430) = yin(430) + dyij*yin(421)
                                      zin(430) = zin(430) + dzij*zin(421)

                                      ! i3 = i4 =  421
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  340

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  340

                                      ! do ni = 1,    3

                                      xin(340) = xin(358) + dxij*xin(331)
                                      yin(340) = yin(358) + dyij*yin(331)
                                      zin(340) = zin(358) + dzij*zin(331)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  367

                                      ! ni =    2

                                      xin(367) = xin(385) + dxij*xin(358)
                                      yin(367) = yin(385) + dyij*yin(358)
                                      zin(367) = zin(385) + dzij*zin(358)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  394

                                      ! ni =    3

                                      xin(394) = xin(412) + dxij*xin(385)
                                      yin(394) = yin(412) + dyij*yin(385)
                                      zin(394) = zin(412) + dzij*zin(385)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  421

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  349

                                      ! nj =    2

                                      ! i4 = i3 =  349

                                      ! do ni = 1,    3

                                      xin(349) = xin(367) + dxij*xin(340)
                                      yin(349) = yin(367) + dyij*yin(340)
                                      zin(349) = zin(367) + dzij*zin(340)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  376

                                      ! ni =    2

                                      xin(376) = xin(394) + dxij*xin(367)
                                      yin(376) = yin(394) + dyij*yin(367)
                                      zin(376) = zin(394) + dzij*zin(367)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  403

                                      ! ni =    3

                                      xin(403) = xin(421) + dxij*xin(394)
                                      yin(403) = yin(421) + dyij*yin(394)
                                      zin(403) = zin(421) + dzij*zin(394)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  430

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  358

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  431

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  422

                                      xin(431) = xin(431) + dxij*xin(422)
                                      yin(431) = yin(431) + dyij*yin(422)
                                      zin(431) = zin(431) + dzij*zin(422)

                                      ! i3 = i4 =  422
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  413

                                      xin(422) = xin(422) + dxij*xin(413)
                                      yin(422) = yin(422) + dyij*yin(413)
                                      zin(422) = zin(422) + dzij*zin(413)

                                      ! i3 = i4 =  413
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  431

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  422

                                      xin(431) = xin(431) + dxij*xin(422)
                                      yin(431) = yin(431) + dyij*yin(422)
                                      zin(431) = zin(431) + dzij*zin(422)

                                      ! i3 = i4 =  422
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  341

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  341

                                      ! do ni = 1,    3

                                      xin(341) = xin(359) + dxij*xin(332)
                                      yin(341) = yin(359) + dyij*yin(332)
                                      zin(341) = zin(359) + dzij*zin(332)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  368

                                      ! ni =    2

                                      xin(368) = xin(386) + dxij*xin(359)
                                      yin(368) = yin(386) + dyij*yin(359)
                                      zin(368) = zin(386) + dzij*zin(359)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  395

                                      ! ni =    3

                                      xin(395) = xin(413) + dxij*xin(386)
                                      yin(395) = yin(413) + dyij*yin(386)
                                      zin(395) = zin(413) + dzij*zin(386)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  422

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  350

                                      ! nj =    2

                                      ! i4 = i3 =  350

                                      ! do ni = 1,    3

                                      xin(350) = xin(368) + dxij*xin(341)
                                      yin(350) = yin(368) + dyij*yin(341)
                                      zin(350) = zin(368) + dzij*zin(341)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  377

                                      ! ni =    2

                                      xin(377) = xin(395) + dxij*xin(368)
                                      yin(377) = yin(395) + dyij*yin(368)
                                      zin(377) = zin(395) + dzij*zin(368)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  404

                                      ! ni =    3

                                      xin(404) = xin(422) + dxij*xin(395)
                                      yin(404) = yin(422) + dyij*yin(395)
                                      zin(404) = zin(422) + dzij*zin(395)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  431

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  359

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  432

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  423

                                      xin(432) = xin(432) + dxij*xin(423)
                                      yin(432) = yin(432) + dyij*yin(423)
                                      zin(432) = zin(432) + dzij*zin(423)

                                      ! i3 = i4 =  423
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  414

                                      xin(423) = xin(423) + dxij*xin(414)
                                      yin(423) = yin(423) + dyij*yin(414)
                                      zin(423) = zin(423) + dzij*zin(414)

                                      ! i3 = i4 =  414
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  432

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  423

                                      xin(432) = xin(432) + dxij*xin(423)
                                      yin(432) = yin(432) + dyij*yin(423)
                                      zin(432) = zin(432) + dzij*zin(423)

                                      ! i3 = i4 =  423
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  342

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  342

                                      ! do ni = 1,    3

                                      xin(342) = xin(360) + dxij*xin(333)
                                      yin(342) = yin(360) + dyij*yin(333)
                                      zin(342) = zin(360) + dzij*zin(333)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  369

                                      ! ni =    2

                                      xin(369) = xin(387) + dxij*xin(360)
                                      yin(369) = yin(387) + dyij*yin(360)
                                      zin(369) = zin(387) + dzij*zin(360)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  396

                                      ! ni =    3

                                      xin(396) = xin(414) + dxij*xin(387)
                                      yin(396) = yin(414) + dyij*yin(387)
                                      zin(396) = zin(414) + dzij*zin(387)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  423

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  351

                                      ! nj =    2

                                      ! i4 = i3 =  351

                                      ! do ni = 1,    3

                                      xin(351) = xin(369) + dxij*xin(342)
                                      yin(351) = yin(369) + dyij*yin(342)
                                      zin(351) = zin(369) + dzij*zin(342)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  378

                                      ! ni =    2

                                      xin(378) = xin(396) + dxij*xin(369)
                                      yin(378) = yin(396) + dyij*yin(369)
                                      zin(378) = zin(396) + dzij*zin(369)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  405

                                      ! ni =    3

                                      xin(405) = xin(423) + dxij*xin(396)
                                      yin(405) = yin(423) + dyij*yin(396)
                                      zin(405) = zin(423) + dzij*zin(396)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  432

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  360

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  325

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  333

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  332

                                      xin(333) = xin(333) + dxkl*xin(332)
                                      yin(333) = yin(333) + dykl*yin(332)
                                      zin(333) = zin(333) + dzkl*zin(332)

                                      ! i3 = i4 =  332
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  331

                                      xin(332) = xin(332) + dxkl*xin(331)
                                      yin(332) = yin(332) + dykl*yin(331)
                                      zin(332) = zin(332) + dzkl*zin(331)

                                      ! i3 = i4 =  331
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  333

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  332

                                      xin(333) = xin(333) + dxkl*xin(332)
                                      yin(333) = yin(333) + dykl*yin(332)
                                      zin(333) = zin(333) + dzkl*zin(332)

                                      ! i3 = i4 =  332
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  326

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  326

                                      ! do nk = 1,    2

                                      xin(326) = xin(328) + dxkl*xin(325)
                                      yin(326) = yin(328) + dykl*yin(325)
                                      zin(326) = zin(328) + dzkl*zin(325)
                                      ! i4 = i4 + lang+1 =  329

                                      ! nk =    2

                                      xin(329) = xin(331) + dxkl*xin(328)
                                      yin(329) = yin(331) + dykl*yin(328)
                                      zin(329) = zin(331) + dzkl*zin(328)
                                      ! i4 = i4 + lang+1 =  332

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  327

                                      ! nl =    2

                                      ! i4 = i3 =  327

                                      ! do nk = 1,    2

                                      xin(327) = xin(329) + dxkl*xin(326)
                                      yin(327) = yin(329) + dykl*yin(326)
                                      zin(327) = zin(329) + dzkl*zin(326)
                                      ! i4 = i4 + lang+1 =  330

                                      ! nk =    2

                                      xin(330) = xin(332) + dxkl*xin(329)
                                      yin(330) = yin(332) + dykl*yin(329)
                                      zin(330) = zin(332) + dzkl*zin(329)
                                      ! i4 = i4 + lang+1 =  333

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  328

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  334

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  342

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  341

                                      xin(342) = xin(342) + dxkl*xin(341)
                                      yin(342) = yin(342) + dykl*yin(341)
                                      zin(342) = zin(342) + dzkl*zin(341)

                                      ! i3 = i4 =  341
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  340

                                      xin(341) = xin(341) + dxkl*xin(340)
                                      yin(341) = yin(341) + dykl*yin(340)
                                      zin(341) = zin(341) + dzkl*zin(340)

                                      ! i3 = i4 =  340
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  342

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  341

                                      xin(342) = xin(342) + dxkl*xin(341)
                                      yin(342) = yin(342) + dykl*yin(341)
                                      zin(342) = zin(342) + dzkl*zin(341)

                                      ! i3 = i4 =  341
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  335

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  335

                                      ! do nk = 1,    2

                                      xin(335) = xin(337) + dxkl*xin(334)
                                      yin(335) = yin(337) + dykl*yin(334)
                                      zin(335) = zin(337) + dzkl*zin(334)
                                      ! i4 = i4 + lang+1 =  338

                                      ! nk =    2

                                      xin(338) = xin(340) + dxkl*xin(337)
                                      yin(338) = yin(340) + dykl*yin(337)
                                      zin(338) = zin(340) + dzkl*zin(337)
                                      ! i4 = i4 + lang+1 =  341

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  336

                                      ! nl =    2

                                      ! i4 = i3 =  336

                                      ! do nk = 1,    2

                                      xin(336) = xin(338) + dxkl*xin(335)
                                      yin(336) = yin(338) + dykl*yin(335)
                                      zin(336) = zin(338) + dzkl*zin(335)
                                      ! i4 = i4 + lang+1 =  339

                                      ! nk =    2

                                      xin(339) = xin(341) + dxkl*xin(338)
                                      yin(339) = yin(341) + dykl*yin(338)
                                      zin(339) = zin(341) + dzkl*zin(338)
                                      ! i4 = i4 + lang+1 =  342

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  337

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  343

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  351

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  350

                                      xin(351) = xin(351) + dxkl*xin(350)
                                      yin(351) = yin(351) + dykl*yin(350)
                                      zin(351) = zin(351) + dzkl*zin(350)

                                      ! i3 = i4 =  350
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  349

                                      xin(350) = xin(350) + dxkl*xin(349)
                                      yin(350) = yin(350) + dykl*yin(349)
                                      zin(350) = zin(350) + dzkl*zin(349)

                                      ! i3 = i4 =  349
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  351

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  350

                                      xin(351) = xin(351) + dxkl*xin(350)
                                      yin(351) = yin(351) + dykl*yin(350)
                                      zin(351) = zin(351) + dzkl*zin(350)

                                      ! i3 = i4 =  350
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  344

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  344

                                      ! do nk = 1,    2

                                      xin(344) = xin(346) + dxkl*xin(343)
                                      yin(344) = yin(346) + dykl*yin(343)
                                      zin(344) = zin(346) + dzkl*zin(343)
                                      ! i4 = i4 + lang+1 =  347

                                      ! nk =    2

                                      xin(347) = xin(349) + dxkl*xin(346)
                                      yin(347) = yin(349) + dykl*yin(346)
                                      zin(347) = zin(349) + dzkl*zin(346)
                                      ! i4 = i4 + lang+1 =  350

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  345

                                      ! nl =    2

                                      ! i4 = i3 =  345

                                      ! do nk = 1,    2

                                      xin(345) = xin(347) + dxkl*xin(344)
                                      yin(345) = yin(347) + dykl*yin(344)
                                      zin(345) = zin(347) + dzkl*zin(344)
                                      ! i4 = i4 + lang+1 =  348

                                      ! nk =    2

                                      xin(348) = xin(350) + dxkl*xin(347)
                                      yin(348) = yin(350) + dykl*yin(347)
                                      zin(348) = zin(350) + dzkl*zin(347)
                                      ! i4 = i4 + lang+1 =  351

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  346

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  352

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  352

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  358

                                      xin(359) = xin(359) + dxkl*xin(358)
                                      yin(359) = yin(359) + dykl*yin(358)
                                      zin(359) = zin(359) + dzkl*zin(358)

                                      ! i3 = i4 =  358
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  360

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  359

                                      xin(360) = xin(360) + dxkl*xin(359)
                                      yin(360) = yin(360) + dykl*yin(359)
                                      zin(360) = zin(360) + dzkl*zin(359)

                                      ! i3 = i4 =  359
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  353

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  353

                                      ! do nk = 1,    2

                                      xin(353) = xin(355) + dxkl*xin(352)
                                      yin(353) = yin(355) + dykl*yin(352)
                                      zin(353) = zin(355) + dzkl*zin(352)
                                      ! i4 = i4 + lang+1 =  356

                                      ! nk =    2

                                      xin(356) = xin(358) + dxkl*xin(355)
                                      yin(356) = yin(358) + dykl*yin(355)
                                      zin(356) = zin(358) + dzkl*zin(355)
                                      ! i4 = i4 + lang+1 =  359

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  354

                                      ! nl =    2

                                      ! i4 = i3 =  354

                                      ! do nk = 1,    2

                                      xin(354) = xin(356) + dxkl*xin(353)
                                      yin(354) = yin(356) + dykl*yin(353)
                                      zin(354) = zin(356) + dzkl*zin(353)
                                      ! i4 = i4 + lang+1 =  357

                                      ! nk =    2

                                      xin(357) = xin(359) + dxkl*xin(356)
                                      yin(357) = yin(359) + dykl*yin(356)
                                      zin(357) = zin(359) + dzkl*zin(356)
                                      ! i4 = i4 + lang+1 =  360

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  355

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  361

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  369

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  368

                                      xin(369) = xin(369) + dxkl*xin(368)
                                      yin(369) = yin(369) + dykl*yin(368)
                                      zin(369) = zin(369) + dzkl*zin(368)

                                      ! i3 = i4 =  368
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  367

                                      xin(368) = xin(368) + dxkl*xin(367)
                                      yin(368) = yin(368) + dykl*yin(367)
                                      zin(368) = zin(368) + dzkl*zin(367)

                                      ! i3 = i4 =  367
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  369

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  368

                                      xin(369) = xin(369) + dxkl*xin(368)
                                      yin(369) = yin(369) + dykl*yin(368)
                                      zin(369) = zin(369) + dzkl*zin(368)

                                      ! i3 = i4 =  368
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  362

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  362

                                      ! do nk = 1,    2

                                      xin(362) = xin(364) + dxkl*xin(361)
                                      yin(362) = yin(364) + dykl*yin(361)
                                      zin(362) = zin(364) + dzkl*zin(361)
                                      ! i4 = i4 + lang+1 =  365

                                      ! nk =    2

                                      xin(365) = xin(367) + dxkl*xin(364)
                                      yin(365) = yin(367) + dykl*yin(364)
                                      zin(365) = zin(367) + dzkl*zin(364)
                                      ! i4 = i4 + lang+1 =  368

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  363

                                      ! nl =    2

                                      ! i4 = i3 =  363

                                      ! do nk = 1,    2

                                      xin(363) = xin(365) + dxkl*xin(362)
                                      yin(363) = yin(365) + dykl*yin(362)
                                      zin(363) = zin(365) + dzkl*zin(362)
                                      ! i4 = i4 + lang+1 =  366

                                      ! nk =    2

                                      xin(366) = xin(368) + dxkl*xin(365)
                                      yin(366) = yin(368) + dykl*yin(365)
                                      zin(366) = zin(368) + dzkl*zin(365)
                                      ! i4 = i4 + lang+1 =  369

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  364

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  370

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  378

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  377

                                      xin(378) = xin(378) + dxkl*xin(377)
                                      yin(378) = yin(378) + dykl*yin(377)
                                      zin(378) = zin(378) + dzkl*zin(377)

                                      ! i3 = i4 =  377
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  376

                                      xin(377) = xin(377) + dxkl*xin(376)
                                      yin(377) = yin(377) + dykl*yin(376)
                                      zin(377) = zin(377) + dzkl*zin(376)

                                      ! i3 = i4 =  376
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  378

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  377

                                      xin(378) = xin(378) + dxkl*xin(377)
                                      yin(378) = yin(378) + dykl*yin(377)
                                      zin(378) = zin(378) + dzkl*zin(377)

                                      ! i3 = i4 =  377
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  371

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  371

                                      ! do nk = 1,    2

                                      xin(371) = xin(373) + dxkl*xin(370)
                                      yin(371) = yin(373) + dykl*yin(370)
                                      zin(371) = zin(373) + dzkl*zin(370)
                                      ! i4 = i4 + lang+1 =  374

                                      ! nk =    2

                                      xin(374) = xin(376) + dxkl*xin(373)
                                      yin(374) = yin(376) + dykl*yin(373)
                                      zin(374) = zin(376) + dzkl*zin(373)
                                      ! i4 = i4 + lang+1 =  377

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  372

                                      ! nl =    2

                                      ! i4 = i3 =  372

                                      ! do nk = 1,    2

                                      xin(372) = xin(374) + dxkl*xin(371)
                                      yin(372) = yin(374) + dykl*yin(371)
                                      zin(372) = zin(374) + dzkl*zin(371)
                                      ! i4 = i4 + lang+1 =  375

                                      ! nk =    2

                                      xin(375) = xin(377) + dxkl*xin(374)
                                      yin(375) = yin(377) + dykl*yin(374)
                                      zin(375) = zin(377) + dzkl*zin(374)
                                      ! i4 = i4 + lang+1 =  378

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  373

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  379

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  379

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  387

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  386

                                      xin(387) = xin(387) + dxkl*xin(386)
                                      yin(387) = yin(387) + dykl*yin(386)
                                      zin(387) = zin(387) + dzkl*zin(386)

                                      ! i3 = i4 =  386
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  385

                                      xin(386) = xin(386) + dxkl*xin(385)
                                      yin(386) = yin(386) + dykl*yin(385)
                                      zin(386) = zin(386) + dzkl*zin(385)

                                      ! i3 = i4 =  385
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  387

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  386

                                      xin(387) = xin(387) + dxkl*xin(386)
                                      yin(387) = yin(387) + dykl*yin(386)
                                      zin(387) = zin(387) + dzkl*zin(386)

                                      ! i3 = i4 =  386
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  380

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  380

                                      ! do nk = 1,    2

                                      xin(380) = xin(382) + dxkl*xin(379)
                                      yin(380) = yin(382) + dykl*yin(379)
                                      zin(380) = zin(382) + dzkl*zin(379)
                                      ! i4 = i4 + lang+1 =  383

                                      ! nk =    2

                                      xin(383) = xin(385) + dxkl*xin(382)
                                      yin(383) = yin(385) + dykl*yin(382)
                                      zin(383) = zin(385) + dzkl*zin(382)
                                      ! i4 = i4 + lang+1 =  386

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  381

                                      ! nl =    2

                                      ! i4 = i3 =  381

                                      ! do nk = 1,    2

                                      xin(381) = xin(383) + dxkl*xin(380)
                                      yin(381) = yin(383) + dykl*yin(380)
                                      zin(381) = zin(383) + dzkl*zin(380)
                                      ! i4 = i4 + lang+1 =  384

                                      ! nk =    2

                                      xin(384) = xin(386) + dxkl*xin(383)
                                      yin(384) = yin(386) + dykl*yin(383)
                                      zin(384) = zin(386) + dzkl*zin(383)
                                      ! i4 = i4 + lang+1 =  387

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  382

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  388

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  396

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  395

                                      xin(396) = xin(396) + dxkl*xin(395)
                                      yin(396) = yin(396) + dykl*yin(395)
                                      zin(396) = zin(396) + dzkl*zin(395)

                                      ! i3 = i4 =  395
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  394

                                      xin(395) = xin(395) + dxkl*xin(394)
                                      yin(395) = yin(395) + dykl*yin(394)
                                      zin(395) = zin(395) + dzkl*zin(394)

                                      ! i3 = i4 =  394
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  396

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  395

                                      xin(396) = xin(396) + dxkl*xin(395)
                                      yin(396) = yin(396) + dykl*yin(395)
                                      zin(396) = zin(396) + dzkl*zin(395)

                                      ! i3 = i4 =  395
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  389

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  389

                                      ! do nk = 1,    2

                                      xin(389) = xin(391) + dxkl*xin(388)
                                      yin(389) = yin(391) + dykl*yin(388)
                                      zin(389) = zin(391) + dzkl*zin(388)
                                      ! i4 = i4 + lang+1 =  392

                                      ! nk =    2

                                      xin(392) = xin(394) + dxkl*xin(391)
                                      yin(392) = yin(394) + dykl*yin(391)
                                      zin(392) = zin(394) + dzkl*zin(391)
                                      ! i4 = i4 + lang+1 =  395

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  390

                                      ! nl =    2

                                      ! i4 = i3 =  390

                                      ! do nk = 1,    2

                                      xin(390) = xin(392) + dxkl*xin(389)
                                      yin(390) = yin(392) + dykl*yin(389)
                                      zin(390) = zin(392) + dzkl*zin(389)
                                      ! i4 = i4 + lang+1 =  393

                                      ! nk =    2

                                      xin(393) = xin(395) + dxkl*xin(392)
                                      yin(393) = yin(395) + dykl*yin(392)
                                      zin(393) = zin(395) + dzkl*zin(392)
                                      ! i4 = i4 + lang+1 =  396

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  391

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  397

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  405

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  404

                                      xin(405) = xin(405) + dxkl*xin(404)
                                      yin(405) = yin(405) + dykl*yin(404)
                                      zin(405) = zin(405) + dzkl*zin(404)

                                      ! i3 = i4 =  404
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  403

                                      xin(404) = xin(404) + dxkl*xin(403)
                                      yin(404) = yin(404) + dykl*yin(403)
                                      zin(404) = zin(404) + dzkl*zin(403)

                                      ! i3 = i4 =  403
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  405

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  404

                                      xin(405) = xin(405) + dxkl*xin(404)
                                      yin(405) = yin(405) + dykl*yin(404)
                                      zin(405) = zin(405) + dzkl*zin(404)

                                      ! i3 = i4 =  404
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  398

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  398

                                      ! do nk = 1,    2

                                      xin(398) = xin(400) + dxkl*xin(397)
                                      yin(398) = yin(400) + dykl*yin(397)
                                      zin(398) = zin(400) + dzkl*zin(397)
                                      ! i4 = i4 + lang+1 =  401

                                      ! nk =    2

                                      xin(401) = xin(403) + dxkl*xin(400)
                                      yin(401) = yin(403) + dykl*yin(400)
                                      zin(401) = zin(403) + dzkl*zin(400)
                                      ! i4 = i4 + lang+1 =  404

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  399

                                      ! nl =    2

                                      ! i4 = i3 =  399

                                      ! do nk = 1,    2

                                      xin(399) = xin(401) + dxkl*xin(398)
                                      yin(399) = yin(401) + dykl*yin(398)
                                      zin(399) = zin(401) + dzkl*zin(398)
                                      ! i4 = i4 + lang+1 =  402

                                      ! nk =    2

                                      xin(402) = xin(404) + dxkl*xin(401)
                                      yin(402) = yin(404) + dykl*yin(401)
                                      zin(402) = zin(404) + dzkl*zin(401)
                                      ! i4 = i4 + lang+1 =  405

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  400

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  406

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  406

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  414

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  413

                                      xin(414) = xin(414) + dxkl*xin(413)
                                      yin(414) = yin(414) + dykl*yin(413)
                                      zin(414) = zin(414) + dzkl*zin(413)

                                      ! i3 = i4 =  413
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  412

                                      xin(413) = xin(413) + dxkl*xin(412)
                                      yin(413) = yin(413) + dykl*yin(412)
                                      zin(413) = zin(413) + dzkl*zin(412)

                                      ! i3 = i4 =  412
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  414

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  413

                                      xin(414) = xin(414) + dxkl*xin(413)
                                      yin(414) = yin(414) + dykl*yin(413)
                                      zin(414) = zin(414) + dzkl*zin(413)

                                      ! i3 = i4 =  413
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  407

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  407

                                      ! do nk = 1,    2

                                      xin(407) = xin(409) + dxkl*xin(406)
                                      yin(407) = yin(409) + dykl*yin(406)
                                      zin(407) = zin(409) + dzkl*zin(406)
                                      ! i4 = i4 + lang+1 =  410

                                      ! nk =    2

                                      xin(410) = xin(412) + dxkl*xin(409)
                                      yin(410) = yin(412) + dykl*yin(409)
                                      zin(410) = zin(412) + dzkl*zin(409)
                                      ! i4 = i4 + lang+1 =  413

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  408

                                      ! nl =    2

                                      ! i4 = i3 =  408

                                      ! do nk = 1,    2

                                      xin(408) = xin(410) + dxkl*xin(407)
                                      yin(408) = yin(410) + dykl*yin(407)
                                      zin(408) = zin(410) + dzkl*zin(407)
                                      ! i4 = i4 + lang+1 =  411

                                      ! nk =    2

                                      xin(411) = xin(413) + dxkl*xin(410)
                                      yin(411) = yin(413) + dykl*yin(410)
                                      zin(411) = zin(413) + dzkl*zin(410)
                                      ! i4 = i4 + lang+1 =  414

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  409

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  415

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  423

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  422

                                      xin(423) = xin(423) + dxkl*xin(422)
                                      yin(423) = yin(423) + dykl*yin(422)
                                      zin(423) = zin(423) + dzkl*zin(422)

                                      ! i3 = i4 =  422
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  421

                                      xin(422) = xin(422) + dxkl*xin(421)
                                      yin(422) = yin(422) + dykl*yin(421)
                                      zin(422) = zin(422) + dzkl*zin(421)

                                      ! i3 = i4 =  421
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  423

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  422

                                      xin(423) = xin(423) + dxkl*xin(422)
                                      yin(423) = yin(423) + dykl*yin(422)
                                      zin(423) = zin(423) + dzkl*zin(422)

                                      ! i3 = i4 =  422
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  416

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  416

                                      ! do nk = 1,    2

                                      xin(416) = xin(418) + dxkl*xin(415)
                                      yin(416) = yin(418) + dykl*yin(415)
                                      zin(416) = zin(418) + dzkl*zin(415)
                                      ! i4 = i4 + lang+1 =  419

                                      ! nk =    2

                                      xin(419) = xin(421) + dxkl*xin(418)
                                      yin(419) = yin(421) + dykl*yin(418)
                                      zin(419) = zin(421) + dzkl*zin(418)
                                      ! i4 = i4 + lang+1 =  422

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  417

                                      ! nl =    2

                                      ! i4 = i3 =  417

                                      ! do nk = 1,    2

                                      xin(417) = xin(419) + dxkl*xin(416)
                                      yin(417) = yin(419) + dykl*yin(416)
                                      zin(417) = zin(419) + dzkl*zin(416)
                                      ! i4 = i4 + lang+1 =  420

                                      ! nk =    2

                                      xin(420) = xin(422) + dxkl*xin(419)
                                      yin(420) = yin(422) + dykl*yin(419)
                                      zin(420) = zin(422) + dzkl*zin(419)
                                      ! i4 = i4 + lang+1 =  423

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  418

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  424

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  432

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  431

                                      xin(432) = xin(432) + dxkl*xin(431)
                                      yin(432) = yin(432) + dykl*yin(431)
                                      zin(432) = zin(432) + dzkl*zin(431)

                                      ! i3 = i4 =  431
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  430

                                      xin(431) = xin(431) + dxkl*xin(430)
                                      yin(431) = yin(431) + dykl*yin(430)
                                      zin(431) = zin(431) + dzkl*zin(430)

                                      ! i3 = i4 =  430
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  432

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  431

                                      xin(432) = xin(432) + dxkl*xin(431)
                                      yin(432) = yin(432) + dykl*yin(431)
                                      zin(432) = zin(432) + dzkl*zin(431)

                                      ! i3 = i4 =  431
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  425

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  425

                                      ! do nk = 1,    2

                                      xin(425) = xin(427) + dxkl*xin(424)
                                      yin(425) = yin(427) + dykl*yin(424)
                                      zin(425) = zin(427) + dzkl*zin(424)
                                      ! i4 = i4 + lang+1 =  428

                                      ! nk =    2

                                      xin(428) = xin(430) + dxkl*xin(427)
                                      yin(428) = yin(430) + dykl*yin(427)
                                      zin(428) = zin(430) + dzkl*zin(427)
                                      ! i4 = i4 + lang+1 =  431

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  426

                                      ! nl =    2

                                      ! i4 = i3 =  426

                                      ! do nk = 1,    2

                                      xin(426) = xin(428) + dxkl*xin(425)
                                      yin(426) = yin(428) + dykl*yin(425)
                                      zin(426) = zin(428) + dzkl*zin(425)
                                      ! i4 = i4 + lang+1 =  429

                                      ! nk =    2

                                      xin(429) = xin(431) + dxkl*xin(428)
                                      yin(429) = yin(431) + dykl*yin(428)
                                      zin(429) = zin(431) + dzkl*zin(428)
                                      ! i4 = i4 + lang+1 =  432

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  427

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  433

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  433

                                      ! end do

                                      ! *** Now root =    5

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  432

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

                                      ! i1 = in(1) =  433

                                      xin(433) = 1.0_dp
                                      yin(433) = 1.0_dp
                                      zin(433) = f00

                                      ! i2 = in(2) =  460
                                      ! k2 = kn(2) =    3
                                      cp10 = b00

                                      ! ----- I(1,0) -----

                                      xin(460) = xc00
                                      yin(460) = yc00
                                      zin(460) = zc00*f00

                                      ! ----- I(0,1) -----

                                      ! i3 = i1+k2 =  436

                                      xin(436) = xcp00
                                      yin(436) = ycp00
                                      zin(436) = zcp00*f00

                                      ! ----- I(1,1) -----

                                      ! i3 = i2+k2 =  463
                                      ! i2 =  460

                                      xin(463) = xcp00*xin(460) + cp10
                                      yin(463) = ycp00*yin(460) + cp10
                                      zin(463) = zcp00*zin(460) + cp10*f00

                                      ! ----- I(N,0) -----

                                      c10 = 0.0_dp

                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  460

                                      ! do n = 2,   5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  487
                                      ! i3 =  433
                                      ! i4 =  460

                                      xin(487) = c10*xin(433) + xc00*xin(460)
                                      yin(487) = c10*yin(433) + yc00*yin(460)
                                      zin(487) = c10*zin(433) + zc00*zin(460)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  490
                                      ! i5 =  487
                                      ! i4 =  460

                                      xin(490) = xcp00*xin(487) + cp10*xin(460)
                                      yin(490) = ycp00*yin(487) + cp10*yin(460)
                                      zin(490) = zcp00*zin(487) + cp10*zin(460)

                                      ! ------------------

                                      ! i3 = i4 =  460
                                      ! i4 = i5 =  487

                                      ! n =    3

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  514
                                      ! i3 =  460
                                      ! i4 =  487

                                      xin(514) = c10*xin(460) + xc00*xin(487)
                                      yin(514) = c10*yin(460) + yc00*yin(487)
                                      zin(514) = c10*zin(460) + zc00*zin(487)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  517
                                      ! i5 =  514
                                      ! i4 =  487

                                      xin(517) = xcp00*xin(514) + cp10*xin(487)
                                      yin(517) = ycp00*yin(514) + cp10*yin(487)
                                      zin(517) = zcp00*zin(514) + cp10*zin(487)

                                      ! ------------------

                                      ! i3 = i4 =  487
                                      ! i4 = i5 =  514

                                      ! n =    4

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  523
                                      ! i3 =  487
                                      ! i4 =  514

                                      xin(523) = c10*xin(487) + xc00*xin(514)
                                      yin(523) = c10*yin(487) + yc00*yin(514)
                                      zin(523) = c10*zin(487) + zc00*zin(514)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  526
                                      ! i5 =  523
                                      ! i4 =  514

                                      xin(526) = xcp00*xin(523) + cp10*xin(514)
                                      yin(526) = ycp00*yin(523) + cp10*yin(514)
                                      zin(526) = zcp00*zin(523) + cp10*zin(514)

                                      ! ------------------

                                      ! i3 = i4 =  514
                                      ! i4 = i5 =  523

                                      ! n =    5

                                      c10 = c10 + b10

                                      ! i5 = in(n+1) =  532
                                      ! i3 =  514
                                      ! i4 =  523

                                      xin(532) = c10*xin(514) + xc00*xin(523)
                                      yin(532) = c10*yin(514) + yc00*yin(523)
                                      zin(532) = c10*zin(514) + zc00*zin(523)

                                      ! ----- I(N,1) -----

                                      cp10 = cp10 + b00

                                      ! i3 = i5 + k2 =  535
                                      ! i5 =  532
                                      ! i4 =  523

                                      xin(535) = xcp00*xin(532) + cp10*xin(523)
                                      yin(535) = ycp00*yin(532) + cp10*yin(523)
                                      zin(535) = zcp00*zin(532) + cp10*zin(523)

                                      ! ------------------

                                      ! i3 = i4 =  523
                                      ! i4 = i5 =  532

                                      ! n =    6

                                      ! end do

                                      ! ----- I(0,M) -----

                                      cp01 = 0.0_dp
                                      c01 = b00

                                      ! i3 = i1 =  433
                                      ! i4 = i1+k2 =  436

                                      ! do n = 2,    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  439
                                      ! i3 =  433
                                      ! i4 =  436

                                      xin(439) = cp01*xin(433) + xcp00*xin(436)
                                      yin(439) = cp01*yin(433) + ycp00*yin(436)
                                      zin(439) = cp01*zin(433) + zcp00*zin(436)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  466

                                      xin(466) = xc00*xin(439) + c01*xin(436)
                                      yin(466) = yc00*yin(439) + c01*yin(436)
                                      zin(466) = zc00*zin(439) + c01*zin(436)

                                      ! ------------------

                                      ! i3 = i4 =  436
                                      ! i4 = i5 =  439

                                      ! n =    3

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  440
                                      ! i3 =  436
                                      ! i4 =  439

                                      xin(440) = cp01*xin(436) + xcp00*xin(439)
                                      yin(440) = cp01*yin(436) + ycp00*yin(439)
                                      zin(440) = cp01*zin(436) + zcp00*zin(439)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  467

                                      xin(467) = xc00*xin(440) + c01*xin(439)
                                      yin(467) = yc00*yin(440) + c01*yin(439)
                                      zin(467) = zc00*zin(440) + c01*zin(439)

                                      ! ------------------

                                      ! i3 = i4 =  439
                                      ! i4 = i5 =  440

                                      ! n =    4

                                      cp01 = cp01 + bp01

                                      ! i5 = i1+kn(n+1) =  441
                                      ! i3 =  439
                                      ! i4 =  440

                                      xin(441) = cp01*xin(439) + xcp00*xin(440)
                                      yin(441) = cp01*yin(439) + ycp00*yin(440)
                                      zin(441) = cp01*zin(439) + zcp00*zin(440)

                                      ! ----- I(1,M) -----

                                      c01 = c01 + b00

                                      ! i3 = i2+kn(n+1) =  468

                                      xin(468) = xc00*xin(441) + c01*xin(440)
                                      yin(468) = yc00*yin(441) + c01*yin(440)
                                      zin(468) = zc00*zin(441) + c01*zin(440)

                                      ! ------------------

                                      ! i3 = i4 =  440
                                      ! i4 = i5 =  441

                                      ! n =    5

                                      ! end do

                                      ! ----- I(N,M) -----

                                      c01 = b00
                                      ! k3 = k2 =    3

                                      ! do n = 2,    4

                                      ! k4 = kn(n+1) =    6
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  460

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  487

                                      xin(493) = c10*xin(439) + xc00*xin(466) + c01*xin(463)
                                      yin(493) = c10*yin(439) + yc00*yin(466) + c01*yin(463)
                                      zin(493) = c10*zin(439) + zc00*zin(466) + c01*zin(463)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  460
                                      ! i4 = i5 =  487

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  514

                                      xin(520) = c10*xin(466) + xc00*xin(493) + c01*xin(490)
                                      yin(520) = c10*yin(466) + yc00*yin(493) + c01*yin(490)
                                      zin(520) = c10*zin(466) + zc00*zin(493) + c01*zin(490)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  487
                                      ! i4 = i5 =  514

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  523

                                      xin(529) = c10*xin(493) + xc00*xin(520) + c01*xin(517)
                                      yin(529) = c10*yin(493) + yc00*yin(520) + c01*yin(517)
                                      zin(529) = c10*zin(493) + zc00*zin(520) + c01*zin(517)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  514
                                      ! i4 = i5 =  523

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  532

                                      xin(538) = c10*xin(520) + xc00*xin(529) + c01*xin(526)
                                      yin(538) = c10*yin(520) + yc00*yin(529) + c01*yin(526)
                                      zin(538) = c10*zin(520) + zc00*zin(529) + c01*zin(526)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  523
                                      ! i4 = i5 =  532

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   6

                                      ! n =    3

                                      ! k4 = kn(n+1) =    7
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  460

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  487

                                      xin(494) = c10*xin(440) + xc00*xin(467) + c01*xin(466)
                                      yin(494) = c10*yin(440) + yc00*yin(467) + c01*yin(466)
                                      zin(494) = c10*zin(440) + zc00*zin(467) + c01*zin(466)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  460
                                      ! i4 = i5 =  487

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  514

                                      xin(521) = c10*xin(467) + xc00*xin(494) + c01*xin(493)
                                      yin(521) = c10*yin(467) + yc00*yin(494) + c01*yin(493)
                                      zin(521) = c10*zin(467) + zc00*zin(494) + c01*zin(493)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  487
                                      ! i4 = i5 =  514

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  523

                                      xin(530) = c10*xin(494) + xc00*xin(521) + c01*xin(520)
                                      yin(530) = c10*yin(494) + yc00*yin(521) + c01*yin(520)
                                      zin(530) = c10*zin(494) + zc00*zin(521) + c01*zin(520)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  514
                                      ! i4 = i5 =  523

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  532

                                      xin(539) = c10*xin(521) + xc00*xin(530) + c01*xin(529)
                                      yin(539) = c10*yin(521) + yc00*yin(530) + c01*yin(529)
                                      zin(539) = c10*zin(521) + zc00*zin(530) + c01*zin(529)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  523
                                      ! i4 = i5 =  532

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   7

                                      ! n =    4

                                      ! k4 = kn(n+1) =    8
                                      ! i3 = i1 =  433
                                      ! i4 = i2 =  460

                                      c01 = c01 + b00
                                      c10 = b10

                                      ! do nn = 2,    5

                                      ! i5 = in(nn+1) =  487

                                      xin(495) = c10*xin(441) + xc00*xin(468) + c01*xin(467)
                                      yin(495) = c10*yin(441) + yc00*yin(468) + c01*yin(467)
                                      zin(495) = c10*zin(441) + zc00*zin(468) + c01*zin(467)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  460
                                      ! i4 = i5 =  487

                                      ! nn =    3

                                      ! i5 = in(nn+1) =  514

                                      xin(522) = c10*xin(468) + xc00*xin(495) + c01*xin(494)
                                      yin(522) = c10*yin(468) + yc00*yin(495) + c01*yin(494)
                                      zin(522) = c10*zin(468) + zc00*zin(495) + c01*zin(494)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  487
                                      ! i4 = i5 =  514

                                      ! nn =    4

                                      ! i5 = in(nn+1) =  523

                                      xin(531) = c10*xin(495) + xc00*xin(522) + c01*xin(521)
                                      yin(531) = c10*yin(495) + yc00*yin(522) + c01*yin(521)
                                      zin(531) = c10*zin(495) + zc00*zin(522) + c01*zin(521)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  514
                                      ! i4 = i5 =  523

                                      ! nn =    5

                                      ! i5 = in(nn+1) =  532

                                      xin(540) = c10*xin(522) + xc00*xin(531) + c01*xin(530)
                                      yin(540) = c10*yin(522) + yc00*yin(531) + c01*yin(530)
                                      zin(540) = c10*zin(522) + zc00*zin(531) + c01*zin(530)

                                      c10 = c10 + b10

                                      ! i3 = i4 =  523
                                      ! i4 = i5 =  532

                                      ! nn =    6

                                      ! end do

                                      ! k3 = k4   8

                                      ! n =    5

                                      ! end do

                                      ! ----- I(NI,NJ,M) -----

                                      ! nm = 0
                                      ! i5 = in(iang+jang+1) =  532

                                      ! do while nm.le.(kang+lang)

                                      ! min = iang

                                      ! km = kn(nm+1) =    0

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  532

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  523

                                      xin(532) = xin(532) + dxij*xin(523)
                                      yin(532) = yin(532) + dyij*yin(523)
                                      zin(532) = zin(532) + dzij*zin(523)

                                      ! i3 = i4 =  523
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  514

                                      xin(523) = xin(523) + dxij*xin(514)
                                      yin(523) = yin(523) + dyij*yin(514)
                                      zin(523) = zin(523) + dzij*zin(514)

                                      ! i3 = i4 =  514
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  532

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  523

                                      xin(532) = xin(532) + dxij*xin(523)
                                      yin(532) = yin(532) + dyij*yin(523)
                                      zin(532) = zin(532) + dzij*zin(523)

                                      ! i3 = i4 =  523
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  442

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  442

                                      ! do ni = 1,    3

                                      xin(442) = xin(460) + dxij*xin(433)
                                      yin(442) = yin(460) + dyij*yin(433)
                                      zin(442) = zin(460) + dzij*zin(433)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  469

                                      ! ni =    2

                                      xin(469) = xin(487) + dxij*xin(460)
                                      yin(469) = yin(487) + dyij*yin(460)
                                      zin(469) = zin(487) + dzij*zin(460)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  496

                                      ! ni =    3

                                      xin(496) = xin(514) + dxij*xin(487)
                                      yin(496) = yin(514) + dyij*yin(487)
                                      zin(496) = zin(514) + dzij*zin(487)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  523

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  451

                                      ! nj =    2

                                      ! i4 = i3 =  451

                                      ! do ni = 1,    3

                                      xin(451) = xin(469) + dxij*xin(442)
                                      yin(451) = yin(469) + dyij*yin(442)
                                      zin(451) = zin(469) + dzij*zin(442)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  478

                                      ! ni =    2

                                      xin(478) = xin(496) + dxij*xin(469)
                                      yin(478) = yin(496) + dyij*yin(469)
                                      zin(478) = zin(496) + dzij*zin(469)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  505

                                      ! ni =    3

                                      xin(505) = xin(523) + dxij*xin(496)
                                      yin(505) = yin(523) + dyij*yin(496)
                                      zin(505) = zin(523) + dzij*zin(496)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  532

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  460

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    1

                                      ! min = iang

                                      ! km = kn(nm+1) =    3

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  535

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  526

                                      xin(535) = xin(535) + dxij*xin(526)
                                      yin(535) = yin(535) + dyij*yin(526)
                                      zin(535) = zin(535) + dzij*zin(526)

                                      ! i3 = i4 =  526
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  517

                                      xin(526) = xin(526) + dxij*xin(517)
                                      yin(526) = yin(526) + dyij*yin(517)
                                      zin(526) = zin(526) + dzij*zin(517)

                                      ! i3 = i4 =  517
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  535

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  526

                                      xin(535) = xin(535) + dxij*xin(526)
                                      yin(535) = yin(535) + dyij*yin(526)
                                      zin(535) = zin(535) + dzij*zin(526)

                                      ! i3 = i4 =  526
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  445

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  445

                                      ! do ni = 1,    3

                                      xin(445) = xin(463) + dxij*xin(436)
                                      yin(445) = yin(463) + dyij*yin(436)
                                      zin(445) = zin(463) + dzij*zin(436)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  472

                                      ! ni =    2

                                      xin(472) = xin(490) + dxij*xin(463)
                                      yin(472) = yin(490) + dyij*yin(463)
                                      zin(472) = zin(490) + dzij*zin(463)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  499

                                      ! ni =    3

                                      xin(499) = xin(517) + dxij*xin(490)
                                      yin(499) = yin(517) + dyij*yin(490)
                                      zin(499) = zin(517) + dzij*zin(490)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  526

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  454

                                      ! nj =    2

                                      ! i4 = i3 =  454

                                      ! do ni = 1,    3

                                      xin(454) = xin(472) + dxij*xin(445)
                                      yin(454) = yin(472) + dyij*yin(445)
                                      zin(454) = zin(472) + dzij*zin(445)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  481

                                      ! ni =    2

                                      xin(481) = xin(499) + dxij*xin(472)
                                      yin(481) = yin(499) + dyij*yin(472)
                                      zin(481) = zin(499) + dzij*zin(472)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  508

                                      ! ni =    3

                                      xin(508) = xin(526) + dxij*xin(499)
                                      yin(508) = yin(526) + dyij*yin(499)
                                      zin(508) = zin(526) + dzij*zin(499)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  535

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  463

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    2

                                      ! min = iang

                                      ! km = kn(nm+1) =    6

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  538

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  529

                                      xin(538) = xin(538) + dxij*xin(529)
                                      yin(538) = yin(538) + dyij*yin(529)
                                      zin(538) = zin(538) + dzij*zin(529)

                                      ! i3 = i4 =  529
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  520

                                      xin(529) = xin(529) + dxij*xin(520)
                                      yin(529) = yin(529) + dyij*yin(520)
                                      zin(529) = zin(529) + dzij*zin(520)

                                      ! i3 = i4 =  520
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  538

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  529

                                      xin(538) = xin(538) + dxij*xin(529)
                                      yin(538) = yin(538) + dyij*yin(529)
                                      zin(538) = zin(538) + dzij*zin(529)

                                      ! i3 = i4 =  529
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  448

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  448

                                      ! do ni = 1,    3

                                      xin(448) = xin(466) + dxij*xin(439)
                                      yin(448) = yin(466) + dyij*yin(439)
                                      zin(448) = zin(466) + dzij*zin(439)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  475

                                      ! ni =    2

                                      xin(475) = xin(493) + dxij*xin(466)
                                      yin(475) = yin(493) + dyij*yin(466)
                                      zin(475) = zin(493) + dzij*zin(466)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  502

                                      ! ni =    3

                                      xin(502) = xin(520) + dxij*xin(493)
                                      yin(502) = yin(520) + dyij*yin(493)
                                      zin(502) = zin(520) + dzij*zin(493)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  529

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  457

                                      ! nj =    2

                                      ! i4 = i3 =  457

                                      ! do ni = 1,    3

                                      xin(457) = xin(475) + dxij*xin(448)
                                      yin(457) = yin(475) + dyij*yin(448)
                                      zin(457) = zin(475) + dzij*zin(448)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  484

                                      ! ni =    2

                                      xin(484) = xin(502) + dxij*xin(475)
                                      yin(484) = yin(502) + dyij*yin(475)
                                      zin(484) = zin(502) + dzij*zin(475)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  511

                                      ! ni =    3

                                      xin(511) = xin(529) + dxij*xin(502)
                                      yin(511) = yin(529) + dyij*yin(502)
                                      zin(511) = zin(529) + dzij*zin(502)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  538

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  466

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    3

                                      ! min = iang

                                      ! km = kn(nm+1) =    7

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  539

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  530

                                      xin(539) = xin(539) + dxij*xin(530)
                                      yin(539) = yin(539) + dyij*yin(530)
                                      zin(539) = zin(539) + dzij*zin(530)

                                      ! i3 = i4 =  530
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  521

                                      xin(530) = xin(530) + dxij*xin(521)
                                      yin(530) = yin(530) + dyij*yin(521)
                                      zin(530) = zin(530) + dzij*zin(521)

                                      ! i3 = i4 =  521
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  539

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  530

                                      xin(539) = xin(539) + dxij*xin(530)
                                      yin(539) = yin(539) + dyij*yin(530)
                                      zin(539) = zin(539) + dzij*zin(530)

                                      ! i3 = i4 =  530
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  449

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  449

                                      ! do ni = 1,    3

                                      xin(449) = xin(467) + dxij*xin(440)
                                      yin(449) = yin(467) + dyij*yin(440)
                                      zin(449) = zin(467) + dzij*zin(440)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  476

                                      ! ni =    2

                                      xin(476) = xin(494) + dxij*xin(467)
                                      yin(476) = yin(494) + dyij*yin(467)
                                      zin(476) = zin(494) + dzij*zin(467)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  503

                                      ! ni =    3

                                      xin(503) = xin(521) + dxij*xin(494)
                                      yin(503) = yin(521) + dyij*yin(494)
                                      zin(503) = zin(521) + dzij*zin(494)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  530

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  458

                                      ! nj =    2

                                      ! i4 = i3 =  458

                                      ! do ni = 1,    3

                                      xin(458) = xin(476) + dxij*xin(449)
                                      yin(458) = yin(476) + dyij*yin(449)
                                      zin(458) = zin(476) + dzij*zin(449)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  485

                                      ! ni =    2

                                      xin(485) = xin(503) + dxij*xin(476)
                                      yin(485) = yin(503) + dyij*yin(476)
                                      zin(485) = zin(503) + dzij*zin(476)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  512

                                      ! ni =    3

                                      xin(512) = xin(530) + dxij*xin(503)
                                      yin(512) = yin(530) + dyij*yin(503)
                                      zin(512) = zin(530) + dzij*zin(503)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  539

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  467

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    4

                                      ! min = iang

                                      ! km = kn(nm+1) =    8

                                      ! do while min.lt.(iang+jang)

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  540

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  531

                                      xin(540) = xin(540) + dxij*xin(531)
                                      yin(540) = yin(540) + dyij*yin(531)
                                      zin(540) = zin(540) + dzij*zin(531)

                                      ! i3 = i4 =  531
                                      ! nn = nn-1 =    4

                                      ! i4 = in(nn)+km =  522

                                      xin(531) = xin(531) + dxij*xin(522)
                                      yin(531) = yin(531) + dyij*yin(522)
                                      zin(531) = zin(531) + dzij*zin(522)

                                      ! i3 = i4 =  522
                                      ! nn = nn-1 =    3

                                      ! end do

                                      ! min = min + 1

                                      ! nn = (iang+jang) =    5

                                      ! i3 = i5 + km =  540

                                      ! do while nn.gt.min

                                      ! i4 = in(nn)+km =  531

                                      xin(540) = xin(540) + dxij*xin(531)
                                      yin(540) = yin(540) + dyij*yin(531)
                                      zin(540) = zin(540) + dzij*zin(531)

                                      ! i3 = i4 =  531
                                      ! nn = nn-1 =    4

                                      ! end do

                                      ! min = min + 1

                                      ! end do

                                      ! i3 = km + i1 + (kang+1)*(lang+1) =  450

                                      ! do nj = 1,    2

                                      ! i4 = i3 =  450

                                      ! do ni = 1,    3

                                      xin(450) = xin(468) + dxij*xin(441)
                                      yin(450) = yin(468) + dyij*yin(441)
                                      zin(450) = zin(468) + dzij*zin(441)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  477

                                      ! ni =    2

                                      xin(477) = xin(495) + dxij*xin(468)
                                      yin(477) = yin(495) + dyij*yin(468)
                                      zin(477) = zin(495) + dzij*zin(468)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  504

                                      ! ni =    3

                                      xin(504) = xin(522) + dxij*xin(495)
                                      yin(504) = yin(522) + dyij*yin(495)
                                      zin(504) = zin(522) + dzij*zin(495)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  531

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  459

                                      ! nj =    2

                                      ! i4 = i3 =  459

                                      ! do ni = 1,    3

                                      xin(459) = xin(477) + dxij*xin(450)
                                      yin(459) = yin(477) + dyij*yin(450)
                                      zin(459) = zin(477) + dzij*zin(450)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  486

                                      ! ni =    2

                                      xin(486) = xin(504) + dxij*xin(477)
                                      yin(486) = yin(504) + dyij*yin(477)
                                      zin(486) = zin(504) + dzij*zin(477)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  513

                                      ! ni =    3

                                      xin(513) = xin(531) + dxij*xin(504)
                                      yin(513) = yin(531) + dyij*yin(504)
                                      zin(513) = zin(531) + dzij*zin(504)

                                      ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  540

                                      ! ni =    4

                                      ! end do

                                      ! i3 = i3 + (kang+1)*(lang+1) =  468

                                      ! nj =    3

                                      ! end do

                                      ! nm = nm + 1 =    5

                                      ! end do

                                      ! ----- I(NI,NJ,NK,NL) -----

                                      ! i5 = kn(kang+lang+1) =    8

                                      ! iaa = i1 =  433

                                      ! ni = 0

                                      ! do while ni.le.iang

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  441

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  440

                                      xin(441) = xin(441) + dxkl*xin(440)
                                      yin(441) = yin(441) + dykl*yin(440)
                                      zin(441) = zin(441) + dzkl*zin(440)

                                      ! i3 = i4 =  440
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  439

                                      xin(440) = xin(440) + dxkl*xin(439)
                                      yin(440) = yin(440) + dykl*yin(439)
                                      zin(440) = zin(440) + dzkl*zin(439)

                                      ! i3 = i4 =  439
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  441

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  440

                                      xin(441) = xin(441) + dxkl*xin(440)
                                      yin(441) = yin(441) + dykl*yin(440)
                                      zin(441) = zin(441) + dzkl*zin(440)

                                      ! i3 = i4 =  440
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  434

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  434

                                      ! do nk = 1,    2

                                      xin(434) = xin(436) + dxkl*xin(433)
                                      yin(434) = yin(436) + dykl*yin(433)
                                      zin(434) = zin(436) + dzkl*zin(433)
                                      ! i4 = i4 + lang+1 =  437

                                      ! nk =    2

                                      xin(437) = xin(439) + dxkl*xin(436)
                                      yin(437) = yin(439) + dykl*yin(436)
                                      zin(437) = zin(439) + dzkl*zin(436)
                                      ! i4 = i4 + lang+1 =  440

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  435

                                      ! nl =    2

                                      ! i4 = i3 =  435

                                      ! do nk = 1,    2

                                      xin(435) = xin(437) + dxkl*xin(434)
                                      yin(435) = yin(437) + dykl*yin(434)
                                      zin(435) = zin(437) + dzkl*zin(434)
                                      ! i4 = i4 + lang+1 =  438

                                      ! nk =    2

                                      xin(438) = xin(440) + dxkl*xin(437)
                                      yin(438) = yin(440) + dykl*yin(437)
                                      zin(438) = zin(440) + dzkl*zin(437)
                                      ! i4 = i4 + lang+1 =  441

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  436

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  442

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  450

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  449

                                      xin(450) = xin(450) + dxkl*xin(449)
                                      yin(450) = yin(450) + dykl*yin(449)
                                      zin(450) = zin(450) + dzkl*zin(449)

                                      ! i3 = i4 =  449
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  448

                                      xin(449) = xin(449) + dxkl*xin(448)
                                      yin(449) = yin(449) + dykl*yin(448)
                                      zin(449) = zin(449) + dzkl*zin(448)

                                      ! i3 = i4 =  448
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  450

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  449

                                      xin(450) = xin(450) + dxkl*xin(449)
                                      yin(450) = yin(450) + dykl*yin(449)
                                      zin(450) = zin(450) + dzkl*zin(449)

                                      ! i3 = i4 =  449
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  443

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  443

                                      ! do nk = 1,    2

                                      xin(443) = xin(445) + dxkl*xin(442)
                                      yin(443) = yin(445) + dykl*yin(442)
                                      zin(443) = zin(445) + dzkl*zin(442)
                                      ! i4 = i4 + lang+1 =  446

                                      ! nk =    2

                                      xin(446) = xin(448) + dxkl*xin(445)
                                      yin(446) = yin(448) + dykl*yin(445)
                                      zin(446) = zin(448) + dzkl*zin(445)
                                      ! i4 = i4 + lang+1 =  449

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  444

                                      ! nl =    2

                                      ! i4 = i3 =  444

                                      ! do nk = 1,    2

                                      xin(444) = xin(446) + dxkl*xin(443)
                                      yin(444) = yin(446) + dykl*yin(443)
                                      zin(444) = zin(446) + dzkl*zin(443)
                                      ! i4 = i4 + lang+1 =  447

                                      ! nk =    2

                                      xin(447) = xin(449) + dxkl*xin(446)
                                      yin(447) = yin(449) + dykl*yin(446)
                                      zin(447) = zin(449) + dzkl*zin(446)
                                      ! i4 = i4 + lang+1 =  450

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  445

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  451

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  459

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  458

                                      xin(459) = xin(459) + dxkl*xin(458)
                                      yin(459) = yin(459) + dykl*yin(458)
                                      zin(459) = zin(459) + dzkl*zin(458)

                                      ! i3 = i4 =  458
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  457

                                      xin(458) = xin(458) + dxkl*xin(457)
                                      yin(458) = yin(458) + dykl*yin(457)
                                      zin(458) = zin(458) + dzkl*zin(457)

                                      ! i3 = i4 =  457
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  459

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  458

                                      xin(459) = xin(459) + dxkl*xin(458)
                                      yin(459) = yin(459) + dykl*yin(458)
                                      zin(459) = zin(459) + dzkl*zin(458)

                                      ! i3 = i4 =  458
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  452

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  452

                                      ! do nk = 1,    2

                                      xin(452) = xin(454) + dxkl*xin(451)
                                      yin(452) = yin(454) + dykl*yin(451)
                                      zin(452) = zin(454) + dzkl*zin(451)
                                      ! i4 = i4 + lang+1 =  455

                                      ! nk =    2

                                      xin(455) = xin(457) + dxkl*xin(454)
                                      yin(455) = yin(457) + dykl*yin(454)
                                      zin(455) = zin(457) + dzkl*zin(454)
                                      ! i4 = i4 + lang+1 =  458

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  453

                                      ! nl =    2

                                      ! i4 = i3 =  453

                                      ! do nk = 1,    2

                                      xin(453) = xin(455) + dxkl*xin(452)
                                      yin(453) = yin(455) + dykl*yin(452)
                                      zin(453) = zin(455) + dzkl*zin(452)
                                      ! i4 = i4 + lang+1 =  456

                                      ! nk =    2

                                      xin(456) = xin(458) + dxkl*xin(455)
                                      yin(456) = yin(458) + dykl*yin(455)
                                      zin(456) = zin(458) + dzkl*zin(455)
                                      ! i4 = i4 + lang+1 =  459

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  454

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  460

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    1

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  460

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  468

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  467

                                      xin(468) = xin(468) + dxkl*xin(467)
                                      yin(468) = yin(468) + dykl*yin(467)
                                      zin(468) = zin(468) + dzkl*zin(467)

                                      ! i3 = i4 =  467
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  466

                                      xin(467) = xin(467) + dxkl*xin(466)
                                      yin(467) = yin(467) + dykl*yin(466)
                                      zin(467) = zin(467) + dzkl*zin(466)

                                      ! i3 = i4 =  466
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  468

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  467

                                      xin(468) = xin(468) + dxkl*xin(467)
                                      yin(468) = yin(468) + dykl*yin(467)
                                      zin(468) = zin(468) + dzkl*zin(467)

                                      ! i3 = i4 =  467
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  461

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  461

                                      ! do nk = 1,    2

                                      xin(461) = xin(463) + dxkl*xin(460)
                                      yin(461) = yin(463) + dykl*yin(460)
                                      zin(461) = zin(463) + dzkl*zin(460)
                                      ! i4 = i4 + lang+1 =  464

                                      ! nk =    2

                                      xin(464) = xin(466) + dxkl*xin(463)
                                      yin(464) = yin(466) + dykl*yin(463)
                                      zin(464) = zin(466) + dzkl*zin(463)
                                      ! i4 = i4 + lang+1 =  467

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  462

                                      ! nl =    2

                                      ! i4 = i3 =  462

                                      ! do nk = 1,    2

                                      xin(462) = xin(464) + dxkl*xin(461)
                                      yin(462) = yin(464) + dykl*yin(461)
                                      zin(462) = zin(464) + dzkl*zin(461)
                                      ! i4 = i4 + lang+1 =  465

                                      ! nk =    2

                                      xin(465) = xin(467) + dxkl*xin(464)
                                      yin(465) = yin(467) + dykl*yin(464)
                                      zin(465) = zin(467) + dzkl*zin(464)
                                      ! i4 = i4 + lang+1 =  468

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  463

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  469

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  477

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  476

                                      xin(477) = xin(477) + dxkl*xin(476)
                                      yin(477) = yin(477) + dykl*yin(476)
                                      zin(477) = zin(477) + dzkl*zin(476)

                                      ! i3 = i4 =  476
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  475

                                      xin(476) = xin(476) + dxkl*xin(475)
                                      yin(476) = yin(476) + dykl*yin(475)
                                      zin(476) = zin(476) + dzkl*zin(475)

                                      ! i3 = i4 =  475
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  477

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  476

                                      xin(477) = xin(477) + dxkl*xin(476)
                                      yin(477) = yin(477) + dykl*yin(476)
                                      zin(477) = zin(477) + dzkl*zin(476)

                                      ! i3 = i4 =  476
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  470

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  470

                                      ! do nk = 1,    2

                                      xin(470) = xin(472) + dxkl*xin(469)
                                      yin(470) = yin(472) + dykl*yin(469)
                                      zin(470) = zin(472) + dzkl*zin(469)
                                      ! i4 = i4 + lang+1 =  473

                                      ! nk =    2

                                      xin(473) = xin(475) + dxkl*xin(472)
                                      yin(473) = yin(475) + dykl*yin(472)
                                      zin(473) = zin(475) + dzkl*zin(472)
                                      ! i4 = i4 + lang+1 =  476

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  471

                                      ! nl =    2

                                      ! i4 = i3 =  471

                                      ! do nk = 1,    2

                                      xin(471) = xin(473) + dxkl*xin(470)
                                      yin(471) = yin(473) + dykl*yin(470)
                                      zin(471) = zin(473) + dzkl*zin(470)
                                      ! i4 = i4 + lang+1 =  474

                                      ! nk =    2

                                      xin(474) = xin(476) + dxkl*xin(473)
                                      yin(474) = yin(476) + dykl*yin(473)
                                      zin(474) = zin(476) + dzkl*zin(473)
                                      ! i4 = i4 + lang+1 =  477

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  472

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  478

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  486

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  485

                                      xin(486) = xin(486) + dxkl*xin(485)
                                      yin(486) = yin(486) + dykl*yin(485)
                                      zin(486) = zin(486) + dzkl*zin(485)

                                      ! i3 = i4 =  485
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  484

                                      xin(485) = xin(485) + dxkl*xin(484)
                                      yin(485) = yin(485) + dykl*yin(484)
                                      zin(485) = zin(485) + dzkl*zin(484)

                                      ! i3 = i4 =  484
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  486

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  485

                                      xin(486) = xin(486) + dxkl*xin(485)
                                      yin(486) = yin(486) + dykl*yin(485)
                                      zin(486) = zin(486) + dzkl*zin(485)

                                      ! i3 = i4 =  485
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  479

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  479

                                      ! do nk = 1,    2

                                      xin(479) = xin(481) + dxkl*xin(478)
                                      yin(479) = yin(481) + dykl*yin(478)
                                      zin(479) = zin(481) + dzkl*zin(478)
                                      ! i4 = i4 + lang+1 =  482

                                      ! nk =    2

                                      xin(482) = xin(484) + dxkl*xin(481)
                                      yin(482) = yin(484) + dykl*yin(481)
                                      zin(482) = zin(484) + dzkl*zin(481)
                                      ! i4 = i4 + lang+1 =  485

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  480

                                      ! nl =    2

                                      ! i4 = i3 =  480

                                      ! do nk = 1,    2

                                      xin(480) = xin(482) + dxkl*xin(479)
                                      yin(480) = yin(482) + dykl*yin(479)
                                      zin(480) = zin(482) + dzkl*zin(479)
                                      ! i4 = i4 + lang+1 =  483

                                      ! nk =    2

                                      xin(483) = xin(485) + dxkl*xin(482)
                                      yin(483) = yin(485) + dykl*yin(482)
                                      zin(483) = zin(485) + dzkl*zin(482)
                                      ! i4 = i4 + lang+1 =  486

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  481

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  487

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    2

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  487

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  495

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  494

                                      xin(495) = xin(495) + dxkl*xin(494)
                                      yin(495) = yin(495) + dykl*yin(494)
                                      zin(495) = zin(495) + dzkl*zin(494)

                                      ! i3 = i4 =  494
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  493

                                      xin(494) = xin(494) + dxkl*xin(493)
                                      yin(494) = yin(494) + dykl*yin(493)
                                      zin(494) = zin(494) + dzkl*zin(493)

                                      ! i3 = i4 =  493
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  495

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  494

                                      xin(495) = xin(495) + dxkl*xin(494)
                                      yin(495) = yin(495) + dykl*yin(494)
                                      zin(495) = zin(495) + dzkl*zin(494)

                                      ! i3 = i4 =  494
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  488

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  488

                                      ! do nk = 1,    2

                                      xin(488) = xin(490) + dxkl*xin(487)
                                      yin(488) = yin(490) + dykl*yin(487)
                                      zin(488) = zin(490) + dzkl*zin(487)
                                      ! i4 = i4 + lang+1 =  491

                                      ! nk =    2

                                      xin(491) = xin(493) + dxkl*xin(490)
                                      yin(491) = yin(493) + dykl*yin(490)
                                      zin(491) = zin(493) + dzkl*zin(490)
                                      ! i4 = i4 + lang+1 =  494

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  489

                                      ! nl =    2

                                      ! i4 = i3 =  489

                                      ! do nk = 1,    2

                                      xin(489) = xin(491) + dxkl*xin(488)
                                      yin(489) = yin(491) + dykl*yin(488)
                                      zin(489) = zin(491) + dzkl*zin(488)
                                      ! i4 = i4 + lang+1 =  492

                                      ! nk =    2

                                      xin(492) = xin(494) + dxkl*xin(491)
                                      yin(492) = yin(494) + dykl*yin(491)
                                      zin(492) = zin(494) + dzkl*zin(491)
                                      ! i4 = i4 + lang+1 =  495

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  490

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  496

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  504

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  503

                                      xin(504) = xin(504) + dxkl*xin(503)
                                      yin(504) = yin(504) + dykl*yin(503)
                                      zin(504) = zin(504) + dzkl*zin(503)

                                      ! i3 = i4 =  503
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  502

                                      xin(503) = xin(503) + dxkl*xin(502)
                                      yin(503) = yin(503) + dykl*yin(502)
                                      zin(503) = zin(503) + dzkl*zin(502)

                                      ! i3 = i4 =  502
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  504

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  503

                                      xin(504) = xin(504) + dxkl*xin(503)
                                      yin(504) = yin(504) + dykl*yin(503)
                                      zin(504) = zin(504) + dzkl*zin(503)

                                      ! i3 = i4 =  503
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  497

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  497

                                      ! do nk = 1,    2

                                      xin(497) = xin(499) + dxkl*xin(496)
                                      yin(497) = yin(499) + dykl*yin(496)
                                      zin(497) = zin(499) + dzkl*zin(496)
                                      ! i4 = i4 + lang+1 =  500

                                      ! nk =    2

                                      xin(500) = xin(502) + dxkl*xin(499)
                                      yin(500) = yin(502) + dykl*yin(499)
                                      zin(500) = zin(502) + dzkl*zin(499)
                                      ! i4 = i4 + lang+1 =  503

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  498

                                      ! nl =    2

                                      ! i4 = i3 =  498

                                      ! do nk = 1,    2

                                      xin(498) = xin(500) + dxkl*xin(497)
                                      yin(498) = yin(500) + dykl*yin(497)
                                      zin(498) = zin(500) + dzkl*zin(497)
                                      ! i4 = i4 + lang+1 =  501

                                      ! nk =    2

                                      xin(501) = xin(503) + dxkl*xin(500)
                                      yin(501) = yin(503) + dykl*yin(500)
                                      zin(501) = zin(503) + dzkl*zin(500)
                                      ! i4 = i4 + lang+1 =  504

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  499

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  505

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  513

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  512

                                      xin(513) = xin(513) + dxkl*xin(512)
                                      yin(513) = yin(513) + dykl*yin(512)
                                      zin(513) = zin(513) + dzkl*zin(512)

                                      ! i3 = i4 =  512
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  511

                                      xin(512) = xin(512) + dxkl*xin(511)
                                      yin(512) = yin(512) + dykl*yin(511)
                                      zin(512) = zin(512) + dzkl*zin(511)

                                      ! i3 = i4 =  511
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  513

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  512

                                      xin(513) = xin(513) + dxkl*xin(512)
                                      yin(513) = yin(513) + dykl*yin(512)
                                      zin(513) = zin(513) + dzkl*zin(512)

                                      ! i3 = i4 =  512
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  506

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  506

                                      ! do nk = 1,    2

                                      xin(506) = xin(508) + dxkl*xin(505)
                                      yin(506) = yin(508) + dykl*yin(505)
                                      zin(506) = zin(508) + dzkl*zin(505)
                                      ! i4 = i4 + lang+1 =  509

                                      ! nk =    2

                                      xin(509) = xin(511) + dxkl*xin(508)
                                      yin(509) = yin(511) + dykl*yin(508)
                                      zin(509) = zin(511) + dzkl*zin(508)
                                      ! i4 = i4 + lang+1 =  512

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  507

                                      ! nl =    2

                                      ! i4 = i3 =  507

                                      ! do nk = 1,    2

                                      xin(507) = xin(509) + dxkl*xin(506)
                                      yin(507) = yin(509) + dykl*yin(506)
                                      zin(507) = zin(509) + dzkl*zin(506)
                                      ! i4 = i4 + lang+1 =  510

                                      ! nk =    2

                                      xin(510) = xin(512) + dxkl*xin(509)
                                      yin(510) = yin(512) + dykl*yin(509)
                                      zin(510) = zin(512) + dzkl*zin(509)
                                      ! i4 = i4 + lang+1 =  513

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  508

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  514

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    3

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  514

                                      ! nj = 0

                                      ! ib = iaa

                                      ! do while nj.le.jang

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  522

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  521

                                      xin(522) = xin(522) + dxkl*xin(521)
                                      yin(522) = yin(522) + dykl*yin(521)
                                      zin(522) = zin(522) + dzkl*zin(521)

                                      ! i3 = i4 =  521
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  520

                                      xin(521) = xin(521) + dxkl*xin(520)
                                      yin(521) = yin(521) + dykl*yin(520)
                                      zin(521) = zin(521) + dzkl*zin(520)

                                      ! i3 = i4 =  520
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  522

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  521

                                      xin(522) = xin(522) + dxkl*xin(521)
                                      yin(522) = yin(522) + dykl*yin(521)
                                      zin(522) = zin(522) + dzkl*zin(521)

                                      ! i3 = i4 =  521
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  515

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  515

                                      ! do nk = 1,    2

                                      xin(515) = xin(517) + dxkl*xin(514)
                                      yin(515) = yin(517) + dykl*yin(514)
                                      zin(515) = zin(517) + dzkl*zin(514)
                                      ! i4 = i4 + lang+1 =  518

                                      ! nk =    2

                                      xin(518) = xin(520) + dxkl*xin(517)
                                      yin(518) = yin(520) + dykl*yin(517)
                                      zin(518) = zin(520) + dzkl*zin(517)
                                      ! i4 = i4 + lang+1 =  521

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  516

                                      ! nl =    2

                                      ! i4 = i3 =  516

                                      ! do nk = 1,    2

                                      xin(516) = xin(518) + dxkl*xin(515)
                                      yin(516) = yin(518) + dykl*yin(515)
                                      zin(516) = zin(518) + dzkl*zin(515)
                                      ! i4 = i4 + lang+1 =  519

                                      ! nk =    2

                                      xin(519) = xin(521) + dxkl*xin(518)
                                      yin(519) = yin(521) + dykl*yin(518)
                                      zin(519) = zin(521) + dzkl*zin(518)
                                      ! i4 = i4 + lang+1 =  522

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  517

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  523

                                      ! nj = nj + 1 =    1

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  531

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  530

                                      xin(531) = xin(531) + dxkl*xin(530)
                                      yin(531) = yin(531) + dykl*yin(530)
                                      zin(531) = zin(531) + dzkl*zin(530)

                                      ! i3 = i4 =  530
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  529

                                      xin(530) = xin(530) + dxkl*xin(529)
                                      yin(530) = yin(530) + dykl*yin(529)
                                      zin(530) = zin(530) + dzkl*zin(529)

                                      ! i3 = i4 =  529
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  531

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  530

                                      xin(531) = xin(531) + dxkl*xin(530)
                                      yin(531) = yin(531) + dykl*yin(530)
                                      zin(531) = zin(531) + dzkl*zin(530)

                                      ! i3 = i4 =  530
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  524

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  524

                                      ! do nk = 1,    2

                                      xin(524) = xin(526) + dxkl*xin(523)
                                      yin(524) = yin(526) + dykl*yin(523)
                                      zin(524) = zin(526) + dzkl*zin(523)
                                      ! i4 = i4 + lang+1 =  527

                                      ! nk =    2

                                      xin(527) = xin(529) + dxkl*xin(526)
                                      yin(527) = yin(529) + dykl*yin(526)
                                      zin(527) = zin(529) + dzkl*zin(526)
                                      ! i4 = i4 + lang+1 =  530

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  525

                                      ! nl =    2

                                      ! i4 = i3 =  525

                                      ! do nk = 1,    2

                                      xin(525) = xin(527) + dxkl*xin(524)
                                      yin(525) = yin(527) + dykl*yin(524)
                                      zin(525) = zin(527) + dzkl*zin(524)
                                      ! i4 = i4 + lang+1 =  528

                                      ! nk =    2

                                      xin(528) = xin(530) + dxkl*xin(527)
                                      yin(528) = yin(530) + dykl*yin(527)
                                      zin(528) = zin(530) + dzkl*zin(527)
                                      ! i4 = i4 + lang+1 =  531

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  526

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  532

                                      ! nj = nj + 1 =    2

                                      ! min = kang

                                      ! do while min.lt.(kang+lang)

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  540

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  539

                                      xin(540) = xin(540) + dxkl*xin(539)
                                      yin(540) = yin(540) + dykl*yin(539)
                                      zin(540) = zin(540) + dzkl*zin(539)

                                      ! i3 = i4 =  539
                                      ! nm = nm -1 =    3

                                      ! i4 = ib+kn(nm) =  538

                                      xin(539) = xin(539) + dxkl*xin(538)
                                      yin(539) = yin(539) + dykl*yin(538)
                                      zin(539) = zin(539) + dzkl*zin(538)

                                      ! i3 = i4 =  538
                                      ! nm = nm -1 =    2

                                      ! end do

                                      ! min = min + 1 =    3

                                      ! nm = (kang+lang) =    4

                                      ! i3 = ib+i5 =  540

                                      ! do while nm.gt.min

                                      ! i4 = ib+kn(nm) =  539

                                      xin(540) = xin(540) + dxkl*xin(539)
                                      yin(540) = yin(540) + dykl*yin(539)
                                      zin(540) = zin(540) + dzkl*zin(539)

                                      ! i3 = i4 =  539
                                      ! nm = nm -1 =    3

                                      ! end do

                                      ! min = min + 1 =    4

                                      ! end do

                                      ! i3 = ib + 1 =  533

                                      ! do nl = 1,    2

                                      ! i4 = i3 =  533

                                      ! do nk = 1,    2

                                      xin(533) = xin(535) + dxkl*xin(532)
                                      yin(533) = yin(535) + dykl*yin(532)
                                      zin(533) = zin(535) + dzkl*zin(532)
                                      ! i4 = i4 + lang+1 =  536

                                      ! nk =    2

                                      xin(536) = xin(538) + dxkl*xin(535)
                                      yin(536) = yin(538) + dykl*yin(535)
                                      zin(536) = zin(538) + dzkl*zin(535)
                                      ! i4 = i4 + lang+1 =  539

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  534

                                      ! nl =    2

                                      ! i4 = i3 =  534

                                      ! do nk = 1,    2

                                      xin(534) = xin(536) + dxkl*xin(533)
                                      yin(534) = yin(536) + dykl*yin(533)
                                      zin(534) = zin(536) + dzkl*zin(533)
                                      ! i4 = i4 + lang+1 =  537

                                      ! nk =    2

                                      xin(537) = xin(539) + dxkl*xin(536)
                                      yin(537) = yin(539) + dykl*yin(536)
                                      zin(537) = zin(539) + dzkl*zin(536)
                                      ! i4 = i4 + lang+1 =  540

                                      ! nk =    3

                                      ! end do

                                      ! i3 = i3 + 1 =  535

                                      ! nl =    3

                                      ! end do

                                      ! ib = ib + (kang+1)*(lang+1) =  541

                                      ! nj = nj + 1 =    3

                                      ! end do

                                      ! ni = ni + 1 =    4

                                      ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  541

                                      ! end do

                                      ! *** Now root =    6

                                      ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  540

                                      !                     --- END XYZINT ---

                                      !                       --- FORMS ---
                                      ! Form final integrals adding the 2D auxiliaries over all roots

                                      j = 1

                                      do n = 1, 2160! loop over all integrals

                                        l = n - 36*(j - 1) ! index for the ket cartesian pair

                                        mx = ijx(j) + klx(l)
                                        my = ijy(j) + kly(l)
                                        mz = ijz(j) + klz(l)

                                        eri_value(n) = eri_value(n) + d23bra(j)*d22ket(l)* &
                                                       (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                        + xin(mx + 108)*yin(my + 108)*zin(mz + 108) & ! root  2
                                                        + xin(mx + 216)*yin(my + 216)*zin(mz + 216) & ! root  3
                                                        + xin(mx + 324)*yin(my + 324)*zin(mz + 324) & ! root  4
                                                        + xin(mx + 432)*yin(my + 432)*zin(mz + 432)) ! root  5

                                        j = int(n/36) + 1 ! index for the next bra cartesian pair

                                      end do

                                      !                     --- END FORMS ---

                                    end do ! ij primitve loop

                                  end do ! kl primitve loop

                                  !                     --- DIRFCK_RHF ---
                                  !          Compute Fock matrix elements from 2EIs

                                  maxl = 6
                                  kandl = ksh .eq. lsh

                                  loci = res%atom_loc(ish) - 1
                                  locj = res%atom_loc(jsh) - 1
                                  lock = res%atom_loc(ksh) - 1
                                  locl = res%atom_loc(lsh) - 1

                                  do i = 1, 10 ! # of cartesians in i

                                    ii1 = i + loci
                                    ip = (i - 1)*216 ! Stride between functions in i

                                    do j = 1, 6 ! # of cartesians in j

                                      maxl2 = maxl

                                      jj1 = j + locj
                                      i2 = ii1
                                      j2 = jj1
                                      if (ii1 .lt. jj1) then ! Sort <ij|
                                        i2 = jj1
                                        j2 = ii1
                                      end if

                                      ijp = (j - 1)*36 + ip ! Add stride between functions in j

                                      do k = 1, 6 ! # of cartesians in k

                                        if (kandl) maxl2 = k

                                        kk1 = k + lock

                                        ijkp = (k - 1)*6 + ijp ! Add stride between functions in k

                                        do l = 1, maxl2

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
                              deallocate (n22ket)
                              deallocate (xint22ket)

                              end subroutine int3222
                              end submodule
