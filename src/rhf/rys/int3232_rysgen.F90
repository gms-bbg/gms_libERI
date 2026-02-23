! The total angular momentum of this class is:          10
! The algorithm chosen is: Rys quadrature
submodule(rys_kernels) int3232_impl
contains
  module subroutine int3232(df_pair, density, fock, res)
    use omp_lib
    use liberi_types, only: shell_pair_t, int32, int64, dp, eri_resources_t
    !use mdi_api
    use parameters
    use boys
    implicit none
    type(shell_pair_t), intent(in) :: df_pair
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)
    type(eri_resources_t), intent(in) :: res

    ! Variables for the class
    integer(kind=int64), allocatable :: n23bra(:)
    real(dp), allocatable :: xint23bra(:)
    integer(kind=int64) :: ndfbra
    real(dp) :: scutdfbra, test
    integer(kind=int64) :: ij, kl, nquarts, nchunk, nquart_start, nquart_end, iijj, kkll
    integer(kind=int64) :: ii, jj, kk, ll, ijkl, ish_tmp, jsh_tmp, ksh_tmp, lsh_tmp
    integer(kind=int64) :: ij_tmp, kl_tmp, ish, jsh, ksh, lsh, i, j, k, l, m, mm, iii
    integer(kind=int64) :: loci, locj, lock, locl, ip, ijp, ijkp, ijklp, jk, jl, il, ik
    integer(kind=int64) :: ii1, i2, ii2, jj1, j2, jj2, kk1, k2, kk2, ll1, l2, maxl, maxl2, nij, nkl, itmp
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
    real(dp) :: roots(6), wghts(6)
    real(dp) :: factr, factw, sum0, sum1, sum2, t, dpp, dg, dr, ds, dc, df, db
    integer(kind=int64) :: mml, mmii, iim1
    real(dp) :: rgrid(40), wgrid(40), p0(40), p1(40), p2(40)
    real(dp) :: rts(6), wts(6), alpha(6), beta(6), wrk(6)
    real(dp) :: xin(864), yin(864), zin(864)
    real(dp) :: eri_value(3600)
    real(dp) :: d23bra(60), d23ket(60)
    integer(kind=int64) :: ix(10), jx(6), kx(10), lx(6)
    integer(kind=int64) :: iy(10), jy(6), ky(10), ly(6)
    integer(kind=int64) :: iz(10), jz(6), kz(10), lz(6)
    integer(kind=int64) :: in(6), in1(6), kn(6)
    integer(kind=int64) :: ijx(60), ijy(60), ijz(60)
    integer(kind=int64) :: klx(60), kly(60), klz(60)
    logical :: same

    ! Fill arrays containing the 2D auxuliary integral indexing
    ! These are probably not necessary anymore since we pre-fill the
    ! ij- and kl-xyz arrays. Leave them for now for clarity.
    in1(1) = 1
    in1(2) = 37
    in1(3) = 73
    in1(4) = 109
    in1(5) = 121
    in1(6) = 133

    kn(1) = 0
    kn(2) = 3
    kn(3) = 6
    kn(4) = 9
    kn(5) = 10
    kn(6) = 11

    ! Fill arrays for accessing of 2D auxiliary integrals

    ! x-arrays

    lx(1) = 2
    lx(2) = 0
    lx(3) = 0
    lx(4) = 1
    lx(5) = 1
    lx(6) = 0

    kx(1) = 9
    kx(2) = 0
    kx(3) = 0
    kx(4) = 6
    kx(5) = 6
    kx(6) = 3
    kx(7) = 0
    kx(8) = 3
    kx(9) = 0
    kx(10) = 3

    jx(1) = 24
    jx(2) = 0
    jx(3) = 0
    jx(4) = 12
    jx(5) = 12
    jx(6) = 0

    ix(1) = 109
    ix(2) = 1
    ix(3) = 1
    ix(4) = 73
    ix(5) = 73
    ix(6) = 37
    ix(7) = 1
    ix(8) = 37
    ix(9) = 1
    ix(10) = 37

    ! y-arrays

    ly(1) = 0
    ly(2) = 2
    ly(3) = 0
    ly(4) = 1
    ly(5) = 0
    ly(6) = 1

    ky(1) = 0
    ky(2) = 9
    ky(3) = 0
    ky(4) = 3
    ky(5) = 0
    ky(6) = 6
    ky(7) = 6
    ky(8) = 0
    ky(9) = 3
    ky(10) = 3

    jy(1) = 0
    jy(2) = 24
    jy(3) = 0
    jy(4) = 12
    jy(5) = 0
    jy(6) = 12

    iy(1) = 1
    iy(2) = 109
    iy(3) = 1
    iy(4) = 37
    iy(5) = 1
    iy(6) = 73
    iy(7) = 73
    iy(8) = 1
    iy(9) = 37
    iy(10) = 37

    ! z-arrays

    lz(1) = 0
    lz(2) = 0
    lz(3) = 2
    lz(4) = 0
    lz(5) = 1
    lz(6) = 1

    kz(1) = 0
    kz(2) = 0
    kz(3) = 9
    kz(4) = 0
    kz(5) = 3
    kz(6) = 0
    kz(7) = 3
    kz(8) = 6
    kz(9) = 6
    kz(10) = 3

    jz(1) = 0
    jz(2) = 0
    jz(3) = 24
    jz(4) = 0
    jz(5) = 12
    jz(6) = 12

    iz(1) = 1
    iz(2) = 1
    iz(3) = 109
    iz(4) = 1
    iz(5) = 37
    iz(6) = 1
    iz(7) = 37
    iz(8) = 73
    iz(9) = 73
    iz(10) = 37

    ! ij-xyz arrays to form final integrals from 2D auxiliaries

    ijx(1) = 133
    ijx(2) = 109
    ijx(3) = 109
    ijx(4) = 121
    ijx(5) = 121
    ijx(6) = 109
    ijx(7) = 25
    ijx(8) = 1
    ijx(9) = 1
    ijx(10) = 13
    ijx(11) = 13
    ijx(12) = 1
    ijx(13) = 25
    ijx(14) = 1
    ijx(15) = 1
    ijx(16) = 13
    ijx(17) = 13
    ijx(18) = 1
    ijx(19) = 97
    ijx(20) = 73
    ijx(21) = 73
    ijx(22) = 85
    ijx(23) = 85
    ijx(24) = 73
    ijx(25) = 97
    ijx(26) = 73
    ijx(27) = 73
    ijx(28) = 85
    ijx(29) = 85
    ijx(30) = 73
    ijx(31) = 61
    ijx(32) = 37
    ijx(33) = 37
    ijx(34) = 49
    ijx(35) = 49
    ijx(36) = 37
    ijx(37) = 25
    ijx(38) = 1
    ijx(39) = 1
    ijx(40) = 13
    ijx(41) = 13
    ijx(42) = 1
    ijx(43) = 61
    ijx(44) = 37
    ijx(45) = 37
    ijx(46) = 49
    ijx(47) = 49
    ijx(48) = 37
    ijx(49) = 25
    ijx(50) = 1
    ijx(51) = 1
    ijx(52) = 13
    ijx(53) = 13
    ijx(54) = 1
    ijx(55) = 61
    ijx(56) = 37
    ijx(57) = 37
    ijx(58) = 49
    ijx(59) = 49
    ijx(60) = 37

    ijy(1) = 1
    ijy(2) = 25
    ijy(3) = 1
    ijy(4) = 13
    ijy(5) = 1
    ijy(6) = 13
    ijy(7) = 109
    ijy(8) = 133
    ijy(9) = 109
    ijy(10) = 121
    ijy(11) = 109
    ijy(12) = 121
    ijy(13) = 1
    ijy(14) = 25
    ijy(15) = 1
    ijy(16) = 13
    ijy(17) = 1
    ijy(18) = 13
    ijy(19) = 37
    ijy(20) = 61
    ijy(21) = 37
    ijy(22) = 49
    ijy(23) = 37
    ijy(24) = 49
    ijy(25) = 1
    ijy(26) = 25
    ijy(27) = 1
    ijy(28) = 13
    ijy(29) = 1
    ijy(30) = 13
    ijy(31) = 73
    ijy(32) = 97
    ijy(33) = 73
    ijy(34) = 85
    ijy(35) = 73
    ijy(36) = 85
    ijy(37) = 73
    ijy(38) = 97
    ijy(39) = 73
    ijy(40) = 85
    ijy(41) = 73
    ijy(42) = 85
    ijy(43) = 1
    ijy(44) = 25
    ijy(45) = 1
    ijy(46) = 13
    ijy(47) = 1
    ijy(48) = 13
    ijy(49) = 37
    ijy(50) = 61
    ijy(51) = 37
    ijy(52) = 49
    ijy(53) = 37
    ijy(54) = 49
    ijy(55) = 37
    ijy(56) = 61
    ijy(57) = 37
    ijy(58) = 49
    ijy(59) = 37
    ijy(60) = 49

    ijz(1) = 1
    ijz(2) = 1
    ijz(3) = 25
    ijz(4) = 1
    ijz(5) = 13
    ijz(6) = 13
    ijz(7) = 1
    ijz(8) = 1
    ijz(9) = 25
    ijz(10) = 1
    ijz(11) = 13
    ijz(12) = 13
    ijz(13) = 109
    ijz(14) = 109
    ijz(15) = 133
    ijz(16) = 109
    ijz(17) = 121
    ijz(18) = 121
    ijz(19) = 1
    ijz(20) = 1
    ijz(21) = 25
    ijz(22) = 1
    ijz(23) = 13
    ijz(24) = 13
    ijz(25) = 37
    ijz(26) = 37
    ijz(27) = 61
    ijz(28) = 37
    ijz(29) = 49
    ijz(30) = 49
    ijz(31) = 1
    ijz(32) = 1
    ijz(33) = 25
    ijz(34) = 1
    ijz(35) = 13
    ijz(36) = 13
    ijz(37) = 37
    ijz(38) = 37
    ijz(39) = 61
    ijz(40) = 37
    ijz(41) = 49
    ijz(42) = 49
    ijz(43) = 73
    ijz(44) = 73
    ijz(45) = 97
    ijz(46) = 73
    ijz(47) = 85
    ijz(48) = 85
    ijz(49) = 73
    ijz(50) = 73
    ijz(51) = 97
    ijz(52) = 73
    ijz(53) = 85
    ijz(54) = 85
    ijz(55) = 37
    ijz(56) = 37
    ijz(57) = 61
    ijz(58) = 37
    ijz(59) = 49
    ijz(60) = 49

    ! kl-xyz arrays to form final integrals from 2D auxiliaries

    klx(1) = 11
    klx(2) = 9
    klx(3) = 9
    klx(4) = 10
    klx(5) = 10
    klx(6) = 9
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
    klx(19) = 8
    klx(20) = 6
    klx(21) = 6
    klx(22) = 7
    klx(23) = 7
    klx(24) = 6
    klx(25) = 8
    klx(26) = 6
    klx(27) = 6
    klx(28) = 7
    klx(29) = 7
    klx(30) = 6
    klx(31) = 5
    klx(32) = 3
    klx(33) = 3
    klx(34) = 4
    klx(35) = 4
    klx(36) = 3
    klx(37) = 2
    klx(38) = 0
    klx(39) = 0
    klx(40) = 1
    klx(41) = 1
    klx(42) = 0
    klx(43) = 5
    klx(44) = 3
    klx(45) = 3
    klx(46) = 4
    klx(47) = 4
    klx(48) = 3
    klx(49) = 2
    klx(50) = 0
    klx(51) = 0
    klx(52) = 1
    klx(53) = 1
    klx(54) = 0
    klx(55) = 5
    klx(56) = 3
    klx(57) = 3
    klx(58) = 4
    klx(59) = 4
    klx(60) = 3

    kly(1) = 0
    kly(2) = 2
    kly(3) = 0
    kly(4) = 1
    kly(5) = 0
    kly(6) = 1
    kly(7) = 9
    kly(8) = 11
    kly(9) = 9
    kly(10) = 10
    kly(11) = 9
    kly(12) = 10
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
    kly(31) = 6
    kly(32) = 8
    kly(33) = 6
    kly(34) = 7
    kly(35) = 6
    kly(36) = 7
    kly(37) = 6
    kly(38) = 8
    kly(39) = 6
    kly(40) = 7
    kly(41) = 6
    kly(42) = 7
    kly(43) = 0
    kly(44) = 2
    kly(45) = 0
    kly(46) = 1
    kly(47) = 0
    kly(48) = 1
    kly(49) = 3
    kly(50) = 5
    kly(51) = 3
    kly(52) = 4
    kly(53) = 3
    kly(54) = 4
    kly(55) = 3
    kly(56) = 5
    kly(57) = 3
    kly(58) = 4
    kly(59) = 3
    kly(60) = 4

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
    klz(13) = 9
    klz(14) = 9
    klz(15) = 11
    klz(16) = 9
    klz(17) = 10
    klz(18) = 10
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
    klz(31) = 0
    klz(32) = 0
    klz(33) = 2
    klz(34) = 0
    klz(35) = 1
    klz(36) = 1
    klz(37) = 3
    klz(38) = 3
    klz(39) = 5
    klz(40) = 3
    klz(41) = 4
    klz(42) = 4
    klz(43) = 6
    klz(44) = 6
    klz(45) = 8
    klz(46) = 6
    klz(47) = 7
    klz(48) = 7
    klz(49) = 6
    klz(50) = 6
    klz(51) = 8
    klz(52) = 6
    klz(53) = 7
    klz(54) = 7
    klz(55) = 3
    klz(56) = 3
    klz(57) = 5
    klz(58) = 3
    klz(59) = 4
    klz(60) = 4

    allocate (n23bra(res%n_d_shl*res%n_f_shl))
    allocate (xint23bra(res%n_d_shl*res%n_f_shl))

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

    ! --multi-gpu--work
    nchunk = ndfbra/res%n_size
    nquart_start = nchunk*res%n_rank + 1
    nquart_end = nquart_start + nchunk - 1
    if (res%n_rank .EQ. res%n_size - 1) nquart_end = ndfbra

    ! Mappings to GPU

 !$omp target teams distribute parallel do collapse(2) default(none) &
 !$omp shared(res, density, fock, nquart_start, nquart_end, ndfbra, xint23bra, n23bra, df_pair) &
 !$omp shared(ijx, ijy, ijz, klx, kly, klz, in1, kn) &
 !$omp private(test,ij,kl,ish_tmp,jsh_tmp,ksh_tmp,lsh_tmp) &
 !$omp private(ish,jsh,ksh,lsh,ijtop,kltop,rab,rcd,dxij,dyij,dzij,dxkl,dykl,dzkl) &
 !$omp private(eri_value,ket_loop,k,t_expon_cd,t_expon_c,t_expon_d,t_inverse_expon_cd) &
 !$omp private(brrk,akxk,akyk,akzk,bbrrk,xb,yb,zb,bxbk,bybk,bzbk,bxbi,bybi,bzbi) &
 !$omp private(d23ket,bra_loop,i,t_expon_ab,t_expon_a,t_expon_b,t_inverse_expon_ab) &
 !$omp private(dum,t_expon_abcd_inverse,expe,rho,xa,ya,za,axak,ayak,azak,axai) &
 !$omp private(ayai,azai,c1x,c2x,c3x,c4x,c1y,c2y,c3y,c4y,c1z,c2z,c3z,c4z,d23bra,xx) &
 !$omp private(factr,factw,rts,wts,rgrid,wgrid,sum0,sum1,sum2,m,alpha,beta,p0,p1) &
 !$omp private(p2,kk,t,wrk,l,jj,dpp,dg,dr,ds,dc,mml,ii,mmii,df,db,iim1,roots,wghts) &
 !$omp private(mm,n,u2,f00,iii,duminv,dm2inv,bp01,b00,b10,xcp00,xc00,ycp00,yc00) &
 !$omp private(in,zcp00,zc00,i1,i3,cp10,i4,i5,c10,cp01,c01,k3,k4,nn,nm,km,iaa,ib,nj,ni,nl,nk) &
 !$omp private(j,dij,nx,ny,nz,mx,my,mz,loci,locj,lock,locl,ii1,ip,jj1,i2,j2,ijp) &
 !$omp private(xin,yin,zin,kk1,ijkp,ijklp,buff,ll1,k2,l2,ll,ii2,jj2,kk2,ik,il) &
 !$omp private(jk,jl,maxl,maxl2,nij,nkl,itmp,same)
              do iijj = nquart_start, nquart_end
                do kkll = 1, ndfbra

                  if (kkll .gt. iijj) cycle

                  test = xint23bra(iijj)*xint23bra(kkll)
                  if (test .gt. cutoff_schwarz) then

                    ij = n23bra(iijj)
                    kl = n23bra(kkll)

                    ish_tmp = mod(ij - 1, res%n_f_shl) + 1
                    jsh_tmp = (ij - 1)/res%n_f_shl + 1
                    ksh_tmp = mod(kl - 1, res%n_f_shl) + 1
                    lsh_tmp = (kl - 1)/res%n_f_shl + 1

                    ish = res%i_f_shl(ish_tmp)
                    jsh = res%i_d_shl(jsh_tmp)
                    ksh = res%i_f_shl(ksh_tmp)
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

                      t_expon_cd = df_pair%t_expon_ab(df_pair%pair_loc(kl) + ket_loop)
                      t_expon_c = df_pair%expon_b(df_pair%pair_loc(kl) + ket_loop)
                      t_expon_d = df_pair%expon_a(df_pair%pair_loc(kl) + ket_loop)
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

                      d23ket(1) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(2) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(3) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(4) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(5) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(6) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(7) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(8) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(9) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(10) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(11) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(12) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(13) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(14) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(15) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2
                      d23ket(16) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(17) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(18) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt3
                      d23ket(19) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(20) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(21) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(22) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(23) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(24) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(25) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(26) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(27) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(28) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(29) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(30) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(31) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(32) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(33) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(34) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(35) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(36) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(37) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(38) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(39) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(40) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(41) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(42) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(43) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(44) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(45) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(46) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(47) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(48) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(49) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(50) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(51) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5
                      d23ket(52) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(53) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(54) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(55) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(56) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(57) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3
                      d23ket(58) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3*sqrt3
                      d23ket(59) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3*sqrt3
                      d23ket(60) = df_pair%d_coeff_alt(df_pair%pair_loc(kl) + ket_loop)*twopi_5_2*sqrt5*sqrt3*sqrt3

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

                        if (xx .ge. 60.0D+00) then ! Asymptotic form

                          factr = 1.0_dp/xx
                          factw = sqrt(factr)

                          rts(1) = factr*0.9874701406848084D-01
                          rts(2) = factr*0.8983028345696176D+00
                          rts(3) = factr*0.2552589802668170D+01
                          rts(4) = factr*0.5196152530054469D+01
                          rts(5) = factr*0.9124248037531169D+01
                          rts(6) = factr*0.1512995978110809D+02

                          wts(1) = factw*0.5701352362624800D+00
                          wts(2) = factw*0.2604923102641603D+00
                          wts(3) = factw*0.5160798561588373D-01
                          wts(4) = factw*0.3905390584629056D-02
                          wts(5) = factw*0.8573687043587825D-04
                          wts(6) = factw*0.2658551684356290D-06

                        else ! "regular" evaluation

                          rgrid(1) = 0.7764167660645311D-06
                          rgrid(2) = 0.2150066216486457D-04
                          rgrid(3) = 0.1292774686851154D-03
                          rgrid(4) = 0.4427485262711727D-03
                          rgrid(5) = 0.1128529682849472D-02
                          rgrid(6) = 0.2396160899229440D-02
                          rgrid(7) = 0.4491713694776070D-02
                          rgrid(8) = 0.7690217393316887D-02
                          rgrid(9) = 0.1228709604735500D-01
                          rgrid(10) = 0.1858883348816642D-01
                          rgrid(11) = 0.2690310419233024D-01
                          rgrid(12) = 0.3752862210285295D-01
                          rgrid(13) = 0.5074496784251226D-01
                          rgrid(14) = 0.6680265670059675D-01
                          rgrid(15) = 0.8591370530679703D-01
                          rgrid(16) = 0.1082429441270552D+00
                          rgrid(17) = 0.1339003060774142D+00
                          rgrid(18) = 0.1629342990513548D+00
                          rgrid(19) = 0.1953268425285067D+00
                          rgrid(20) = 0.2309896163367905D+00
                          rgrid(21) = 0.2697620338428412D+00
                          rgrid(22) = 0.3114109132037622D+00
                          rgrid(23) = 0.3556318797527262D+00
                          rgrid(24) = 0.4020524910846677D+00
                          rgrid(25) = 0.4502370349528138D+00
                          rgrid(26) = 0.4996929096784019D+00
                          rgrid(27) = 0.5498784583867753D+00
                          rgrid(28) = 0.6002120929376404D+00
                          rgrid(29) = 0.6500825117708330D+00
                          rgrid(30) = 0.6988597888065096D+00
                          rgrid(31) = 0.7459070886780932D+00
                          rgrid(32) = 0.7905927474738742D+00
                          rgrid(33) = 0.8323024482266286D+00
                          rgrid(34) = 0.8704512169070354D+00
                          rgrid(35) = 0.9044949678681030D+00
                          rgrid(36) = 0.9339413379615261D+00
                          rgrid(37) = 0.9583595677400624D+00
                          rgrid(38) = 0.9773892274524588D+00
                          rgrid(39) = 0.9907477393616217D+00
                          rgrid(40) = 0.9982384861273251D+00

                          wgrid(1) = 0.2260638549266816D-02*exp(-xx*0.7764167660645311D-06)
                          wgrid(2) = 0.5249142265576177D-02*exp(-xx*0.2150066216486457D-04)
                          wgrid(3) = 0.8210529190954175D-02*exp(-xx*0.1292774686851154D-03)
                          wgrid(4) = 0.1112292459708341D-01*exp(-xx*0.4427485262711727D-03)
                          wgrid(5) = 0.1396850349001161D-01*exp(-xx*0.1128529682849472D-02)
                          wgrid(6) = 0.1673009764127369D-01*exp(-xx*0.2396160899229440D-02)
                          wgrid(7) = 0.1939108398723607D-01*exp(-xx*0.4491713694776070D-02)
                          wgrid(8) = 0.2193545409283397D-01*exp(-xx*0.7690217393316887D-02)
                          wgrid(9) = 0.2434790381753858D-01*exp(-xx*0.1228709604735500D-01)
                          wgrid(10) = 0.2661392349196827D-01*exp(-xx*0.1858883348816642D-01)
                          wgrid(11) = 0.2871988454969593D-01*exp(-xx*0.2690310419233024D-01)
                          wgrid(12) = 0.3065312124646436D-01*exp(-xx*0.3752862210285295D-01)
                          wgrid(13) = 0.3240200672830060D-01*exp(-xx*0.5074496784251226D-01)
                          wgrid(14) = 0.3395602290761694D-01*exp(-xx*0.6680265670059675D-01)
                          wgrid(15) = 0.3530582369564367D-01*exp(-xx*0.8591370530679703D-01)
                          wgrid(16) = 0.3644329119790184D-01*exp(-xx*0.1082429441270552D+00)
                          wgrid(17) = 0.3736158452898425D-01*exp(-xx*0.1339003060774142D+00)
                          wgrid(18) = 0.3805518095031275D-01*exp(-xx*0.1629342990513548D+00)
                          wgrid(19) = 0.3851990908212372D-01*exp(-xx*0.1953268425285067D+00)
                          wgrid(20) = 0.3875297398921261D-01*exp(-xx*0.2309896163367905D+00)
                          wgrid(21) = 0.3875297398921173D-01*exp(-xx*0.2697620338428412D+00)
                          wgrid(22) = 0.3851990908212474D-01*exp(-xx*0.3114109132037622D+00)
                          wgrid(23) = 0.3805518095031294D-01*exp(-xx*0.3556318797527262D+00)
                          wgrid(24) = 0.3736158452898516D-01*exp(-xx*0.4020524910846677D+00)
                          wgrid(25) = 0.3644329119790151D-01*exp(-xx*0.4502370349528138D+00)
                          wgrid(26) = 0.3530582369564279D-01*exp(-xx*0.4996929096784019D+00)
                          wgrid(27) = 0.3395602290761690D-01*exp(-xx*0.5498784583867753D+00)
                          wgrid(28) = 0.3240200672830062D-01*exp(-xx*0.6002120929376404D+00)
                          wgrid(29) = 0.3065312124646425D-01*exp(-xx*0.6500825117708330D+00)
                          wgrid(30) = 0.2871988454969628D-01*exp(-xx*0.6988597888065096D+00)
                          wgrid(31) = 0.2661392349196852D-01*exp(-xx*0.7459070886780932D+00)
                          wgrid(32) = 0.2434790381753621D-01*exp(-xx*0.7905927474738742D+00)
                          wgrid(33) = 0.2193545409283680D-01*exp(-xx*0.8323024482266286D+00)
                          wgrid(34) = 0.1939108398723554D-01*exp(-xx*0.8704512169070354D+00)
                          wgrid(35) = 0.1673009764127440D-01*exp(-xx*0.9044949678681030D+00)
                          wgrid(36) = 0.1396850349001189D-01*exp(-xx*0.9339413379615261D+00)
                          wgrid(37) = 0.1112292459708395D-01*exp(-xx*0.9583595677400624D+00)
                          wgrid(38) = 0.8210529190953284D-02*exp(-xx*0.9773892274524588D+00)
                          wgrid(39) = 0.5249142265576468D-02*exp(-xx*0.9907477393616217D+00)
                          wgrid(40) = 0.2260638549266104D-02*exp(-xx*0.9982384861273251D+00)

                          ! Call to RYSDS

                          sum0 = 0.0D+00
                          sum1 = 0.0D+00

                          do m = 1, 40
                            sum0 = sum0 + wgrid(m)
                            sum1 = sum1 + wgrid(m)*rgrid(m)
                          end do

                          alpha(1) = sum1/sum0
                          beta(1) = sum0

                          do m = 1, 40
                            p1(m) = 0.0D+00
                            p2(m) = 1.0D+00
                          end do

                          do kk = 1, 5

                            sum1 = 0.0D+00
                            sum2 = 0.0D+00

                            do 30 m = 1, 40

                              if (wgrid(m) .eq. 0.0d+00) goto 30

                              p0(m) = p1(m)
                              p1(m) = p2(m)
                              p2(m) = (rgrid(m) - alpha(kk))*p1(m) - beta(kk)*p0(m)
                              t = wgrid(m)*p2(m)*p2(m)
                              sum1 = sum1 + t
                              sum2 = sum2 + t*rgrid(m)

30                            continue

                              alpha(kk + 1) = sum2/sum1
                              beta(kk + 1) = sum1/sum0
                              sum0 = sum1

                            end do

                            ! End of RYSDS

                            ! Call to RYSGW

                            rts(1) = alpha(1)
                            wts(1) = 1.0D+00
                            wrk(6) = 0.0D+00
                            do 100 kk = 2, 6

                              rts(kk) = alpha(kk)
                              wrk(kk - 1) = sqrt(beta(kk))
                              wts(kk) = 0.0D+00

100                           continue

                              do 240 l = 1, 6

                                jj = 0

105                             do 110 m = l, 6
                                  if (m .eq. 6) go to 120
                                  if (abs(wrk(m)) .le. (1.0D-14)*(abs(rts(m)) + abs(rts(m + 1)))) go to 120
110                               continue

120                               dpp = rts(l)
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
150                                 ds = df/dg
                                    dr = sqrt(ds*ds + 1.0D+00)
                                    wrk(mmii + 1) = dg*dr
                                    dc = 1.0D+00/dr
                                    ds = ds*dc
160                                 dg = rts(mmii + 1) - dpp
                                    dr = (rts(mmii) - dg)*ds + 2.0D+00*dc*db
                                    dpp = ds*dr
                                    rts(mmii + 1) = dg + dpp
                                    dg = dc*dr - db
                                    df = wts(mmii + 1)
                                    wts(mmii + 1) = ds*wts(mmii) + dc*df
                                    wts(mmii) = dc*wts(mmii) - ds*df

200                                 continue

                                    rts(l) = rts(l) - dpp
                                    wrk(l) = dg
                                    wrk(m) = 0.0D+00
                                    go to 105

240                                 continue

                                    do 300 ii = 2, 6

                                      iim1 = ii - 1
                                      kk = iim1
                                      dpp = rts(iim1)

                                      do 260 jj = ii, 6
                                        if (rts(jj) .ge. dpp) go to 260
                                        kk = jj
                                        dpp = rts(jj)
260                                     continue

                                        if (kk .eq. iim1) go to 300

                                        rts(kk) = rts(iim1)
                                        rts(iim1) = dpp
                                        dpp = wts(iim1)
                                        wts(iim1) = wts(kk)
                                        wts(kk) = dpp

300                                     continue

                                        do 310 kk = 1, 6
                                          wts(kk) = beta(1)*wts(kk)*wts(kk)
310                                       continue

                                          end if

                                          do kk = 1, 6
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

                                          ! i2 = in(2) =   37
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(37) = xc00
                                          yin(37) = yc00
                                          zin(37) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =    4

                                          xin(4) = xcp00
                                          yin(4) = ycp00
                                          zin(4) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =   40
                                          ! i2 =   37

                                          xin(40) = xcp00*xin(37) + cp10
                                          yin(40) = ycp00*yin(37) + cp10
                                          zin(40) = zcp00*zin(37) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   37

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =   73
                                          ! i3 =    1
                                          ! i4 =   37

                                          xin(73) = c10*xin(1) + xc00*xin(37)
                                          yin(73) = c10*yin(1) + yc00*yin(37)
                                          zin(73) = c10*zin(1) + zc00*zin(37)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =   76
                                          ! i5 =   73
                                          ! i4 =   37

                                          xin(76) = xcp00*xin(73) + cp10*xin(37)
                                          yin(76) = ycp00*yin(73) + cp10*yin(37)
                                          zin(76) = zcp00*zin(73) + cp10*zin(37)

                                          ! ------------------

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   73

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  109
                                          ! i3 =   37
                                          ! i4 =   73

                                          xin(109) = c10*xin(37) + xc00*xin(73)
                                          yin(109) = c10*yin(37) + yc00*yin(73)
                                          zin(109) = c10*zin(37) + zc00*zin(73)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  112
                                          ! i5 =  109
                                          ! i4 =   73

                                          xin(112) = xcp00*xin(109) + cp10*xin(73)
                                          yin(112) = ycp00*yin(109) + cp10*yin(73)
                                          zin(112) = zcp00*zin(109) + cp10*zin(73)

                                          ! ------------------

                                          ! i3 = i4 =   73
                                          ! i4 = i5 =  109

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  121
                                          ! i3 =   73
                                          ! i4 =  109

                                          xin(121) = c10*xin(73) + xc00*xin(109)
                                          yin(121) = c10*yin(73) + yc00*yin(109)
                                          zin(121) = c10*zin(73) + zc00*zin(109)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  124
                                          ! i5 =  121
                                          ! i4 =  109

                                          xin(124) = xcp00*xin(121) + cp10*xin(109)
                                          yin(124) = ycp00*yin(121) + cp10*yin(109)
                                          zin(124) = zcp00*zin(121) + cp10*zin(109)

                                          ! ------------------

                                          ! i3 = i4 =  109
                                          ! i4 = i5 =  121

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  133
                                          ! i3 =  109
                                          ! i4 =  121

                                          xin(133) = c10*xin(109) + xc00*xin(121)
                                          yin(133) = c10*yin(109) + yc00*yin(121)
                                          zin(133) = c10*zin(109) + zc00*zin(121)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  136
                                          ! i5 =  133
                                          ! i4 =  121

                                          xin(136) = xcp00*xin(133) + cp10*xin(121)
                                          yin(136) = ycp00*yin(133) + cp10*yin(121)
                                          zin(136) = zcp00*zin(133) + cp10*zin(121)

                                          ! ------------------

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =    1
                                          ! i4 = i1+k2 =    4

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =    7
                                          ! i3 =    1
                                          ! i4 =    4

                                          xin(7) = cp01*xin(1) + xcp00*xin(4)
                                          yin(7) = cp01*yin(1) + ycp00*yin(4)
                                          zin(7) = cp01*zin(1) + zcp00*zin(4)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   43

                                          xin(43) = xc00*xin(7) + c01*xin(4)
                                          yin(43) = yc00*yin(7) + c01*yin(4)
                                          zin(43) = zc00*zin(7) + c01*zin(4)

                                          ! ------------------

                                          ! i3 = i4 =    4
                                          ! i4 = i5 =    7

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   10
                                          ! i3 =    4
                                          ! i4 =    7

                                          xin(10) = cp01*xin(4) + xcp00*xin(7)
                                          yin(10) = cp01*yin(4) + ycp00*yin(7)
                                          zin(10) = cp01*zin(4) + zcp00*zin(7)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   46

                                          xin(46) = xc00*xin(10) + c01*xin(7)
                                          yin(46) = yc00*yin(10) + c01*yin(7)
                                          zin(46) = zc00*zin(10) + c01*zin(7)

                                          ! ------------------

                                          ! i3 = i4 =    7
                                          ! i4 = i5 =   10

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   11
                                          ! i3 =    7
                                          ! i4 =   10

                                          xin(11) = cp01*xin(7) + xcp00*xin(10)
                                          yin(11) = cp01*yin(7) + ycp00*yin(10)
                                          zin(11) = cp01*zin(7) + zcp00*zin(10)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   47

                                          xin(47) = xc00*xin(11) + c01*xin(10)
                                          yin(47) = yc00*yin(11) + c01*yin(10)
                                          zin(47) = zc00*zin(11) + c01*zin(10)

                                          ! ------------------

                                          ! i3 = i4 =   10
                                          ! i4 = i5 =   11

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =   12
                                          ! i3 =   10
                                          ! i4 =   11

                                          xin(12) = cp01*xin(10) + xcp00*xin(11)
                                          yin(12) = cp01*yin(10) + ycp00*yin(11)
                                          zin(12) = cp01*zin(10) + zcp00*zin(11)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =   48

                                          xin(48) = xc00*xin(12) + c01*xin(11)
                                          yin(48) = yc00*yin(12) + c01*yin(11)
                                          zin(48) = zc00*zin(12) + c01*zin(11)

                                          ! ------------------

                                          ! i3 = i4 =   11
                                          ! i4 = i5 =   12

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =   73

                                          xin(79) = c10*xin(7) + xc00*xin(43) + c01*xin(40)
                                          yin(79) = c10*yin(7) + yc00*yin(43) + c01*yin(40)
                                          zin(79) = c10*zin(7) + zc00*zin(43) + c01*zin(40)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   73

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  109

                                          xin(115) = c10*xin(43) + xc00*xin(79) + c01*xin(76)
                                          yin(115) = c10*yin(43) + yc00*yin(79) + c01*yin(76)
                                          zin(115) = c10*zin(43) + zc00*zin(79) + c01*zin(76)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   73
                                          ! i4 = i5 =  109

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  121

                                          xin(127) = c10*xin(79) + xc00*xin(115) + c01*xin(112)
                                          yin(127) = c10*yin(79) + yc00*yin(115) + c01*yin(112)
                                          zin(127) = c10*zin(79) + zc00*zin(115) + c01*zin(112)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  109
                                          ! i4 = i5 =  121

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  133

                                          xin(139) = c10*xin(115) + xc00*xin(127) + c01*xin(124)
                                          yin(139) = c10*yin(115) + yc00*yin(127) + c01*yin(124)
                                          zin(139) = c10*zin(115) + zc00*zin(127) + c01*zin(124)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =   73

                                          xin(82) = c10*xin(10) + xc00*xin(46) + c01*xin(43)
                                          yin(82) = c10*yin(10) + yc00*yin(46) + c01*yin(43)
                                          zin(82) = c10*zin(10) + zc00*zin(46) + c01*zin(43)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   73

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  109

                                          xin(118) = c10*xin(46) + xc00*xin(82) + c01*xin(79)
                                          yin(118) = c10*yin(46) + yc00*yin(82) + c01*yin(79)
                                          zin(118) = c10*zin(46) + zc00*zin(82) + c01*zin(79)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   73
                                          ! i4 = i5 =  109

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  121

                                          xin(130) = c10*xin(82) + xc00*xin(118) + c01*xin(115)
                                          yin(130) = c10*yin(82) + yc00*yin(118) + c01*yin(115)
                                          zin(130) = c10*zin(82) + zc00*zin(118) + c01*zin(115)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  109
                                          ! i4 = i5 =  121

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  133

                                          xin(142) = c10*xin(118) + xc00*xin(130) + c01*xin(127)
                                          yin(142) = c10*yin(118) + yc00*yin(130) + c01*yin(127)
                                          zin(142) = c10*zin(118) + zc00*zin(130) + c01*zin(127)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =   73

                                          xin(83) = c10*xin(11) + xc00*xin(47) + c01*xin(46)
                                          yin(83) = c10*yin(11) + yc00*yin(47) + c01*yin(46)
                                          zin(83) = c10*zin(11) + zc00*zin(47) + c01*zin(46)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   73

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  109

                                          xin(119) = c10*xin(47) + xc00*xin(83) + c01*xin(82)
                                          yin(119) = c10*yin(47) + yc00*yin(83) + c01*yin(82)
                                          zin(119) = c10*zin(47) + zc00*zin(83) + c01*zin(82)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   73
                                          ! i4 = i5 =  109

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  121

                                          xin(131) = c10*xin(83) + xc00*xin(119) + c01*xin(118)
                                          yin(131) = c10*yin(83) + yc00*yin(119) + c01*yin(118)
                                          zin(131) = c10*zin(83) + zc00*zin(119) + c01*zin(118)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  109
                                          ! i4 = i5 =  121

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  133

                                          xin(143) = c10*xin(119) + xc00*xin(131) + c01*xin(130)
                                          yin(143) = c10*yin(119) + yc00*yin(131) + c01*yin(130)
                                          zin(143) = c10*zin(119) + zc00*zin(131) + c01*zin(130)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =    1
                                          ! i4 = i2 =   37

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =   73

                                          xin(84) = c10*xin(12) + xc00*xin(48) + c01*xin(47)
                                          yin(84) = c10*yin(12) + yc00*yin(48) + c01*yin(47)
                                          zin(84) = c10*zin(12) + zc00*zin(48) + c01*zin(47)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   37
                                          ! i4 = i5 =   73

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  109

                                          xin(120) = c10*xin(48) + xc00*xin(84) + c01*xin(83)
                                          yin(120) = c10*yin(48) + yc00*yin(84) + c01*yin(83)
                                          zin(120) = c10*zin(48) + zc00*zin(84) + c01*zin(83)

                                          c10 = c10 + b10

                                          ! i3 = i4 =   73
                                          ! i4 = i5 =  109

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  121

                                          xin(132) = c10*xin(84) + xc00*xin(120) + c01*xin(119)
                                          yin(132) = c10*yin(84) + yc00*yin(120) + c01*yin(119)
                                          zin(132) = c10*zin(84) + zc00*zin(120) + c01*zin(119)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  109
                                          ! i4 = i5 =  121

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  133

                                          xin(144) = c10*xin(120) + xc00*xin(132) + c01*xin(131)
                                          yin(144) = c10*yin(120) + yc00*yin(132) + c01*yin(131)
                                          zin(144) = c10*zin(120) + zc00*zin(132) + c01*zin(131)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  121
                                          ! i4 = i5 =  133

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  133

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  133

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  121

                                          xin(133) = xin(133) + dxij*xin(121)
                                          yin(133) = yin(133) + dyij*yin(121)
                                          zin(133) = zin(133) + dzij*zin(121)

                                          ! i3 = i4 =  121
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  109

                                          xin(121) = xin(121) + dxij*xin(109)
                                          yin(121) = yin(121) + dyij*yin(109)
                                          zin(121) = zin(121) + dzij*zin(109)

                                          ! i3 = i4 =  109
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  133

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  121

                                          xin(133) = xin(133) + dxij*xin(121)
                                          yin(133) = yin(133) + dyij*yin(121)
                                          zin(133) = zin(133) + dzij*zin(121)

                                          ! i3 = i4 =  121
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   13

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   13

                                          ! do ni = 1,    3

                                          xin(13) = xin(37) + dxij*xin(1)
                                          yin(13) = yin(37) + dyij*yin(1)
                                          zin(13) = zin(37) + dzij*zin(1)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   49

                                          ! ni =    2

                                          xin(49) = xin(73) + dxij*xin(37)
                                          yin(49) = yin(73) + dyij*yin(37)
                                          zin(49) = zin(73) + dzij*zin(37)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   85

                                          ! ni =    3

                                          xin(85) = xin(109) + dxij*xin(73)
                                          yin(85) = yin(109) + dyij*yin(73)
                                          zin(85) = zin(109) + dzij*zin(73)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  121

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   25

                                          ! nj =    2

                                          ! i4 = i3 =   25

                                          ! do ni = 1,    3

                                          xin(25) = xin(49) + dxij*xin(13)
                                          yin(25) = yin(49) + dyij*yin(13)
                                          zin(25) = zin(49) + dzij*zin(13)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   61

                                          ! ni =    2

                                          xin(61) = xin(85) + dxij*xin(49)
                                          yin(61) = yin(85) + dyij*yin(49)
                                          zin(61) = zin(85) + dzij*zin(49)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   97

                                          ! ni =    3

                                          xin(97) = xin(121) + dxij*xin(85)
                                          yin(97) = yin(121) + dyij*yin(85)
                                          zin(97) = zin(121) + dzij*zin(85)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  133

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   37

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  136

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  124

                                          xin(136) = xin(136) + dxij*xin(124)
                                          yin(136) = yin(136) + dyij*yin(124)
                                          zin(136) = zin(136) + dzij*zin(124)

                                          ! i3 = i4 =  124
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  112

                                          xin(124) = xin(124) + dxij*xin(112)
                                          yin(124) = yin(124) + dyij*yin(112)
                                          zin(124) = zin(124) + dzij*zin(112)

                                          ! i3 = i4 =  112
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  136

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  124

                                          xin(136) = xin(136) + dxij*xin(124)
                                          yin(136) = yin(136) + dyij*yin(124)
                                          zin(136) = zin(136) + dzij*zin(124)

                                          ! i3 = i4 =  124
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   16

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   16

                                          ! do ni = 1,    3

                                          xin(16) = xin(40) + dxij*xin(4)
                                          yin(16) = yin(40) + dyij*yin(4)
                                          zin(16) = zin(40) + dzij*zin(4)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   52

                                          ! ni =    2

                                          xin(52) = xin(76) + dxij*xin(40)
                                          yin(52) = yin(76) + dyij*yin(40)
                                          zin(52) = zin(76) + dzij*zin(40)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   88

                                          ! ni =    3

                                          xin(88) = xin(112) + dxij*xin(76)
                                          yin(88) = yin(112) + dyij*yin(76)
                                          zin(88) = zin(112) + dzij*zin(76)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  124

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   28

                                          ! nj =    2

                                          ! i4 = i3 =   28

                                          ! do ni = 1,    3

                                          xin(28) = xin(52) + dxij*xin(16)
                                          yin(28) = yin(52) + dyij*yin(16)
                                          zin(28) = zin(52) + dzij*zin(16)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   64

                                          ! ni =    2

                                          xin(64) = xin(88) + dxij*xin(52)
                                          yin(64) = yin(88) + dyij*yin(52)
                                          zin(64) = zin(88) + dzij*zin(52)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  100

                                          ! ni =    3

                                          xin(100) = xin(124) + dxij*xin(88)
                                          yin(100) = yin(124) + dyij*yin(88)
                                          zin(100) = zin(124) + dzij*zin(88)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  136

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   40

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  139

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  127

                                          xin(139) = xin(139) + dxij*xin(127)
                                          yin(139) = yin(139) + dyij*yin(127)
                                          zin(139) = zin(139) + dzij*zin(127)

                                          ! i3 = i4 =  127
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  115

                                          xin(127) = xin(127) + dxij*xin(115)
                                          yin(127) = yin(127) + dyij*yin(115)
                                          zin(127) = zin(127) + dzij*zin(115)

                                          ! i3 = i4 =  115
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  139

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  127

                                          xin(139) = xin(139) + dxij*xin(127)
                                          yin(139) = yin(139) + dyij*yin(127)
                                          zin(139) = zin(139) + dzij*zin(127)

                                          ! i3 = i4 =  127
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   19

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   19

                                          ! do ni = 1,    3

                                          xin(19) = xin(43) + dxij*xin(7)
                                          yin(19) = yin(43) + dyij*yin(7)
                                          zin(19) = zin(43) + dzij*zin(7)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   55

                                          ! ni =    2

                                          xin(55) = xin(79) + dxij*xin(43)
                                          yin(55) = yin(79) + dyij*yin(43)
                                          zin(55) = zin(79) + dzij*zin(43)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   91

                                          ! ni =    3

                                          xin(91) = xin(115) + dxij*xin(79)
                                          yin(91) = yin(115) + dyij*yin(79)
                                          zin(91) = zin(115) + dzij*zin(79)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  127

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   31

                                          ! nj =    2

                                          ! i4 = i3 =   31

                                          ! do ni = 1,    3

                                          xin(31) = xin(55) + dxij*xin(19)
                                          yin(31) = yin(55) + dyij*yin(19)
                                          zin(31) = zin(55) + dzij*zin(19)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   67

                                          ! ni =    2

                                          xin(67) = xin(91) + dxij*xin(55)
                                          yin(67) = yin(91) + dyij*yin(55)
                                          zin(67) = zin(91) + dzij*zin(55)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  103

                                          ! ni =    3

                                          xin(103) = xin(127) + dxij*xin(91)
                                          yin(103) = yin(127) + dyij*yin(91)
                                          zin(103) = zin(127) + dzij*zin(91)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  139

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   43

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  142

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  130

                                          xin(142) = xin(142) + dxij*xin(130)
                                          yin(142) = yin(142) + dyij*yin(130)
                                          zin(142) = zin(142) + dzij*zin(130)

                                          ! i3 = i4 =  130
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  118

                                          xin(130) = xin(130) + dxij*xin(118)
                                          yin(130) = yin(130) + dyij*yin(118)
                                          zin(130) = zin(130) + dzij*zin(118)

                                          ! i3 = i4 =  118
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  142

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  130

                                          xin(142) = xin(142) + dxij*xin(130)
                                          yin(142) = yin(142) + dyij*yin(130)
                                          zin(142) = zin(142) + dzij*zin(130)

                                          ! i3 = i4 =  130
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   22

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   22

                                          ! do ni = 1,    3

                                          xin(22) = xin(46) + dxij*xin(10)
                                          yin(22) = yin(46) + dyij*yin(10)
                                          zin(22) = zin(46) + dzij*zin(10)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   58

                                          ! ni =    2

                                          xin(58) = xin(82) + dxij*xin(46)
                                          yin(58) = yin(82) + dyij*yin(46)
                                          zin(58) = zin(82) + dzij*zin(46)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   94

                                          ! ni =    3

                                          xin(94) = xin(118) + dxij*xin(82)
                                          yin(94) = yin(118) + dyij*yin(82)
                                          zin(94) = zin(118) + dzij*zin(82)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  130

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   34

                                          ! nj =    2

                                          ! i4 = i3 =   34

                                          ! do ni = 1,    3

                                          xin(34) = xin(58) + dxij*xin(22)
                                          yin(34) = yin(58) + dyij*yin(22)
                                          zin(34) = zin(58) + dzij*zin(22)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   70

                                          ! ni =    2

                                          xin(70) = xin(94) + dxij*xin(58)
                                          yin(70) = yin(94) + dyij*yin(58)
                                          zin(70) = zin(94) + dzij*zin(58)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  106

                                          ! ni =    3

                                          xin(106) = xin(130) + dxij*xin(94)
                                          yin(106) = yin(130) + dyij*yin(94)
                                          zin(106) = zin(130) + dzij*zin(94)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  142

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   46

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  143

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  131

                                          xin(143) = xin(143) + dxij*xin(131)
                                          yin(143) = yin(143) + dyij*yin(131)
                                          zin(143) = zin(143) + dzij*zin(131)

                                          ! i3 = i4 =  131
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  119

                                          xin(131) = xin(131) + dxij*xin(119)
                                          yin(131) = yin(131) + dyij*yin(119)
                                          zin(131) = zin(131) + dzij*zin(119)

                                          ! i3 = i4 =  119
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  143

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  131

                                          xin(143) = xin(143) + dxij*xin(131)
                                          yin(143) = yin(143) + dyij*yin(131)
                                          zin(143) = zin(143) + dzij*zin(131)

                                          ! i3 = i4 =  131
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   23

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   23

                                          ! do ni = 1,    3

                                          xin(23) = xin(47) + dxij*xin(11)
                                          yin(23) = yin(47) + dyij*yin(11)
                                          zin(23) = zin(47) + dzij*zin(11)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   59

                                          ! ni =    2

                                          xin(59) = xin(83) + dxij*xin(47)
                                          yin(59) = yin(83) + dyij*yin(47)
                                          zin(59) = zin(83) + dzij*zin(47)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   95

                                          ! ni =    3

                                          xin(95) = xin(119) + dxij*xin(83)
                                          yin(95) = yin(119) + dyij*yin(83)
                                          zin(95) = zin(119) + dzij*zin(83)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  131

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   35

                                          ! nj =    2

                                          ! i4 = i3 =   35

                                          ! do ni = 1,    3

                                          xin(35) = xin(59) + dxij*xin(23)
                                          yin(35) = yin(59) + dyij*yin(23)
                                          zin(35) = zin(59) + dzij*zin(23)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   71

                                          ! ni =    2

                                          xin(71) = xin(95) + dxij*xin(59)
                                          yin(71) = yin(95) + dyij*yin(59)
                                          zin(71) = zin(95) + dzij*zin(59)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  107

                                          ! ni =    3

                                          xin(107) = xin(131) + dxij*xin(95)
                                          yin(107) = yin(131) + dyij*yin(95)
                                          zin(107) = zin(131) + dzij*zin(95)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  143

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   47

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  144

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  132

                                          xin(144) = xin(144) + dxij*xin(132)
                                          yin(144) = yin(144) + dyij*yin(132)
                                          zin(144) = zin(144) + dzij*zin(132)

                                          ! i3 = i4 =  132
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  120

                                          xin(132) = xin(132) + dxij*xin(120)
                                          yin(132) = yin(132) + dyij*yin(120)
                                          zin(132) = zin(132) + dzij*zin(120)

                                          ! i3 = i4 =  120
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  144

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  132

                                          xin(144) = xin(144) + dxij*xin(132)
                                          yin(144) = yin(144) + dyij*yin(132)
                                          zin(144) = zin(144) + dzij*zin(132)

                                          ! i3 = i4 =  132
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =   24

                                          ! do nj = 1,    2

                                          ! i4 = i3 =   24

                                          ! do ni = 1,    3

                                          xin(24) = xin(48) + dxij*xin(12)
                                          yin(24) = yin(48) + dyij*yin(12)
                                          zin(24) = zin(48) + dzij*zin(12)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   60

                                          ! ni =    2

                                          xin(60) = xin(84) + dxij*xin(48)
                                          yin(60) = yin(84) + dyij*yin(48)
                                          zin(60) = zin(84) + dzij*zin(48)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   96

                                          ! ni =    3

                                          xin(96) = xin(120) + dxij*xin(84)
                                          yin(96) = yin(120) + dyij*yin(84)
                                          zin(96) = zin(120) + dzij*zin(84)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  132

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   36

                                          ! nj =    2

                                          ! i4 = i3 =   36

                                          ! do ni = 1,    3

                                          xin(36) = xin(60) + dxij*xin(24)
                                          yin(36) = yin(60) + dyij*yin(24)
                                          zin(36) = zin(60) + dzij*zin(24)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =   72

                                          ! ni =    2

                                          xin(72) = xin(96) + dxij*xin(60)
                                          yin(72) = yin(96) + dyij*yin(60)
                                          zin(72) = zin(96) + dzij*zin(60)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  108

                                          ! ni =    3

                                          xin(108) = xin(132) + dxij*xin(96)
                                          yin(108) = yin(132) + dyij*yin(96)
                                          zin(108) = zin(132) + dzij*zin(96)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  144

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =   48

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =    1

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   12

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   11

                                          xin(12) = xin(12) + dxkl*xin(11)
                                          yin(12) = yin(12) + dykl*yin(11)
                                          zin(12) = zin(12) + dzkl*zin(11)

                                          ! i3 = i4 =   11
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   10

                                          xin(11) = xin(11) + dxkl*xin(10)
                                          yin(11) = yin(11) + dykl*yin(10)
                                          zin(11) = zin(11) + dzkl*zin(10)

                                          ! i3 = i4 =   10
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   12

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   11

                                          xin(12) = xin(12) + dxkl*xin(11)
                                          yin(12) = yin(12) + dykl*yin(11)
                                          zin(12) = zin(12) + dzkl*zin(11)

                                          ! i3 = i4 =   11
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =    2

                                          ! do nl = 1,    2

                                          ! i4 = i3 =    2

                                          ! do nk = 1,    3

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

                                          xin(8) = xin(10) + dxkl*xin(7)
                                          yin(8) = yin(10) + dykl*yin(7)
                                          zin(8) = zin(10) + dzkl*zin(7)
                                          ! i4 = i4 + lang+1 =   11

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =    3

                                          ! nl =    2

                                          ! i4 = i3 =    3

                                          ! do nk = 1,    3

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

                                          xin(9) = xin(11) + dxkl*xin(8)
                                          yin(9) = yin(11) + dykl*yin(8)
                                          zin(9) = zin(11) + dzkl*zin(8)
                                          ! i4 = i4 + lang+1 =   12

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =    4

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   13

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   24

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   23

                                          xin(24) = xin(24) + dxkl*xin(23)
                                          yin(24) = yin(24) + dykl*yin(23)
                                          zin(24) = zin(24) + dzkl*zin(23)

                                          ! i3 = i4 =   23
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   22

                                          xin(23) = xin(23) + dxkl*xin(22)
                                          yin(23) = yin(23) + dykl*yin(22)
                                          zin(23) = zin(23) + dzkl*zin(22)

                                          ! i3 = i4 =   22
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   24

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   23

                                          xin(24) = xin(24) + dxkl*xin(23)
                                          yin(24) = yin(24) + dykl*yin(23)
                                          zin(24) = zin(24) + dzkl*zin(23)

                                          ! i3 = i4 =   23
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   14

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   14

                                          ! do nk = 1,    3

                                          xin(14) = xin(16) + dxkl*xin(13)
                                          yin(14) = yin(16) + dykl*yin(13)
                                          zin(14) = zin(16) + dzkl*zin(13)
                                          ! i4 = i4 + lang+1 =   17

                                          ! nk =    2

                                          xin(17) = xin(19) + dxkl*xin(16)
                                          yin(17) = yin(19) + dykl*yin(16)
                                          zin(17) = zin(19) + dzkl*zin(16)
                                          ! i4 = i4 + lang+1 =   20

                                          ! nk =    3

                                          xin(20) = xin(22) + dxkl*xin(19)
                                          yin(20) = yin(22) + dykl*yin(19)
                                          zin(20) = zin(22) + dzkl*zin(19)
                                          ! i4 = i4 + lang+1 =   23

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   15

                                          ! nl =    2

                                          ! i4 = i3 =   15

                                          ! do nk = 1,    3

                                          xin(15) = xin(17) + dxkl*xin(14)
                                          yin(15) = yin(17) + dykl*yin(14)
                                          zin(15) = zin(17) + dzkl*zin(14)
                                          ! i4 = i4 + lang+1 =   18

                                          ! nk =    2

                                          xin(18) = xin(20) + dxkl*xin(17)
                                          yin(18) = yin(20) + dykl*yin(17)
                                          zin(18) = zin(20) + dzkl*zin(17)
                                          ! i4 = i4 + lang+1 =   21

                                          ! nk =    3

                                          xin(21) = xin(23) + dxkl*xin(20)
                                          yin(21) = yin(23) + dykl*yin(20)
                                          zin(21) = zin(23) + dzkl*zin(20)
                                          ! i4 = i4 + lang+1 =   24

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   16

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   25

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   36

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   35

                                          xin(36) = xin(36) + dxkl*xin(35)
                                          yin(36) = yin(36) + dykl*yin(35)
                                          zin(36) = zin(36) + dzkl*zin(35)

                                          ! i3 = i4 =   35
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   34

                                          xin(35) = xin(35) + dxkl*xin(34)
                                          yin(35) = yin(35) + dykl*yin(34)
                                          zin(35) = zin(35) + dzkl*zin(34)

                                          ! i3 = i4 =   34
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   36

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   35

                                          xin(36) = xin(36) + dxkl*xin(35)
                                          yin(36) = yin(36) + dykl*yin(35)
                                          zin(36) = zin(36) + dzkl*zin(35)

                                          ! i3 = i4 =   35
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   26

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   26

                                          ! do nk = 1,    3

                                          xin(26) = xin(28) + dxkl*xin(25)
                                          yin(26) = yin(28) + dykl*yin(25)
                                          zin(26) = zin(28) + dzkl*zin(25)
                                          ! i4 = i4 + lang+1 =   29

                                          ! nk =    2

                                          xin(29) = xin(31) + dxkl*xin(28)
                                          yin(29) = yin(31) + dykl*yin(28)
                                          zin(29) = zin(31) + dzkl*zin(28)
                                          ! i4 = i4 + lang+1 =   32

                                          ! nk =    3

                                          xin(32) = xin(34) + dxkl*xin(31)
                                          yin(32) = yin(34) + dykl*yin(31)
                                          zin(32) = zin(34) + dzkl*zin(31)
                                          ! i4 = i4 + lang+1 =   35

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   27

                                          ! nl =    2

                                          ! i4 = i3 =   27

                                          ! do nk = 1,    3

                                          xin(27) = xin(29) + dxkl*xin(26)
                                          yin(27) = yin(29) + dykl*yin(26)
                                          zin(27) = zin(29) + dzkl*zin(26)
                                          ! i4 = i4 + lang+1 =   30

                                          ! nk =    2

                                          xin(30) = xin(32) + dxkl*xin(29)
                                          yin(30) = yin(32) + dykl*yin(29)
                                          zin(30) = zin(32) + dzkl*zin(29)
                                          ! i4 = i4 + lang+1 =   33

                                          ! nk =    3

                                          xin(33) = xin(35) + dxkl*xin(32)
                                          yin(33) = yin(35) + dykl*yin(32)
                                          zin(33) = zin(35) + dzkl*zin(32)
                                          ! i4 = i4 + lang+1 =   36

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   28

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   37

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   37

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   48

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   47

                                          xin(48) = xin(48) + dxkl*xin(47)
                                          yin(48) = yin(48) + dykl*yin(47)
                                          zin(48) = zin(48) + dzkl*zin(47)

                                          ! i3 = i4 =   47
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   46

                                          xin(47) = xin(47) + dxkl*xin(46)
                                          yin(47) = yin(47) + dykl*yin(46)
                                          zin(47) = zin(47) + dzkl*zin(46)

                                          ! i3 = i4 =   46
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   48

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   47

                                          xin(48) = xin(48) + dxkl*xin(47)
                                          yin(48) = yin(48) + dykl*yin(47)
                                          zin(48) = zin(48) + dzkl*zin(47)

                                          ! i3 = i4 =   47
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   38

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   38

                                          ! do nk = 1,    3

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

                                          xin(44) = xin(46) + dxkl*xin(43)
                                          yin(44) = yin(46) + dykl*yin(43)
                                          zin(44) = zin(46) + dzkl*zin(43)
                                          ! i4 = i4 + lang+1 =   47

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   39

                                          ! nl =    2

                                          ! i4 = i3 =   39

                                          ! do nk = 1,    3

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

                                          xin(45) = xin(47) + dxkl*xin(44)
                                          yin(45) = yin(47) + dykl*yin(44)
                                          zin(45) = zin(47) + dzkl*zin(44)
                                          ! i4 = i4 + lang+1 =   48

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   40

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   49

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   60

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   59

                                          xin(60) = xin(60) + dxkl*xin(59)
                                          yin(60) = yin(60) + dykl*yin(59)
                                          zin(60) = zin(60) + dzkl*zin(59)

                                          ! i3 = i4 =   59
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   58

                                          xin(59) = xin(59) + dxkl*xin(58)
                                          yin(59) = yin(59) + dykl*yin(58)
                                          zin(59) = zin(59) + dzkl*zin(58)

                                          ! i3 = i4 =   58
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   60

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   59

                                          xin(60) = xin(60) + dxkl*xin(59)
                                          yin(60) = yin(60) + dykl*yin(59)
                                          zin(60) = zin(60) + dzkl*zin(59)

                                          ! i3 = i4 =   59
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   50

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   50

                                          ! do nk = 1,    3

                                          xin(50) = xin(52) + dxkl*xin(49)
                                          yin(50) = yin(52) + dykl*yin(49)
                                          zin(50) = zin(52) + dzkl*zin(49)
                                          ! i4 = i4 + lang+1 =   53

                                          ! nk =    2

                                          xin(53) = xin(55) + dxkl*xin(52)
                                          yin(53) = yin(55) + dykl*yin(52)
                                          zin(53) = zin(55) + dzkl*zin(52)
                                          ! i4 = i4 + lang+1 =   56

                                          ! nk =    3

                                          xin(56) = xin(58) + dxkl*xin(55)
                                          yin(56) = yin(58) + dykl*yin(55)
                                          zin(56) = zin(58) + dzkl*zin(55)
                                          ! i4 = i4 + lang+1 =   59

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   51

                                          ! nl =    2

                                          ! i4 = i3 =   51

                                          ! do nk = 1,    3

                                          xin(51) = xin(53) + dxkl*xin(50)
                                          yin(51) = yin(53) + dykl*yin(50)
                                          zin(51) = zin(53) + dzkl*zin(50)
                                          ! i4 = i4 + lang+1 =   54

                                          ! nk =    2

                                          xin(54) = xin(56) + dxkl*xin(53)
                                          yin(54) = yin(56) + dykl*yin(53)
                                          zin(54) = zin(56) + dzkl*zin(53)
                                          ! i4 = i4 + lang+1 =   57

                                          ! nk =    3

                                          xin(57) = xin(59) + dxkl*xin(56)
                                          yin(57) = yin(59) + dykl*yin(56)
                                          zin(57) = zin(59) + dzkl*zin(56)
                                          ! i4 = i4 + lang+1 =   60

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   52

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   61

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   72

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   71

                                          xin(72) = xin(72) + dxkl*xin(71)
                                          yin(72) = yin(72) + dykl*yin(71)
                                          zin(72) = zin(72) + dzkl*zin(71)

                                          ! i3 = i4 =   71
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   70

                                          xin(71) = xin(71) + dxkl*xin(70)
                                          yin(71) = yin(71) + dykl*yin(70)
                                          zin(71) = zin(71) + dzkl*zin(70)

                                          ! i3 = i4 =   70
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   72

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   71

                                          xin(72) = xin(72) + dxkl*xin(71)
                                          yin(72) = yin(72) + dykl*yin(71)
                                          zin(72) = zin(72) + dzkl*zin(71)

                                          ! i3 = i4 =   71
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   62

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   62

                                          ! do nk = 1,    3

                                          xin(62) = xin(64) + dxkl*xin(61)
                                          yin(62) = yin(64) + dykl*yin(61)
                                          zin(62) = zin(64) + dzkl*zin(61)
                                          ! i4 = i4 + lang+1 =   65

                                          ! nk =    2

                                          xin(65) = xin(67) + dxkl*xin(64)
                                          yin(65) = yin(67) + dykl*yin(64)
                                          zin(65) = zin(67) + dzkl*zin(64)
                                          ! i4 = i4 + lang+1 =   68

                                          ! nk =    3

                                          xin(68) = xin(70) + dxkl*xin(67)
                                          yin(68) = yin(70) + dykl*yin(67)
                                          zin(68) = zin(70) + dzkl*zin(67)
                                          ! i4 = i4 + lang+1 =   71

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   63

                                          ! nl =    2

                                          ! i4 = i3 =   63

                                          ! do nk = 1,    3

                                          xin(63) = xin(65) + dxkl*xin(62)
                                          yin(63) = yin(65) + dykl*yin(62)
                                          zin(63) = zin(65) + dzkl*zin(62)
                                          ! i4 = i4 + lang+1 =   66

                                          ! nk =    2

                                          xin(66) = xin(68) + dxkl*xin(65)
                                          yin(66) = yin(68) + dykl*yin(65)
                                          zin(66) = zin(68) + dzkl*zin(65)
                                          ! i4 = i4 + lang+1 =   69

                                          ! nk =    3

                                          xin(69) = xin(71) + dxkl*xin(68)
                                          yin(69) = yin(71) + dykl*yin(68)
                                          zin(69) = zin(71) + dzkl*zin(68)
                                          ! i4 = i4 + lang+1 =   72

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   64

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   73

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =   73

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   84

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   83

                                          xin(84) = xin(84) + dxkl*xin(83)
                                          yin(84) = yin(84) + dykl*yin(83)
                                          zin(84) = zin(84) + dzkl*zin(83)

                                          ! i3 = i4 =   83
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   82

                                          xin(83) = xin(83) + dxkl*xin(82)
                                          yin(83) = yin(83) + dykl*yin(82)
                                          zin(83) = zin(83) + dzkl*zin(82)

                                          ! i3 = i4 =   82
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   84

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   83

                                          xin(84) = xin(84) + dxkl*xin(83)
                                          yin(84) = yin(84) + dykl*yin(83)
                                          zin(84) = zin(84) + dzkl*zin(83)

                                          ! i3 = i4 =   83
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   74

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   74

                                          ! do nk = 1,    3

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

                                          xin(80) = xin(82) + dxkl*xin(79)
                                          yin(80) = yin(82) + dykl*yin(79)
                                          zin(80) = zin(82) + dzkl*zin(79)
                                          ! i4 = i4 + lang+1 =   83

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   75

                                          ! nl =    2

                                          ! i4 = i3 =   75

                                          ! do nk = 1,    3

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

                                          xin(81) = xin(83) + dxkl*xin(80)
                                          yin(81) = yin(83) + dykl*yin(80)
                                          zin(81) = zin(83) + dzkl*zin(80)
                                          ! i4 = i4 + lang+1 =   84

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   76

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   85

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   96

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   95

                                          xin(96) = xin(96) + dxkl*xin(95)
                                          yin(96) = yin(96) + dykl*yin(95)
                                          zin(96) = zin(96) + dzkl*zin(95)

                                          ! i3 = i4 =   95
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =   94

                                          xin(95) = xin(95) + dxkl*xin(94)
                                          yin(95) = yin(95) + dykl*yin(94)
                                          zin(95) = zin(95) + dzkl*zin(94)

                                          ! i3 = i4 =   94
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =   96

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =   95

                                          xin(96) = xin(96) + dxkl*xin(95)
                                          yin(96) = yin(96) + dykl*yin(95)
                                          zin(96) = zin(96) + dzkl*zin(95)

                                          ! i3 = i4 =   95
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   86

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   86

                                          ! do nk = 1,    3

                                          xin(86) = xin(88) + dxkl*xin(85)
                                          yin(86) = yin(88) + dykl*yin(85)
                                          zin(86) = zin(88) + dzkl*zin(85)
                                          ! i4 = i4 + lang+1 =   89

                                          ! nk =    2

                                          xin(89) = xin(91) + dxkl*xin(88)
                                          yin(89) = yin(91) + dykl*yin(88)
                                          zin(89) = zin(91) + dzkl*zin(88)
                                          ! i4 = i4 + lang+1 =   92

                                          ! nk =    3

                                          xin(92) = xin(94) + dxkl*xin(91)
                                          yin(92) = yin(94) + dykl*yin(91)
                                          zin(92) = zin(94) + dzkl*zin(91)
                                          ! i4 = i4 + lang+1 =   95

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   87

                                          ! nl =    2

                                          ! i4 = i3 =   87

                                          ! do nk = 1,    3

                                          xin(87) = xin(89) + dxkl*xin(86)
                                          yin(87) = yin(89) + dykl*yin(86)
                                          zin(87) = zin(89) + dzkl*zin(86)
                                          ! i4 = i4 + lang+1 =   90

                                          ! nk =    2

                                          xin(90) = xin(92) + dxkl*xin(89)
                                          yin(90) = yin(92) + dykl*yin(89)
                                          zin(90) = zin(92) + dzkl*zin(89)
                                          ! i4 = i4 + lang+1 =   93

                                          ! nk =    3

                                          xin(93) = xin(95) + dxkl*xin(92)
                                          yin(93) = yin(95) + dykl*yin(92)
                                          zin(93) = zin(95) + dzkl*zin(92)
                                          ! i4 = i4 + lang+1 =   96

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   88

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =   97

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  108

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  107

                                          xin(108) = xin(108) + dxkl*xin(107)
                                          yin(108) = yin(108) + dykl*yin(107)
                                          zin(108) = zin(108) + dzkl*zin(107)

                                          ! i3 = i4 =  107
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  106

                                          xin(107) = xin(107) + dxkl*xin(106)
                                          yin(107) = yin(107) + dykl*yin(106)
                                          zin(107) = zin(107) + dzkl*zin(106)

                                          ! i3 = i4 =  106
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  108

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  107

                                          xin(108) = xin(108) + dxkl*xin(107)
                                          yin(108) = yin(108) + dykl*yin(107)
                                          zin(108) = zin(108) + dzkl*zin(107)

                                          ! i3 = i4 =  107
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =   98

                                          ! do nl = 1,    2

                                          ! i4 = i3 =   98

                                          ! do nk = 1,    3

                                          xin(98) = xin(100) + dxkl*xin(97)
                                          yin(98) = yin(100) + dykl*yin(97)
                                          zin(98) = zin(100) + dzkl*zin(97)
                                          ! i4 = i4 + lang+1 =  101

                                          ! nk =    2

                                          xin(101) = xin(103) + dxkl*xin(100)
                                          yin(101) = yin(103) + dykl*yin(100)
                                          zin(101) = zin(103) + dzkl*zin(100)
                                          ! i4 = i4 + lang+1 =  104

                                          ! nk =    3

                                          xin(104) = xin(106) + dxkl*xin(103)
                                          yin(104) = yin(106) + dykl*yin(103)
                                          zin(104) = zin(106) + dzkl*zin(103)
                                          ! i4 = i4 + lang+1 =  107

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =   99

                                          ! nl =    2

                                          ! i4 = i3 =   99

                                          ! do nk = 1,    3

                                          xin(99) = xin(101) + dxkl*xin(98)
                                          yin(99) = yin(101) + dykl*yin(98)
                                          zin(99) = zin(101) + dzkl*zin(98)
                                          ! i4 = i4 + lang+1 =  102

                                          ! nk =    2

                                          xin(102) = xin(104) + dxkl*xin(101)
                                          yin(102) = yin(104) + dykl*yin(101)
                                          zin(102) = zin(104) + dzkl*zin(101)
                                          ! i4 = i4 + lang+1 =  105

                                          ! nk =    3

                                          xin(105) = xin(107) + dxkl*xin(104)
                                          yin(105) = yin(107) + dykl*yin(104)
                                          zin(105) = zin(107) + dzkl*zin(104)
                                          ! i4 = i4 + lang+1 =  108

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  100

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  109

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  109

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  120

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  119

                                          xin(120) = xin(120) + dxkl*xin(119)
                                          yin(120) = yin(120) + dykl*yin(119)
                                          zin(120) = zin(120) + dzkl*zin(119)

                                          ! i3 = i4 =  119
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  118

                                          xin(119) = xin(119) + dxkl*xin(118)
                                          yin(119) = yin(119) + dykl*yin(118)
                                          zin(119) = zin(119) + dzkl*zin(118)

                                          ! i3 = i4 =  118
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  120

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  119

                                          xin(120) = xin(120) + dxkl*xin(119)
                                          yin(120) = yin(120) + dykl*yin(119)
                                          zin(120) = zin(120) + dzkl*zin(119)

                                          ! i3 = i4 =  119
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  110

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  110

                                          ! do nk = 1,    3

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

                                          xin(116) = xin(118) + dxkl*xin(115)
                                          yin(116) = yin(118) + dykl*yin(115)
                                          zin(116) = zin(118) + dzkl*zin(115)
                                          ! i4 = i4 + lang+1 =  119

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  111

                                          ! nl =    2

                                          ! i4 = i3 =  111

                                          ! do nk = 1,    3

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

                                          xin(117) = xin(119) + dxkl*xin(116)
                                          yin(117) = yin(119) + dykl*yin(116)
                                          zin(117) = zin(119) + dzkl*zin(116)
                                          ! i4 = i4 + lang+1 =  120

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  112

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  121

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  132

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  131

                                          xin(132) = xin(132) + dxkl*xin(131)
                                          yin(132) = yin(132) + dykl*yin(131)
                                          zin(132) = zin(132) + dzkl*zin(131)

                                          ! i3 = i4 =  131
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  130

                                          xin(131) = xin(131) + dxkl*xin(130)
                                          yin(131) = yin(131) + dykl*yin(130)
                                          zin(131) = zin(131) + dzkl*zin(130)

                                          ! i3 = i4 =  130
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  132

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  131

                                          xin(132) = xin(132) + dxkl*xin(131)
                                          yin(132) = yin(132) + dykl*yin(131)
                                          zin(132) = zin(132) + dzkl*zin(131)

                                          ! i3 = i4 =  131
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  122

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  122

                                          ! do nk = 1,    3

                                          xin(122) = xin(124) + dxkl*xin(121)
                                          yin(122) = yin(124) + dykl*yin(121)
                                          zin(122) = zin(124) + dzkl*zin(121)
                                          ! i4 = i4 + lang+1 =  125

                                          ! nk =    2

                                          xin(125) = xin(127) + dxkl*xin(124)
                                          yin(125) = yin(127) + dykl*yin(124)
                                          zin(125) = zin(127) + dzkl*zin(124)
                                          ! i4 = i4 + lang+1 =  128

                                          ! nk =    3

                                          xin(128) = xin(130) + dxkl*xin(127)
                                          yin(128) = yin(130) + dykl*yin(127)
                                          zin(128) = zin(130) + dzkl*zin(127)
                                          ! i4 = i4 + lang+1 =  131

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  123

                                          ! nl =    2

                                          ! i4 = i3 =  123

                                          ! do nk = 1,    3

                                          xin(123) = xin(125) + dxkl*xin(122)
                                          yin(123) = yin(125) + dykl*yin(122)
                                          zin(123) = zin(125) + dzkl*zin(122)
                                          ! i4 = i4 + lang+1 =  126

                                          ! nk =    2

                                          xin(126) = xin(128) + dxkl*xin(125)
                                          yin(126) = yin(128) + dykl*yin(125)
                                          zin(126) = zin(128) + dzkl*zin(125)
                                          ! i4 = i4 + lang+1 =  129

                                          ! nk =    3

                                          xin(129) = xin(131) + dxkl*xin(128)
                                          yin(129) = yin(131) + dykl*yin(128)
                                          zin(129) = zin(131) + dzkl*zin(128)
                                          ! i4 = i4 + lang+1 =  132

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  124

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  133

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  144

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  143

                                          xin(144) = xin(144) + dxkl*xin(143)
                                          yin(144) = yin(144) + dykl*yin(143)
                                          zin(144) = zin(144) + dzkl*zin(143)

                                          ! i3 = i4 =  143
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  142

                                          xin(143) = xin(143) + dxkl*xin(142)
                                          yin(143) = yin(143) + dykl*yin(142)
                                          zin(143) = zin(143) + dzkl*zin(142)

                                          ! i3 = i4 =  142
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  144

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  143

                                          xin(144) = xin(144) + dxkl*xin(143)
                                          yin(144) = yin(144) + dykl*yin(143)
                                          zin(144) = zin(144) + dzkl*zin(143)

                                          ! i3 = i4 =  143
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  134

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  134

                                          ! do nk = 1,    3

                                          xin(134) = xin(136) + dxkl*xin(133)
                                          yin(134) = yin(136) + dykl*yin(133)
                                          zin(134) = zin(136) + dzkl*zin(133)
                                          ! i4 = i4 + lang+1 =  137

                                          ! nk =    2

                                          xin(137) = xin(139) + dxkl*xin(136)
                                          yin(137) = yin(139) + dykl*yin(136)
                                          zin(137) = zin(139) + dzkl*zin(136)
                                          ! i4 = i4 + lang+1 =  140

                                          ! nk =    3

                                          xin(140) = xin(142) + dxkl*xin(139)
                                          yin(140) = yin(142) + dykl*yin(139)
                                          zin(140) = zin(142) + dzkl*zin(139)
                                          ! i4 = i4 + lang+1 =  143

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  135

                                          ! nl =    2

                                          ! i4 = i3 =  135

                                          ! do nk = 1,    3

                                          xin(135) = xin(137) + dxkl*xin(134)
                                          yin(135) = yin(137) + dykl*yin(134)
                                          zin(135) = zin(137) + dzkl*zin(134)
                                          ! i4 = i4 + lang+1 =  138

                                          ! nk =    2

                                          xin(138) = xin(140) + dxkl*xin(137)
                                          yin(138) = yin(140) + dykl*yin(137)
                                          zin(138) = zin(140) + dzkl*zin(137)
                                          ! i4 = i4 + lang+1 =  141

                                          ! nk =    3

                                          xin(141) = xin(143) + dxkl*xin(140)
                                          yin(141) = yin(143) + dykl*yin(140)
                                          zin(141) = zin(143) + dzkl*zin(140)
                                          ! i4 = i4 + lang+1 =  144

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  136

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  145

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  145

                                          ! end do

                                          ! *** Now root =    2

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  144

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

                                          ! i1 = in(1) =  145

                                          xin(145) = 1.0_dp
                                          yin(145) = 1.0_dp
                                          zin(145) = f00

                                          ! i2 = in(2) =  181
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(181) = xc00
                                          yin(181) = yc00
                                          zin(181) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  148

                                          xin(148) = xcp00
                                          yin(148) = ycp00
                                          zin(148) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  184
                                          ! i2 =  181

                                          xin(184) = xcp00*xin(181) + cp10
                                          yin(184) = ycp00*yin(181) + cp10
                                          zin(184) = zcp00*zin(181) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  145
                                          ! i4 = i2 =  181

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  217
                                          ! i3 =  145
                                          ! i4 =  181

                                          xin(217) = c10*xin(145) + xc00*xin(181)
                                          yin(217) = c10*yin(145) + yc00*yin(181)
                                          zin(217) = c10*zin(145) + zc00*zin(181)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  220
                                          ! i5 =  217
                                          ! i4 =  181

                                          xin(220) = xcp00*xin(217) + cp10*xin(181)
                                          yin(220) = ycp00*yin(217) + cp10*yin(181)
                                          zin(220) = zcp00*zin(217) + cp10*zin(181)

                                          ! ------------------

                                          ! i3 = i4 =  181
                                          ! i4 = i5 =  217

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  253
                                          ! i3 =  181
                                          ! i4 =  217

                                          xin(253) = c10*xin(181) + xc00*xin(217)
                                          yin(253) = c10*yin(181) + yc00*yin(217)
                                          zin(253) = c10*zin(181) + zc00*zin(217)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  256
                                          ! i5 =  253
                                          ! i4 =  217

                                          xin(256) = xcp00*xin(253) + cp10*xin(217)
                                          yin(256) = ycp00*yin(253) + cp10*yin(217)
                                          zin(256) = zcp00*zin(253) + cp10*zin(217)

                                          ! ------------------

                                          ! i3 = i4 =  217
                                          ! i4 = i5 =  253

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  265
                                          ! i3 =  217
                                          ! i4 =  253

                                          xin(265) = c10*xin(217) + xc00*xin(253)
                                          yin(265) = c10*yin(217) + yc00*yin(253)
                                          zin(265) = c10*zin(217) + zc00*zin(253)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  268
                                          ! i5 =  265
                                          ! i4 =  253

                                          xin(268) = xcp00*xin(265) + cp10*xin(253)
                                          yin(268) = ycp00*yin(265) + cp10*yin(253)
                                          zin(268) = zcp00*zin(265) + cp10*zin(253)

                                          ! ------------------

                                          ! i3 = i4 =  253
                                          ! i4 = i5 =  265

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  277
                                          ! i3 =  253
                                          ! i4 =  265

                                          xin(277) = c10*xin(253) + xc00*xin(265)
                                          yin(277) = c10*yin(253) + yc00*yin(265)
                                          zin(277) = c10*zin(253) + zc00*zin(265)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  280
                                          ! i5 =  277
                                          ! i4 =  265

                                          xin(280) = xcp00*xin(277) + cp10*xin(265)
                                          yin(280) = ycp00*yin(277) + cp10*yin(265)
                                          zin(280) = zcp00*zin(277) + cp10*zin(265)

                                          ! ------------------

                                          ! i3 = i4 =  265
                                          ! i4 = i5 =  277

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =  145
                                          ! i4 = i1+k2 =  148

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  151
                                          ! i3 =  145
                                          ! i4 =  148

                                          xin(151) = cp01*xin(145) + xcp00*xin(148)
                                          yin(151) = cp01*yin(145) + ycp00*yin(148)
                                          zin(151) = cp01*zin(145) + zcp00*zin(148)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  187

                                          xin(187) = xc00*xin(151) + c01*xin(148)
                                          yin(187) = yc00*yin(151) + c01*yin(148)
                                          zin(187) = zc00*zin(151) + c01*zin(148)

                                          ! ------------------

                                          ! i3 = i4 =  148
                                          ! i4 = i5 =  151

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  154
                                          ! i3 =  148
                                          ! i4 =  151

                                          xin(154) = cp01*xin(148) + xcp00*xin(151)
                                          yin(154) = cp01*yin(148) + ycp00*yin(151)
                                          zin(154) = cp01*zin(148) + zcp00*zin(151)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  190

                                          xin(190) = xc00*xin(154) + c01*xin(151)
                                          yin(190) = yc00*yin(154) + c01*yin(151)
                                          zin(190) = zc00*zin(154) + c01*zin(151)

                                          ! ------------------

                                          ! i3 = i4 =  151
                                          ! i4 = i5 =  154

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  155
                                          ! i3 =  151
                                          ! i4 =  154

                                          xin(155) = cp01*xin(151) + xcp00*xin(154)
                                          yin(155) = cp01*yin(151) + ycp00*yin(154)
                                          zin(155) = cp01*zin(151) + zcp00*zin(154)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  191

                                          xin(191) = xc00*xin(155) + c01*xin(154)
                                          yin(191) = yc00*yin(155) + c01*yin(154)
                                          zin(191) = zc00*zin(155) + c01*zin(154)

                                          ! ------------------

                                          ! i3 = i4 =  154
                                          ! i4 = i5 =  155

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  156
                                          ! i3 =  154
                                          ! i4 =  155

                                          xin(156) = cp01*xin(154) + xcp00*xin(155)
                                          yin(156) = cp01*yin(154) + ycp00*yin(155)
                                          zin(156) = cp01*zin(154) + zcp00*zin(155)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  192

                                          xin(192) = xc00*xin(156) + c01*xin(155)
                                          yin(192) = yc00*yin(156) + c01*yin(155)
                                          zin(192) = zc00*zin(156) + c01*zin(155)

                                          ! ------------------

                                          ! i3 = i4 =  155
                                          ! i4 = i5 =  156

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =  145
                                          ! i4 = i2 =  181

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  217

                                          xin(223) = c10*xin(151) + xc00*xin(187) + c01*xin(184)
                                          yin(223) = c10*yin(151) + yc00*yin(187) + c01*yin(184)
                                          zin(223) = c10*zin(151) + zc00*zin(187) + c01*zin(184)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  181
                                          ! i4 = i5 =  217

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  253

                                          xin(259) = c10*xin(187) + xc00*xin(223) + c01*xin(220)
                                          yin(259) = c10*yin(187) + yc00*yin(223) + c01*yin(220)
                                          zin(259) = c10*zin(187) + zc00*zin(223) + c01*zin(220)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  217
                                          ! i4 = i5 =  253

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  265

                                          xin(271) = c10*xin(223) + xc00*xin(259) + c01*xin(256)
                                          yin(271) = c10*yin(223) + yc00*yin(259) + c01*yin(256)
                                          zin(271) = c10*zin(223) + zc00*zin(259) + c01*zin(256)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  253
                                          ! i4 = i5 =  265

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  277

                                          xin(283) = c10*xin(259) + xc00*xin(271) + c01*xin(268)
                                          yin(283) = c10*yin(259) + yc00*yin(271) + c01*yin(268)
                                          zin(283) = c10*zin(259) + zc00*zin(271) + c01*zin(268)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  265
                                          ! i4 = i5 =  277

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =  145
                                          ! i4 = i2 =  181

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  217

                                          xin(226) = c10*xin(154) + xc00*xin(190) + c01*xin(187)
                                          yin(226) = c10*yin(154) + yc00*yin(190) + c01*yin(187)
                                          zin(226) = c10*zin(154) + zc00*zin(190) + c01*zin(187)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  181
                                          ! i4 = i5 =  217

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  253

                                          xin(262) = c10*xin(190) + xc00*xin(226) + c01*xin(223)
                                          yin(262) = c10*yin(190) + yc00*yin(226) + c01*yin(223)
                                          zin(262) = c10*zin(190) + zc00*zin(226) + c01*zin(223)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  217
                                          ! i4 = i5 =  253

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  265

                                          xin(274) = c10*xin(226) + xc00*xin(262) + c01*xin(259)
                                          yin(274) = c10*yin(226) + yc00*yin(262) + c01*yin(259)
                                          zin(274) = c10*zin(226) + zc00*zin(262) + c01*zin(259)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  253
                                          ! i4 = i5 =  265

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  277

                                          xin(286) = c10*xin(262) + xc00*xin(274) + c01*xin(271)
                                          yin(286) = c10*yin(262) + yc00*yin(274) + c01*yin(271)
                                          zin(286) = c10*zin(262) + zc00*zin(274) + c01*zin(271)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  265
                                          ! i4 = i5 =  277

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =  145
                                          ! i4 = i2 =  181

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  217

                                          xin(227) = c10*xin(155) + xc00*xin(191) + c01*xin(190)
                                          yin(227) = c10*yin(155) + yc00*yin(191) + c01*yin(190)
                                          zin(227) = c10*zin(155) + zc00*zin(191) + c01*zin(190)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  181
                                          ! i4 = i5 =  217

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  253

                                          xin(263) = c10*xin(191) + xc00*xin(227) + c01*xin(226)
                                          yin(263) = c10*yin(191) + yc00*yin(227) + c01*yin(226)
                                          zin(263) = c10*zin(191) + zc00*zin(227) + c01*zin(226)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  217
                                          ! i4 = i5 =  253

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  265

                                          xin(275) = c10*xin(227) + xc00*xin(263) + c01*xin(262)
                                          yin(275) = c10*yin(227) + yc00*yin(263) + c01*yin(262)
                                          zin(275) = c10*zin(227) + zc00*zin(263) + c01*zin(262)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  253
                                          ! i4 = i5 =  265

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  277

                                          xin(287) = c10*xin(263) + xc00*xin(275) + c01*xin(274)
                                          yin(287) = c10*yin(263) + yc00*yin(275) + c01*yin(274)
                                          zin(287) = c10*zin(263) + zc00*zin(275) + c01*zin(274)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  265
                                          ! i4 = i5 =  277

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =  145
                                          ! i4 = i2 =  181

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  217

                                          xin(228) = c10*xin(156) + xc00*xin(192) + c01*xin(191)
                                          yin(228) = c10*yin(156) + yc00*yin(192) + c01*yin(191)
                                          zin(228) = c10*zin(156) + zc00*zin(192) + c01*zin(191)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  181
                                          ! i4 = i5 =  217

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  253

                                          xin(264) = c10*xin(192) + xc00*xin(228) + c01*xin(227)
                                          yin(264) = c10*yin(192) + yc00*yin(228) + c01*yin(227)
                                          zin(264) = c10*zin(192) + zc00*zin(228) + c01*zin(227)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  217
                                          ! i4 = i5 =  253

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  265

                                          xin(276) = c10*xin(228) + xc00*xin(264) + c01*xin(263)
                                          yin(276) = c10*yin(228) + yc00*yin(264) + c01*yin(263)
                                          zin(276) = c10*zin(228) + zc00*zin(264) + c01*zin(263)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  253
                                          ! i4 = i5 =  265

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  277

                                          xin(288) = c10*xin(264) + xc00*xin(276) + c01*xin(275)
                                          yin(288) = c10*yin(264) + yc00*yin(276) + c01*yin(275)
                                          zin(288) = c10*zin(264) + zc00*zin(276) + c01*zin(275)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  265
                                          ! i4 = i5 =  277

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  277

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  277

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  265

                                          xin(277) = xin(277) + dxij*xin(265)
                                          yin(277) = yin(277) + dyij*yin(265)
                                          zin(277) = zin(277) + dzij*zin(265)

                                          ! i3 = i4 =  265
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  253

                                          xin(265) = xin(265) + dxij*xin(253)
                                          yin(265) = yin(265) + dyij*yin(253)
                                          zin(265) = zin(265) + dzij*zin(253)

                                          ! i3 = i4 =  253
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  277

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  265

                                          xin(277) = xin(277) + dxij*xin(265)
                                          yin(277) = yin(277) + dyij*yin(265)
                                          zin(277) = zin(277) + dzij*zin(265)

                                          ! i3 = i4 =  265
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  157

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  157

                                          ! do ni = 1,    3

                                          xin(157) = xin(181) + dxij*xin(145)
                                          yin(157) = yin(181) + dyij*yin(145)
                                          zin(157) = zin(181) + dzij*zin(145)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  193

                                          ! ni =    2

                                          xin(193) = xin(217) + dxij*xin(181)
                                          yin(193) = yin(217) + dyij*yin(181)
                                          zin(193) = zin(217) + dzij*zin(181)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  229

                                          ! ni =    3

                                          xin(229) = xin(253) + dxij*xin(217)
                                          yin(229) = yin(253) + dyij*yin(217)
                                          zin(229) = zin(253) + dzij*zin(217)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  265

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  169

                                          ! nj =    2

                                          ! i4 = i3 =  169

                                          ! do ni = 1,    3

                                          xin(169) = xin(193) + dxij*xin(157)
                                          yin(169) = yin(193) + dyij*yin(157)
                                          zin(169) = zin(193) + dzij*zin(157)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  205

                                          ! ni =    2

                                          xin(205) = xin(229) + dxij*xin(193)
                                          yin(205) = yin(229) + dyij*yin(193)
                                          zin(205) = zin(229) + dzij*zin(193)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  241

                                          ! ni =    3

                                          xin(241) = xin(265) + dxij*xin(229)
                                          yin(241) = yin(265) + dyij*yin(229)
                                          zin(241) = zin(265) + dzij*zin(229)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  277

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  181

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  280

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  268

                                          xin(280) = xin(280) + dxij*xin(268)
                                          yin(280) = yin(280) + dyij*yin(268)
                                          zin(280) = zin(280) + dzij*zin(268)

                                          ! i3 = i4 =  268
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  256

                                          xin(268) = xin(268) + dxij*xin(256)
                                          yin(268) = yin(268) + dyij*yin(256)
                                          zin(268) = zin(268) + dzij*zin(256)

                                          ! i3 = i4 =  256
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  280

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  268

                                          xin(280) = xin(280) + dxij*xin(268)
                                          yin(280) = yin(280) + dyij*yin(268)
                                          zin(280) = zin(280) + dzij*zin(268)

                                          ! i3 = i4 =  268
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  160

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  160

                                          ! do ni = 1,    3

                                          xin(160) = xin(184) + dxij*xin(148)
                                          yin(160) = yin(184) + dyij*yin(148)
                                          zin(160) = zin(184) + dzij*zin(148)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  196

                                          ! ni =    2

                                          xin(196) = xin(220) + dxij*xin(184)
                                          yin(196) = yin(220) + dyij*yin(184)
                                          zin(196) = zin(220) + dzij*zin(184)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  232

                                          ! ni =    3

                                          xin(232) = xin(256) + dxij*xin(220)
                                          yin(232) = yin(256) + dyij*yin(220)
                                          zin(232) = zin(256) + dzij*zin(220)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  268

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  172

                                          ! nj =    2

                                          ! i4 = i3 =  172

                                          ! do ni = 1,    3

                                          xin(172) = xin(196) + dxij*xin(160)
                                          yin(172) = yin(196) + dyij*yin(160)
                                          zin(172) = zin(196) + dzij*zin(160)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  208

                                          ! ni =    2

                                          xin(208) = xin(232) + dxij*xin(196)
                                          yin(208) = yin(232) + dyij*yin(196)
                                          zin(208) = zin(232) + dzij*zin(196)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  244

                                          ! ni =    3

                                          xin(244) = xin(268) + dxij*xin(232)
                                          yin(244) = yin(268) + dyij*yin(232)
                                          zin(244) = zin(268) + dzij*zin(232)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  280

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  184

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  283

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  271

                                          xin(283) = xin(283) + dxij*xin(271)
                                          yin(283) = yin(283) + dyij*yin(271)
                                          zin(283) = zin(283) + dzij*zin(271)

                                          ! i3 = i4 =  271
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  259

                                          xin(271) = xin(271) + dxij*xin(259)
                                          yin(271) = yin(271) + dyij*yin(259)
                                          zin(271) = zin(271) + dzij*zin(259)

                                          ! i3 = i4 =  259
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  283

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  271

                                          xin(283) = xin(283) + dxij*xin(271)
                                          yin(283) = yin(283) + dyij*yin(271)
                                          zin(283) = zin(283) + dzij*zin(271)

                                          ! i3 = i4 =  271
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  163

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  163

                                          ! do ni = 1,    3

                                          xin(163) = xin(187) + dxij*xin(151)
                                          yin(163) = yin(187) + dyij*yin(151)
                                          zin(163) = zin(187) + dzij*zin(151)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  199

                                          ! ni =    2

                                          xin(199) = xin(223) + dxij*xin(187)
                                          yin(199) = yin(223) + dyij*yin(187)
                                          zin(199) = zin(223) + dzij*zin(187)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  235

                                          ! ni =    3

                                          xin(235) = xin(259) + dxij*xin(223)
                                          yin(235) = yin(259) + dyij*yin(223)
                                          zin(235) = zin(259) + dzij*zin(223)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  271

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  175

                                          ! nj =    2

                                          ! i4 = i3 =  175

                                          ! do ni = 1,    3

                                          xin(175) = xin(199) + dxij*xin(163)
                                          yin(175) = yin(199) + dyij*yin(163)
                                          zin(175) = zin(199) + dzij*zin(163)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  211

                                          ! ni =    2

                                          xin(211) = xin(235) + dxij*xin(199)
                                          yin(211) = yin(235) + dyij*yin(199)
                                          zin(211) = zin(235) + dzij*zin(199)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  247

                                          ! ni =    3

                                          xin(247) = xin(271) + dxij*xin(235)
                                          yin(247) = yin(271) + dyij*yin(235)
                                          zin(247) = zin(271) + dzij*zin(235)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  283

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  187

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  286

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  274

                                          xin(286) = xin(286) + dxij*xin(274)
                                          yin(286) = yin(286) + dyij*yin(274)
                                          zin(286) = zin(286) + dzij*zin(274)

                                          ! i3 = i4 =  274
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  262

                                          xin(274) = xin(274) + dxij*xin(262)
                                          yin(274) = yin(274) + dyij*yin(262)
                                          zin(274) = zin(274) + dzij*zin(262)

                                          ! i3 = i4 =  262
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  286

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  274

                                          xin(286) = xin(286) + dxij*xin(274)
                                          yin(286) = yin(286) + dyij*yin(274)
                                          zin(286) = zin(286) + dzij*zin(274)

                                          ! i3 = i4 =  274
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  166

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  166

                                          ! do ni = 1,    3

                                          xin(166) = xin(190) + dxij*xin(154)
                                          yin(166) = yin(190) + dyij*yin(154)
                                          zin(166) = zin(190) + dzij*zin(154)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  202

                                          ! ni =    2

                                          xin(202) = xin(226) + dxij*xin(190)
                                          yin(202) = yin(226) + dyij*yin(190)
                                          zin(202) = zin(226) + dzij*zin(190)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  238

                                          ! ni =    3

                                          xin(238) = xin(262) + dxij*xin(226)
                                          yin(238) = yin(262) + dyij*yin(226)
                                          zin(238) = zin(262) + dzij*zin(226)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  274

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  178

                                          ! nj =    2

                                          ! i4 = i3 =  178

                                          ! do ni = 1,    3

                                          xin(178) = xin(202) + dxij*xin(166)
                                          yin(178) = yin(202) + dyij*yin(166)
                                          zin(178) = zin(202) + dzij*zin(166)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  214

                                          ! ni =    2

                                          xin(214) = xin(238) + dxij*xin(202)
                                          yin(214) = yin(238) + dyij*yin(202)
                                          zin(214) = zin(238) + dzij*zin(202)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  250

                                          ! ni =    3

                                          xin(250) = xin(274) + dxij*xin(238)
                                          yin(250) = yin(274) + dyij*yin(238)
                                          zin(250) = zin(274) + dzij*zin(238)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  286

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  190

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  287

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  275

                                          xin(287) = xin(287) + dxij*xin(275)
                                          yin(287) = yin(287) + dyij*yin(275)
                                          zin(287) = zin(287) + dzij*zin(275)

                                          ! i3 = i4 =  275
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  263

                                          xin(275) = xin(275) + dxij*xin(263)
                                          yin(275) = yin(275) + dyij*yin(263)
                                          zin(275) = zin(275) + dzij*zin(263)

                                          ! i3 = i4 =  263
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  287

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  275

                                          xin(287) = xin(287) + dxij*xin(275)
                                          yin(287) = yin(287) + dyij*yin(275)
                                          zin(287) = zin(287) + dzij*zin(275)

                                          ! i3 = i4 =  275
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  167

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  167

                                          ! do ni = 1,    3

                                          xin(167) = xin(191) + dxij*xin(155)
                                          yin(167) = yin(191) + dyij*yin(155)
                                          zin(167) = zin(191) + dzij*zin(155)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  203

                                          ! ni =    2

                                          xin(203) = xin(227) + dxij*xin(191)
                                          yin(203) = yin(227) + dyij*yin(191)
                                          zin(203) = zin(227) + dzij*zin(191)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  239

                                          ! ni =    3

                                          xin(239) = xin(263) + dxij*xin(227)
                                          yin(239) = yin(263) + dyij*yin(227)
                                          zin(239) = zin(263) + dzij*zin(227)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  275

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  179

                                          ! nj =    2

                                          ! i4 = i3 =  179

                                          ! do ni = 1,    3

                                          xin(179) = xin(203) + dxij*xin(167)
                                          yin(179) = yin(203) + dyij*yin(167)
                                          zin(179) = zin(203) + dzij*zin(167)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  215

                                          ! ni =    2

                                          xin(215) = xin(239) + dxij*xin(203)
                                          yin(215) = yin(239) + dyij*yin(203)
                                          zin(215) = zin(239) + dzij*zin(203)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  251

                                          ! ni =    3

                                          xin(251) = xin(275) + dxij*xin(239)
                                          yin(251) = yin(275) + dyij*yin(239)
                                          zin(251) = zin(275) + dzij*zin(239)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  287

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  191

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  288

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  276

                                          xin(288) = xin(288) + dxij*xin(276)
                                          yin(288) = yin(288) + dyij*yin(276)
                                          zin(288) = zin(288) + dzij*zin(276)

                                          ! i3 = i4 =  276
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  264

                                          xin(276) = xin(276) + dxij*xin(264)
                                          yin(276) = yin(276) + dyij*yin(264)
                                          zin(276) = zin(276) + dzij*zin(264)

                                          ! i3 = i4 =  264
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  288

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  276

                                          xin(288) = xin(288) + dxij*xin(276)
                                          yin(288) = yin(288) + dyij*yin(276)
                                          zin(288) = zin(288) + dzij*zin(276)

                                          ! i3 = i4 =  276
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  168

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  168

                                          ! do ni = 1,    3

                                          xin(168) = xin(192) + dxij*xin(156)
                                          yin(168) = yin(192) + dyij*yin(156)
                                          zin(168) = zin(192) + dzij*zin(156)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  204

                                          ! ni =    2

                                          xin(204) = xin(228) + dxij*xin(192)
                                          yin(204) = yin(228) + dyij*yin(192)
                                          zin(204) = zin(228) + dzij*zin(192)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  240

                                          ! ni =    3

                                          xin(240) = xin(264) + dxij*xin(228)
                                          yin(240) = yin(264) + dyij*yin(228)
                                          zin(240) = zin(264) + dzij*zin(228)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  276

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  180

                                          ! nj =    2

                                          ! i4 = i3 =  180

                                          ! do ni = 1,    3

                                          xin(180) = xin(204) + dxij*xin(168)
                                          yin(180) = yin(204) + dyij*yin(168)
                                          zin(180) = zin(204) + dzij*zin(168)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  216

                                          ! ni =    2

                                          xin(216) = xin(240) + dxij*xin(204)
                                          yin(216) = yin(240) + dyij*yin(204)
                                          zin(216) = zin(240) + dzij*zin(204)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  252

                                          ! ni =    3

                                          xin(252) = xin(276) + dxij*xin(240)
                                          yin(252) = yin(276) + dyij*yin(240)
                                          zin(252) = zin(276) + dzij*zin(240)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  288

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  192

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =  145

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  156

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  155

                                          xin(156) = xin(156) + dxkl*xin(155)
                                          yin(156) = yin(156) + dykl*yin(155)
                                          zin(156) = zin(156) + dzkl*zin(155)

                                          ! i3 = i4 =  155
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  154

                                          xin(155) = xin(155) + dxkl*xin(154)
                                          yin(155) = yin(155) + dykl*yin(154)
                                          zin(155) = zin(155) + dzkl*zin(154)

                                          ! i3 = i4 =  154
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  156

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  155

                                          xin(156) = xin(156) + dxkl*xin(155)
                                          yin(156) = yin(156) + dykl*yin(155)
                                          zin(156) = zin(156) + dzkl*zin(155)

                                          ! i3 = i4 =  155
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  146

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  146

                                          ! do nk = 1,    3

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

                                          xin(152) = xin(154) + dxkl*xin(151)
                                          yin(152) = yin(154) + dykl*yin(151)
                                          zin(152) = zin(154) + dzkl*zin(151)
                                          ! i4 = i4 + lang+1 =  155

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  147

                                          ! nl =    2

                                          ! i4 = i3 =  147

                                          ! do nk = 1,    3

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

                                          xin(153) = xin(155) + dxkl*xin(152)
                                          yin(153) = yin(155) + dykl*yin(152)
                                          zin(153) = zin(155) + dzkl*zin(152)
                                          ! i4 = i4 + lang+1 =  156

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  148

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  157

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  168

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  167

                                          xin(168) = xin(168) + dxkl*xin(167)
                                          yin(168) = yin(168) + dykl*yin(167)
                                          zin(168) = zin(168) + dzkl*zin(167)

                                          ! i3 = i4 =  167
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  166

                                          xin(167) = xin(167) + dxkl*xin(166)
                                          yin(167) = yin(167) + dykl*yin(166)
                                          zin(167) = zin(167) + dzkl*zin(166)

                                          ! i3 = i4 =  166
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  168

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  167

                                          xin(168) = xin(168) + dxkl*xin(167)
                                          yin(168) = yin(168) + dykl*yin(167)
                                          zin(168) = zin(168) + dzkl*zin(167)

                                          ! i3 = i4 =  167
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  158

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  158

                                          ! do nk = 1,    3

                                          xin(158) = xin(160) + dxkl*xin(157)
                                          yin(158) = yin(160) + dykl*yin(157)
                                          zin(158) = zin(160) + dzkl*zin(157)
                                          ! i4 = i4 + lang+1 =  161

                                          ! nk =    2

                                          xin(161) = xin(163) + dxkl*xin(160)
                                          yin(161) = yin(163) + dykl*yin(160)
                                          zin(161) = zin(163) + dzkl*zin(160)
                                          ! i4 = i4 + lang+1 =  164

                                          ! nk =    3

                                          xin(164) = xin(166) + dxkl*xin(163)
                                          yin(164) = yin(166) + dykl*yin(163)
                                          zin(164) = zin(166) + dzkl*zin(163)
                                          ! i4 = i4 + lang+1 =  167

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  159

                                          ! nl =    2

                                          ! i4 = i3 =  159

                                          ! do nk = 1,    3

                                          xin(159) = xin(161) + dxkl*xin(158)
                                          yin(159) = yin(161) + dykl*yin(158)
                                          zin(159) = zin(161) + dzkl*zin(158)
                                          ! i4 = i4 + lang+1 =  162

                                          ! nk =    2

                                          xin(162) = xin(164) + dxkl*xin(161)
                                          yin(162) = yin(164) + dykl*yin(161)
                                          zin(162) = zin(164) + dzkl*zin(161)
                                          ! i4 = i4 + lang+1 =  165

                                          ! nk =    3

                                          xin(165) = xin(167) + dxkl*xin(164)
                                          yin(165) = yin(167) + dykl*yin(164)
                                          zin(165) = zin(167) + dzkl*zin(164)
                                          ! i4 = i4 + lang+1 =  168

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  160

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  169

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  180

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  179

                                          xin(180) = xin(180) + dxkl*xin(179)
                                          yin(180) = yin(180) + dykl*yin(179)
                                          zin(180) = zin(180) + dzkl*zin(179)

                                          ! i3 = i4 =  179
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  178

                                          xin(179) = xin(179) + dxkl*xin(178)
                                          yin(179) = yin(179) + dykl*yin(178)
                                          zin(179) = zin(179) + dzkl*zin(178)

                                          ! i3 = i4 =  178
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  180

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  179

                                          xin(180) = xin(180) + dxkl*xin(179)
                                          yin(180) = yin(180) + dykl*yin(179)
                                          zin(180) = zin(180) + dzkl*zin(179)

                                          ! i3 = i4 =  179
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  170

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  170

                                          ! do nk = 1,    3

                                          xin(170) = xin(172) + dxkl*xin(169)
                                          yin(170) = yin(172) + dykl*yin(169)
                                          zin(170) = zin(172) + dzkl*zin(169)
                                          ! i4 = i4 + lang+1 =  173

                                          ! nk =    2

                                          xin(173) = xin(175) + dxkl*xin(172)
                                          yin(173) = yin(175) + dykl*yin(172)
                                          zin(173) = zin(175) + dzkl*zin(172)
                                          ! i4 = i4 + lang+1 =  176

                                          ! nk =    3

                                          xin(176) = xin(178) + dxkl*xin(175)
                                          yin(176) = yin(178) + dykl*yin(175)
                                          zin(176) = zin(178) + dzkl*zin(175)
                                          ! i4 = i4 + lang+1 =  179

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  171

                                          ! nl =    2

                                          ! i4 = i3 =  171

                                          ! do nk = 1,    3

                                          xin(171) = xin(173) + dxkl*xin(170)
                                          yin(171) = yin(173) + dykl*yin(170)
                                          zin(171) = zin(173) + dzkl*zin(170)
                                          ! i4 = i4 + lang+1 =  174

                                          ! nk =    2

                                          xin(174) = xin(176) + dxkl*xin(173)
                                          yin(174) = yin(176) + dykl*yin(173)
                                          zin(174) = zin(176) + dzkl*zin(173)
                                          ! i4 = i4 + lang+1 =  177

                                          ! nk =    3

                                          xin(177) = xin(179) + dxkl*xin(176)
                                          yin(177) = yin(179) + dykl*yin(176)
                                          zin(177) = zin(179) + dzkl*zin(176)
                                          ! i4 = i4 + lang+1 =  180

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  172

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  181

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  181

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  192

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  191

                                          xin(192) = xin(192) + dxkl*xin(191)
                                          yin(192) = yin(192) + dykl*yin(191)
                                          zin(192) = zin(192) + dzkl*zin(191)

                                          ! i3 = i4 =  191
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  190

                                          xin(191) = xin(191) + dxkl*xin(190)
                                          yin(191) = yin(191) + dykl*yin(190)
                                          zin(191) = zin(191) + dzkl*zin(190)

                                          ! i3 = i4 =  190
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  192

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  191

                                          xin(192) = xin(192) + dxkl*xin(191)
                                          yin(192) = yin(192) + dykl*yin(191)
                                          zin(192) = zin(192) + dzkl*zin(191)

                                          ! i3 = i4 =  191
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  182

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  182

                                          ! do nk = 1,    3

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

                                          xin(188) = xin(190) + dxkl*xin(187)
                                          yin(188) = yin(190) + dykl*yin(187)
                                          zin(188) = zin(190) + dzkl*zin(187)
                                          ! i4 = i4 + lang+1 =  191

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  183

                                          ! nl =    2

                                          ! i4 = i3 =  183

                                          ! do nk = 1,    3

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

                                          xin(189) = xin(191) + dxkl*xin(188)
                                          yin(189) = yin(191) + dykl*yin(188)
                                          zin(189) = zin(191) + dzkl*zin(188)
                                          ! i4 = i4 + lang+1 =  192

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  184

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  193

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  204

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  203

                                          xin(204) = xin(204) + dxkl*xin(203)
                                          yin(204) = yin(204) + dykl*yin(203)
                                          zin(204) = zin(204) + dzkl*zin(203)

                                          ! i3 = i4 =  203
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  202

                                          xin(203) = xin(203) + dxkl*xin(202)
                                          yin(203) = yin(203) + dykl*yin(202)
                                          zin(203) = zin(203) + dzkl*zin(202)

                                          ! i3 = i4 =  202
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  204

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  203

                                          xin(204) = xin(204) + dxkl*xin(203)
                                          yin(204) = yin(204) + dykl*yin(203)
                                          zin(204) = zin(204) + dzkl*zin(203)

                                          ! i3 = i4 =  203
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  194

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  194

                                          ! do nk = 1,    3

                                          xin(194) = xin(196) + dxkl*xin(193)
                                          yin(194) = yin(196) + dykl*yin(193)
                                          zin(194) = zin(196) + dzkl*zin(193)
                                          ! i4 = i4 + lang+1 =  197

                                          ! nk =    2

                                          xin(197) = xin(199) + dxkl*xin(196)
                                          yin(197) = yin(199) + dykl*yin(196)
                                          zin(197) = zin(199) + dzkl*zin(196)
                                          ! i4 = i4 + lang+1 =  200

                                          ! nk =    3

                                          xin(200) = xin(202) + dxkl*xin(199)
                                          yin(200) = yin(202) + dykl*yin(199)
                                          zin(200) = zin(202) + dzkl*zin(199)
                                          ! i4 = i4 + lang+1 =  203

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  195

                                          ! nl =    2

                                          ! i4 = i3 =  195

                                          ! do nk = 1,    3

                                          xin(195) = xin(197) + dxkl*xin(194)
                                          yin(195) = yin(197) + dykl*yin(194)
                                          zin(195) = zin(197) + dzkl*zin(194)
                                          ! i4 = i4 + lang+1 =  198

                                          ! nk =    2

                                          xin(198) = xin(200) + dxkl*xin(197)
                                          yin(198) = yin(200) + dykl*yin(197)
                                          zin(198) = zin(200) + dzkl*zin(197)
                                          ! i4 = i4 + lang+1 =  201

                                          ! nk =    3

                                          xin(201) = xin(203) + dxkl*xin(200)
                                          yin(201) = yin(203) + dykl*yin(200)
                                          zin(201) = zin(203) + dzkl*zin(200)
                                          ! i4 = i4 + lang+1 =  204

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  196

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  205

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  216

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  215

                                          xin(216) = xin(216) + dxkl*xin(215)
                                          yin(216) = yin(216) + dykl*yin(215)
                                          zin(216) = zin(216) + dzkl*zin(215)

                                          ! i3 = i4 =  215
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  214

                                          xin(215) = xin(215) + dxkl*xin(214)
                                          yin(215) = yin(215) + dykl*yin(214)
                                          zin(215) = zin(215) + dzkl*zin(214)

                                          ! i3 = i4 =  214
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  216

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  215

                                          xin(216) = xin(216) + dxkl*xin(215)
                                          yin(216) = yin(216) + dykl*yin(215)
                                          zin(216) = zin(216) + dzkl*zin(215)

                                          ! i3 = i4 =  215
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  206

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  206

                                          ! do nk = 1,    3

                                          xin(206) = xin(208) + dxkl*xin(205)
                                          yin(206) = yin(208) + dykl*yin(205)
                                          zin(206) = zin(208) + dzkl*zin(205)
                                          ! i4 = i4 + lang+1 =  209

                                          ! nk =    2

                                          xin(209) = xin(211) + dxkl*xin(208)
                                          yin(209) = yin(211) + dykl*yin(208)
                                          zin(209) = zin(211) + dzkl*zin(208)
                                          ! i4 = i4 + lang+1 =  212

                                          ! nk =    3

                                          xin(212) = xin(214) + dxkl*xin(211)
                                          yin(212) = yin(214) + dykl*yin(211)
                                          zin(212) = zin(214) + dzkl*zin(211)
                                          ! i4 = i4 + lang+1 =  215

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  207

                                          ! nl =    2

                                          ! i4 = i3 =  207

                                          ! do nk = 1,    3

                                          xin(207) = xin(209) + dxkl*xin(206)
                                          yin(207) = yin(209) + dykl*yin(206)
                                          zin(207) = zin(209) + dzkl*zin(206)
                                          ! i4 = i4 + lang+1 =  210

                                          ! nk =    2

                                          xin(210) = xin(212) + dxkl*xin(209)
                                          yin(210) = yin(212) + dykl*yin(209)
                                          zin(210) = zin(212) + dzkl*zin(209)
                                          ! i4 = i4 + lang+1 =  213

                                          ! nk =    3

                                          xin(213) = xin(215) + dxkl*xin(212)
                                          yin(213) = yin(215) + dykl*yin(212)
                                          zin(213) = zin(215) + dzkl*zin(212)
                                          ! i4 = i4 + lang+1 =  216

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  208

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  217

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  217

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  228

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  227

                                          xin(228) = xin(228) + dxkl*xin(227)
                                          yin(228) = yin(228) + dykl*yin(227)
                                          zin(228) = zin(228) + dzkl*zin(227)

                                          ! i3 = i4 =  227
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  226

                                          xin(227) = xin(227) + dxkl*xin(226)
                                          yin(227) = yin(227) + dykl*yin(226)
                                          zin(227) = zin(227) + dzkl*zin(226)

                                          ! i3 = i4 =  226
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  228

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  227

                                          xin(228) = xin(228) + dxkl*xin(227)
                                          yin(228) = yin(228) + dykl*yin(227)
                                          zin(228) = zin(228) + dzkl*zin(227)

                                          ! i3 = i4 =  227
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  218

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  218

                                          ! do nk = 1,    3

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

                                          xin(224) = xin(226) + dxkl*xin(223)
                                          yin(224) = yin(226) + dykl*yin(223)
                                          zin(224) = zin(226) + dzkl*zin(223)
                                          ! i4 = i4 + lang+1 =  227

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  219

                                          ! nl =    2

                                          ! i4 = i3 =  219

                                          ! do nk = 1,    3

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

                                          xin(225) = xin(227) + dxkl*xin(224)
                                          yin(225) = yin(227) + dykl*yin(224)
                                          zin(225) = zin(227) + dzkl*zin(224)
                                          ! i4 = i4 + lang+1 =  228

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  220

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  229

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  240

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  239

                                          xin(240) = xin(240) + dxkl*xin(239)
                                          yin(240) = yin(240) + dykl*yin(239)
                                          zin(240) = zin(240) + dzkl*zin(239)

                                          ! i3 = i4 =  239
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  238

                                          xin(239) = xin(239) + dxkl*xin(238)
                                          yin(239) = yin(239) + dykl*yin(238)
                                          zin(239) = zin(239) + dzkl*zin(238)

                                          ! i3 = i4 =  238
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  240

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  239

                                          xin(240) = xin(240) + dxkl*xin(239)
                                          yin(240) = yin(240) + dykl*yin(239)
                                          zin(240) = zin(240) + dzkl*zin(239)

                                          ! i3 = i4 =  239
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  230

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  230

                                          ! do nk = 1,    3

                                          xin(230) = xin(232) + dxkl*xin(229)
                                          yin(230) = yin(232) + dykl*yin(229)
                                          zin(230) = zin(232) + dzkl*zin(229)
                                          ! i4 = i4 + lang+1 =  233

                                          ! nk =    2

                                          xin(233) = xin(235) + dxkl*xin(232)
                                          yin(233) = yin(235) + dykl*yin(232)
                                          zin(233) = zin(235) + dzkl*zin(232)
                                          ! i4 = i4 + lang+1 =  236

                                          ! nk =    3

                                          xin(236) = xin(238) + dxkl*xin(235)
                                          yin(236) = yin(238) + dykl*yin(235)
                                          zin(236) = zin(238) + dzkl*zin(235)
                                          ! i4 = i4 + lang+1 =  239

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  231

                                          ! nl =    2

                                          ! i4 = i3 =  231

                                          ! do nk = 1,    3

                                          xin(231) = xin(233) + dxkl*xin(230)
                                          yin(231) = yin(233) + dykl*yin(230)
                                          zin(231) = zin(233) + dzkl*zin(230)
                                          ! i4 = i4 + lang+1 =  234

                                          ! nk =    2

                                          xin(234) = xin(236) + dxkl*xin(233)
                                          yin(234) = yin(236) + dykl*yin(233)
                                          zin(234) = zin(236) + dzkl*zin(233)
                                          ! i4 = i4 + lang+1 =  237

                                          ! nk =    3

                                          xin(237) = xin(239) + dxkl*xin(236)
                                          yin(237) = yin(239) + dykl*yin(236)
                                          zin(237) = zin(239) + dzkl*zin(236)
                                          ! i4 = i4 + lang+1 =  240

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  232

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  241

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  252

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  251

                                          xin(252) = xin(252) + dxkl*xin(251)
                                          yin(252) = yin(252) + dykl*yin(251)
                                          zin(252) = zin(252) + dzkl*zin(251)

                                          ! i3 = i4 =  251
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  250

                                          xin(251) = xin(251) + dxkl*xin(250)
                                          yin(251) = yin(251) + dykl*yin(250)
                                          zin(251) = zin(251) + dzkl*zin(250)

                                          ! i3 = i4 =  250
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  252

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  251

                                          xin(252) = xin(252) + dxkl*xin(251)
                                          yin(252) = yin(252) + dykl*yin(251)
                                          zin(252) = zin(252) + dzkl*zin(251)

                                          ! i3 = i4 =  251
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  242

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  242

                                          ! do nk = 1,    3

                                          xin(242) = xin(244) + dxkl*xin(241)
                                          yin(242) = yin(244) + dykl*yin(241)
                                          zin(242) = zin(244) + dzkl*zin(241)
                                          ! i4 = i4 + lang+1 =  245

                                          ! nk =    2

                                          xin(245) = xin(247) + dxkl*xin(244)
                                          yin(245) = yin(247) + dykl*yin(244)
                                          zin(245) = zin(247) + dzkl*zin(244)
                                          ! i4 = i4 + lang+1 =  248

                                          ! nk =    3

                                          xin(248) = xin(250) + dxkl*xin(247)
                                          yin(248) = yin(250) + dykl*yin(247)
                                          zin(248) = zin(250) + dzkl*zin(247)
                                          ! i4 = i4 + lang+1 =  251

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  243

                                          ! nl =    2

                                          ! i4 = i3 =  243

                                          ! do nk = 1,    3

                                          xin(243) = xin(245) + dxkl*xin(242)
                                          yin(243) = yin(245) + dykl*yin(242)
                                          zin(243) = zin(245) + dzkl*zin(242)
                                          ! i4 = i4 + lang+1 =  246

                                          ! nk =    2

                                          xin(246) = xin(248) + dxkl*xin(245)
                                          yin(246) = yin(248) + dykl*yin(245)
                                          zin(246) = zin(248) + dzkl*zin(245)
                                          ! i4 = i4 + lang+1 =  249

                                          ! nk =    3

                                          xin(249) = xin(251) + dxkl*xin(248)
                                          yin(249) = yin(251) + dykl*yin(248)
                                          zin(249) = zin(251) + dzkl*zin(248)
                                          ! i4 = i4 + lang+1 =  252

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  244

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  253

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  253

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  264

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  263

                                          xin(264) = xin(264) + dxkl*xin(263)
                                          yin(264) = yin(264) + dykl*yin(263)
                                          zin(264) = zin(264) + dzkl*zin(263)

                                          ! i3 = i4 =  263
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  262

                                          xin(263) = xin(263) + dxkl*xin(262)
                                          yin(263) = yin(263) + dykl*yin(262)
                                          zin(263) = zin(263) + dzkl*zin(262)

                                          ! i3 = i4 =  262
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  264

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  263

                                          xin(264) = xin(264) + dxkl*xin(263)
                                          yin(264) = yin(264) + dykl*yin(263)
                                          zin(264) = zin(264) + dzkl*zin(263)

                                          ! i3 = i4 =  263
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  254

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  254

                                          ! do nk = 1,    3

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

                                          xin(260) = xin(262) + dxkl*xin(259)
                                          yin(260) = yin(262) + dykl*yin(259)
                                          zin(260) = zin(262) + dzkl*zin(259)
                                          ! i4 = i4 + lang+1 =  263

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  255

                                          ! nl =    2

                                          ! i4 = i3 =  255

                                          ! do nk = 1,    3

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

                                          xin(261) = xin(263) + dxkl*xin(260)
                                          yin(261) = yin(263) + dykl*yin(260)
                                          zin(261) = zin(263) + dzkl*zin(260)
                                          ! i4 = i4 + lang+1 =  264

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  256

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  265

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  276

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  275

                                          xin(276) = xin(276) + dxkl*xin(275)
                                          yin(276) = yin(276) + dykl*yin(275)
                                          zin(276) = zin(276) + dzkl*zin(275)

                                          ! i3 = i4 =  275
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  274

                                          xin(275) = xin(275) + dxkl*xin(274)
                                          yin(275) = yin(275) + dykl*yin(274)
                                          zin(275) = zin(275) + dzkl*zin(274)

                                          ! i3 = i4 =  274
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  276

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  275

                                          xin(276) = xin(276) + dxkl*xin(275)
                                          yin(276) = yin(276) + dykl*yin(275)
                                          zin(276) = zin(276) + dzkl*zin(275)

                                          ! i3 = i4 =  275
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  266

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  266

                                          ! do nk = 1,    3

                                          xin(266) = xin(268) + dxkl*xin(265)
                                          yin(266) = yin(268) + dykl*yin(265)
                                          zin(266) = zin(268) + dzkl*zin(265)
                                          ! i4 = i4 + lang+1 =  269

                                          ! nk =    2

                                          xin(269) = xin(271) + dxkl*xin(268)
                                          yin(269) = yin(271) + dykl*yin(268)
                                          zin(269) = zin(271) + dzkl*zin(268)
                                          ! i4 = i4 + lang+1 =  272

                                          ! nk =    3

                                          xin(272) = xin(274) + dxkl*xin(271)
                                          yin(272) = yin(274) + dykl*yin(271)
                                          zin(272) = zin(274) + dzkl*zin(271)
                                          ! i4 = i4 + lang+1 =  275

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  267

                                          ! nl =    2

                                          ! i4 = i3 =  267

                                          ! do nk = 1,    3

                                          xin(267) = xin(269) + dxkl*xin(266)
                                          yin(267) = yin(269) + dykl*yin(266)
                                          zin(267) = zin(269) + dzkl*zin(266)
                                          ! i4 = i4 + lang+1 =  270

                                          ! nk =    2

                                          xin(270) = xin(272) + dxkl*xin(269)
                                          yin(270) = yin(272) + dykl*yin(269)
                                          zin(270) = zin(272) + dzkl*zin(269)
                                          ! i4 = i4 + lang+1 =  273

                                          ! nk =    3

                                          xin(273) = xin(275) + dxkl*xin(272)
                                          yin(273) = yin(275) + dykl*yin(272)
                                          zin(273) = zin(275) + dzkl*zin(272)
                                          ! i4 = i4 + lang+1 =  276

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  268

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  277

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  288

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  287

                                          xin(288) = xin(288) + dxkl*xin(287)
                                          yin(288) = yin(288) + dykl*yin(287)
                                          zin(288) = zin(288) + dzkl*zin(287)

                                          ! i3 = i4 =  287
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  286

                                          xin(287) = xin(287) + dxkl*xin(286)
                                          yin(287) = yin(287) + dykl*yin(286)
                                          zin(287) = zin(287) + dzkl*zin(286)

                                          ! i3 = i4 =  286
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  288

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  287

                                          xin(288) = xin(288) + dxkl*xin(287)
                                          yin(288) = yin(288) + dykl*yin(287)
                                          zin(288) = zin(288) + dzkl*zin(287)

                                          ! i3 = i4 =  287
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  278

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  278

                                          ! do nk = 1,    3

                                          xin(278) = xin(280) + dxkl*xin(277)
                                          yin(278) = yin(280) + dykl*yin(277)
                                          zin(278) = zin(280) + dzkl*zin(277)
                                          ! i4 = i4 + lang+1 =  281

                                          ! nk =    2

                                          xin(281) = xin(283) + dxkl*xin(280)
                                          yin(281) = yin(283) + dykl*yin(280)
                                          zin(281) = zin(283) + dzkl*zin(280)
                                          ! i4 = i4 + lang+1 =  284

                                          ! nk =    3

                                          xin(284) = xin(286) + dxkl*xin(283)
                                          yin(284) = yin(286) + dykl*yin(283)
                                          zin(284) = zin(286) + dzkl*zin(283)
                                          ! i4 = i4 + lang+1 =  287

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  279

                                          ! nl =    2

                                          ! i4 = i3 =  279

                                          ! do nk = 1,    3

                                          xin(279) = xin(281) + dxkl*xin(278)
                                          yin(279) = yin(281) + dykl*yin(278)
                                          zin(279) = zin(281) + dzkl*zin(278)
                                          ! i4 = i4 + lang+1 =  282

                                          ! nk =    2

                                          xin(282) = xin(284) + dxkl*xin(281)
                                          yin(282) = yin(284) + dykl*yin(281)
                                          zin(282) = zin(284) + dzkl*zin(281)
                                          ! i4 = i4 + lang+1 =  285

                                          ! nk =    3

                                          xin(285) = xin(287) + dxkl*xin(284)
                                          yin(285) = yin(287) + dykl*yin(284)
                                          zin(285) = zin(287) + dzkl*zin(284)
                                          ! i4 = i4 + lang+1 =  288

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  280

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  289

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  289

                                          ! end do

                                          ! *** Now root =    3

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  288

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

                                          ! i1 = in(1) =  289

                                          xin(289) = 1.0_dp
                                          yin(289) = 1.0_dp
                                          zin(289) = f00

                                          ! i2 = in(2) =  325
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(325) = xc00
                                          yin(325) = yc00
                                          zin(325) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  292

                                          xin(292) = xcp00
                                          yin(292) = ycp00
                                          zin(292) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  328
                                          ! i2 =  325

                                          xin(328) = xcp00*xin(325) + cp10
                                          yin(328) = ycp00*yin(325) + cp10
                                          zin(328) = zcp00*zin(325) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  289
                                          ! i4 = i2 =  325

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  361
                                          ! i3 =  289
                                          ! i4 =  325

                                          xin(361) = c10*xin(289) + xc00*xin(325)
                                          yin(361) = c10*yin(289) + yc00*yin(325)
                                          zin(361) = c10*zin(289) + zc00*zin(325)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  364
                                          ! i5 =  361
                                          ! i4 =  325

                                          xin(364) = xcp00*xin(361) + cp10*xin(325)
                                          yin(364) = ycp00*yin(361) + cp10*yin(325)
                                          zin(364) = zcp00*zin(361) + cp10*zin(325)

                                          ! ------------------

                                          ! i3 = i4 =  325
                                          ! i4 = i5 =  361

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  397
                                          ! i3 =  325
                                          ! i4 =  361

                                          xin(397) = c10*xin(325) + xc00*xin(361)
                                          yin(397) = c10*yin(325) + yc00*yin(361)
                                          zin(397) = c10*zin(325) + zc00*zin(361)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  400
                                          ! i5 =  397
                                          ! i4 =  361

                                          xin(400) = xcp00*xin(397) + cp10*xin(361)
                                          yin(400) = ycp00*yin(397) + cp10*yin(361)
                                          zin(400) = zcp00*zin(397) + cp10*zin(361)

                                          ! ------------------

                                          ! i3 = i4 =  361
                                          ! i4 = i5 =  397

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  409
                                          ! i3 =  361
                                          ! i4 =  397

                                          xin(409) = c10*xin(361) + xc00*xin(397)
                                          yin(409) = c10*yin(361) + yc00*yin(397)
                                          zin(409) = c10*zin(361) + zc00*zin(397)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  412
                                          ! i5 =  409
                                          ! i4 =  397

                                          xin(412) = xcp00*xin(409) + cp10*xin(397)
                                          yin(412) = ycp00*yin(409) + cp10*yin(397)
                                          zin(412) = zcp00*zin(409) + cp10*zin(397)

                                          ! ------------------

                                          ! i3 = i4 =  397
                                          ! i4 = i5 =  409

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  421
                                          ! i3 =  397
                                          ! i4 =  409

                                          xin(421) = c10*xin(397) + xc00*xin(409)
                                          yin(421) = c10*yin(397) + yc00*yin(409)
                                          zin(421) = c10*zin(397) + zc00*zin(409)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  424
                                          ! i5 =  421
                                          ! i4 =  409

                                          xin(424) = xcp00*xin(421) + cp10*xin(409)
                                          yin(424) = ycp00*yin(421) + cp10*yin(409)
                                          zin(424) = zcp00*zin(421) + cp10*zin(409)

                                          ! ------------------

                                          ! i3 = i4 =  409
                                          ! i4 = i5 =  421

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =  289
                                          ! i4 = i1+k2 =  292

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  295
                                          ! i3 =  289
                                          ! i4 =  292

                                          xin(295) = cp01*xin(289) + xcp00*xin(292)
                                          yin(295) = cp01*yin(289) + ycp00*yin(292)
                                          zin(295) = cp01*zin(289) + zcp00*zin(292)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  331

                                          xin(331) = xc00*xin(295) + c01*xin(292)
                                          yin(331) = yc00*yin(295) + c01*yin(292)
                                          zin(331) = zc00*zin(295) + c01*zin(292)

                                          ! ------------------

                                          ! i3 = i4 =  292
                                          ! i4 = i5 =  295

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  298
                                          ! i3 =  292
                                          ! i4 =  295

                                          xin(298) = cp01*xin(292) + xcp00*xin(295)
                                          yin(298) = cp01*yin(292) + ycp00*yin(295)
                                          zin(298) = cp01*zin(292) + zcp00*zin(295)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  334

                                          xin(334) = xc00*xin(298) + c01*xin(295)
                                          yin(334) = yc00*yin(298) + c01*yin(295)
                                          zin(334) = zc00*zin(298) + c01*zin(295)

                                          ! ------------------

                                          ! i3 = i4 =  295
                                          ! i4 = i5 =  298

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  299
                                          ! i3 =  295
                                          ! i4 =  298

                                          xin(299) = cp01*xin(295) + xcp00*xin(298)
                                          yin(299) = cp01*yin(295) + ycp00*yin(298)
                                          zin(299) = cp01*zin(295) + zcp00*zin(298)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  335

                                          xin(335) = xc00*xin(299) + c01*xin(298)
                                          yin(335) = yc00*yin(299) + c01*yin(298)
                                          zin(335) = zc00*zin(299) + c01*zin(298)

                                          ! ------------------

                                          ! i3 = i4 =  298
                                          ! i4 = i5 =  299

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  300
                                          ! i3 =  298
                                          ! i4 =  299

                                          xin(300) = cp01*xin(298) + xcp00*xin(299)
                                          yin(300) = cp01*yin(298) + ycp00*yin(299)
                                          zin(300) = cp01*zin(298) + zcp00*zin(299)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  336

                                          xin(336) = xc00*xin(300) + c01*xin(299)
                                          yin(336) = yc00*yin(300) + c01*yin(299)
                                          zin(336) = zc00*zin(300) + c01*zin(299)

                                          ! ------------------

                                          ! i3 = i4 =  299
                                          ! i4 = i5 =  300

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =  289
                                          ! i4 = i2 =  325

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  361

                                          xin(367) = c10*xin(295) + xc00*xin(331) + c01*xin(328)
                                          yin(367) = c10*yin(295) + yc00*yin(331) + c01*yin(328)
                                          zin(367) = c10*zin(295) + zc00*zin(331) + c01*zin(328)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  325
                                          ! i4 = i5 =  361

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  397

                                          xin(403) = c10*xin(331) + xc00*xin(367) + c01*xin(364)
                                          yin(403) = c10*yin(331) + yc00*yin(367) + c01*yin(364)
                                          zin(403) = c10*zin(331) + zc00*zin(367) + c01*zin(364)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  361
                                          ! i4 = i5 =  397

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  409

                                          xin(415) = c10*xin(367) + xc00*xin(403) + c01*xin(400)
                                          yin(415) = c10*yin(367) + yc00*yin(403) + c01*yin(400)
                                          zin(415) = c10*zin(367) + zc00*zin(403) + c01*zin(400)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  397
                                          ! i4 = i5 =  409

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  421

                                          xin(427) = c10*xin(403) + xc00*xin(415) + c01*xin(412)
                                          yin(427) = c10*yin(403) + yc00*yin(415) + c01*yin(412)
                                          zin(427) = c10*zin(403) + zc00*zin(415) + c01*zin(412)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  409
                                          ! i4 = i5 =  421

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =  289
                                          ! i4 = i2 =  325

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  361

                                          xin(370) = c10*xin(298) + xc00*xin(334) + c01*xin(331)
                                          yin(370) = c10*yin(298) + yc00*yin(334) + c01*yin(331)
                                          zin(370) = c10*zin(298) + zc00*zin(334) + c01*zin(331)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  325
                                          ! i4 = i5 =  361

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  397

                                          xin(406) = c10*xin(334) + xc00*xin(370) + c01*xin(367)
                                          yin(406) = c10*yin(334) + yc00*yin(370) + c01*yin(367)
                                          zin(406) = c10*zin(334) + zc00*zin(370) + c01*zin(367)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  361
                                          ! i4 = i5 =  397

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  409

                                          xin(418) = c10*xin(370) + xc00*xin(406) + c01*xin(403)
                                          yin(418) = c10*yin(370) + yc00*yin(406) + c01*yin(403)
                                          zin(418) = c10*zin(370) + zc00*zin(406) + c01*zin(403)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  397
                                          ! i4 = i5 =  409

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  421

                                          xin(430) = c10*xin(406) + xc00*xin(418) + c01*xin(415)
                                          yin(430) = c10*yin(406) + yc00*yin(418) + c01*yin(415)
                                          zin(430) = c10*zin(406) + zc00*zin(418) + c01*zin(415)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  409
                                          ! i4 = i5 =  421

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =  289
                                          ! i4 = i2 =  325

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  361

                                          xin(371) = c10*xin(299) + xc00*xin(335) + c01*xin(334)
                                          yin(371) = c10*yin(299) + yc00*yin(335) + c01*yin(334)
                                          zin(371) = c10*zin(299) + zc00*zin(335) + c01*zin(334)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  325
                                          ! i4 = i5 =  361

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  397

                                          xin(407) = c10*xin(335) + xc00*xin(371) + c01*xin(370)
                                          yin(407) = c10*yin(335) + yc00*yin(371) + c01*yin(370)
                                          zin(407) = c10*zin(335) + zc00*zin(371) + c01*zin(370)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  361
                                          ! i4 = i5 =  397

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  409

                                          xin(419) = c10*xin(371) + xc00*xin(407) + c01*xin(406)
                                          yin(419) = c10*yin(371) + yc00*yin(407) + c01*yin(406)
                                          zin(419) = c10*zin(371) + zc00*zin(407) + c01*zin(406)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  397
                                          ! i4 = i5 =  409

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  421

                                          xin(431) = c10*xin(407) + xc00*xin(419) + c01*xin(418)
                                          yin(431) = c10*yin(407) + yc00*yin(419) + c01*yin(418)
                                          zin(431) = c10*zin(407) + zc00*zin(419) + c01*zin(418)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  409
                                          ! i4 = i5 =  421

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =  289
                                          ! i4 = i2 =  325

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  361

                                          xin(372) = c10*xin(300) + xc00*xin(336) + c01*xin(335)
                                          yin(372) = c10*yin(300) + yc00*yin(336) + c01*yin(335)
                                          zin(372) = c10*zin(300) + zc00*zin(336) + c01*zin(335)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  325
                                          ! i4 = i5 =  361

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  397

                                          xin(408) = c10*xin(336) + xc00*xin(372) + c01*xin(371)
                                          yin(408) = c10*yin(336) + yc00*yin(372) + c01*yin(371)
                                          zin(408) = c10*zin(336) + zc00*zin(372) + c01*zin(371)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  361
                                          ! i4 = i5 =  397

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  409

                                          xin(420) = c10*xin(372) + xc00*xin(408) + c01*xin(407)
                                          yin(420) = c10*yin(372) + yc00*yin(408) + c01*yin(407)
                                          zin(420) = c10*zin(372) + zc00*zin(408) + c01*zin(407)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  397
                                          ! i4 = i5 =  409

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  421

                                          xin(432) = c10*xin(408) + xc00*xin(420) + c01*xin(419)
                                          yin(432) = c10*yin(408) + yc00*yin(420) + c01*yin(419)
                                          zin(432) = c10*zin(408) + zc00*zin(420) + c01*zin(419)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  409
                                          ! i4 = i5 =  421

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  421

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  421

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  409

                                          xin(421) = xin(421) + dxij*xin(409)
                                          yin(421) = yin(421) + dyij*yin(409)
                                          zin(421) = zin(421) + dzij*zin(409)

                                          ! i3 = i4 =  409
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  397

                                          xin(409) = xin(409) + dxij*xin(397)
                                          yin(409) = yin(409) + dyij*yin(397)
                                          zin(409) = zin(409) + dzij*zin(397)

                                          ! i3 = i4 =  397
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  421

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  409

                                          xin(421) = xin(421) + dxij*xin(409)
                                          yin(421) = yin(421) + dyij*yin(409)
                                          zin(421) = zin(421) + dzij*zin(409)

                                          ! i3 = i4 =  409
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  301

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  301

                                          ! do ni = 1,    3

                                          xin(301) = xin(325) + dxij*xin(289)
                                          yin(301) = yin(325) + dyij*yin(289)
                                          zin(301) = zin(325) + dzij*zin(289)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  337

                                          ! ni =    2

                                          xin(337) = xin(361) + dxij*xin(325)
                                          yin(337) = yin(361) + dyij*yin(325)
                                          zin(337) = zin(361) + dzij*zin(325)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  373

                                          ! ni =    3

                                          xin(373) = xin(397) + dxij*xin(361)
                                          yin(373) = yin(397) + dyij*yin(361)
                                          zin(373) = zin(397) + dzij*zin(361)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  409

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  313

                                          ! nj =    2

                                          ! i4 = i3 =  313

                                          ! do ni = 1,    3

                                          xin(313) = xin(337) + dxij*xin(301)
                                          yin(313) = yin(337) + dyij*yin(301)
                                          zin(313) = zin(337) + dzij*zin(301)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  349

                                          ! ni =    2

                                          xin(349) = xin(373) + dxij*xin(337)
                                          yin(349) = yin(373) + dyij*yin(337)
                                          zin(349) = zin(373) + dzij*zin(337)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  385

                                          ! ni =    3

                                          xin(385) = xin(409) + dxij*xin(373)
                                          yin(385) = yin(409) + dyij*yin(373)
                                          zin(385) = zin(409) + dzij*zin(373)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  421

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  325

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  424

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  412

                                          xin(424) = xin(424) + dxij*xin(412)
                                          yin(424) = yin(424) + dyij*yin(412)
                                          zin(424) = zin(424) + dzij*zin(412)

                                          ! i3 = i4 =  412
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  400

                                          xin(412) = xin(412) + dxij*xin(400)
                                          yin(412) = yin(412) + dyij*yin(400)
                                          zin(412) = zin(412) + dzij*zin(400)

                                          ! i3 = i4 =  400
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  424

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  412

                                          xin(424) = xin(424) + dxij*xin(412)
                                          yin(424) = yin(424) + dyij*yin(412)
                                          zin(424) = zin(424) + dzij*zin(412)

                                          ! i3 = i4 =  412
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  304

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  304

                                          ! do ni = 1,    3

                                          xin(304) = xin(328) + dxij*xin(292)
                                          yin(304) = yin(328) + dyij*yin(292)
                                          zin(304) = zin(328) + dzij*zin(292)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  340

                                          ! ni =    2

                                          xin(340) = xin(364) + dxij*xin(328)
                                          yin(340) = yin(364) + dyij*yin(328)
                                          zin(340) = zin(364) + dzij*zin(328)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  376

                                          ! ni =    3

                                          xin(376) = xin(400) + dxij*xin(364)
                                          yin(376) = yin(400) + dyij*yin(364)
                                          zin(376) = zin(400) + dzij*zin(364)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  412

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  316

                                          ! nj =    2

                                          ! i4 = i3 =  316

                                          ! do ni = 1,    3

                                          xin(316) = xin(340) + dxij*xin(304)
                                          yin(316) = yin(340) + dyij*yin(304)
                                          zin(316) = zin(340) + dzij*zin(304)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  352

                                          ! ni =    2

                                          xin(352) = xin(376) + dxij*xin(340)
                                          yin(352) = yin(376) + dyij*yin(340)
                                          zin(352) = zin(376) + dzij*zin(340)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  388

                                          ! ni =    3

                                          xin(388) = xin(412) + dxij*xin(376)
                                          yin(388) = yin(412) + dyij*yin(376)
                                          zin(388) = zin(412) + dzij*zin(376)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  424

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  328

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  427

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  415

                                          xin(427) = xin(427) + dxij*xin(415)
                                          yin(427) = yin(427) + dyij*yin(415)
                                          zin(427) = zin(427) + dzij*zin(415)

                                          ! i3 = i4 =  415
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  403

                                          xin(415) = xin(415) + dxij*xin(403)
                                          yin(415) = yin(415) + dyij*yin(403)
                                          zin(415) = zin(415) + dzij*zin(403)

                                          ! i3 = i4 =  403
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  427

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  415

                                          xin(427) = xin(427) + dxij*xin(415)
                                          yin(427) = yin(427) + dyij*yin(415)
                                          zin(427) = zin(427) + dzij*zin(415)

                                          ! i3 = i4 =  415
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  307

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  307

                                          ! do ni = 1,    3

                                          xin(307) = xin(331) + dxij*xin(295)
                                          yin(307) = yin(331) + dyij*yin(295)
                                          zin(307) = zin(331) + dzij*zin(295)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  343

                                          ! ni =    2

                                          xin(343) = xin(367) + dxij*xin(331)
                                          yin(343) = yin(367) + dyij*yin(331)
                                          zin(343) = zin(367) + dzij*zin(331)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  379

                                          ! ni =    3

                                          xin(379) = xin(403) + dxij*xin(367)
                                          yin(379) = yin(403) + dyij*yin(367)
                                          zin(379) = zin(403) + dzij*zin(367)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  415

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  319

                                          ! nj =    2

                                          ! i4 = i3 =  319

                                          ! do ni = 1,    3

                                          xin(319) = xin(343) + dxij*xin(307)
                                          yin(319) = yin(343) + dyij*yin(307)
                                          zin(319) = zin(343) + dzij*zin(307)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  355

                                          ! ni =    2

                                          xin(355) = xin(379) + dxij*xin(343)
                                          yin(355) = yin(379) + dyij*yin(343)
                                          zin(355) = zin(379) + dzij*zin(343)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  391

                                          ! ni =    3

                                          xin(391) = xin(415) + dxij*xin(379)
                                          yin(391) = yin(415) + dyij*yin(379)
                                          zin(391) = zin(415) + dzij*zin(379)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  427

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  331

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  430

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  418

                                          xin(430) = xin(430) + dxij*xin(418)
                                          yin(430) = yin(430) + dyij*yin(418)
                                          zin(430) = zin(430) + dzij*zin(418)

                                          ! i3 = i4 =  418
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  406

                                          xin(418) = xin(418) + dxij*xin(406)
                                          yin(418) = yin(418) + dyij*yin(406)
                                          zin(418) = zin(418) + dzij*zin(406)

                                          ! i3 = i4 =  406
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  430

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  418

                                          xin(430) = xin(430) + dxij*xin(418)
                                          yin(430) = yin(430) + dyij*yin(418)
                                          zin(430) = zin(430) + dzij*zin(418)

                                          ! i3 = i4 =  418
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  310

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  310

                                          ! do ni = 1,    3

                                          xin(310) = xin(334) + dxij*xin(298)
                                          yin(310) = yin(334) + dyij*yin(298)
                                          zin(310) = zin(334) + dzij*zin(298)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  346

                                          ! ni =    2

                                          xin(346) = xin(370) + dxij*xin(334)
                                          yin(346) = yin(370) + dyij*yin(334)
                                          zin(346) = zin(370) + dzij*zin(334)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  382

                                          ! ni =    3

                                          xin(382) = xin(406) + dxij*xin(370)
                                          yin(382) = yin(406) + dyij*yin(370)
                                          zin(382) = zin(406) + dzij*zin(370)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  418

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  322

                                          ! nj =    2

                                          ! i4 = i3 =  322

                                          ! do ni = 1,    3

                                          xin(322) = xin(346) + dxij*xin(310)
                                          yin(322) = yin(346) + dyij*yin(310)
                                          zin(322) = zin(346) + dzij*zin(310)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  358

                                          ! ni =    2

                                          xin(358) = xin(382) + dxij*xin(346)
                                          yin(358) = yin(382) + dyij*yin(346)
                                          zin(358) = zin(382) + dzij*zin(346)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  394

                                          ! ni =    3

                                          xin(394) = xin(418) + dxij*xin(382)
                                          yin(394) = yin(418) + dyij*yin(382)
                                          zin(394) = zin(418) + dzij*zin(382)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  430

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  334

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  431

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  419

                                          xin(431) = xin(431) + dxij*xin(419)
                                          yin(431) = yin(431) + dyij*yin(419)
                                          zin(431) = zin(431) + dzij*zin(419)

                                          ! i3 = i4 =  419
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  407

                                          xin(419) = xin(419) + dxij*xin(407)
                                          yin(419) = yin(419) + dyij*yin(407)
                                          zin(419) = zin(419) + dzij*zin(407)

                                          ! i3 = i4 =  407
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  431

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  419

                                          xin(431) = xin(431) + dxij*xin(419)
                                          yin(431) = yin(431) + dyij*yin(419)
                                          zin(431) = zin(431) + dzij*zin(419)

                                          ! i3 = i4 =  419
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  311

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  311

                                          ! do ni = 1,    3

                                          xin(311) = xin(335) + dxij*xin(299)
                                          yin(311) = yin(335) + dyij*yin(299)
                                          zin(311) = zin(335) + dzij*zin(299)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  347

                                          ! ni =    2

                                          xin(347) = xin(371) + dxij*xin(335)
                                          yin(347) = yin(371) + dyij*yin(335)
                                          zin(347) = zin(371) + dzij*zin(335)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  383

                                          ! ni =    3

                                          xin(383) = xin(407) + dxij*xin(371)
                                          yin(383) = yin(407) + dyij*yin(371)
                                          zin(383) = zin(407) + dzij*zin(371)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  419

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  323

                                          ! nj =    2

                                          ! i4 = i3 =  323

                                          ! do ni = 1,    3

                                          xin(323) = xin(347) + dxij*xin(311)
                                          yin(323) = yin(347) + dyij*yin(311)
                                          zin(323) = zin(347) + dzij*zin(311)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  359

                                          ! ni =    2

                                          xin(359) = xin(383) + dxij*xin(347)
                                          yin(359) = yin(383) + dyij*yin(347)
                                          zin(359) = zin(383) + dzij*zin(347)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  395

                                          ! ni =    3

                                          xin(395) = xin(419) + dxij*xin(383)
                                          yin(395) = yin(419) + dyij*yin(383)
                                          zin(395) = zin(419) + dzij*zin(383)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  431

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  335

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  432

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  420

                                          xin(432) = xin(432) + dxij*xin(420)
                                          yin(432) = yin(432) + dyij*yin(420)
                                          zin(432) = zin(432) + dzij*zin(420)

                                          ! i3 = i4 =  420
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  408

                                          xin(420) = xin(420) + dxij*xin(408)
                                          yin(420) = yin(420) + dyij*yin(408)
                                          zin(420) = zin(420) + dzij*zin(408)

                                          ! i3 = i4 =  408
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  432

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  420

                                          xin(432) = xin(432) + dxij*xin(420)
                                          yin(432) = yin(432) + dyij*yin(420)
                                          zin(432) = zin(432) + dzij*zin(420)

                                          ! i3 = i4 =  420
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  312

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  312

                                          ! do ni = 1,    3

                                          xin(312) = xin(336) + dxij*xin(300)
                                          yin(312) = yin(336) + dyij*yin(300)
                                          zin(312) = zin(336) + dzij*zin(300)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  348

                                          ! ni =    2

                                          xin(348) = xin(372) + dxij*xin(336)
                                          yin(348) = yin(372) + dyij*yin(336)
                                          zin(348) = zin(372) + dzij*zin(336)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  384

                                          ! ni =    3

                                          xin(384) = xin(408) + dxij*xin(372)
                                          yin(384) = yin(408) + dyij*yin(372)
                                          zin(384) = zin(408) + dzij*zin(372)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  420

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  324

                                          ! nj =    2

                                          ! i4 = i3 =  324

                                          ! do ni = 1,    3

                                          xin(324) = xin(348) + dxij*xin(312)
                                          yin(324) = yin(348) + dyij*yin(312)
                                          zin(324) = zin(348) + dzij*zin(312)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  360

                                          ! ni =    2

                                          xin(360) = xin(384) + dxij*xin(348)
                                          yin(360) = yin(384) + dyij*yin(348)
                                          zin(360) = zin(384) + dzij*zin(348)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  396

                                          ! ni =    3

                                          xin(396) = xin(420) + dxij*xin(384)
                                          yin(396) = yin(420) + dyij*yin(384)
                                          zin(396) = zin(420) + dzij*zin(384)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  432

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  336

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =  289

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  300

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  299

                                          xin(300) = xin(300) + dxkl*xin(299)
                                          yin(300) = yin(300) + dykl*yin(299)
                                          zin(300) = zin(300) + dzkl*zin(299)

                                          ! i3 = i4 =  299
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  298

                                          xin(299) = xin(299) + dxkl*xin(298)
                                          yin(299) = yin(299) + dykl*yin(298)
                                          zin(299) = zin(299) + dzkl*zin(298)

                                          ! i3 = i4 =  298
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  300

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  299

                                          xin(300) = xin(300) + dxkl*xin(299)
                                          yin(300) = yin(300) + dykl*yin(299)
                                          zin(300) = zin(300) + dzkl*zin(299)

                                          ! i3 = i4 =  299
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  290

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  290

                                          ! do nk = 1,    3

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

                                          xin(296) = xin(298) + dxkl*xin(295)
                                          yin(296) = yin(298) + dykl*yin(295)
                                          zin(296) = zin(298) + dzkl*zin(295)
                                          ! i4 = i4 + lang+1 =  299

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  291

                                          ! nl =    2

                                          ! i4 = i3 =  291

                                          ! do nk = 1,    3

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

                                          xin(297) = xin(299) + dxkl*xin(296)
                                          yin(297) = yin(299) + dykl*yin(296)
                                          zin(297) = zin(299) + dzkl*zin(296)
                                          ! i4 = i4 + lang+1 =  300

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  292

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  301

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  312

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  311

                                          xin(312) = xin(312) + dxkl*xin(311)
                                          yin(312) = yin(312) + dykl*yin(311)
                                          zin(312) = zin(312) + dzkl*zin(311)

                                          ! i3 = i4 =  311
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  310

                                          xin(311) = xin(311) + dxkl*xin(310)
                                          yin(311) = yin(311) + dykl*yin(310)
                                          zin(311) = zin(311) + dzkl*zin(310)

                                          ! i3 = i4 =  310
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  312

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  311

                                          xin(312) = xin(312) + dxkl*xin(311)
                                          yin(312) = yin(312) + dykl*yin(311)
                                          zin(312) = zin(312) + dzkl*zin(311)

                                          ! i3 = i4 =  311
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  302

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  302

                                          ! do nk = 1,    3

                                          xin(302) = xin(304) + dxkl*xin(301)
                                          yin(302) = yin(304) + dykl*yin(301)
                                          zin(302) = zin(304) + dzkl*zin(301)
                                          ! i4 = i4 + lang+1 =  305

                                          ! nk =    2

                                          xin(305) = xin(307) + dxkl*xin(304)
                                          yin(305) = yin(307) + dykl*yin(304)
                                          zin(305) = zin(307) + dzkl*zin(304)
                                          ! i4 = i4 + lang+1 =  308

                                          ! nk =    3

                                          xin(308) = xin(310) + dxkl*xin(307)
                                          yin(308) = yin(310) + dykl*yin(307)
                                          zin(308) = zin(310) + dzkl*zin(307)
                                          ! i4 = i4 + lang+1 =  311

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  303

                                          ! nl =    2

                                          ! i4 = i3 =  303

                                          ! do nk = 1,    3

                                          xin(303) = xin(305) + dxkl*xin(302)
                                          yin(303) = yin(305) + dykl*yin(302)
                                          zin(303) = zin(305) + dzkl*zin(302)
                                          ! i4 = i4 + lang+1 =  306

                                          ! nk =    2

                                          xin(306) = xin(308) + dxkl*xin(305)
                                          yin(306) = yin(308) + dykl*yin(305)
                                          zin(306) = zin(308) + dzkl*zin(305)
                                          ! i4 = i4 + lang+1 =  309

                                          ! nk =    3

                                          xin(309) = xin(311) + dxkl*xin(308)
                                          yin(309) = yin(311) + dykl*yin(308)
                                          zin(309) = zin(311) + dzkl*zin(308)
                                          ! i4 = i4 + lang+1 =  312

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  304

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  313

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  324

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  323

                                          xin(324) = xin(324) + dxkl*xin(323)
                                          yin(324) = yin(324) + dykl*yin(323)
                                          zin(324) = zin(324) + dzkl*zin(323)

                                          ! i3 = i4 =  323
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  322

                                          xin(323) = xin(323) + dxkl*xin(322)
                                          yin(323) = yin(323) + dykl*yin(322)
                                          zin(323) = zin(323) + dzkl*zin(322)

                                          ! i3 = i4 =  322
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  324

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  323

                                          xin(324) = xin(324) + dxkl*xin(323)
                                          yin(324) = yin(324) + dykl*yin(323)
                                          zin(324) = zin(324) + dzkl*zin(323)

                                          ! i3 = i4 =  323
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  314

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  314

                                          ! do nk = 1,    3

                                          xin(314) = xin(316) + dxkl*xin(313)
                                          yin(314) = yin(316) + dykl*yin(313)
                                          zin(314) = zin(316) + dzkl*zin(313)
                                          ! i4 = i4 + lang+1 =  317

                                          ! nk =    2

                                          xin(317) = xin(319) + dxkl*xin(316)
                                          yin(317) = yin(319) + dykl*yin(316)
                                          zin(317) = zin(319) + dzkl*zin(316)
                                          ! i4 = i4 + lang+1 =  320

                                          ! nk =    3

                                          xin(320) = xin(322) + dxkl*xin(319)
                                          yin(320) = yin(322) + dykl*yin(319)
                                          zin(320) = zin(322) + dzkl*zin(319)
                                          ! i4 = i4 + lang+1 =  323

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  315

                                          ! nl =    2

                                          ! i4 = i3 =  315

                                          ! do nk = 1,    3

                                          xin(315) = xin(317) + dxkl*xin(314)
                                          yin(315) = yin(317) + dykl*yin(314)
                                          zin(315) = zin(317) + dzkl*zin(314)
                                          ! i4 = i4 + lang+1 =  318

                                          ! nk =    2

                                          xin(318) = xin(320) + dxkl*xin(317)
                                          yin(318) = yin(320) + dykl*yin(317)
                                          zin(318) = zin(320) + dzkl*zin(317)
                                          ! i4 = i4 + lang+1 =  321

                                          ! nk =    3

                                          xin(321) = xin(323) + dxkl*xin(320)
                                          yin(321) = yin(323) + dykl*yin(320)
                                          zin(321) = zin(323) + dzkl*zin(320)
                                          ! i4 = i4 + lang+1 =  324

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  316

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  325

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  325

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  336

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  335

                                          xin(336) = xin(336) + dxkl*xin(335)
                                          yin(336) = yin(336) + dykl*yin(335)
                                          zin(336) = zin(336) + dzkl*zin(335)

                                          ! i3 = i4 =  335
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  334

                                          xin(335) = xin(335) + dxkl*xin(334)
                                          yin(335) = yin(335) + dykl*yin(334)
                                          zin(335) = zin(335) + dzkl*zin(334)

                                          ! i3 = i4 =  334
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  336

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  335

                                          xin(336) = xin(336) + dxkl*xin(335)
                                          yin(336) = yin(336) + dykl*yin(335)
                                          zin(336) = zin(336) + dzkl*zin(335)

                                          ! i3 = i4 =  335
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  326

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  326

                                          ! do nk = 1,    3

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

                                          xin(332) = xin(334) + dxkl*xin(331)
                                          yin(332) = yin(334) + dykl*yin(331)
                                          zin(332) = zin(334) + dzkl*zin(331)
                                          ! i4 = i4 + lang+1 =  335

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  327

                                          ! nl =    2

                                          ! i4 = i3 =  327

                                          ! do nk = 1,    3

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

                                          xin(333) = xin(335) + dxkl*xin(332)
                                          yin(333) = yin(335) + dykl*yin(332)
                                          zin(333) = zin(335) + dzkl*zin(332)
                                          ! i4 = i4 + lang+1 =  336

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  328

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  337

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  348

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  347

                                          xin(348) = xin(348) + dxkl*xin(347)
                                          yin(348) = yin(348) + dykl*yin(347)
                                          zin(348) = zin(348) + dzkl*zin(347)

                                          ! i3 = i4 =  347
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  346

                                          xin(347) = xin(347) + dxkl*xin(346)
                                          yin(347) = yin(347) + dykl*yin(346)
                                          zin(347) = zin(347) + dzkl*zin(346)

                                          ! i3 = i4 =  346
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  348

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  347

                                          xin(348) = xin(348) + dxkl*xin(347)
                                          yin(348) = yin(348) + dykl*yin(347)
                                          zin(348) = zin(348) + dzkl*zin(347)

                                          ! i3 = i4 =  347
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  338

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  338

                                          ! do nk = 1,    3

                                          xin(338) = xin(340) + dxkl*xin(337)
                                          yin(338) = yin(340) + dykl*yin(337)
                                          zin(338) = zin(340) + dzkl*zin(337)
                                          ! i4 = i4 + lang+1 =  341

                                          ! nk =    2

                                          xin(341) = xin(343) + dxkl*xin(340)
                                          yin(341) = yin(343) + dykl*yin(340)
                                          zin(341) = zin(343) + dzkl*zin(340)
                                          ! i4 = i4 + lang+1 =  344

                                          ! nk =    3

                                          xin(344) = xin(346) + dxkl*xin(343)
                                          yin(344) = yin(346) + dykl*yin(343)
                                          zin(344) = zin(346) + dzkl*zin(343)
                                          ! i4 = i4 + lang+1 =  347

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  339

                                          ! nl =    2

                                          ! i4 = i3 =  339

                                          ! do nk = 1,    3

                                          xin(339) = xin(341) + dxkl*xin(338)
                                          yin(339) = yin(341) + dykl*yin(338)
                                          zin(339) = zin(341) + dzkl*zin(338)
                                          ! i4 = i4 + lang+1 =  342

                                          ! nk =    2

                                          xin(342) = xin(344) + dxkl*xin(341)
                                          yin(342) = yin(344) + dykl*yin(341)
                                          zin(342) = zin(344) + dzkl*zin(341)
                                          ! i4 = i4 + lang+1 =  345

                                          ! nk =    3

                                          xin(345) = xin(347) + dxkl*xin(344)
                                          yin(345) = yin(347) + dykl*yin(344)
                                          zin(345) = zin(347) + dzkl*zin(344)
                                          ! i4 = i4 + lang+1 =  348

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  340

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  349

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  360

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  359

                                          xin(360) = xin(360) + dxkl*xin(359)
                                          yin(360) = yin(360) + dykl*yin(359)
                                          zin(360) = zin(360) + dzkl*zin(359)

                                          ! i3 = i4 =  359
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  358

                                          xin(359) = xin(359) + dxkl*xin(358)
                                          yin(359) = yin(359) + dykl*yin(358)
                                          zin(359) = zin(359) + dzkl*zin(358)

                                          ! i3 = i4 =  358
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  360

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  359

                                          xin(360) = xin(360) + dxkl*xin(359)
                                          yin(360) = yin(360) + dykl*yin(359)
                                          zin(360) = zin(360) + dzkl*zin(359)

                                          ! i3 = i4 =  359
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  350

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  350

                                          ! do nk = 1,    3

                                          xin(350) = xin(352) + dxkl*xin(349)
                                          yin(350) = yin(352) + dykl*yin(349)
                                          zin(350) = zin(352) + dzkl*zin(349)
                                          ! i4 = i4 + lang+1 =  353

                                          ! nk =    2

                                          xin(353) = xin(355) + dxkl*xin(352)
                                          yin(353) = yin(355) + dykl*yin(352)
                                          zin(353) = zin(355) + dzkl*zin(352)
                                          ! i4 = i4 + lang+1 =  356

                                          ! nk =    3

                                          xin(356) = xin(358) + dxkl*xin(355)
                                          yin(356) = yin(358) + dykl*yin(355)
                                          zin(356) = zin(358) + dzkl*zin(355)
                                          ! i4 = i4 + lang+1 =  359

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  351

                                          ! nl =    2

                                          ! i4 = i3 =  351

                                          ! do nk = 1,    3

                                          xin(351) = xin(353) + dxkl*xin(350)
                                          yin(351) = yin(353) + dykl*yin(350)
                                          zin(351) = zin(353) + dzkl*zin(350)
                                          ! i4 = i4 + lang+1 =  354

                                          ! nk =    2

                                          xin(354) = xin(356) + dxkl*xin(353)
                                          yin(354) = yin(356) + dykl*yin(353)
                                          zin(354) = zin(356) + dzkl*zin(353)
                                          ! i4 = i4 + lang+1 =  357

                                          ! nk =    3

                                          xin(357) = xin(359) + dxkl*xin(356)
                                          yin(357) = yin(359) + dykl*yin(356)
                                          zin(357) = zin(359) + dzkl*zin(356)
                                          ! i4 = i4 + lang+1 =  360

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  352

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  361

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  361

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  372

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  371

                                          xin(372) = xin(372) + dxkl*xin(371)
                                          yin(372) = yin(372) + dykl*yin(371)
                                          zin(372) = zin(372) + dzkl*zin(371)

                                          ! i3 = i4 =  371
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  370

                                          xin(371) = xin(371) + dxkl*xin(370)
                                          yin(371) = yin(371) + dykl*yin(370)
                                          zin(371) = zin(371) + dzkl*zin(370)

                                          ! i3 = i4 =  370
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  372

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  371

                                          xin(372) = xin(372) + dxkl*xin(371)
                                          yin(372) = yin(372) + dykl*yin(371)
                                          zin(372) = zin(372) + dzkl*zin(371)

                                          ! i3 = i4 =  371
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  362

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  362

                                          ! do nk = 1,    3

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

                                          xin(368) = xin(370) + dxkl*xin(367)
                                          yin(368) = yin(370) + dykl*yin(367)
                                          zin(368) = zin(370) + dzkl*zin(367)
                                          ! i4 = i4 + lang+1 =  371

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  363

                                          ! nl =    2

                                          ! i4 = i3 =  363

                                          ! do nk = 1,    3

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

                                          xin(369) = xin(371) + dxkl*xin(368)
                                          yin(369) = yin(371) + dykl*yin(368)
                                          zin(369) = zin(371) + dzkl*zin(368)
                                          ! i4 = i4 + lang+1 =  372

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  364

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  373

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  384

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  383

                                          xin(384) = xin(384) + dxkl*xin(383)
                                          yin(384) = yin(384) + dykl*yin(383)
                                          zin(384) = zin(384) + dzkl*zin(383)

                                          ! i3 = i4 =  383
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  382

                                          xin(383) = xin(383) + dxkl*xin(382)
                                          yin(383) = yin(383) + dykl*yin(382)
                                          zin(383) = zin(383) + dzkl*zin(382)

                                          ! i3 = i4 =  382
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  384

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  383

                                          xin(384) = xin(384) + dxkl*xin(383)
                                          yin(384) = yin(384) + dykl*yin(383)
                                          zin(384) = zin(384) + dzkl*zin(383)

                                          ! i3 = i4 =  383
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  374

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  374

                                          ! do nk = 1,    3

                                          xin(374) = xin(376) + dxkl*xin(373)
                                          yin(374) = yin(376) + dykl*yin(373)
                                          zin(374) = zin(376) + dzkl*zin(373)
                                          ! i4 = i4 + lang+1 =  377

                                          ! nk =    2

                                          xin(377) = xin(379) + dxkl*xin(376)
                                          yin(377) = yin(379) + dykl*yin(376)
                                          zin(377) = zin(379) + dzkl*zin(376)
                                          ! i4 = i4 + lang+1 =  380

                                          ! nk =    3

                                          xin(380) = xin(382) + dxkl*xin(379)
                                          yin(380) = yin(382) + dykl*yin(379)
                                          zin(380) = zin(382) + dzkl*zin(379)
                                          ! i4 = i4 + lang+1 =  383

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  375

                                          ! nl =    2

                                          ! i4 = i3 =  375

                                          ! do nk = 1,    3

                                          xin(375) = xin(377) + dxkl*xin(374)
                                          yin(375) = yin(377) + dykl*yin(374)
                                          zin(375) = zin(377) + dzkl*zin(374)
                                          ! i4 = i4 + lang+1 =  378

                                          ! nk =    2

                                          xin(378) = xin(380) + dxkl*xin(377)
                                          yin(378) = yin(380) + dykl*yin(377)
                                          zin(378) = zin(380) + dzkl*zin(377)
                                          ! i4 = i4 + lang+1 =  381

                                          ! nk =    3

                                          xin(381) = xin(383) + dxkl*xin(380)
                                          yin(381) = yin(383) + dykl*yin(380)
                                          zin(381) = zin(383) + dzkl*zin(380)
                                          ! i4 = i4 + lang+1 =  384

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  376

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  385

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  396

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  395

                                          xin(396) = xin(396) + dxkl*xin(395)
                                          yin(396) = yin(396) + dykl*yin(395)
                                          zin(396) = zin(396) + dzkl*zin(395)

                                          ! i3 = i4 =  395
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  394

                                          xin(395) = xin(395) + dxkl*xin(394)
                                          yin(395) = yin(395) + dykl*yin(394)
                                          zin(395) = zin(395) + dzkl*zin(394)

                                          ! i3 = i4 =  394
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  396

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  395

                                          xin(396) = xin(396) + dxkl*xin(395)
                                          yin(396) = yin(396) + dykl*yin(395)
                                          zin(396) = zin(396) + dzkl*zin(395)

                                          ! i3 = i4 =  395
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  386

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  386

                                          ! do nk = 1,    3

                                          xin(386) = xin(388) + dxkl*xin(385)
                                          yin(386) = yin(388) + dykl*yin(385)
                                          zin(386) = zin(388) + dzkl*zin(385)
                                          ! i4 = i4 + lang+1 =  389

                                          ! nk =    2

                                          xin(389) = xin(391) + dxkl*xin(388)
                                          yin(389) = yin(391) + dykl*yin(388)
                                          zin(389) = zin(391) + dzkl*zin(388)
                                          ! i4 = i4 + lang+1 =  392

                                          ! nk =    3

                                          xin(392) = xin(394) + dxkl*xin(391)
                                          yin(392) = yin(394) + dykl*yin(391)
                                          zin(392) = zin(394) + dzkl*zin(391)
                                          ! i4 = i4 + lang+1 =  395

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  387

                                          ! nl =    2

                                          ! i4 = i3 =  387

                                          ! do nk = 1,    3

                                          xin(387) = xin(389) + dxkl*xin(386)
                                          yin(387) = yin(389) + dykl*yin(386)
                                          zin(387) = zin(389) + dzkl*zin(386)
                                          ! i4 = i4 + lang+1 =  390

                                          ! nk =    2

                                          xin(390) = xin(392) + dxkl*xin(389)
                                          yin(390) = yin(392) + dykl*yin(389)
                                          zin(390) = zin(392) + dzkl*zin(389)
                                          ! i4 = i4 + lang+1 =  393

                                          ! nk =    3

                                          xin(393) = xin(395) + dxkl*xin(392)
                                          yin(393) = yin(395) + dykl*yin(392)
                                          zin(393) = zin(395) + dzkl*zin(392)
                                          ! i4 = i4 + lang+1 =  396

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  388

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  397

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  397

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  408

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  407

                                          xin(408) = xin(408) + dxkl*xin(407)
                                          yin(408) = yin(408) + dykl*yin(407)
                                          zin(408) = zin(408) + dzkl*zin(407)

                                          ! i3 = i4 =  407
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  406

                                          xin(407) = xin(407) + dxkl*xin(406)
                                          yin(407) = yin(407) + dykl*yin(406)
                                          zin(407) = zin(407) + dzkl*zin(406)

                                          ! i3 = i4 =  406
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  408

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  407

                                          xin(408) = xin(408) + dxkl*xin(407)
                                          yin(408) = yin(408) + dykl*yin(407)
                                          zin(408) = zin(408) + dzkl*zin(407)

                                          ! i3 = i4 =  407
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  398

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  398

                                          ! do nk = 1,    3

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

                                          xin(404) = xin(406) + dxkl*xin(403)
                                          yin(404) = yin(406) + dykl*yin(403)
                                          zin(404) = zin(406) + dzkl*zin(403)
                                          ! i4 = i4 + lang+1 =  407

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  399

                                          ! nl =    2

                                          ! i4 = i3 =  399

                                          ! do nk = 1,    3

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

                                          xin(405) = xin(407) + dxkl*xin(404)
                                          yin(405) = yin(407) + dykl*yin(404)
                                          zin(405) = zin(407) + dzkl*zin(404)
                                          ! i4 = i4 + lang+1 =  408

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  400

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  409

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  420

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  419

                                          xin(420) = xin(420) + dxkl*xin(419)
                                          yin(420) = yin(420) + dykl*yin(419)
                                          zin(420) = zin(420) + dzkl*zin(419)

                                          ! i3 = i4 =  419
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  418

                                          xin(419) = xin(419) + dxkl*xin(418)
                                          yin(419) = yin(419) + dykl*yin(418)
                                          zin(419) = zin(419) + dzkl*zin(418)

                                          ! i3 = i4 =  418
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  420

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  419

                                          xin(420) = xin(420) + dxkl*xin(419)
                                          yin(420) = yin(420) + dykl*yin(419)
                                          zin(420) = zin(420) + dzkl*zin(419)

                                          ! i3 = i4 =  419
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  410

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  410

                                          ! do nk = 1,    3

                                          xin(410) = xin(412) + dxkl*xin(409)
                                          yin(410) = yin(412) + dykl*yin(409)
                                          zin(410) = zin(412) + dzkl*zin(409)
                                          ! i4 = i4 + lang+1 =  413

                                          ! nk =    2

                                          xin(413) = xin(415) + dxkl*xin(412)
                                          yin(413) = yin(415) + dykl*yin(412)
                                          zin(413) = zin(415) + dzkl*zin(412)
                                          ! i4 = i4 + lang+1 =  416

                                          ! nk =    3

                                          xin(416) = xin(418) + dxkl*xin(415)
                                          yin(416) = yin(418) + dykl*yin(415)
                                          zin(416) = zin(418) + dzkl*zin(415)
                                          ! i4 = i4 + lang+1 =  419

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  411

                                          ! nl =    2

                                          ! i4 = i3 =  411

                                          ! do nk = 1,    3

                                          xin(411) = xin(413) + dxkl*xin(410)
                                          yin(411) = yin(413) + dykl*yin(410)
                                          zin(411) = zin(413) + dzkl*zin(410)
                                          ! i4 = i4 + lang+1 =  414

                                          ! nk =    2

                                          xin(414) = xin(416) + dxkl*xin(413)
                                          yin(414) = yin(416) + dykl*yin(413)
                                          zin(414) = zin(416) + dzkl*zin(413)
                                          ! i4 = i4 + lang+1 =  417

                                          ! nk =    3

                                          xin(417) = xin(419) + dxkl*xin(416)
                                          yin(417) = yin(419) + dykl*yin(416)
                                          zin(417) = zin(419) + dzkl*zin(416)
                                          ! i4 = i4 + lang+1 =  420

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  412

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  421

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  432

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  431

                                          xin(432) = xin(432) + dxkl*xin(431)
                                          yin(432) = yin(432) + dykl*yin(431)
                                          zin(432) = zin(432) + dzkl*zin(431)

                                          ! i3 = i4 =  431
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  430

                                          xin(431) = xin(431) + dxkl*xin(430)
                                          yin(431) = yin(431) + dykl*yin(430)
                                          zin(431) = zin(431) + dzkl*zin(430)

                                          ! i3 = i4 =  430
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  432

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  431

                                          xin(432) = xin(432) + dxkl*xin(431)
                                          yin(432) = yin(432) + dykl*yin(431)
                                          zin(432) = zin(432) + dzkl*zin(431)

                                          ! i3 = i4 =  431
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  422

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  422

                                          ! do nk = 1,    3

                                          xin(422) = xin(424) + dxkl*xin(421)
                                          yin(422) = yin(424) + dykl*yin(421)
                                          zin(422) = zin(424) + dzkl*zin(421)
                                          ! i4 = i4 + lang+1 =  425

                                          ! nk =    2

                                          xin(425) = xin(427) + dxkl*xin(424)
                                          yin(425) = yin(427) + dykl*yin(424)
                                          zin(425) = zin(427) + dzkl*zin(424)
                                          ! i4 = i4 + lang+1 =  428

                                          ! nk =    3

                                          xin(428) = xin(430) + dxkl*xin(427)
                                          yin(428) = yin(430) + dykl*yin(427)
                                          zin(428) = zin(430) + dzkl*zin(427)
                                          ! i4 = i4 + lang+1 =  431

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  423

                                          ! nl =    2

                                          ! i4 = i3 =  423

                                          ! do nk = 1,    3

                                          xin(423) = xin(425) + dxkl*xin(422)
                                          yin(423) = yin(425) + dykl*yin(422)
                                          zin(423) = zin(425) + dzkl*zin(422)
                                          ! i4 = i4 + lang+1 =  426

                                          ! nk =    2

                                          xin(426) = xin(428) + dxkl*xin(425)
                                          yin(426) = yin(428) + dykl*yin(425)
                                          zin(426) = zin(428) + dzkl*zin(425)
                                          ! i4 = i4 + lang+1 =  429

                                          ! nk =    3

                                          xin(429) = xin(431) + dxkl*xin(428)
                                          yin(429) = yin(431) + dykl*yin(428)
                                          zin(429) = zin(431) + dzkl*zin(428)
                                          ! i4 = i4 + lang+1 =  432

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  424

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  433

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  433

                                          ! end do

                                          ! *** Now root =    4

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  432

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

                                          ! i1 = in(1) =  433

                                          xin(433) = 1.0_dp
                                          yin(433) = 1.0_dp
                                          zin(433) = f00

                                          ! i2 = in(2) =  469
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(469) = xc00
                                          yin(469) = yc00
                                          zin(469) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  436

                                          xin(436) = xcp00
                                          yin(436) = ycp00
                                          zin(436) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  472
                                          ! i2 =  469

                                          xin(472) = xcp00*xin(469) + cp10
                                          yin(472) = ycp00*yin(469) + cp10
                                          zin(472) = zcp00*zin(469) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  433
                                          ! i4 = i2 =  469

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  505
                                          ! i3 =  433
                                          ! i4 =  469

                                          xin(505) = c10*xin(433) + xc00*xin(469)
                                          yin(505) = c10*yin(433) + yc00*yin(469)
                                          zin(505) = c10*zin(433) + zc00*zin(469)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  508
                                          ! i5 =  505
                                          ! i4 =  469

                                          xin(508) = xcp00*xin(505) + cp10*xin(469)
                                          yin(508) = ycp00*yin(505) + cp10*yin(469)
                                          zin(508) = zcp00*zin(505) + cp10*zin(469)

                                          ! ------------------

                                          ! i3 = i4 =  469
                                          ! i4 = i5 =  505

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  541
                                          ! i3 =  469
                                          ! i4 =  505

                                          xin(541) = c10*xin(469) + xc00*xin(505)
                                          yin(541) = c10*yin(469) + yc00*yin(505)
                                          zin(541) = c10*zin(469) + zc00*zin(505)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  544
                                          ! i5 =  541
                                          ! i4 =  505

                                          xin(544) = xcp00*xin(541) + cp10*xin(505)
                                          yin(544) = ycp00*yin(541) + cp10*yin(505)
                                          zin(544) = zcp00*zin(541) + cp10*zin(505)

                                          ! ------------------

                                          ! i3 = i4 =  505
                                          ! i4 = i5 =  541

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  553
                                          ! i3 =  505
                                          ! i4 =  541

                                          xin(553) = c10*xin(505) + xc00*xin(541)
                                          yin(553) = c10*yin(505) + yc00*yin(541)
                                          zin(553) = c10*zin(505) + zc00*zin(541)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  556
                                          ! i5 =  553
                                          ! i4 =  541

                                          xin(556) = xcp00*xin(553) + cp10*xin(541)
                                          yin(556) = ycp00*yin(553) + cp10*yin(541)
                                          zin(556) = zcp00*zin(553) + cp10*zin(541)

                                          ! ------------------

                                          ! i3 = i4 =  541
                                          ! i4 = i5 =  553

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  565
                                          ! i3 =  541
                                          ! i4 =  553

                                          xin(565) = c10*xin(541) + xc00*xin(553)
                                          yin(565) = c10*yin(541) + yc00*yin(553)
                                          zin(565) = c10*zin(541) + zc00*zin(553)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  568
                                          ! i5 =  565
                                          ! i4 =  553

                                          xin(568) = xcp00*xin(565) + cp10*xin(553)
                                          yin(568) = ycp00*yin(565) + cp10*yin(553)
                                          zin(568) = zcp00*zin(565) + cp10*zin(553)

                                          ! ------------------

                                          ! i3 = i4 =  553
                                          ! i4 = i5 =  565

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =  433
                                          ! i4 = i1+k2 =  436

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  439
                                          ! i3 =  433
                                          ! i4 =  436

                                          xin(439) = cp01*xin(433) + xcp00*xin(436)
                                          yin(439) = cp01*yin(433) + ycp00*yin(436)
                                          zin(439) = cp01*zin(433) + zcp00*zin(436)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  475

                                          xin(475) = xc00*xin(439) + c01*xin(436)
                                          yin(475) = yc00*yin(439) + c01*yin(436)
                                          zin(475) = zc00*zin(439) + c01*zin(436)

                                          ! ------------------

                                          ! i3 = i4 =  436
                                          ! i4 = i5 =  439

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  442
                                          ! i3 =  436
                                          ! i4 =  439

                                          xin(442) = cp01*xin(436) + xcp00*xin(439)
                                          yin(442) = cp01*yin(436) + ycp00*yin(439)
                                          zin(442) = cp01*zin(436) + zcp00*zin(439)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  478

                                          xin(478) = xc00*xin(442) + c01*xin(439)
                                          yin(478) = yc00*yin(442) + c01*yin(439)
                                          zin(478) = zc00*zin(442) + c01*zin(439)

                                          ! ------------------

                                          ! i3 = i4 =  439
                                          ! i4 = i5 =  442

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  443
                                          ! i3 =  439
                                          ! i4 =  442

                                          xin(443) = cp01*xin(439) + xcp00*xin(442)
                                          yin(443) = cp01*yin(439) + ycp00*yin(442)
                                          zin(443) = cp01*zin(439) + zcp00*zin(442)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  479

                                          xin(479) = xc00*xin(443) + c01*xin(442)
                                          yin(479) = yc00*yin(443) + c01*yin(442)
                                          zin(479) = zc00*zin(443) + c01*zin(442)

                                          ! ------------------

                                          ! i3 = i4 =  442
                                          ! i4 = i5 =  443

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  444
                                          ! i3 =  442
                                          ! i4 =  443

                                          xin(444) = cp01*xin(442) + xcp00*xin(443)
                                          yin(444) = cp01*yin(442) + ycp00*yin(443)
                                          zin(444) = cp01*zin(442) + zcp00*zin(443)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  480

                                          xin(480) = xc00*xin(444) + c01*xin(443)
                                          yin(480) = yc00*yin(444) + c01*yin(443)
                                          zin(480) = zc00*zin(444) + c01*zin(443)

                                          ! ------------------

                                          ! i3 = i4 =  443
                                          ! i4 = i5 =  444

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =  433
                                          ! i4 = i2 =  469

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  505

                                          xin(511) = c10*xin(439) + xc00*xin(475) + c01*xin(472)
                                          yin(511) = c10*yin(439) + yc00*yin(475) + c01*yin(472)
                                          zin(511) = c10*zin(439) + zc00*zin(475) + c01*zin(472)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  469
                                          ! i4 = i5 =  505

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  541

                                          xin(547) = c10*xin(475) + xc00*xin(511) + c01*xin(508)
                                          yin(547) = c10*yin(475) + yc00*yin(511) + c01*yin(508)
                                          zin(547) = c10*zin(475) + zc00*zin(511) + c01*zin(508)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  505
                                          ! i4 = i5 =  541

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  553

                                          xin(559) = c10*xin(511) + xc00*xin(547) + c01*xin(544)
                                          yin(559) = c10*yin(511) + yc00*yin(547) + c01*yin(544)
                                          zin(559) = c10*zin(511) + zc00*zin(547) + c01*zin(544)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  541
                                          ! i4 = i5 =  553

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  565

                                          xin(571) = c10*xin(547) + xc00*xin(559) + c01*xin(556)
                                          yin(571) = c10*yin(547) + yc00*yin(559) + c01*yin(556)
                                          zin(571) = c10*zin(547) + zc00*zin(559) + c01*zin(556)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  553
                                          ! i4 = i5 =  565

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =  433
                                          ! i4 = i2 =  469

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  505

                                          xin(514) = c10*xin(442) + xc00*xin(478) + c01*xin(475)
                                          yin(514) = c10*yin(442) + yc00*yin(478) + c01*yin(475)
                                          zin(514) = c10*zin(442) + zc00*zin(478) + c01*zin(475)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  469
                                          ! i4 = i5 =  505

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  541

                                          xin(550) = c10*xin(478) + xc00*xin(514) + c01*xin(511)
                                          yin(550) = c10*yin(478) + yc00*yin(514) + c01*yin(511)
                                          zin(550) = c10*zin(478) + zc00*zin(514) + c01*zin(511)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  505
                                          ! i4 = i5 =  541

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  553

                                          xin(562) = c10*xin(514) + xc00*xin(550) + c01*xin(547)
                                          yin(562) = c10*yin(514) + yc00*yin(550) + c01*yin(547)
                                          zin(562) = c10*zin(514) + zc00*zin(550) + c01*zin(547)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  541
                                          ! i4 = i5 =  553

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  565

                                          xin(574) = c10*xin(550) + xc00*xin(562) + c01*xin(559)
                                          yin(574) = c10*yin(550) + yc00*yin(562) + c01*yin(559)
                                          zin(574) = c10*zin(550) + zc00*zin(562) + c01*zin(559)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  553
                                          ! i4 = i5 =  565

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =  433
                                          ! i4 = i2 =  469

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  505

                                          xin(515) = c10*xin(443) + xc00*xin(479) + c01*xin(478)
                                          yin(515) = c10*yin(443) + yc00*yin(479) + c01*yin(478)
                                          zin(515) = c10*zin(443) + zc00*zin(479) + c01*zin(478)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  469
                                          ! i4 = i5 =  505

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  541

                                          xin(551) = c10*xin(479) + xc00*xin(515) + c01*xin(514)
                                          yin(551) = c10*yin(479) + yc00*yin(515) + c01*yin(514)
                                          zin(551) = c10*zin(479) + zc00*zin(515) + c01*zin(514)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  505
                                          ! i4 = i5 =  541

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  553

                                          xin(563) = c10*xin(515) + xc00*xin(551) + c01*xin(550)
                                          yin(563) = c10*yin(515) + yc00*yin(551) + c01*yin(550)
                                          zin(563) = c10*zin(515) + zc00*zin(551) + c01*zin(550)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  541
                                          ! i4 = i5 =  553

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  565

                                          xin(575) = c10*xin(551) + xc00*xin(563) + c01*xin(562)
                                          yin(575) = c10*yin(551) + yc00*yin(563) + c01*yin(562)
                                          zin(575) = c10*zin(551) + zc00*zin(563) + c01*zin(562)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  553
                                          ! i4 = i5 =  565

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =  433
                                          ! i4 = i2 =  469

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  505

                                          xin(516) = c10*xin(444) + xc00*xin(480) + c01*xin(479)
                                          yin(516) = c10*yin(444) + yc00*yin(480) + c01*yin(479)
                                          zin(516) = c10*zin(444) + zc00*zin(480) + c01*zin(479)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  469
                                          ! i4 = i5 =  505

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  541

                                          xin(552) = c10*xin(480) + xc00*xin(516) + c01*xin(515)
                                          yin(552) = c10*yin(480) + yc00*yin(516) + c01*yin(515)
                                          zin(552) = c10*zin(480) + zc00*zin(516) + c01*zin(515)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  505
                                          ! i4 = i5 =  541

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  553

                                          xin(564) = c10*xin(516) + xc00*xin(552) + c01*xin(551)
                                          yin(564) = c10*yin(516) + yc00*yin(552) + c01*yin(551)
                                          zin(564) = c10*zin(516) + zc00*zin(552) + c01*zin(551)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  541
                                          ! i4 = i5 =  553

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  565

                                          xin(576) = c10*xin(552) + xc00*xin(564) + c01*xin(563)
                                          yin(576) = c10*yin(552) + yc00*yin(564) + c01*yin(563)
                                          zin(576) = c10*zin(552) + zc00*zin(564) + c01*zin(563)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  553
                                          ! i4 = i5 =  565

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  565

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  565

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  553

                                          xin(565) = xin(565) + dxij*xin(553)
                                          yin(565) = yin(565) + dyij*yin(553)
                                          zin(565) = zin(565) + dzij*zin(553)

                                          ! i3 = i4 =  553
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  541

                                          xin(553) = xin(553) + dxij*xin(541)
                                          yin(553) = yin(553) + dyij*yin(541)
                                          zin(553) = zin(553) + dzij*zin(541)

                                          ! i3 = i4 =  541
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  565

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  553

                                          xin(565) = xin(565) + dxij*xin(553)
                                          yin(565) = yin(565) + dyij*yin(553)
                                          zin(565) = zin(565) + dzij*zin(553)

                                          ! i3 = i4 =  553
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  445

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  445

                                          ! do ni = 1,    3

                                          xin(445) = xin(469) + dxij*xin(433)
                                          yin(445) = yin(469) + dyij*yin(433)
                                          zin(445) = zin(469) + dzij*zin(433)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  481

                                          ! ni =    2

                                          xin(481) = xin(505) + dxij*xin(469)
                                          yin(481) = yin(505) + dyij*yin(469)
                                          zin(481) = zin(505) + dzij*zin(469)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  517

                                          ! ni =    3

                                          xin(517) = xin(541) + dxij*xin(505)
                                          yin(517) = yin(541) + dyij*yin(505)
                                          zin(517) = zin(541) + dzij*zin(505)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  553

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  457

                                          ! nj =    2

                                          ! i4 = i3 =  457

                                          ! do ni = 1,    3

                                          xin(457) = xin(481) + dxij*xin(445)
                                          yin(457) = yin(481) + dyij*yin(445)
                                          zin(457) = zin(481) + dzij*zin(445)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  493

                                          ! ni =    2

                                          xin(493) = xin(517) + dxij*xin(481)
                                          yin(493) = yin(517) + dyij*yin(481)
                                          zin(493) = zin(517) + dzij*zin(481)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  529

                                          ! ni =    3

                                          xin(529) = xin(553) + dxij*xin(517)
                                          yin(529) = yin(553) + dyij*yin(517)
                                          zin(529) = zin(553) + dzij*zin(517)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  565

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  469

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  568

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  556

                                          xin(568) = xin(568) + dxij*xin(556)
                                          yin(568) = yin(568) + dyij*yin(556)
                                          zin(568) = zin(568) + dzij*zin(556)

                                          ! i3 = i4 =  556
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  544

                                          xin(556) = xin(556) + dxij*xin(544)
                                          yin(556) = yin(556) + dyij*yin(544)
                                          zin(556) = zin(556) + dzij*zin(544)

                                          ! i3 = i4 =  544
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  568

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  556

                                          xin(568) = xin(568) + dxij*xin(556)
                                          yin(568) = yin(568) + dyij*yin(556)
                                          zin(568) = zin(568) + dzij*zin(556)

                                          ! i3 = i4 =  556
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  448

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  448

                                          ! do ni = 1,    3

                                          xin(448) = xin(472) + dxij*xin(436)
                                          yin(448) = yin(472) + dyij*yin(436)
                                          zin(448) = zin(472) + dzij*zin(436)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  484

                                          ! ni =    2

                                          xin(484) = xin(508) + dxij*xin(472)
                                          yin(484) = yin(508) + dyij*yin(472)
                                          zin(484) = zin(508) + dzij*zin(472)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  520

                                          ! ni =    3

                                          xin(520) = xin(544) + dxij*xin(508)
                                          yin(520) = yin(544) + dyij*yin(508)
                                          zin(520) = zin(544) + dzij*zin(508)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  556

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  460

                                          ! nj =    2

                                          ! i4 = i3 =  460

                                          ! do ni = 1,    3

                                          xin(460) = xin(484) + dxij*xin(448)
                                          yin(460) = yin(484) + dyij*yin(448)
                                          zin(460) = zin(484) + dzij*zin(448)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  496

                                          ! ni =    2

                                          xin(496) = xin(520) + dxij*xin(484)
                                          yin(496) = yin(520) + dyij*yin(484)
                                          zin(496) = zin(520) + dzij*zin(484)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  532

                                          ! ni =    3

                                          xin(532) = xin(556) + dxij*xin(520)
                                          yin(532) = yin(556) + dyij*yin(520)
                                          zin(532) = zin(556) + dzij*zin(520)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  568

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  472

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  571

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  559

                                          xin(571) = xin(571) + dxij*xin(559)
                                          yin(571) = yin(571) + dyij*yin(559)
                                          zin(571) = zin(571) + dzij*zin(559)

                                          ! i3 = i4 =  559
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  547

                                          xin(559) = xin(559) + dxij*xin(547)
                                          yin(559) = yin(559) + dyij*yin(547)
                                          zin(559) = zin(559) + dzij*zin(547)

                                          ! i3 = i4 =  547
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  571

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  559

                                          xin(571) = xin(571) + dxij*xin(559)
                                          yin(571) = yin(571) + dyij*yin(559)
                                          zin(571) = zin(571) + dzij*zin(559)

                                          ! i3 = i4 =  559
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  451

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  451

                                          ! do ni = 1,    3

                                          xin(451) = xin(475) + dxij*xin(439)
                                          yin(451) = yin(475) + dyij*yin(439)
                                          zin(451) = zin(475) + dzij*zin(439)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  487

                                          ! ni =    2

                                          xin(487) = xin(511) + dxij*xin(475)
                                          yin(487) = yin(511) + dyij*yin(475)
                                          zin(487) = zin(511) + dzij*zin(475)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  523

                                          ! ni =    3

                                          xin(523) = xin(547) + dxij*xin(511)
                                          yin(523) = yin(547) + dyij*yin(511)
                                          zin(523) = zin(547) + dzij*zin(511)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  559

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  463

                                          ! nj =    2

                                          ! i4 = i3 =  463

                                          ! do ni = 1,    3

                                          xin(463) = xin(487) + dxij*xin(451)
                                          yin(463) = yin(487) + dyij*yin(451)
                                          zin(463) = zin(487) + dzij*zin(451)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  499

                                          ! ni =    2

                                          xin(499) = xin(523) + dxij*xin(487)
                                          yin(499) = yin(523) + dyij*yin(487)
                                          zin(499) = zin(523) + dzij*zin(487)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  535

                                          ! ni =    3

                                          xin(535) = xin(559) + dxij*xin(523)
                                          yin(535) = yin(559) + dyij*yin(523)
                                          zin(535) = zin(559) + dzij*zin(523)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  571

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  475

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  574

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  562

                                          xin(574) = xin(574) + dxij*xin(562)
                                          yin(574) = yin(574) + dyij*yin(562)
                                          zin(574) = zin(574) + dzij*zin(562)

                                          ! i3 = i4 =  562
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  550

                                          xin(562) = xin(562) + dxij*xin(550)
                                          yin(562) = yin(562) + dyij*yin(550)
                                          zin(562) = zin(562) + dzij*zin(550)

                                          ! i3 = i4 =  550
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  574

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  562

                                          xin(574) = xin(574) + dxij*xin(562)
                                          yin(574) = yin(574) + dyij*yin(562)
                                          zin(574) = zin(574) + dzij*zin(562)

                                          ! i3 = i4 =  562
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  454

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  454

                                          ! do ni = 1,    3

                                          xin(454) = xin(478) + dxij*xin(442)
                                          yin(454) = yin(478) + dyij*yin(442)
                                          zin(454) = zin(478) + dzij*zin(442)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  490

                                          ! ni =    2

                                          xin(490) = xin(514) + dxij*xin(478)
                                          yin(490) = yin(514) + dyij*yin(478)
                                          zin(490) = zin(514) + dzij*zin(478)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  526

                                          ! ni =    3

                                          xin(526) = xin(550) + dxij*xin(514)
                                          yin(526) = yin(550) + dyij*yin(514)
                                          zin(526) = zin(550) + dzij*zin(514)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  562

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  466

                                          ! nj =    2

                                          ! i4 = i3 =  466

                                          ! do ni = 1,    3

                                          xin(466) = xin(490) + dxij*xin(454)
                                          yin(466) = yin(490) + dyij*yin(454)
                                          zin(466) = zin(490) + dzij*zin(454)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  502

                                          ! ni =    2

                                          xin(502) = xin(526) + dxij*xin(490)
                                          yin(502) = yin(526) + dyij*yin(490)
                                          zin(502) = zin(526) + dzij*zin(490)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  538

                                          ! ni =    3

                                          xin(538) = xin(562) + dxij*xin(526)
                                          yin(538) = yin(562) + dyij*yin(526)
                                          zin(538) = zin(562) + dzij*zin(526)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  574

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  478

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  575

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  563

                                          xin(575) = xin(575) + dxij*xin(563)
                                          yin(575) = yin(575) + dyij*yin(563)
                                          zin(575) = zin(575) + dzij*zin(563)

                                          ! i3 = i4 =  563
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  551

                                          xin(563) = xin(563) + dxij*xin(551)
                                          yin(563) = yin(563) + dyij*yin(551)
                                          zin(563) = zin(563) + dzij*zin(551)

                                          ! i3 = i4 =  551
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  575

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  563

                                          xin(575) = xin(575) + dxij*xin(563)
                                          yin(575) = yin(575) + dyij*yin(563)
                                          zin(575) = zin(575) + dzij*zin(563)

                                          ! i3 = i4 =  563
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  455

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  455

                                          ! do ni = 1,    3

                                          xin(455) = xin(479) + dxij*xin(443)
                                          yin(455) = yin(479) + dyij*yin(443)
                                          zin(455) = zin(479) + dzij*zin(443)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  491

                                          ! ni =    2

                                          xin(491) = xin(515) + dxij*xin(479)
                                          yin(491) = yin(515) + dyij*yin(479)
                                          zin(491) = zin(515) + dzij*zin(479)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  527

                                          ! ni =    3

                                          xin(527) = xin(551) + dxij*xin(515)
                                          yin(527) = yin(551) + dyij*yin(515)
                                          zin(527) = zin(551) + dzij*zin(515)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  563

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  467

                                          ! nj =    2

                                          ! i4 = i3 =  467

                                          ! do ni = 1,    3

                                          xin(467) = xin(491) + dxij*xin(455)
                                          yin(467) = yin(491) + dyij*yin(455)
                                          zin(467) = zin(491) + dzij*zin(455)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  503

                                          ! ni =    2

                                          xin(503) = xin(527) + dxij*xin(491)
                                          yin(503) = yin(527) + dyij*yin(491)
                                          zin(503) = zin(527) + dzij*zin(491)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  539

                                          ! ni =    3

                                          xin(539) = xin(563) + dxij*xin(527)
                                          yin(539) = yin(563) + dyij*yin(527)
                                          zin(539) = zin(563) + dzij*zin(527)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  575

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  479

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  576

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  564

                                          xin(576) = xin(576) + dxij*xin(564)
                                          yin(576) = yin(576) + dyij*yin(564)
                                          zin(576) = zin(576) + dzij*zin(564)

                                          ! i3 = i4 =  564
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  552

                                          xin(564) = xin(564) + dxij*xin(552)
                                          yin(564) = yin(564) + dyij*yin(552)
                                          zin(564) = zin(564) + dzij*zin(552)

                                          ! i3 = i4 =  552
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  576

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  564

                                          xin(576) = xin(576) + dxij*xin(564)
                                          yin(576) = yin(576) + dyij*yin(564)
                                          zin(576) = zin(576) + dzij*zin(564)

                                          ! i3 = i4 =  564
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  456

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  456

                                          ! do ni = 1,    3

                                          xin(456) = xin(480) + dxij*xin(444)
                                          yin(456) = yin(480) + dyij*yin(444)
                                          zin(456) = zin(480) + dzij*zin(444)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  492

                                          ! ni =    2

                                          xin(492) = xin(516) + dxij*xin(480)
                                          yin(492) = yin(516) + dyij*yin(480)
                                          zin(492) = zin(516) + dzij*zin(480)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  528

                                          ! ni =    3

                                          xin(528) = xin(552) + dxij*xin(516)
                                          yin(528) = yin(552) + dyij*yin(516)
                                          zin(528) = zin(552) + dzij*zin(516)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  564

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  468

                                          ! nj =    2

                                          ! i4 = i3 =  468

                                          ! do ni = 1,    3

                                          xin(468) = xin(492) + dxij*xin(456)
                                          yin(468) = yin(492) + dyij*yin(456)
                                          zin(468) = zin(492) + dzij*zin(456)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  504

                                          ! ni =    2

                                          xin(504) = xin(528) + dxij*xin(492)
                                          yin(504) = yin(528) + dyij*yin(492)
                                          zin(504) = zin(528) + dzij*zin(492)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  540

                                          ! ni =    3

                                          xin(540) = xin(564) + dxij*xin(528)
                                          yin(540) = yin(564) + dyij*yin(528)
                                          zin(540) = zin(564) + dzij*zin(528)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  576

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  480

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =  433

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  444

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  443

                                          xin(444) = xin(444) + dxkl*xin(443)
                                          yin(444) = yin(444) + dykl*yin(443)
                                          zin(444) = zin(444) + dzkl*zin(443)

                                          ! i3 = i4 =  443
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  442

                                          xin(443) = xin(443) + dxkl*xin(442)
                                          yin(443) = yin(443) + dykl*yin(442)
                                          zin(443) = zin(443) + dzkl*zin(442)

                                          ! i3 = i4 =  442
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  444

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  443

                                          xin(444) = xin(444) + dxkl*xin(443)
                                          yin(444) = yin(444) + dykl*yin(443)
                                          zin(444) = zin(444) + dzkl*zin(443)

                                          ! i3 = i4 =  443
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  434

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  434

                                          ! do nk = 1,    3

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

                                          xin(440) = xin(442) + dxkl*xin(439)
                                          yin(440) = yin(442) + dykl*yin(439)
                                          zin(440) = zin(442) + dzkl*zin(439)
                                          ! i4 = i4 + lang+1 =  443

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  435

                                          ! nl =    2

                                          ! i4 = i3 =  435

                                          ! do nk = 1,    3

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

                                          xin(441) = xin(443) + dxkl*xin(440)
                                          yin(441) = yin(443) + dykl*yin(440)
                                          zin(441) = zin(443) + dzkl*zin(440)
                                          ! i4 = i4 + lang+1 =  444

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  436

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  445

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  456

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  455

                                          xin(456) = xin(456) + dxkl*xin(455)
                                          yin(456) = yin(456) + dykl*yin(455)
                                          zin(456) = zin(456) + dzkl*zin(455)

                                          ! i3 = i4 =  455
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  454

                                          xin(455) = xin(455) + dxkl*xin(454)
                                          yin(455) = yin(455) + dykl*yin(454)
                                          zin(455) = zin(455) + dzkl*zin(454)

                                          ! i3 = i4 =  454
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  456

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  455

                                          xin(456) = xin(456) + dxkl*xin(455)
                                          yin(456) = yin(456) + dykl*yin(455)
                                          zin(456) = zin(456) + dzkl*zin(455)

                                          ! i3 = i4 =  455
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  446

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  446

                                          ! do nk = 1,    3

                                          xin(446) = xin(448) + dxkl*xin(445)
                                          yin(446) = yin(448) + dykl*yin(445)
                                          zin(446) = zin(448) + dzkl*zin(445)
                                          ! i4 = i4 + lang+1 =  449

                                          ! nk =    2

                                          xin(449) = xin(451) + dxkl*xin(448)
                                          yin(449) = yin(451) + dykl*yin(448)
                                          zin(449) = zin(451) + dzkl*zin(448)
                                          ! i4 = i4 + lang+1 =  452

                                          ! nk =    3

                                          xin(452) = xin(454) + dxkl*xin(451)
                                          yin(452) = yin(454) + dykl*yin(451)
                                          zin(452) = zin(454) + dzkl*zin(451)
                                          ! i4 = i4 + lang+1 =  455

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  447

                                          ! nl =    2

                                          ! i4 = i3 =  447

                                          ! do nk = 1,    3

                                          xin(447) = xin(449) + dxkl*xin(446)
                                          yin(447) = yin(449) + dykl*yin(446)
                                          zin(447) = zin(449) + dzkl*zin(446)
                                          ! i4 = i4 + lang+1 =  450

                                          ! nk =    2

                                          xin(450) = xin(452) + dxkl*xin(449)
                                          yin(450) = yin(452) + dykl*yin(449)
                                          zin(450) = zin(452) + dzkl*zin(449)
                                          ! i4 = i4 + lang+1 =  453

                                          ! nk =    3

                                          xin(453) = xin(455) + dxkl*xin(452)
                                          yin(453) = yin(455) + dykl*yin(452)
                                          zin(453) = zin(455) + dzkl*zin(452)
                                          ! i4 = i4 + lang+1 =  456

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  448

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  457

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  468

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  467

                                          xin(468) = xin(468) + dxkl*xin(467)
                                          yin(468) = yin(468) + dykl*yin(467)
                                          zin(468) = zin(468) + dzkl*zin(467)

                                          ! i3 = i4 =  467
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  466

                                          xin(467) = xin(467) + dxkl*xin(466)
                                          yin(467) = yin(467) + dykl*yin(466)
                                          zin(467) = zin(467) + dzkl*zin(466)

                                          ! i3 = i4 =  466
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  468

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  467

                                          xin(468) = xin(468) + dxkl*xin(467)
                                          yin(468) = yin(468) + dykl*yin(467)
                                          zin(468) = zin(468) + dzkl*zin(467)

                                          ! i3 = i4 =  467
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  458

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  458

                                          ! do nk = 1,    3

                                          xin(458) = xin(460) + dxkl*xin(457)
                                          yin(458) = yin(460) + dykl*yin(457)
                                          zin(458) = zin(460) + dzkl*zin(457)
                                          ! i4 = i4 + lang+1 =  461

                                          ! nk =    2

                                          xin(461) = xin(463) + dxkl*xin(460)
                                          yin(461) = yin(463) + dykl*yin(460)
                                          zin(461) = zin(463) + dzkl*zin(460)
                                          ! i4 = i4 + lang+1 =  464

                                          ! nk =    3

                                          xin(464) = xin(466) + dxkl*xin(463)
                                          yin(464) = yin(466) + dykl*yin(463)
                                          zin(464) = zin(466) + dzkl*zin(463)
                                          ! i4 = i4 + lang+1 =  467

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  459

                                          ! nl =    2

                                          ! i4 = i3 =  459

                                          ! do nk = 1,    3

                                          xin(459) = xin(461) + dxkl*xin(458)
                                          yin(459) = yin(461) + dykl*yin(458)
                                          zin(459) = zin(461) + dzkl*zin(458)
                                          ! i4 = i4 + lang+1 =  462

                                          ! nk =    2

                                          xin(462) = xin(464) + dxkl*xin(461)
                                          yin(462) = yin(464) + dykl*yin(461)
                                          zin(462) = zin(464) + dzkl*zin(461)
                                          ! i4 = i4 + lang+1 =  465

                                          ! nk =    3

                                          xin(465) = xin(467) + dxkl*xin(464)
                                          yin(465) = yin(467) + dykl*yin(464)
                                          zin(465) = zin(467) + dzkl*zin(464)
                                          ! i4 = i4 + lang+1 =  468

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  460

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  469

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  469

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  480

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  479

                                          xin(480) = xin(480) + dxkl*xin(479)
                                          yin(480) = yin(480) + dykl*yin(479)
                                          zin(480) = zin(480) + dzkl*zin(479)

                                          ! i3 = i4 =  479
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  478

                                          xin(479) = xin(479) + dxkl*xin(478)
                                          yin(479) = yin(479) + dykl*yin(478)
                                          zin(479) = zin(479) + dzkl*zin(478)

                                          ! i3 = i4 =  478
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  480

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  479

                                          xin(480) = xin(480) + dxkl*xin(479)
                                          yin(480) = yin(480) + dykl*yin(479)
                                          zin(480) = zin(480) + dzkl*zin(479)

                                          ! i3 = i4 =  479
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  470

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  470

                                          ! do nk = 1,    3

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

                                          xin(476) = xin(478) + dxkl*xin(475)
                                          yin(476) = yin(478) + dykl*yin(475)
                                          zin(476) = zin(478) + dzkl*zin(475)
                                          ! i4 = i4 + lang+1 =  479

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  471

                                          ! nl =    2

                                          ! i4 = i3 =  471

                                          ! do nk = 1,    3

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

                                          xin(477) = xin(479) + dxkl*xin(476)
                                          yin(477) = yin(479) + dykl*yin(476)
                                          zin(477) = zin(479) + dzkl*zin(476)
                                          ! i4 = i4 + lang+1 =  480

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  472

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  481

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  492

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  491

                                          xin(492) = xin(492) + dxkl*xin(491)
                                          yin(492) = yin(492) + dykl*yin(491)
                                          zin(492) = zin(492) + dzkl*zin(491)

                                          ! i3 = i4 =  491
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  490

                                          xin(491) = xin(491) + dxkl*xin(490)
                                          yin(491) = yin(491) + dykl*yin(490)
                                          zin(491) = zin(491) + dzkl*zin(490)

                                          ! i3 = i4 =  490
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  492

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  491

                                          xin(492) = xin(492) + dxkl*xin(491)
                                          yin(492) = yin(492) + dykl*yin(491)
                                          zin(492) = zin(492) + dzkl*zin(491)

                                          ! i3 = i4 =  491
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  482

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  482

                                          ! do nk = 1,    3

                                          xin(482) = xin(484) + dxkl*xin(481)
                                          yin(482) = yin(484) + dykl*yin(481)
                                          zin(482) = zin(484) + dzkl*zin(481)
                                          ! i4 = i4 + lang+1 =  485

                                          ! nk =    2

                                          xin(485) = xin(487) + dxkl*xin(484)
                                          yin(485) = yin(487) + dykl*yin(484)
                                          zin(485) = zin(487) + dzkl*zin(484)
                                          ! i4 = i4 + lang+1 =  488

                                          ! nk =    3

                                          xin(488) = xin(490) + dxkl*xin(487)
                                          yin(488) = yin(490) + dykl*yin(487)
                                          zin(488) = zin(490) + dzkl*zin(487)
                                          ! i4 = i4 + lang+1 =  491

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  483

                                          ! nl =    2

                                          ! i4 = i3 =  483

                                          ! do nk = 1,    3

                                          xin(483) = xin(485) + dxkl*xin(482)
                                          yin(483) = yin(485) + dykl*yin(482)
                                          zin(483) = zin(485) + dzkl*zin(482)
                                          ! i4 = i4 + lang+1 =  486

                                          ! nk =    2

                                          xin(486) = xin(488) + dxkl*xin(485)
                                          yin(486) = yin(488) + dykl*yin(485)
                                          zin(486) = zin(488) + dzkl*zin(485)
                                          ! i4 = i4 + lang+1 =  489

                                          ! nk =    3

                                          xin(489) = xin(491) + dxkl*xin(488)
                                          yin(489) = yin(491) + dykl*yin(488)
                                          zin(489) = zin(491) + dzkl*zin(488)
                                          ! i4 = i4 + lang+1 =  492

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  484

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  493

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  504

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  503

                                          xin(504) = xin(504) + dxkl*xin(503)
                                          yin(504) = yin(504) + dykl*yin(503)
                                          zin(504) = zin(504) + dzkl*zin(503)

                                          ! i3 = i4 =  503
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  502

                                          xin(503) = xin(503) + dxkl*xin(502)
                                          yin(503) = yin(503) + dykl*yin(502)
                                          zin(503) = zin(503) + dzkl*zin(502)

                                          ! i3 = i4 =  502
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  504

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  503

                                          xin(504) = xin(504) + dxkl*xin(503)
                                          yin(504) = yin(504) + dykl*yin(503)
                                          zin(504) = zin(504) + dzkl*zin(503)

                                          ! i3 = i4 =  503
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  494

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  494

                                          ! do nk = 1,    3

                                          xin(494) = xin(496) + dxkl*xin(493)
                                          yin(494) = yin(496) + dykl*yin(493)
                                          zin(494) = zin(496) + dzkl*zin(493)
                                          ! i4 = i4 + lang+1 =  497

                                          ! nk =    2

                                          xin(497) = xin(499) + dxkl*xin(496)
                                          yin(497) = yin(499) + dykl*yin(496)
                                          zin(497) = zin(499) + dzkl*zin(496)
                                          ! i4 = i4 + lang+1 =  500

                                          ! nk =    3

                                          xin(500) = xin(502) + dxkl*xin(499)
                                          yin(500) = yin(502) + dykl*yin(499)
                                          zin(500) = zin(502) + dzkl*zin(499)
                                          ! i4 = i4 + lang+1 =  503

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  495

                                          ! nl =    2

                                          ! i4 = i3 =  495

                                          ! do nk = 1,    3

                                          xin(495) = xin(497) + dxkl*xin(494)
                                          yin(495) = yin(497) + dykl*yin(494)
                                          zin(495) = zin(497) + dzkl*zin(494)
                                          ! i4 = i4 + lang+1 =  498

                                          ! nk =    2

                                          xin(498) = xin(500) + dxkl*xin(497)
                                          yin(498) = yin(500) + dykl*yin(497)
                                          zin(498) = zin(500) + dzkl*zin(497)
                                          ! i4 = i4 + lang+1 =  501

                                          ! nk =    3

                                          xin(501) = xin(503) + dxkl*xin(500)
                                          yin(501) = yin(503) + dykl*yin(500)
                                          zin(501) = zin(503) + dzkl*zin(500)
                                          ! i4 = i4 + lang+1 =  504

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  496

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  505

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  505

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  516

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  515

                                          xin(516) = xin(516) + dxkl*xin(515)
                                          yin(516) = yin(516) + dykl*yin(515)
                                          zin(516) = zin(516) + dzkl*zin(515)

                                          ! i3 = i4 =  515
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  514

                                          xin(515) = xin(515) + dxkl*xin(514)
                                          yin(515) = yin(515) + dykl*yin(514)
                                          zin(515) = zin(515) + dzkl*zin(514)

                                          ! i3 = i4 =  514
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  516

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  515

                                          xin(516) = xin(516) + dxkl*xin(515)
                                          yin(516) = yin(516) + dykl*yin(515)
                                          zin(516) = zin(516) + dzkl*zin(515)

                                          ! i3 = i4 =  515
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  506

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  506

                                          ! do nk = 1,    3

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

                                          xin(512) = xin(514) + dxkl*xin(511)
                                          yin(512) = yin(514) + dykl*yin(511)
                                          zin(512) = zin(514) + dzkl*zin(511)
                                          ! i4 = i4 + lang+1 =  515

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  507

                                          ! nl =    2

                                          ! i4 = i3 =  507

                                          ! do nk = 1,    3

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

                                          xin(513) = xin(515) + dxkl*xin(512)
                                          yin(513) = yin(515) + dykl*yin(512)
                                          zin(513) = zin(515) + dzkl*zin(512)
                                          ! i4 = i4 + lang+1 =  516

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  508

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  517

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  528

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  527

                                          xin(528) = xin(528) + dxkl*xin(527)
                                          yin(528) = yin(528) + dykl*yin(527)
                                          zin(528) = zin(528) + dzkl*zin(527)

                                          ! i3 = i4 =  527
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  526

                                          xin(527) = xin(527) + dxkl*xin(526)
                                          yin(527) = yin(527) + dykl*yin(526)
                                          zin(527) = zin(527) + dzkl*zin(526)

                                          ! i3 = i4 =  526
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  528

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  527

                                          xin(528) = xin(528) + dxkl*xin(527)
                                          yin(528) = yin(528) + dykl*yin(527)
                                          zin(528) = zin(528) + dzkl*zin(527)

                                          ! i3 = i4 =  527
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  518

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  518

                                          ! do nk = 1,    3

                                          xin(518) = xin(520) + dxkl*xin(517)
                                          yin(518) = yin(520) + dykl*yin(517)
                                          zin(518) = zin(520) + dzkl*zin(517)
                                          ! i4 = i4 + lang+1 =  521

                                          ! nk =    2

                                          xin(521) = xin(523) + dxkl*xin(520)
                                          yin(521) = yin(523) + dykl*yin(520)
                                          zin(521) = zin(523) + dzkl*zin(520)
                                          ! i4 = i4 + lang+1 =  524

                                          ! nk =    3

                                          xin(524) = xin(526) + dxkl*xin(523)
                                          yin(524) = yin(526) + dykl*yin(523)
                                          zin(524) = zin(526) + dzkl*zin(523)
                                          ! i4 = i4 + lang+1 =  527

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  519

                                          ! nl =    2

                                          ! i4 = i3 =  519

                                          ! do nk = 1,    3

                                          xin(519) = xin(521) + dxkl*xin(518)
                                          yin(519) = yin(521) + dykl*yin(518)
                                          zin(519) = zin(521) + dzkl*zin(518)
                                          ! i4 = i4 + lang+1 =  522

                                          ! nk =    2

                                          xin(522) = xin(524) + dxkl*xin(521)
                                          yin(522) = yin(524) + dykl*yin(521)
                                          zin(522) = zin(524) + dzkl*zin(521)
                                          ! i4 = i4 + lang+1 =  525

                                          ! nk =    3

                                          xin(525) = xin(527) + dxkl*xin(524)
                                          yin(525) = yin(527) + dykl*yin(524)
                                          zin(525) = zin(527) + dzkl*zin(524)
                                          ! i4 = i4 + lang+1 =  528

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  520

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  529

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  540

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  539

                                          xin(540) = xin(540) + dxkl*xin(539)
                                          yin(540) = yin(540) + dykl*yin(539)
                                          zin(540) = zin(540) + dzkl*zin(539)

                                          ! i3 = i4 =  539
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  538

                                          xin(539) = xin(539) + dxkl*xin(538)
                                          yin(539) = yin(539) + dykl*yin(538)
                                          zin(539) = zin(539) + dzkl*zin(538)

                                          ! i3 = i4 =  538
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  540

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  539

                                          xin(540) = xin(540) + dxkl*xin(539)
                                          yin(540) = yin(540) + dykl*yin(539)
                                          zin(540) = zin(540) + dzkl*zin(539)

                                          ! i3 = i4 =  539
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  530

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  530

                                          ! do nk = 1,    3

                                          xin(530) = xin(532) + dxkl*xin(529)
                                          yin(530) = yin(532) + dykl*yin(529)
                                          zin(530) = zin(532) + dzkl*zin(529)
                                          ! i4 = i4 + lang+1 =  533

                                          ! nk =    2

                                          xin(533) = xin(535) + dxkl*xin(532)
                                          yin(533) = yin(535) + dykl*yin(532)
                                          zin(533) = zin(535) + dzkl*zin(532)
                                          ! i4 = i4 + lang+1 =  536

                                          ! nk =    3

                                          xin(536) = xin(538) + dxkl*xin(535)
                                          yin(536) = yin(538) + dykl*yin(535)
                                          zin(536) = zin(538) + dzkl*zin(535)
                                          ! i4 = i4 + lang+1 =  539

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  531

                                          ! nl =    2

                                          ! i4 = i3 =  531

                                          ! do nk = 1,    3

                                          xin(531) = xin(533) + dxkl*xin(530)
                                          yin(531) = yin(533) + dykl*yin(530)
                                          zin(531) = zin(533) + dzkl*zin(530)
                                          ! i4 = i4 + lang+1 =  534

                                          ! nk =    2

                                          xin(534) = xin(536) + dxkl*xin(533)
                                          yin(534) = yin(536) + dykl*yin(533)
                                          zin(534) = zin(536) + dzkl*zin(533)
                                          ! i4 = i4 + lang+1 =  537

                                          ! nk =    3

                                          xin(537) = xin(539) + dxkl*xin(536)
                                          yin(537) = yin(539) + dykl*yin(536)
                                          zin(537) = zin(539) + dzkl*zin(536)
                                          ! i4 = i4 + lang+1 =  540

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  532

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  541

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  541

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  552

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  551

                                          xin(552) = xin(552) + dxkl*xin(551)
                                          yin(552) = yin(552) + dykl*yin(551)
                                          zin(552) = zin(552) + dzkl*zin(551)

                                          ! i3 = i4 =  551
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  550

                                          xin(551) = xin(551) + dxkl*xin(550)
                                          yin(551) = yin(551) + dykl*yin(550)
                                          zin(551) = zin(551) + dzkl*zin(550)

                                          ! i3 = i4 =  550
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  552

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  551

                                          xin(552) = xin(552) + dxkl*xin(551)
                                          yin(552) = yin(552) + dykl*yin(551)
                                          zin(552) = zin(552) + dzkl*zin(551)

                                          ! i3 = i4 =  551
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  542

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  542

                                          ! do nk = 1,    3

                                          xin(542) = xin(544) + dxkl*xin(541)
                                          yin(542) = yin(544) + dykl*yin(541)
                                          zin(542) = zin(544) + dzkl*zin(541)
                                          ! i4 = i4 + lang+1 =  545

                                          ! nk =    2

                                          xin(545) = xin(547) + dxkl*xin(544)
                                          yin(545) = yin(547) + dykl*yin(544)
                                          zin(545) = zin(547) + dzkl*zin(544)
                                          ! i4 = i4 + lang+1 =  548

                                          ! nk =    3

                                          xin(548) = xin(550) + dxkl*xin(547)
                                          yin(548) = yin(550) + dykl*yin(547)
                                          zin(548) = zin(550) + dzkl*zin(547)
                                          ! i4 = i4 + lang+1 =  551

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  543

                                          ! nl =    2

                                          ! i4 = i3 =  543

                                          ! do nk = 1,    3

                                          xin(543) = xin(545) + dxkl*xin(542)
                                          yin(543) = yin(545) + dykl*yin(542)
                                          zin(543) = zin(545) + dzkl*zin(542)
                                          ! i4 = i4 + lang+1 =  546

                                          ! nk =    2

                                          xin(546) = xin(548) + dxkl*xin(545)
                                          yin(546) = yin(548) + dykl*yin(545)
                                          zin(546) = zin(548) + dzkl*zin(545)
                                          ! i4 = i4 + lang+1 =  549

                                          ! nk =    3

                                          xin(549) = xin(551) + dxkl*xin(548)
                                          yin(549) = yin(551) + dykl*yin(548)
                                          zin(549) = zin(551) + dzkl*zin(548)
                                          ! i4 = i4 + lang+1 =  552

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  544

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  553

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  564

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  563

                                          xin(564) = xin(564) + dxkl*xin(563)
                                          yin(564) = yin(564) + dykl*yin(563)
                                          zin(564) = zin(564) + dzkl*zin(563)

                                          ! i3 = i4 =  563
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  562

                                          xin(563) = xin(563) + dxkl*xin(562)
                                          yin(563) = yin(563) + dykl*yin(562)
                                          zin(563) = zin(563) + dzkl*zin(562)

                                          ! i3 = i4 =  562
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  564

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  563

                                          xin(564) = xin(564) + dxkl*xin(563)
                                          yin(564) = yin(564) + dykl*yin(563)
                                          zin(564) = zin(564) + dzkl*zin(563)

                                          ! i3 = i4 =  563
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  554

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  554

                                          ! do nk = 1,    3

                                          xin(554) = xin(556) + dxkl*xin(553)
                                          yin(554) = yin(556) + dykl*yin(553)
                                          zin(554) = zin(556) + dzkl*zin(553)
                                          ! i4 = i4 + lang+1 =  557

                                          ! nk =    2

                                          xin(557) = xin(559) + dxkl*xin(556)
                                          yin(557) = yin(559) + dykl*yin(556)
                                          zin(557) = zin(559) + dzkl*zin(556)
                                          ! i4 = i4 + lang+1 =  560

                                          ! nk =    3

                                          xin(560) = xin(562) + dxkl*xin(559)
                                          yin(560) = yin(562) + dykl*yin(559)
                                          zin(560) = zin(562) + dzkl*zin(559)
                                          ! i4 = i4 + lang+1 =  563

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  555

                                          ! nl =    2

                                          ! i4 = i3 =  555

                                          ! do nk = 1,    3

                                          xin(555) = xin(557) + dxkl*xin(554)
                                          yin(555) = yin(557) + dykl*yin(554)
                                          zin(555) = zin(557) + dzkl*zin(554)
                                          ! i4 = i4 + lang+1 =  558

                                          ! nk =    2

                                          xin(558) = xin(560) + dxkl*xin(557)
                                          yin(558) = yin(560) + dykl*yin(557)
                                          zin(558) = zin(560) + dzkl*zin(557)
                                          ! i4 = i4 + lang+1 =  561

                                          ! nk =    3

                                          xin(561) = xin(563) + dxkl*xin(560)
                                          yin(561) = yin(563) + dykl*yin(560)
                                          zin(561) = zin(563) + dzkl*zin(560)
                                          ! i4 = i4 + lang+1 =  564

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  556

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  565

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  576

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  575

                                          xin(576) = xin(576) + dxkl*xin(575)
                                          yin(576) = yin(576) + dykl*yin(575)
                                          zin(576) = zin(576) + dzkl*zin(575)

                                          ! i3 = i4 =  575
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  574

                                          xin(575) = xin(575) + dxkl*xin(574)
                                          yin(575) = yin(575) + dykl*yin(574)
                                          zin(575) = zin(575) + dzkl*zin(574)

                                          ! i3 = i4 =  574
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  576

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  575

                                          xin(576) = xin(576) + dxkl*xin(575)
                                          yin(576) = yin(576) + dykl*yin(575)
                                          zin(576) = zin(576) + dzkl*zin(575)

                                          ! i3 = i4 =  575
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  566

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  566

                                          ! do nk = 1,    3

                                          xin(566) = xin(568) + dxkl*xin(565)
                                          yin(566) = yin(568) + dykl*yin(565)
                                          zin(566) = zin(568) + dzkl*zin(565)
                                          ! i4 = i4 + lang+1 =  569

                                          ! nk =    2

                                          xin(569) = xin(571) + dxkl*xin(568)
                                          yin(569) = yin(571) + dykl*yin(568)
                                          zin(569) = zin(571) + dzkl*zin(568)
                                          ! i4 = i4 + lang+1 =  572

                                          ! nk =    3

                                          xin(572) = xin(574) + dxkl*xin(571)
                                          yin(572) = yin(574) + dykl*yin(571)
                                          zin(572) = zin(574) + dzkl*zin(571)
                                          ! i4 = i4 + lang+1 =  575

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  567

                                          ! nl =    2

                                          ! i4 = i3 =  567

                                          ! do nk = 1,    3

                                          xin(567) = xin(569) + dxkl*xin(566)
                                          yin(567) = yin(569) + dykl*yin(566)
                                          zin(567) = zin(569) + dzkl*zin(566)
                                          ! i4 = i4 + lang+1 =  570

                                          ! nk =    2

                                          xin(570) = xin(572) + dxkl*xin(569)
                                          yin(570) = yin(572) + dykl*yin(569)
                                          zin(570) = zin(572) + dzkl*zin(569)
                                          ! i4 = i4 + lang+1 =  573

                                          ! nk =    3

                                          xin(573) = xin(575) + dxkl*xin(572)
                                          yin(573) = yin(575) + dykl*yin(572)
                                          zin(573) = zin(575) + dzkl*zin(572)
                                          ! i4 = i4 + lang+1 =  576

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  568

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  577

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  577

                                          ! end do

                                          ! *** Now root =    5

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  576

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

                                          ! i1 = in(1) =  577

                                          xin(577) = 1.0_dp
                                          yin(577) = 1.0_dp
                                          zin(577) = f00

                                          ! i2 = in(2) =  613
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(613) = xc00
                                          yin(613) = yc00
                                          zin(613) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  580

                                          xin(580) = xcp00
                                          yin(580) = ycp00
                                          zin(580) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  616
                                          ! i2 =  613

                                          xin(616) = xcp00*xin(613) + cp10
                                          yin(616) = ycp00*yin(613) + cp10
                                          zin(616) = zcp00*zin(613) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  577
                                          ! i4 = i2 =  613

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  649
                                          ! i3 =  577
                                          ! i4 =  613

                                          xin(649) = c10*xin(577) + xc00*xin(613)
                                          yin(649) = c10*yin(577) + yc00*yin(613)
                                          zin(649) = c10*zin(577) + zc00*zin(613)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  652
                                          ! i5 =  649
                                          ! i4 =  613

                                          xin(652) = xcp00*xin(649) + cp10*xin(613)
                                          yin(652) = ycp00*yin(649) + cp10*yin(613)
                                          zin(652) = zcp00*zin(649) + cp10*zin(613)

                                          ! ------------------

                                          ! i3 = i4 =  613
                                          ! i4 = i5 =  649

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  685
                                          ! i3 =  613
                                          ! i4 =  649

                                          xin(685) = c10*xin(613) + xc00*xin(649)
                                          yin(685) = c10*yin(613) + yc00*yin(649)
                                          zin(685) = c10*zin(613) + zc00*zin(649)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  688
                                          ! i5 =  685
                                          ! i4 =  649

                                          xin(688) = xcp00*xin(685) + cp10*xin(649)
                                          yin(688) = ycp00*yin(685) + cp10*yin(649)
                                          zin(688) = zcp00*zin(685) + cp10*zin(649)

                                          ! ------------------

                                          ! i3 = i4 =  649
                                          ! i4 = i5 =  685

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  697
                                          ! i3 =  649
                                          ! i4 =  685

                                          xin(697) = c10*xin(649) + xc00*xin(685)
                                          yin(697) = c10*yin(649) + yc00*yin(685)
                                          zin(697) = c10*zin(649) + zc00*zin(685)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  700
                                          ! i5 =  697
                                          ! i4 =  685

                                          xin(700) = xcp00*xin(697) + cp10*xin(685)
                                          yin(700) = ycp00*yin(697) + cp10*yin(685)
                                          zin(700) = zcp00*zin(697) + cp10*zin(685)

                                          ! ------------------

                                          ! i3 = i4 =  685
                                          ! i4 = i5 =  697

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  709
                                          ! i3 =  685
                                          ! i4 =  697

                                          xin(709) = c10*xin(685) + xc00*xin(697)
                                          yin(709) = c10*yin(685) + yc00*yin(697)
                                          zin(709) = c10*zin(685) + zc00*zin(697)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  712
                                          ! i5 =  709
                                          ! i4 =  697

                                          xin(712) = xcp00*xin(709) + cp10*xin(697)
                                          yin(712) = ycp00*yin(709) + cp10*yin(697)
                                          zin(712) = zcp00*zin(709) + cp10*zin(697)

                                          ! ------------------

                                          ! i3 = i4 =  697
                                          ! i4 = i5 =  709

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =  577
                                          ! i4 = i1+k2 =  580

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  583
                                          ! i3 =  577
                                          ! i4 =  580

                                          xin(583) = cp01*xin(577) + xcp00*xin(580)
                                          yin(583) = cp01*yin(577) + ycp00*yin(580)
                                          zin(583) = cp01*zin(577) + zcp00*zin(580)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  619

                                          xin(619) = xc00*xin(583) + c01*xin(580)
                                          yin(619) = yc00*yin(583) + c01*yin(580)
                                          zin(619) = zc00*zin(583) + c01*zin(580)

                                          ! ------------------

                                          ! i3 = i4 =  580
                                          ! i4 = i5 =  583

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  586
                                          ! i3 =  580
                                          ! i4 =  583

                                          xin(586) = cp01*xin(580) + xcp00*xin(583)
                                          yin(586) = cp01*yin(580) + ycp00*yin(583)
                                          zin(586) = cp01*zin(580) + zcp00*zin(583)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  622

                                          xin(622) = xc00*xin(586) + c01*xin(583)
                                          yin(622) = yc00*yin(586) + c01*yin(583)
                                          zin(622) = zc00*zin(586) + c01*zin(583)

                                          ! ------------------

                                          ! i3 = i4 =  583
                                          ! i4 = i5 =  586

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  587
                                          ! i3 =  583
                                          ! i4 =  586

                                          xin(587) = cp01*xin(583) + xcp00*xin(586)
                                          yin(587) = cp01*yin(583) + ycp00*yin(586)
                                          zin(587) = cp01*zin(583) + zcp00*zin(586)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  623

                                          xin(623) = xc00*xin(587) + c01*xin(586)
                                          yin(623) = yc00*yin(587) + c01*yin(586)
                                          zin(623) = zc00*zin(587) + c01*zin(586)

                                          ! ------------------

                                          ! i3 = i4 =  586
                                          ! i4 = i5 =  587

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  588
                                          ! i3 =  586
                                          ! i4 =  587

                                          xin(588) = cp01*xin(586) + xcp00*xin(587)
                                          yin(588) = cp01*yin(586) + ycp00*yin(587)
                                          zin(588) = cp01*zin(586) + zcp00*zin(587)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  624

                                          xin(624) = xc00*xin(588) + c01*xin(587)
                                          yin(624) = yc00*yin(588) + c01*yin(587)
                                          zin(624) = zc00*zin(588) + c01*zin(587)

                                          ! ------------------

                                          ! i3 = i4 =  587
                                          ! i4 = i5 =  588

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =  577
                                          ! i4 = i2 =  613

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  649

                                          xin(655) = c10*xin(583) + xc00*xin(619) + c01*xin(616)
                                          yin(655) = c10*yin(583) + yc00*yin(619) + c01*yin(616)
                                          zin(655) = c10*zin(583) + zc00*zin(619) + c01*zin(616)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  613
                                          ! i4 = i5 =  649

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  685

                                          xin(691) = c10*xin(619) + xc00*xin(655) + c01*xin(652)
                                          yin(691) = c10*yin(619) + yc00*yin(655) + c01*yin(652)
                                          zin(691) = c10*zin(619) + zc00*zin(655) + c01*zin(652)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  649
                                          ! i4 = i5 =  685

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  697

                                          xin(703) = c10*xin(655) + xc00*xin(691) + c01*xin(688)
                                          yin(703) = c10*yin(655) + yc00*yin(691) + c01*yin(688)
                                          zin(703) = c10*zin(655) + zc00*zin(691) + c01*zin(688)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  685
                                          ! i4 = i5 =  697

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  709

                                          xin(715) = c10*xin(691) + xc00*xin(703) + c01*xin(700)
                                          yin(715) = c10*yin(691) + yc00*yin(703) + c01*yin(700)
                                          zin(715) = c10*zin(691) + zc00*zin(703) + c01*zin(700)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  697
                                          ! i4 = i5 =  709

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =  577
                                          ! i4 = i2 =  613

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  649

                                          xin(658) = c10*xin(586) + xc00*xin(622) + c01*xin(619)
                                          yin(658) = c10*yin(586) + yc00*yin(622) + c01*yin(619)
                                          zin(658) = c10*zin(586) + zc00*zin(622) + c01*zin(619)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  613
                                          ! i4 = i5 =  649

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  685

                                          xin(694) = c10*xin(622) + xc00*xin(658) + c01*xin(655)
                                          yin(694) = c10*yin(622) + yc00*yin(658) + c01*yin(655)
                                          zin(694) = c10*zin(622) + zc00*zin(658) + c01*zin(655)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  649
                                          ! i4 = i5 =  685

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  697

                                          xin(706) = c10*xin(658) + xc00*xin(694) + c01*xin(691)
                                          yin(706) = c10*yin(658) + yc00*yin(694) + c01*yin(691)
                                          zin(706) = c10*zin(658) + zc00*zin(694) + c01*zin(691)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  685
                                          ! i4 = i5 =  697

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  709

                                          xin(718) = c10*xin(694) + xc00*xin(706) + c01*xin(703)
                                          yin(718) = c10*yin(694) + yc00*yin(706) + c01*yin(703)
                                          zin(718) = c10*zin(694) + zc00*zin(706) + c01*zin(703)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  697
                                          ! i4 = i5 =  709

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =  577
                                          ! i4 = i2 =  613

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  649

                                          xin(659) = c10*xin(587) + xc00*xin(623) + c01*xin(622)
                                          yin(659) = c10*yin(587) + yc00*yin(623) + c01*yin(622)
                                          zin(659) = c10*zin(587) + zc00*zin(623) + c01*zin(622)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  613
                                          ! i4 = i5 =  649

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  685

                                          xin(695) = c10*xin(623) + xc00*xin(659) + c01*xin(658)
                                          yin(695) = c10*yin(623) + yc00*yin(659) + c01*yin(658)
                                          zin(695) = c10*zin(623) + zc00*zin(659) + c01*zin(658)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  649
                                          ! i4 = i5 =  685

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  697

                                          xin(707) = c10*xin(659) + xc00*xin(695) + c01*xin(694)
                                          yin(707) = c10*yin(659) + yc00*yin(695) + c01*yin(694)
                                          zin(707) = c10*zin(659) + zc00*zin(695) + c01*zin(694)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  685
                                          ! i4 = i5 =  697

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  709

                                          xin(719) = c10*xin(695) + xc00*xin(707) + c01*xin(706)
                                          yin(719) = c10*yin(695) + yc00*yin(707) + c01*yin(706)
                                          zin(719) = c10*zin(695) + zc00*zin(707) + c01*zin(706)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  697
                                          ! i4 = i5 =  709

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =  577
                                          ! i4 = i2 =  613

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  649

                                          xin(660) = c10*xin(588) + xc00*xin(624) + c01*xin(623)
                                          yin(660) = c10*yin(588) + yc00*yin(624) + c01*yin(623)
                                          zin(660) = c10*zin(588) + zc00*zin(624) + c01*zin(623)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  613
                                          ! i4 = i5 =  649

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  685

                                          xin(696) = c10*xin(624) + xc00*xin(660) + c01*xin(659)
                                          yin(696) = c10*yin(624) + yc00*yin(660) + c01*yin(659)
                                          zin(696) = c10*zin(624) + zc00*zin(660) + c01*zin(659)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  649
                                          ! i4 = i5 =  685

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  697

                                          xin(708) = c10*xin(660) + xc00*xin(696) + c01*xin(695)
                                          yin(708) = c10*yin(660) + yc00*yin(696) + c01*yin(695)
                                          zin(708) = c10*zin(660) + zc00*zin(696) + c01*zin(695)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  685
                                          ! i4 = i5 =  697

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  709

                                          xin(720) = c10*xin(696) + xc00*xin(708) + c01*xin(707)
                                          yin(720) = c10*yin(696) + yc00*yin(708) + c01*yin(707)
                                          zin(720) = c10*zin(696) + zc00*zin(708) + c01*zin(707)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  697
                                          ! i4 = i5 =  709

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  709

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  709

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  697

                                          xin(709) = xin(709) + dxij*xin(697)
                                          yin(709) = yin(709) + dyij*yin(697)
                                          zin(709) = zin(709) + dzij*zin(697)

                                          ! i3 = i4 =  697
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  685

                                          xin(697) = xin(697) + dxij*xin(685)
                                          yin(697) = yin(697) + dyij*yin(685)
                                          zin(697) = zin(697) + dzij*zin(685)

                                          ! i3 = i4 =  685
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  709

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  697

                                          xin(709) = xin(709) + dxij*xin(697)
                                          yin(709) = yin(709) + dyij*yin(697)
                                          zin(709) = zin(709) + dzij*zin(697)

                                          ! i3 = i4 =  697
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  589

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  589

                                          ! do ni = 1,    3

                                          xin(589) = xin(613) + dxij*xin(577)
                                          yin(589) = yin(613) + dyij*yin(577)
                                          zin(589) = zin(613) + dzij*zin(577)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  625

                                          ! ni =    2

                                          xin(625) = xin(649) + dxij*xin(613)
                                          yin(625) = yin(649) + dyij*yin(613)
                                          zin(625) = zin(649) + dzij*zin(613)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  661

                                          ! ni =    3

                                          xin(661) = xin(685) + dxij*xin(649)
                                          yin(661) = yin(685) + dyij*yin(649)
                                          zin(661) = zin(685) + dzij*zin(649)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  697

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  601

                                          ! nj =    2

                                          ! i4 = i3 =  601

                                          ! do ni = 1,    3

                                          xin(601) = xin(625) + dxij*xin(589)
                                          yin(601) = yin(625) + dyij*yin(589)
                                          zin(601) = zin(625) + dzij*zin(589)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  637

                                          ! ni =    2

                                          xin(637) = xin(661) + dxij*xin(625)
                                          yin(637) = yin(661) + dyij*yin(625)
                                          zin(637) = zin(661) + dzij*zin(625)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  673

                                          ! ni =    3

                                          xin(673) = xin(697) + dxij*xin(661)
                                          yin(673) = yin(697) + dyij*yin(661)
                                          zin(673) = zin(697) + dzij*zin(661)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  709

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  613

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  712

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  700

                                          xin(712) = xin(712) + dxij*xin(700)
                                          yin(712) = yin(712) + dyij*yin(700)
                                          zin(712) = zin(712) + dzij*zin(700)

                                          ! i3 = i4 =  700
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  688

                                          xin(700) = xin(700) + dxij*xin(688)
                                          yin(700) = yin(700) + dyij*yin(688)
                                          zin(700) = zin(700) + dzij*zin(688)

                                          ! i3 = i4 =  688
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  712

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  700

                                          xin(712) = xin(712) + dxij*xin(700)
                                          yin(712) = yin(712) + dyij*yin(700)
                                          zin(712) = zin(712) + dzij*zin(700)

                                          ! i3 = i4 =  700
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  592

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  592

                                          ! do ni = 1,    3

                                          xin(592) = xin(616) + dxij*xin(580)
                                          yin(592) = yin(616) + dyij*yin(580)
                                          zin(592) = zin(616) + dzij*zin(580)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  628

                                          ! ni =    2

                                          xin(628) = xin(652) + dxij*xin(616)
                                          yin(628) = yin(652) + dyij*yin(616)
                                          zin(628) = zin(652) + dzij*zin(616)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  664

                                          ! ni =    3

                                          xin(664) = xin(688) + dxij*xin(652)
                                          yin(664) = yin(688) + dyij*yin(652)
                                          zin(664) = zin(688) + dzij*zin(652)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  700

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  604

                                          ! nj =    2

                                          ! i4 = i3 =  604

                                          ! do ni = 1,    3

                                          xin(604) = xin(628) + dxij*xin(592)
                                          yin(604) = yin(628) + dyij*yin(592)
                                          zin(604) = zin(628) + dzij*zin(592)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  640

                                          ! ni =    2

                                          xin(640) = xin(664) + dxij*xin(628)
                                          yin(640) = yin(664) + dyij*yin(628)
                                          zin(640) = zin(664) + dzij*zin(628)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  676

                                          ! ni =    3

                                          xin(676) = xin(700) + dxij*xin(664)
                                          yin(676) = yin(700) + dyij*yin(664)
                                          zin(676) = zin(700) + dzij*zin(664)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  712

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  616

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  715

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  703

                                          xin(715) = xin(715) + dxij*xin(703)
                                          yin(715) = yin(715) + dyij*yin(703)
                                          zin(715) = zin(715) + dzij*zin(703)

                                          ! i3 = i4 =  703
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  691

                                          xin(703) = xin(703) + dxij*xin(691)
                                          yin(703) = yin(703) + dyij*yin(691)
                                          zin(703) = zin(703) + dzij*zin(691)

                                          ! i3 = i4 =  691
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  715

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  703

                                          xin(715) = xin(715) + dxij*xin(703)
                                          yin(715) = yin(715) + dyij*yin(703)
                                          zin(715) = zin(715) + dzij*zin(703)

                                          ! i3 = i4 =  703
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  595

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  595

                                          ! do ni = 1,    3

                                          xin(595) = xin(619) + dxij*xin(583)
                                          yin(595) = yin(619) + dyij*yin(583)
                                          zin(595) = zin(619) + dzij*zin(583)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  631

                                          ! ni =    2

                                          xin(631) = xin(655) + dxij*xin(619)
                                          yin(631) = yin(655) + dyij*yin(619)
                                          zin(631) = zin(655) + dzij*zin(619)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  667

                                          ! ni =    3

                                          xin(667) = xin(691) + dxij*xin(655)
                                          yin(667) = yin(691) + dyij*yin(655)
                                          zin(667) = zin(691) + dzij*zin(655)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  703

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  607

                                          ! nj =    2

                                          ! i4 = i3 =  607

                                          ! do ni = 1,    3

                                          xin(607) = xin(631) + dxij*xin(595)
                                          yin(607) = yin(631) + dyij*yin(595)
                                          zin(607) = zin(631) + dzij*zin(595)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  643

                                          ! ni =    2

                                          xin(643) = xin(667) + dxij*xin(631)
                                          yin(643) = yin(667) + dyij*yin(631)
                                          zin(643) = zin(667) + dzij*zin(631)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  679

                                          ! ni =    3

                                          xin(679) = xin(703) + dxij*xin(667)
                                          yin(679) = yin(703) + dyij*yin(667)
                                          zin(679) = zin(703) + dzij*zin(667)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  715

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  619

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  718

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  706

                                          xin(718) = xin(718) + dxij*xin(706)
                                          yin(718) = yin(718) + dyij*yin(706)
                                          zin(718) = zin(718) + dzij*zin(706)

                                          ! i3 = i4 =  706
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  694

                                          xin(706) = xin(706) + dxij*xin(694)
                                          yin(706) = yin(706) + dyij*yin(694)
                                          zin(706) = zin(706) + dzij*zin(694)

                                          ! i3 = i4 =  694
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  718

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  706

                                          xin(718) = xin(718) + dxij*xin(706)
                                          yin(718) = yin(718) + dyij*yin(706)
                                          zin(718) = zin(718) + dzij*zin(706)

                                          ! i3 = i4 =  706
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  598

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  598

                                          ! do ni = 1,    3

                                          xin(598) = xin(622) + dxij*xin(586)
                                          yin(598) = yin(622) + dyij*yin(586)
                                          zin(598) = zin(622) + dzij*zin(586)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  634

                                          ! ni =    2

                                          xin(634) = xin(658) + dxij*xin(622)
                                          yin(634) = yin(658) + dyij*yin(622)
                                          zin(634) = zin(658) + dzij*zin(622)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  670

                                          ! ni =    3

                                          xin(670) = xin(694) + dxij*xin(658)
                                          yin(670) = yin(694) + dyij*yin(658)
                                          zin(670) = zin(694) + dzij*zin(658)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  706

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  610

                                          ! nj =    2

                                          ! i4 = i3 =  610

                                          ! do ni = 1,    3

                                          xin(610) = xin(634) + dxij*xin(598)
                                          yin(610) = yin(634) + dyij*yin(598)
                                          zin(610) = zin(634) + dzij*zin(598)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  646

                                          ! ni =    2

                                          xin(646) = xin(670) + dxij*xin(634)
                                          yin(646) = yin(670) + dyij*yin(634)
                                          zin(646) = zin(670) + dzij*zin(634)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  682

                                          ! ni =    3

                                          xin(682) = xin(706) + dxij*xin(670)
                                          yin(682) = yin(706) + dyij*yin(670)
                                          zin(682) = zin(706) + dzij*zin(670)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  718

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  622

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  719

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  707

                                          xin(719) = xin(719) + dxij*xin(707)
                                          yin(719) = yin(719) + dyij*yin(707)
                                          zin(719) = zin(719) + dzij*zin(707)

                                          ! i3 = i4 =  707
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  695

                                          xin(707) = xin(707) + dxij*xin(695)
                                          yin(707) = yin(707) + dyij*yin(695)
                                          zin(707) = zin(707) + dzij*zin(695)

                                          ! i3 = i4 =  695
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  719

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  707

                                          xin(719) = xin(719) + dxij*xin(707)
                                          yin(719) = yin(719) + dyij*yin(707)
                                          zin(719) = zin(719) + dzij*zin(707)

                                          ! i3 = i4 =  707
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  599

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  599

                                          ! do ni = 1,    3

                                          xin(599) = xin(623) + dxij*xin(587)
                                          yin(599) = yin(623) + dyij*yin(587)
                                          zin(599) = zin(623) + dzij*zin(587)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  635

                                          ! ni =    2

                                          xin(635) = xin(659) + dxij*xin(623)
                                          yin(635) = yin(659) + dyij*yin(623)
                                          zin(635) = zin(659) + dzij*zin(623)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  671

                                          ! ni =    3

                                          xin(671) = xin(695) + dxij*xin(659)
                                          yin(671) = yin(695) + dyij*yin(659)
                                          zin(671) = zin(695) + dzij*zin(659)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  707

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  611

                                          ! nj =    2

                                          ! i4 = i3 =  611

                                          ! do ni = 1,    3

                                          xin(611) = xin(635) + dxij*xin(599)
                                          yin(611) = yin(635) + dyij*yin(599)
                                          zin(611) = zin(635) + dzij*zin(599)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  647

                                          ! ni =    2

                                          xin(647) = xin(671) + dxij*xin(635)
                                          yin(647) = yin(671) + dyij*yin(635)
                                          zin(647) = zin(671) + dzij*zin(635)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  683

                                          ! ni =    3

                                          xin(683) = xin(707) + dxij*xin(671)
                                          yin(683) = yin(707) + dyij*yin(671)
                                          zin(683) = zin(707) + dzij*zin(671)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  719

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  623

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  720

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  708

                                          xin(720) = xin(720) + dxij*xin(708)
                                          yin(720) = yin(720) + dyij*yin(708)
                                          zin(720) = zin(720) + dzij*zin(708)

                                          ! i3 = i4 =  708
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  696

                                          xin(708) = xin(708) + dxij*xin(696)
                                          yin(708) = yin(708) + dyij*yin(696)
                                          zin(708) = zin(708) + dzij*zin(696)

                                          ! i3 = i4 =  696
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  720

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  708

                                          xin(720) = xin(720) + dxij*xin(708)
                                          yin(720) = yin(720) + dyij*yin(708)
                                          zin(720) = zin(720) + dzij*zin(708)

                                          ! i3 = i4 =  708
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  600

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  600

                                          ! do ni = 1,    3

                                          xin(600) = xin(624) + dxij*xin(588)
                                          yin(600) = yin(624) + dyij*yin(588)
                                          zin(600) = zin(624) + dzij*zin(588)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  636

                                          ! ni =    2

                                          xin(636) = xin(660) + dxij*xin(624)
                                          yin(636) = yin(660) + dyij*yin(624)
                                          zin(636) = zin(660) + dzij*zin(624)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  672

                                          ! ni =    3

                                          xin(672) = xin(696) + dxij*xin(660)
                                          yin(672) = yin(696) + dyij*yin(660)
                                          zin(672) = zin(696) + dzij*zin(660)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  708

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  612

                                          ! nj =    2

                                          ! i4 = i3 =  612

                                          ! do ni = 1,    3

                                          xin(612) = xin(636) + dxij*xin(600)
                                          yin(612) = yin(636) + dyij*yin(600)
                                          zin(612) = zin(636) + dzij*zin(600)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  648

                                          ! ni =    2

                                          xin(648) = xin(672) + dxij*xin(636)
                                          yin(648) = yin(672) + dyij*yin(636)
                                          zin(648) = zin(672) + dzij*zin(636)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  684

                                          ! ni =    3

                                          xin(684) = xin(708) + dxij*xin(672)
                                          yin(684) = yin(708) + dyij*yin(672)
                                          zin(684) = zin(708) + dzij*zin(672)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  720

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  624

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =  577

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  588

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  587

                                          xin(588) = xin(588) + dxkl*xin(587)
                                          yin(588) = yin(588) + dykl*yin(587)
                                          zin(588) = zin(588) + dzkl*zin(587)

                                          ! i3 = i4 =  587
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  586

                                          xin(587) = xin(587) + dxkl*xin(586)
                                          yin(587) = yin(587) + dykl*yin(586)
                                          zin(587) = zin(587) + dzkl*zin(586)

                                          ! i3 = i4 =  586
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  588

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  587

                                          xin(588) = xin(588) + dxkl*xin(587)
                                          yin(588) = yin(588) + dykl*yin(587)
                                          zin(588) = zin(588) + dzkl*zin(587)

                                          ! i3 = i4 =  587
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  578

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  578

                                          ! do nk = 1,    3

                                          xin(578) = xin(580) + dxkl*xin(577)
                                          yin(578) = yin(580) + dykl*yin(577)
                                          zin(578) = zin(580) + dzkl*zin(577)
                                          ! i4 = i4 + lang+1 =  581

                                          ! nk =    2

                                          xin(581) = xin(583) + dxkl*xin(580)
                                          yin(581) = yin(583) + dykl*yin(580)
                                          zin(581) = zin(583) + dzkl*zin(580)
                                          ! i4 = i4 + lang+1 =  584

                                          ! nk =    3

                                          xin(584) = xin(586) + dxkl*xin(583)
                                          yin(584) = yin(586) + dykl*yin(583)
                                          zin(584) = zin(586) + dzkl*zin(583)
                                          ! i4 = i4 + lang+1 =  587

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  579

                                          ! nl =    2

                                          ! i4 = i3 =  579

                                          ! do nk = 1,    3

                                          xin(579) = xin(581) + dxkl*xin(578)
                                          yin(579) = yin(581) + dykl*yin(578)
                                          zin(579) = zin(581) + dzkl*zin(578)
                                          ! i4 = i4 + lang+1 =  582

                                          ! nk =    2

                                          xin(582) = xin(584) + dxkl*xin(581)
                                          yin(582) = yin(584) + dykl*yin(581)
                                          zin(582) = zin(584) + dzkl*zin(581)
                                          ! i4 = i4 + lang+1 =  585

                                          ! nk =    3

                                          xin(585) = xin(587) + dxkl*xin(584)
                                          yin(585) = yin(587) + dykl*yin(584)
                                          zin(585) = zin(587) + dzkl*zin(584)
                                          ! i4 = i4 + lang+1 =  588

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  580

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  589

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  600

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  599

                                          xin(600) = xin(600) + dxkl*xin(599)
                                          yin(600) = yin(600) + dykl*yin(599)
                                          zin(600) = zin(600) + dzkl*zin(599)

                                          ! i3 = i4 =  599
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  598

                                          xin(599) = xin(599) + dxkl*xin(598)
                                          yin(599) = yin(599) + dykl*yin(598)
                                          zin(599) = zin(599) + dzkl*zin(598)

                                          ! i3 = i4 =  598
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  600

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  599

                                          xin(600) = xin(600) + dxkl*xin(599)
                                          yin(600) = yin(600) + dykl*yin(599)
                                          zin(600) = zin(600) + dzkl*zin(599)

                                          ! i3 = i4 =  599
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  590

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  590

                                          ! do nk = 1,    3

                                          xin(590) = xin(592) + dxkl*xin(589)
                                          yin(590) = yin(592) + dykl*yin(589)
                                          zin(590) = zin(592) + dzkl*zin(589)
                                          ! i4 = i4 + lang+1 =  593

                                          ! nk =    2

                                          xin(593) = xin(595) + dxkl*xin(592)
                                          yin(593) = yin(595) + dykl*yin(592)
                                          zin(593) = zin(595) + dzkl*zin(592)
                                          ! i4 = i4 + lang+1 =  596

                                          ! nk =    3

                                          xin(596) = xin(598) + dxkl*xin(595)
                                          yin(596) = yin(598) + dykl*yin(595)
                                          zin(596) = zin(598) + dzkl*zin(595)
                                          ! i4 = i4 + lang+1 =  599

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  591

                                          ! nl =    2

                                          ! i4 = i3 =  591

                                          ! do nk = 1,    3

                                          xin(591) = xin(593) + dxkl*xin(590)
                                          yin(591) = yin(593) + dykl*yin(590)
                                          zin(591) = zin(593) + dzkl*zin(590)
                                          ! i4 = i4 + lang+1 =  594

                                          ! nk =    2

                                          xin(594) = xin(596) + dxkl*xin(593)
                                          yin(594) = yin(596) + dykl*yin(593)
                                          zin(594) = zin(596) + dzkl*zin(593)
                                          ! i4 = i4 + lang+1 =  597

                                          ! nk =    3

                                          xin(597) = xin(599) + dxkl*xin(596)
                                          yin(597) = yin(599) + dykl*yin(596)
                                          zin(597) = zin(599) + dzkl*zin(596)
                                          ! i4 = i4 + lang+1 =  600

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  592

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  601

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  612

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  611

                                          xin(612) = xin(612) + dxkl*xin(611)
                                          yin(612) = yin(612) + dykl*yin(611)
                                          zin(612) = zin(612) + dzkl*zin(611)

                                          ! i3 = i4 =  611
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  610

                                          xin(611) = xin(611) + dxkl*xin(610)
                                          yin(611) = yin(611) + dykl*yin(610)
                                          zin(611) = zin(611) + dzkl*zin(610)

                                          ! i3 = i4 =  610
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  612

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  611

                                          xin(612) = xin(612) + dxkl*xin(611)
                                          yin(612) = yin(612) + dykl*yin(611)
                                          zin(612) = zin(612) + dzkl*zin(611)

                                          ! i3 = i4 =  611
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  602

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  602

                                          ! do nk = 1,    3

                                          xin(602) = xin(604) + dxkl*xin(601)
                                          yin(602) = yin(604) + dykl*yin(601)
                                          zin(602) = zin(604) + dzkl*zin(601)
                                          ! i4 = i4 + lang+1 =  605

                                          ! nk =    2

                                          xin(605) = xin(607) + dxkl*xin(604)
                                          yin(605) = yin(607) + dykl*yin(604)
                                          zin(605) = zin(607) + dzkl*zin(604)
                                          ! i4 = i4 + lang+1 =  608

                                          ! nk =    3

                                          xin(608) = xin(610) + dxkl*xin(607)
                                          yin(608) = yin(610) + dykl*yin(607)
                                          zin(608) = zin(610) + dzkl*zin(607)
                                          ! i4 = i4 + lang+1 =  611

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  603

                                          ! nl =    2

                                          ! i4 = i3 =  603

                                          ! do nk = 1,    3

                                          xin(603) = xin(605) + dxkl*xin(602)
                                          yin(603) = yin(605) + dykl*yin(602)
                                          zin(603) = zin(605) + dzkl*zin(602)
                                          ! i4 = i4 + lang+1 =  606

                                          ! nk =    2

                                          xin(606) = xin(608) + dxkl*xin(605)
                                          yin(606) = yin(608) + dykl*yin(605)
                                          zin(606) = zin(608) + dzkl*zin(605)
                                          ! i4 = i4 + lang+1 =  609

                                          ! nk =    3

                                          xin(609) = xin(611) + dxkl*xin(608)
                                          yin(609) = yin(611) + dykl*yin(608)
                                          zin(609) = zin(611) + dzkl*zin(608)
                                          ! i4 = i4 + lang+1 =  612

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  604

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  613

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  613

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  624

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  623

                                          xin(624) = xin(624) + dxkl*xin(623)
                                          yin(624) = yin(624) + dykl*yin(623)
                                          zin(624) = zin(624) + dzkl*zin(623)

                                          ! i3 = i4 =  623
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  622

                                          xin(623) = xin(623) + dxkl*xin(622)
                                          yin(623) = yin(623) + dykl*yin(622)
                                          zin(623) = zin(623) + dzkl*zin(622)

                                          ! i3 = i4 =  622
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  624

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  623

                                          xin(624) = xin(624) + dxkl*xin(623)
                                          yin(624) = yin(624) + dykl*yin(623)
                                          zin(624) = zin(624) + dzkl*zin(623)

                                          ! i3 = i4 =  623
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  614

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  614

                                          ! do nk = 1,    3

                                          xin(614) = xin(616) + dxkl*xin(613)
                                          yin(614) = yin(616) + dykl*yin(613)
                                          zin(614) = zin(616) + dzkl*zin(613)
                                          ! i4 = i4 + lang+1 =  617

                                          ! nk =    2

                                          xin(617) = xin(619) + dxkl*xin(616)
                                          yin(617) = yin(619) + dykl*yin(616)
                                          zin(617) = zin(619) + dzkl*zin(616)
                                          ! i4 = i4 + lang+1 =  620

                                          ! nk =    3

                                          xin(620) = xin(622) + dxkl*xin(619)
                                          yin(620) = yin(622) + dykl*yin(619)
                                          zin(620) = zin(622) + dzkl*zin(619)
                                          ! i4 = i4 + lang+1 =  623

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  615

                                          ! nl =    2

                                          ! i4 = i3 =  615

                                          ! do nk = 1,    3

                                          xin(615) = xin(617) + dxkl*xin(614)
                                          yin(615) = yin(617) + dykl*yin(614)
                                          zin(615) = zin(617) + dzkl*zin(614)
                                          ! i4 = i4 + lang+1 =  618

                                          ! nk =    2

                                          xin(618) = xin(620) + dxkl*xin(617)
                                          yin(618) = yin(620) + dykl*yin(617)
                                          zin(618) = zin(620) + dzkl*zin(617)
                                          ! i4 = i4 + lang+1 =  621

                                          ! nk =    3

                                          xin(621) = xin(623) + dxkl*xin(620)
                                          yin(621) = yin(623) + dykl*yin(620)
                                          zin(621) = zin(623) + dzkl*zin(620)
                                          ! i4 = i4 + lang+1 =  624

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  616

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  625

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  636

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  635

                                          xin(636) = xin(636) + dxkl*xin(635)
                                          yin(636) = yin(636) + dykl*yin(635)
                                          zin(636) = zin(636) + dzkl*zin(635)

                                          ! i3 = i4 =  635
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  634

                                          xin(635) = xin(635) + dxkl*xin(634)
                                          yin(635) = yin(635) + dykl*yin(634)
                                          zin(635) = zin(635) + dzkl*zin(634)

                                          ! i3 = i4 =  634
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  636

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  635

                                          xin(636) = xin(636) + dxkl*xin(635)
                                          yin(636) = yin(636) + dykl*yin(635)
                                          zin(636) = zin(636) + dzkl*zin(635)

                                          ! i3 = i4 =  635
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  626

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  626

                                          ! do nk = 1,    3

                                          xin(626) = xin(628) + dxkl*xin(625)
                                          yin(626) = yin(628) + dykl*yin(625)
                                          zin(626) = zin(628) + dzkl*zin(625)
                                          ! i4 = i4 + lang+1 =  629

                                          ! nk =    2

                                          xin(629) = xin(631) + dxkl*xin(628)
                                          yin(629) = yin(631) + dykl*yin(628)
                                          zin(629) = zin(631) + dzkl*zin(628)
                                          ! i4 = i4 + lang+1 =  632

                                          ! nk =    3

                                          xin(632) = xin(634) + dxkl*xin(631)
                                          yin(632) = yin(634) + dykl*yin(631)
                                          zin(632) = zin(634) + dzkl*zin(631)
                                          ! i4 = i4 + lang+1 =  635

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  627

                                          ! nl =    2

                                          ! i4 = i3 =  627

                                          ! do nk = 1,    3

                                          xin(627) = xin(629) + dxkl*xin(626)
                                          yin(627) = yin(629) + dykl*yin(626)
                                          zin(627) = zin(629) + dzkl*zin(626)
                                          ! i4 = i4 + lang+1 =  630

                                          ! nk =    2

                                          xin(630) = xin(632) + dxkl*xin(629)
                                          yin(630) = yin(632) + dykl*yin(629)
                                          zin(630) = zin(632) + dzkl*zin(629)
                                          ! i4 = i4 + lang+1 =  633

                                          ! nk =    3

                                          xin(633) = xin(635) + dxkl*xin(632)
                                          yin(633) = yin(635) + dykl*yin(632)
                                          zin(633) = zin(635) + dzkl*zin(632)
                                          ! i4 = i4 + lang+1 =  636

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  628

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  637

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  648

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  647

                                          xin(648) = xin(648) + dxkl*xin(647)
                                          yin(648) = yin(648) + dykl*yin(647)
                                          zin(648) = zin(648) + dzkl*zin(647)

                                          ! i3 = i4 =  647
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  646

                                          xin(647) = xin(647) + dxkl*xin(646)
                                          yin(647) = yin(647) + dykl*yin(646)
                                          zin(647) = zin(647) + dzkl*zin(646)

                                          ! i3 = i4 =  646
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  648

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  647

                                          xin(648) = xin(648) + dxkl*xin(647)
                                          yin(648) = yin(648) + dykl*yin(647)
                                          zin(648) = zin(648) + dzkl*zin(647)

                                          ! i3 = i4 =  647
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  638

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  638

                                          ! do nk = 1,    3

                                          xin(638) = xin(640) + dxkl*xin(637)
                                          yin(638) = yin(640) + dykl*yin(637)
                                          zin(638) = zin(640) + dzkl*zin(637)
                                          ! i4 = i4 + lang+1 =  641

                                          ! nk =    2

                                          xin(641) = xin(643) + dxkl*xin(640)
                                          yin(641) = yin(643) + dykl*yin(640)
                                          zin(641) = zin(643) + dzkl*zin(640)
                                          ! i4 = i4 + lang+1 =  644

                                          ! nk =    3

                                          xin(644) = xin(646) + dxkl*xin(643)
                                          yin(644) = yin(646) + dykl*yin(643)
                                          zin(644) = zin(646) + dzkl*zin(643)
                                          ! i4 = i4 + lang+1 =  647

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  639

                                          ! nl =    2

                                          ! i4 = i3 =  639

                                          ! do nk = 1,    3

                                          xin(639) = xin(641) + dxkl*xin(638)
                                          yin(639) = yin(641) + dykl*yin(638)
                                          zin(639) = zin(641) + dzkl*zin(638)
                                          ! i4 = i4 + lang+1 =  642

                                          ! nk =    2

                                          xin(642) = xin(644) + dxkl*xin(641)
                                          yin(642) = yin(644) + dykl*yin(641)
                                          zin(642) = zin(644) + dzkl*zin(641)
                                          ! i4 = i4 + lang+1 =  645

                                          ! nk =    3

                                          xin(645) = xin(647) + dxkl*xin(644)
                                          yin(645) = yin(647) + dykl*yin(644)
                                          zin(645) = zin(647) + dzkl*zin(644)
                                          ! i4 = i4 + lang+1 =  648

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  640

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  649

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  649

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  660

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  659

                                          xin(660) = xin(660) + dxkl*xin(659)
                                          yin(660) = yin(660) + dykl*yin(659)
                                          zin(660) = zin(660) + dzkl*zin(659)

                                          ! i3 = i4 =  659
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  658

                                          xin(659) = xin(659) + dxkl*xin(658)
                                          yin(659) = yin(659) + dykl*yin(658)
                                          zin(659) = zin(659) + dzkl*zin(658)

                                          ! i3 = i4 =  658
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  660

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  659

                                          xin(660) = xin(660) + dxkl*xin(659)
                                          yin(660) = yin(660) + dykl*yin(659)
                                          zin(660) = zin(660) + dzkl*zin(659)

                                          ! i3 = i4 =  659
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  650

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  650

                                          ! do nk = 1,    3

                                          xin(650) = xin(652) + dxkl*xin(649)
                                          yin(650) = yin(652) + dykl*yin(649)
                                          zin(650) = zin(652) + dzkl*zin(649)
                                          ! i4 = i4 + lang+1 =  653

                                          ! nk =    2

                                          xin(653) = xin(655) + dxkl*xin(652)
                                          yin(653) = yin(655) + dykl*yin(652)
                                          zin(653) = zin(655) + dzkl*zin(652)
                                          ! i4 = i4 + lang+1 =  656

                                          ! nk =    3

                                          xin(656) = xin(658) + dxkl*xin(655)
                                          yin(656) = yin(658) + dykl*yin(655)
                                          zin(656) = zin(658) + dzkl*zin(655)
                                          ! i4 = i4 + lang+1 =  659

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  651

                                          ! nl =    2

                                          ! i4 = i3 =  651

                                          ! do nk = 1,    3

                                          xin(651) = xin(653) + dxkl*xin(650)
                                          yin(651) = yin(653) + dykl*yin(650)
                                          zin(651) = zin(653) + dzkl*zin(650)
                                          ! i4 = i4 + lang+1 =  654

                                          ! nk =    2

                                          xin(654) = xin(656) + dxkl*xin(653)
                                          yin(654) = yin(656) + dykl*yin(653)
                                          zin(654) = zin(656) + dzkl*zin(653)
                                          ! i4 = i4 + lang+1 =  657

                                          ! nk =    3

                                          xin(657) = xin(659) + dxkl*xin(656)
                                          yin(657) = yin(659) + dykl*yin(656)
                                          zin(657) = zin(659) + dzkl*zin(656)
                                          ! i4 = i4 + lang+1 =  660

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  652

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  661

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  672

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  671

                                          xin(672) = xin(672) + dxkl*xin(671)
                                          yin(672) = yin(672) + dykl*yin(671)
                                          zin(672) = zin(672) + dzkl*zin(671)

                                          ! i3 = i4 =  671
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  670

                                          xin(671) = xin(671) + dxkl*xin(670)
                                          yin(671) = yin(671) + dykl*yin(670)
                                          zin(671) = zin(671) + dzkl*zin(670)

                                          ! i3 = i4 =  670
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  672

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  671

                                          xin(672) = xin(672) + dxkl*xin(671)
                                          yin(672) = yin(672) + dykl*yin(671)
                                          zin(672) = zin(672) + dzkl*zin(671)

                                          ! i3 = i4 =  671
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  662

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  662

                                          ! do nk = 1,    3

                                          xin(662) = xin(664) + dxkl*xin(661)
                                          yin(662) = yin(664) + dykl*yin(661)
                                          zin(662) = zin(664) + dzkl*zin(661)
                                          ! i4 = i4 + lang+1 =  665

                                          ! nk =    2

                                          xin(665) = xin(667) + dxkl*xin(664)
                                          yin(665) = yin(667) + dykl*yin(664)
                                          zin(665) = zin(667) + dzkl*zin(664)
                                          ! i4 = i4 + lang+1 =  668

                                          ! nk =    3

                                          xin(668) = xin(670) + dxkl*xin(667)
                                          yin(668) = yin(670) + dykl*yin(667)
                                          zin(668) = zin(670) + dzkl*zin(667)
                                          ! i4 = i4 + lang+1 =  671

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  663

                                          ! nl =    2

                                          ! i4 = i3 =  663

                                          ! do nk = 1,    3

                                          xin(663) = xin(665) + dxkl*xin(662)
                                          yin(663) = yin(665) + dykl*yin(662)
                                          zin(663) = zin(665) + dzkl*zin(662)
                                          ! i4 = i4 + lang+1 =  666

                                          ! nk =    2

                                          xin(666) = xin(668) + dxkl*xin(665)
                                          yin(666) = yin(668) + dykl*yin(665)
                                          zin(666) = zin(668) + dzkl*zin(665)
                                          ! i4 = i4 + lang+1 =  669

                                          ! nk =    3

                                          xin(669) = xin(671) + dxkl*xin(668)
                                          yin(669) = yin(671) + dykl*yin(668)
                                          zin(669) = zin(671) + dzkl*zin(668)
                                          ! i4 = i4 + lang+1 =  672

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  664

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  673

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  684

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  683

                                          xin(684) = xin(684) + dxkl*xin(683)
                                          yin(684) = yin(684) + dykl*yin(683)
                                          zin(684) = zin(684) + dzkl*zin(683)

                                          ! i3 = i4 =  683
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  682

                                          xin(683) = xin(683) + dxkl*xin(682)
                                          yin(683) = yin(683) + dykl*yin(682)
                                          zin(683) = zin(683) + dzkl*zin(682)

                                          ! i3 = i4 =  682
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  684

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  683

                                          xin(684) = xin(684) + dxkl*xin(683)
                                          yin(684) = yin(684) + dykl*yin(683)
                                          zin(684) = zin(684) + dzkl*zin(683)

                                          ! i3 = i4 =  683
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  674

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  674

                                          ! do nk = 1,    3

                                          xin(674) = xin(676) + dxkl*xin(673)
                                          yin(674) = yin(676) + dykl*yin(673)
                                          zin(674) = zin(676) + dzkl*zin(673)
                                          ! i4 = i4 + lang+1 =  677

                                          ! nk =    2

                                          xin(677) = xin(679) + dxkl*xin(676)
                                          yin(677) = yin(679) + dykl*yin(676)
                                          zin(677) = zin(679) + dzkl*zin(676)
                                          ! i4 = i4 + lang+1 =  680

                                          ! nk =    3

                                          xin(680) = xin(682) + dxkl*xin(679)
                                          yin(680) = yin(682) + dykl*yin(679)
                                          zin(680) = zin(682) + dzkl*zin(679)
                                          ! i4 = i4 + lang+1 =  683

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  675

                                          ! nl =    2

                                          ! i4 = i3 =  675

                                          ! do nk = 1,    3

                                          xin(675) = xin(677) + dxkl*xin(674)
                                          yin(675) = yin(677) + dykl*yin(674)
                                          zin(675) = zin(677) + dzkl*zin(674)
                                          ! i4 = i4 + lang+1 =  678

                                          ! nk =    2

                                          xin(678) = xin(680) + dxkl*xin(677)
                                          yin(678) = yin(680) + dykl*yin(677)
                                          zin(678) = zin(680) + dzkl*zin(677)
                                          ! i4 = i4 + lang+1 =  681

                                          ! nk =    3

                                          xin(681) = xin(683) + dxkl*xin(680)
                                          yin(681) = yin(683) + dykl*yin(680)
                                          zin(681) = zin(683) + dzkl*zin(680)
                                          ! i4 = i4 + lang+1 =  684

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  676

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  685

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  685

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  696

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  695

                                          xin(696) = xin(696) + dxkl*xin(695)
                                          yin(696) = yin(696) + dykl*yin(695)
                                          zin(696) = zin(696) + dzkl*zin(695)

                                          ! i3 = i4 =  695
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  694

                                          xin(695) = xin(695) + dxkl*xin(694)
                                          yin(695) = yin(695) + dykl*yin(694)
                                          zin(695) = zin(695) + dzkl*zin(694)

                                          ! i3 = i4 =  694
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  696

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  695

                                          xin(696) = xin(696) + dxkl*xin(695)
                                          yin(696) = yin(696) + dykl*yin(695)
                                          zin(696) = zin(696) + dzkl*zin(695)

                                          ! i3 = i4 =  695
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  686

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  686

                                          ! do nk = 1,    3

                                          xin(686) = xin(688) + dxkl*xin(685)
                                          yin(686) = yin(688) + dykl*yin(685)
                                          zin(686) = zin(688) + dzkl*zin(685)
                                          ! i4 = i4 + lang+1 =  689

                                          ! nk =    2

                                          xin(689) = xin(691) + dxkl*xin(688)
                                          yin(689) = yin(691) + dykl*yin(688)
                                          zin(689) = zin(691) + dzkl*zin(688)
                                          ! i4 = i4 + lang+1 =  692

                                          ! nk =    3

                                          xin(692) = xin(694) + dxkl*xin(691)
                                          yin(692) = yin(694) + dykl*yin(691)
                                          zin(692) = zin(694) + dzkl*zin(691)
                                          ! i4 = i4 + lang+1 =  695

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  687

                                          ! nl =    2

                                          ! i4 = i3 =  687

                                          ! do nk = 1,    3

                                          xin(687) = xin(689) + dxkl*xin(686)
                                          yin(687) = yin(689) + dykl*yin(686)
                                          zin(687) = zin(689) + dzkl*zin(686)
                                          ! i4 = i4 + lang+1 =  690

                                          ! nk =    2

                                          xin(690) = xin(692) + dxkl*xin(689)
                                          yin(690) = yin(692) + dykl*yin(689)
                                          zin(690) = zin(692) + dzkl*zin(689)
                                          ! i4 = i4 + lang+1 =  693

                                          ! nk =    3

                                          xin(693) = xin(695) + dxkl*xin(692)
                                          yin(693) = yin(695) + dykl*yin(692)
                                          zin(693) = zin(695) + dzkl*zin(692)
                                          ! i4 = i4 + lang+1 =  696

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  688

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  697

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  708

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  707

                                          xin(708) = xin(708) + dxkl*xin(707)
                                          yin(708) = yin(708) + dykl*yin(707)
                                          zin(708) = zin(708) + dzkl*zin(707)

                                          ! i3 = i4 =  707
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  706

                                          xin(707) = xin(707) + dxkl*xin(706)
                                          yin(707) = yin(707) + dykl*yin(706)
                                          zin(707) = zin(707) + dzkl*zin(706)

                                          ! i3 = i4 =  706
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  708

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  707

                                          xin(708) = xin(708) + dxkl*xin(707)
                                          yin(708) = yin(708) + dykl*yin(707)
                                          zin(708) = zin(708) + dzkl*zin(707)

                                          ! i3 = i4 =  707
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  698

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  698

                                          ! do nk = 1,    3

                                          xin(698) = xin(700) + dxkl*xin(697)
                                          yin(698) = yin(700) + dykl*yin(697)
                                          zin(698) = zin(700) + dzkl*zin(697)
                                          ! i4 = i4 + lang+1 =  701

                                          ! nk =    2

                                          xin(701) = xin(703) + dxkl*xin(700)
                                          yin(701) = yin(703) + dykl*yin(700)
                                          zin(701) = zin(703) + dzkl*zin(700)
                                          ! i4 = i4 + lang+1 =  704

                                          ! nk =    3

                                          xin(704) = xin(706) + dxkl*xin(703)
                                          yin(704) = yin(706) + dykl*yin(703)
                                          zin(704) = zin(706) + dzkl*zin(703)
                                          ! i4 = i4 + lang+1 =  707

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  699

                                          ! nl =    2

                                          ! i4 = i3 =  699

                                          ! do nk = 1,    3

                                          xin(699) = xin(701) + dxkl*xin(698)
                                          yin(699) = yin(701) + dykl*yin(698)
                                          zin(699) = zin(701) + dzkl*zin(698)
                                          ! i4 = i4 + lang+1 =  702

                                          ! nk =    2

                                          xin(702) = xin(704) + dxkl*xin(701)
                                          yin(702) = yin(704) + dykl*yin(701)
                                          zin(702) = zin(704) + dzkl*zin(701)
                                          ! i4 = i4 + lang+1 =  705

                                          ! nk =    3

                                          xin(705) = xin(707) + dxkl*xin(704)
                                          yin(705) = yin(707) + dykl*yin(704)
                                          zin(705) = zin(707) + dzkl*zin(704)
                                          ! i4 = i4 + lang+1 =  708

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  700

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  709

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  720

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  719

                                          xin(720) = xin(720) + dxkl*xin(719)
                                          yin(720) = yin(720) + dykl*yin(719)
                                          zin(720) = zin(720) + dzkl*zin(719)

                                          ! i3 = i4 =  719
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  718

                                          xin(719) = xin(719) + dxkl*xin(718)
                                          yin(719) = yin(719) + dykl*yin(718)
                                          zin(719) = zin(719) + dzkl*zin(718)

                                          ! i3 = i4 =  718
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  720

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  719

                                          xin(720) = xin(720) + dxkl*xin(719)
                                          yin(720) = yin(720) + dykl*yin(719)
                                          zin(720) = zin(720) + dzkl*zin(719)

                                          ! i3 = i4 =  719
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  710

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  710

                                          ! do nk = 1,    3

                                          xin(710) = xin(712) + dxkl*xin(709)
                                          yin(710) = yin(712) + dykl*yin(709)
                                          zin(710) = zin(712) + dzkl*zin(709)
                                          ! i4 = i4 + lang+1 =  713

                                          ! nk =    2

                                          xin(713) = xin(715) + dxkl*xin(712)
                                          yin(713) = yin(715) + dykl*yin(712)
                                          zin(713) = zin(715) + dzkl*zin(712)
                                          ! i4 = i4 + lang+1 =  716

                                          ! nk =    3

                                          xin(716) = xin(718) + dxkl*xin(715)
                                          yin(716) = yin(718) + dykl*yin(715)
                                          zin(716) = zin(718) + dzkl*zin(715)
                                          ! i4 = i4 + lang+1 =  719

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  711

                                          ! nl =    2

                                          ! i4 = i3 =  711

                                          ! do nk = 1,    3

                                          xin(711) = xin(713) + dxkl*xin(710)
                                          yin(711) = yin(713) + dykl*yin(710)
                                          zin(711) = zin(713) + dzkl*zin(710)
                                          ! i4 = i4 + lang+1 =  714

                                          ! nk =    2

                                          xin(714) = xin(716) + dxkl*xin(713)
                                          yin(714) = yin(716) + dykl*yin(713)
                                          zin(714) = zin(716) + dzkl*zin(713)
                                          ! i4 = i4 + lang+1 =  717

                                          ! nk =    3

                                          xin(717) = xin(719) + dxkl*xin(716)
                                          yin(717) = yin(719) + dykl*yin(716)
                                          zin(717) = zin(719) + dzkl*zin(716)
                                          ! i4 = i4 + lang+1 =  720

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  712

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  721

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  721

                                          ! end do

                                          ! *** Now root =    6

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  720

                                          u2 = roots(6)*rho
                                          f00 = expe*wghts(6)

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

                                          ! i1 = in(1) =  721

                                          xin(721) = 1.0_dp
                                          yin(721) = 1.0_dp
                                          zin(721) = f00

                                          ! i2 = in(2) =  757
                                          ! k2 = kn(2) =    3
                                          cp10 = b00

                                          ! ----- I(1,0) -----

                                          xin(757) = xc00
                                          yin(757) = yc00
                                          zin(757) = zc00*f00

                                          ! ----- I(0,1) -----

                                          ! i3 = i1+k2 =  724

                                          xin(724) = xcp00
                                          yin(724) = ycp00
                                          zin(724) = zcp00*f00

                                          ! ----- I(1,1) -----

                                          ! i3 = i2+k2 =  760
                                          ! i2 =  757

                                          xin(760) = xcp00*xin(757) + cp10
                                          yin(760) = ycp00*yin(757) + cp10
                                          zin(760) = zcp00*zin(757) + cp10*f00

                                          ! ----- I(N,0) -----

                                          c10 = 0.0_dp

                                          ! i3 = i1 =  721
                                          ! i4 = i2 =  757

                                          ! do n = 2,   5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  793
                                          ! i3 =  721
                                          ! i4 =  757

                                          xin(793) = c10*xin(721) + xc00*xin(757)
                                          yin(793) = c10*yin(721) + yc00*yin(757)
                                          zin(793) = c10*zin(721) + zc00*zin(757)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  796
                                          ! i5 =  793
                                          ! i4 =  757

                                          xin(796) = xcp00*xin(793) + cp10*xin(757)
                                          yin(796) = ycp00*yin(793) + cp10*yin(757)
                                          zin(796) = zcp00*zin(793) + cp10*zin(757)

                                          ! ------------------

                                          ! i3 = i4 =  757
                                          ! i4 = i5 =  793

                                          ! n =    3

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  829
                                          ! i3 =  757
                                          ! i4 =  793

                                          xin(829) = c10*xin(757) + xc00*xin(793)
                                          yin(829) = c10*yin(757) + yc00*yin(793)
                                          zin(829) = c10*zin(757) + zc00*zin(793)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  832
                                          ! i5 =  829
                                          ! i4 =  793

                                          xin(832) = xcp00*xin(829) + cp10*xin(793)
                                          yin(832) = ycp00*yin(829) + cp10*yin(793)
                                          zin(832) = zcp00*zin(829) + cp10*zin(793)

                                          ! ------------------

                                          ! i3 = i4 =  793
                                          ! i4 = i5 =  829

                                          ! n =    4

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  841
                                          ! i3 =  793
                                          ! i4 =  829

                                          xin(841) = c10*xin(793) + xc00*xin(829)
                                          yin(841) = c10*yin(793) + yc00*yin(829)
                                          zin(841) = c10*zin(793) + zc00*zin(829)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  844
                                          ! i5 =  841
                                          ! i4 =  829

                                          xin(844) = xcp00*xin(841) + cp10*xin(829)
                                          yin(844) = ycp00*yin(841) + cp10*yin(829)
                                          zin(844) = zcp00*zin(841) + cp10*zin(829)

                                          ! ------------------

                                          ! i3 = i4 =  829
                                          ! i4 = i5 =  841

                                          ! n =    5

                                          c10 = c10 + b10

                                          ! i5 = in(n+1) =  853
                                          ! i3 =  829
                                          ! i4 =  841

                                          xin(853) = c10*xin(829) + xc00*xin(841)
                                          yin(853) = c10*yin(829) + yc00*yin(841)
                                          zin(853) = c10*zin(829) + zc00*zin(841)

                                          ! ----- I(N,1) -----

                                          cp10 = cp10 + b00

                                          ! i3 = i5 + k2 =  856
                                          ! i5 =  853
                                          ! i4 =  841

                                          xin(856) = xcp00*xin(853) + cp10*xin(841)
                                          yin(856) = ycp00*yin(853) + cp10*yin(841)
                                          zin(856) = zcp00*zin(853) + cp10*zin(841)

                                          ! ------------------

                                          ! i3 = i4 =  841
                                          ! i4 = i5 =  853

                                          ! n =    6

                                          ! end do

                                          ! ----- I(0,M) -----

                                          cp01 = 0.0_dp
                                          c01 = b00

                                          ! i3 = i1 =  721
                                          ! i4 = i1+k2 =  724

                                          ! do n = 2,    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  727
                                          ! i3 =  721
                                          ! i4 =  724

                                          xin(727) = cp01*xin(721) + xcp00*xin(724)
                                          yin(727) = cp01*yin(721) + ycp00*yin(724)
                                          zin(727) = cp01*zin(721) + zcp00*zin(724)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  763

                                          xin(763) = xc00*xin(727) + c01*xin(724)
                                          yin(763) = yc00*yin(727) + c01*yin(724)
                                          zin(763) = zc00*zin(727) + c01*zin(724)

                                          ! ------------------

                                          ! i3 = i4 =  724
                                          ! i4 = i5 =  727

                                          ! n =    3

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  730
                                          ! i3 =  724
                                          ! i4 =  727

                                          xin(730) = cp01*xin(724) + xcp00*xin(727)
                                          yin(730) = cp01*yin(724) + ycp00*yin(727)
                                          zin(730) = cp01*zin(724) + zcp00*zin(727)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  766

                                          xin(766) = xc00*xin(730) + c01*xin(727)
                                          yin(766) = yc00*yin(730) + c01*yin(727)
                                          zin(766) = zc00*zin(730) + c01*zin(727)

                                          ! ------------------

                                          ! i3 = i4 =  727
                                          ! i4 = i5 =  730

                                          ! n =    4

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  731
                                          ! i3 =  727
                                          ! i4 =  730

                                          xin(731) = cp01*xin(727) + xcp00*xin(730)
                                          yin(731) = cp01*yin(727) + ycp00*yin(730)
                                          zin(731) = cp01*zin(727) + zcp00*zin(730)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  767

                                          xin(767) = xc00*xin(731) + c01*xin(730)
                                          yin(767) = yc00*yin(731) + c01*yin(730)
                                          zin(767) = zc00*zin(731) + c01*zin(730)

                                          ! ------------------

                                          ! i3 = i4 =  730
                                          ! i4 = i5 =  731

                                          ! n =    5

                                          cp01 = cp01 + bp01

                                          ! i5 = i1+kn(n+1) =  732
                                          ! i3 =  730
                                          ! i4 =  731

                                          xin(732) = cp01*xin(730) + xcp00*xin(731)
                                          yin(732) = cp01*yin(730) + ycp00*yin(731)
                                          zin(732) = cp01*zin(730) + zcp00*zin(731)

                                          ! ----- I(1,M) -----

                                          c01 = c01 + b00

                                          ! i3 = i2+kn(n+1) =  768

                                          xin(768) = xc00*xin(732) + c01*xin(731)
                                          yin(768) = yc00*yin(732) + c01*yin(731)
                                          zin(768) = zc00*zin(732) + c01*zin(731)

                                          ! ------------------

                                          ! i3 = i4 =  731
                                          ! i4 = i5 =  732

                                          ! n =    6

                                          ! end do

                                          ! ----- I(N,M) -----

                                          c01 = b00
                                          ! k3 = k2 =    3

                                          ! do n = 2,    5

                                          ! k4 = kn(n+1) =    6
                                          ! i3 = i1 =  721
                                          ! i4 = i2 =  757

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  793

                                          xin(799) = c10*xin(727) + xc00*xin(763) + c01*xin(760)
                                          yin(799) = c10*yin(727) + yc00*yin(763) + c01*yin(760)
                                          zin(799) = c10*zin(727) + zc00*zin(763) + c01*zin(760)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  757
                                          ! i4 = i5 =  793

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  829

                                          xin(835) = c10*xin(763) + xc00*xin(799) + c01*xin(796)
                                          yin(835) = c10*yin(763) + yc00*yin(799) + c01*yin(796)
                                          zin(835) = c10*zin(763) + zc00*zin(799) + c01*zin(796)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  793
                                          ! i4 = i5 =  829

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  841

                                          xin(847) = c10*xin(799) + xc00*xin(835) + c01*xin(832)
                                          yin(847) = c10*yin(799) + yc00*yin(835) + c01*yin(832)
                                          zin(847) = c10*zin(799) + zc00*zin(835) + c01*zin(832)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  829
                                          ! i4 = i5 =  841

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  853

                                          xin(859) = c10*xin(835) + xc00*xin(847) + c01*xin(844)
                                          yin(859) = c10*yin(835) + yc00*yin(847) + c01*yin(844)
                                          zin(859) = c10*zin(835) + zc00*zin(847) + c01*zin(844)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  841
                                          ! i4 = i5 =  853

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   6

                                          ! n =    3

                                          ! k4 = kn(n+1) =    9
                                          ! i3 = i1 =  721
                                          ! i4 = i2 =  757

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  793

                                          xin(802) = c10*xin(730) + xc00*xin(766) + c01*xin(763)
                                          yin(802) = c10*yin(730) + yc00*yin(766) + c01*yin(763)
                                          zin(802) = c10*zin(730) + zc00*zin(766) + c01*zin(763)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  757
                                          ! i4 = i5 =  793

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  829

                                          xin(838) = c10*xin(766) + xc00*xin(802) + c01*xin(799)
                                          yin(838) = c10*yin(766) + yc00*yin(802) + c01*yin(799)
                                          zin(838) = c10*zin(766) + zc00*zin(802) + c01*zin(799)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  793
                                          ! i4 = i5 =  829

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  841

                                          xin(850) = c10*xin(802) + xc00*xin(838) + c01*xin(835)
                                          yin(850) = c10*yin(802) + yc00*yin(838) + c01*yin(835)
                                          zin(850) = c10*zin(802) + zc00*zin(838) + c01*zin(835)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  829
                                          ! i4 = i5 =  841

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  853

                                          xin(862) = c10*xin(838) + xc00*xin(850) + c01*xin(847)
                                          yin(862) = c10*yin(838) + yc00*yin(850) + c01*yin(847)
                                          zin(862) = c10*zin(838) + zc00*zin(850) + c01*zin(847)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  841
                                          ! i4 = i5 =  853

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4   9

                                          ! n =    4

                                          ! k4 = kn(n+1) =   10
                                          ! i3 = i1 =  721
                                          ! i4 = i2 =  757

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  793

                                          xin(803) = c10*xin(731) + xc00*xin(767) + c01*xin(766)
                                          yin(803) = c10*yin(731) + yc00*yin(767) + c01*yin(766)
                                          zin(803) = c10*zin(731) + zc00*zin(767) + c01*zin(766)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  757
                                          ! i4 = i5 =  793

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  829

                                          xin(839) = c10*xin(767) + xc00*xin(803) + c01*xin(802)
                                          yin(839) = c10*yin(767) + yc00*yin(803) + c01*yin(802)
                                          zin(839) = c10*zin(767) + zc00*zin(803) + c01*zin(802)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  793
                                          ! i4 = i5 =  829

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  841

                                          xin(851) = c10*xin(803) + xc00*xin(839) + c01*xin(838)
                                          yin(851) = c10*yin(803) + yc00*yin(839) + c01*yin(838)
                                          zin(851) = c10*zin(803) + zc00*zin(839) + c01*zin(838)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  829
                                          ! i4 = i5 =  841

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  853

                                          xin(863) = c10*xin(839) + xc00*xin(851) + c01*xin(850)
                                          yin(863) = c10*yin(839) + yc00*yin(851) + c01*yin(850)
                                          zin(863) = c10*zin(839) + zc00*zin(851) + c01*zin(850)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  841
                                          ! i4 = i5 =  853

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  10

                                          ! n =    5

                                          ! k4 = kn(n+1) =   11
                                          ! i3 = i1 =  721
                                          ! i4 = i2 =  757

                                          c01 = c01 + b00
                                          c10 = b10

                                          ! do nn = 2,    5

                                          ! i5 = in(nn+1) =  793

                                          xin(804) = c10*xin(732) + xc00*xin(768) + c01*xin(767)
                                          yin(804) = c10*yin(732) + yc00*yin(768) + c01*yin(767)
                                          zin(804) = c10*zin(732) + zc00*zin(768) + c01*zin(767)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  757
                                          ! i4 = i5 =  793

                                          ! nn =    3

                                          ! i5 = in(nn+1) =  829

                                          xin(840) = c10*xin(768) + xc00*xin(804) + c01*xin(803)
                                          yin(840) = c10*yin(768) + yc00*yin(804) + c01*yin(803)
                                          zin(840) = c10*zin(768) + zc00*zin(804) + c01*zin(803)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  793
                                          ! i4 = i5 =  829

                                          ! nn =    4

                                          ! i5 = in(nn+1) =  841

                                          xin(852) = c10*xin(804) + xc00*xin(840) + c01*xin(839)
                                          yin(852) = c10*yin(804) + yc00*yin(840) + c01*yin(839)
                                          zin(852) = c10*zin(804) + zc00*zin(840) + c01*zin(839)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  829
                                          ! i4 = i5 =  841

                                          ! nn =    5

                                          ! i5 = in(nn+1) =  853

                                          xin(864) = c10*xin(840) + xc00*xin(852) + c01*xin(851)
                                          yin(864) = c10*yin(840) + yc00*yin(852) + c01*yin(851)
                                          zin(864) = c10*zin(840) + zc00*zin(852) + c01*zin(851)

                                          c10 = c10 + b10

                                          ! i3 = i4 =  841
                                          ! i4 = i5 =  853

                                          ! nn =    6

                                          ! end do

                                          ! k3 = k4  11

                                          ! n =    6

                                          ! end do

                                          ! ----- I(NI,NJ,M) -----

                                          ! nm = 0
                                          ! i5 = in(iang+jang+1) =  853

                                          ! do while nm.le.(kang+lang)

                                          ! min = iang

                                          ! km = kn(nm+1) =    0

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  853

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  841

                                          xin(853) = xin(853) + dxij*xin(841)
                                          yin(853) = yin(853) + dyij*yin(841)
                                          zin(853) = zin(853) + dzij*zin(841)

                                          ! i3 = i4 =  841
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  829

                                          xin(841) = xin(841) + dxij*xin(829)
                                          yin(841) = yin(841) + dyij*yin(829)
                                          zin(841) = zin(841) + dzij*zin(829)

                                          ! i3 = i4 =  829
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  853

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  841

                                          xin(853) = xin(853) + dxij*xin(841)
                                          yin(853) = yin(853) + dyij*yin(841)
                                          zin(853) = zin(853) + dzij*zin(841)

                                          ! i3 = i4 =  841
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  733

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  733

                                          ! do ni = 1,    3

                                          xin(733) = xin(757) + dxij*xin(721)
                                          yin(733) = yin(757) + dyij*yin(721)
                                          zin(733) = zin(757) + dzij*zin(721)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  769

                                          ! ni =    2

                                          xin(769) = xin(793) + dxij*xin(757)
                                          yin(769) = yin(793) + dyij*yin(757)
                                          zin(769) = zin(793) + dzij*zin(757)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  805

                                          ! ni =    3

                                          xin(805) = xin(829) + dxij*xin(793)
                                          yin(805) = yin(829) + dyij*yin(793)
                                          zin(805) = zin(829) + dzij*zin(793)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  841

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  745

                                          ! nj =    2

                                          ! i4 = i3 =  745

                                          ! do ni = 1,    3

                                          xin(745) = xin(769) + dxij*xin(733)
                                          yin(745) = yin(769) + dyij*yin(733)
                                          zin(745) = zin(769) + dzij*zin(733)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  781

                                          ! ni =    2

                                          xin(781) = xin(805) + dxij*xin(769)
                                          yin(781) = yin(805) + dyij*yin(769)
                                          zin(781) = zin(805) + dzij*zin(769)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  817

                                          ! ni =    3

                                          xin(817) = xin(841) + dxij*xin(805)
                                          yin(817) = yin(841) + dyij*yin(805)
                                          zin(817) = zin(841) + dzij*zin(805)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  853

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  757

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    1

                                          ! min = iang

                                          ! km = kn(nm+1) =    3

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  856

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  844

                                          xin(856) = xin(856) + dxij*xin(844)
                                          yin(856) = yin(856) + dyij*yin(844)
                                          zin(856) = zin(856) + dzij*zin(844)

                                          ! i3 = i4 =  844
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  832

                                          xin(844) = xin(844) + dxij*xin(832)
                                          yin(844) = yin(844) + dyij*yin(832)
                                          zin(844) = zin(844) + dzij*zin(832)

                                          ! i3 = i4 =  832
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  856

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  844

                                          xin(856) = xin(856) + dxij*xin(844)
                                          yin(856) = yin(856) + dyij*yin(844)
                                          zin(856) = zin(856) + dzij*zin(844)

                                          ! i3 = i4 =  844
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  736

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  736

                                          ! do ni = 1,    3

                                          xin(736) = xin(760) + dxij*xin(724)
                                          yin(736) = yin(760) + dyij*yin(724)
                                          zin(736) = zin(760) + dzij*zin(724)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  772

                                          ! ni =    2

                                          xin(772) = xin(796) + dxij*xin(760)
                                          yin(772) = yin(796) + dyij*yin(760)
                                          zin(772) = zin(796) + dzij*zin(760)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  808

                                          ! ni =    3

                                          xin(808) = xin(832) + dxij*xin(796)
                                          yin(808) = yin(832) + dyij*yin(796)
                                          zin(808) = zin(832) + dzij*zin(796)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  844

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  748

                                          ! nj =    2

                                          ! i4 = i3 =  748

                                          ! do ni = 1,    3

                                          xin(748) = xin(772) + dxij*xin(736)
                                          yin(748) = yin(772) + dyij*yin(736)
                                          zin(748) = zin(772) + dzij*zin(736)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  784

                                          ! ni =    2

                                          xin(784) = xin(808) + dxij*xin(772)
                                          yin(784) = yin(808) + dyij*yin(772)
                                          zin(784) = zin(808) + dzij*zin(772)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  820

                                          ! ni =    3

                                          xin(820) = xin(844) + dxij*xin(808)
                                          yin(820) = yin(844) + dyij*yin(808)
                                          zin(820) = zin(844) + dzij*zin(808)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  856

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  760

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    2

                                          ! min = iang

                                          ! km = kn(nm+1) =    6

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  859

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  847

                                          xin(859) = xin(859) + dxij*xin(847)
                                          yin(859) = yin(859) + dyij*yin(847)
                                          zin(859) = zin(859) + dzij*zin(847)

                                          ! i3 = i4 =  847
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  835

                                          xin(847) = xin(847) + dxij*xin(835)
                                          yin(847) = yin(847) + dyij*yin(835)
                                          zin(847) = zin(847) + dzij*zin(835)

                                          ! i3 = i4 =  835
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  859

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  847

                                          xin(859) = xin(859) + dxij*xin(847)
                                          yin(859) = yin(859) + dyij*yin(847)
                                          zin(859) = zin(859) + dzij*zin(847)

                                          ! i3 = i4 =  847
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  739

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  739

                                          ! do ni = 1,    3

                                          xin(739) = xin(763) + dxij*xin(727)
                                          yin(739) = yin(763) + dyij*yin(727)
                                          zin(739) = zin(763) + dzij*zin(727)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  775

                                          ! ni =    2

                                          xin(775) = xin(799) + dxij*xin(763)
                                          yin(775) = yin(799) + dyij*yin(763)
                                          zin(775) = zin(799) + dzij*zin(763)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  811

                                          ! ni =    3

                                          xin(811) = xin(835) + dxij*xin(799)
                                          yin(811) = yin(835) + dyij*yin(799)
                                          zin(811) = zin(835) + dzij*zin(799)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  847

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  751

                                          ! nj =    2

                                          ! i4 = i3 =  751

                                          ! do ni = 1,    3

                                          xin(751) = xin(775) + dxij*xin(739)
                                          yin(751) = yin(775) + dyij*yin(739)
                                          zin(751) = zin(775) + dzij*zin(739)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  787

                                          ! ni =    2

                                          xin(787) = xin(811) + dxij*xin(775)
                                          yin(787) = yin(811) + dyij*yin(775)
                                          zin(787) = zin(811) + dzij*zin(775)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  823

                                          ! ni =    3

                                          xin(823) = xin(847) + dxij*xin(811)
                                          yin(823) = yin(847) + dyij*yin(811)
                                          zin(823) = zin(847) + dzij*zin(811)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  859

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  763

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    3

                                          ! min = iang

                                          ! km = kn(nm+1) =    9

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  862

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  850

                                          xin(862) = xin(862) + dxij*xin(850)
                                          yin(862) = yin(862) + dyij*yin(850)
                                          zin(862) = zin(862) + dzij*zin(850)

                                          ! i3 = i4 =  850
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  838

                                          xin(850) = xin(850) + dxij*xin(838)
                                          yin(850) = yin(850) + dyij*yin(838)
                                          zin(850) = zin(850) + dzij*zin(838)

                                          ! i3 = i4 =  838
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  862

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  850

                                          xin(862) = xin(862) + dxij*xin(850)
                                          yin(862) = yin(862) + dyij*yin(850)
                                          zin(862) = zin(862) + dzij*zin(850)

                                          ! i3 = i4 =  850
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  742

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  742

                                          ! do ni = 1,    3

                                          xin(742) = xin(766) + dxij*xin(730)
                                          yin(742) = yin(766) + dyij*yin(730)
                                          zin(742) = zin(766) + dzij*zin(730)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  778

                                          ! ni =    2

                                          xin(778) = xin(802) + dxij*xin(766)
                                          yin(778) = yin(802) + dyij*yin(766)
                                          zin(778) = zin(802) + dzij*zin(766)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  814

                                          ! ni =    3

                                          xin(814) = xin(838) + dxij*xin(802)
                                          yin(814) = yin(838) + dyij*yin(802)
                                          zin(814) = zin(838) + dzij*zin(802)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  850

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  754

                                          ! nj =    2

                                          ! i4 = i3 =  754

                                          ! do ni = 1,    3

                                          xin(754) = xin(778) + dxij*xin(742)
                                          yin(754) = yin(778) + dyij*yin(742)
                                          zin(754) = zin(778) + dzij*zin(742)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  790

                                          ! ni =    2

                                          xin(790) = xin(814) + dxij*xin(778)
                                          yin(790) = yin(814) + dyij*yin(778)
                                          zin(790) = zin(814) + dzij*zin(778)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  826

                                          ! ni =    3

                                          xin(826) = xin(850) + dxij*xin(814)
                                          yin(826) = yin(850) + dyij*yin(814)
                                          zin(826) = zin(850) + dzij*zin(814)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  862

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  766

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    4

                                          ! min = iang

                                          ! km = kn(nm+1) =   10

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  863

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  851

                                          xin(863) = xin(863) + dxij*xin(851)
                                          yin(863) = yin(863) + dyij*yin(851)
                                          zin(863) = zin(863) + dzij*zin(851)

                                          ! i3 = i4 =  851
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  839

                                          xin(851) = xin(851) + dxij*xin(839)
                                          yin(851) = yin(851) + dyij*yin(839)
                                          zin(851) = zin(851) + dzij*zin(839)

                                          ! i3 = i4 =  839
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  863

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  851

                                          xin(863) = xin(863) + dxij*xin(851)
                                          yin(863) = yin(863) + dyij*yin(851)
                                          zin(863) = zin(863) + dzij*zin(851)

                                          ! i3 = i4 =  851
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  743

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  743

                                          ! do ni = 1,    3

                                          xin(743) = xin(767) + dxij*xin(731)
                                          yin(743) = yin(767) + dyij*yin(731)
                                          zin(743) = zin(767) + dzij*zin(731)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  779

                                          ! ni =    2

                                          xin(779) = xin(803) + dxij*xin(767)
                                          yin(779) = yin(803) + dyij*yin(767)
                                          zin(779) = zin(803) + dzij*zin(767)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  815

                                          ! ni =    3

                                          xin(815) = xin(839) + dxij*xin(803)
                                          yin(815) = yin(839) + dyij*yin(803)
                                          zin(815) = zin(839) + dzij*zin(803)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  851

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  755

                                          ! nj =    2

                                          ! i4 = i3 =  755

                                          ! do ni = 1,    3

                                          xin(755) = xin(779) + dxij*xin(743)
                                          yin(755) = yin(779) + dyij*yin(743)
                                          zin(755) = zin(779) + dzij*zin(743)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  791

                                          ! ni =    2

                                          xin(791) = xin(815) + dxij*xin(779)
                                          yin(791) = yin(815) + dyij*yin(779)
                                          zin(791) = zin(815) + dzij*zin(779)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  827

                                          ! ni =    3

                                          xin(827) = xin(851) + dxij*xin(815)
                                          yin(827) = yin(851) + dyij*yin(815)
                                          zin(827) = zin(851) + dzij*zin(815)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  863

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  767

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    5

                                          ! min = iang

                                          ! km = kn(nm+1) =   11

                                          ! do while min.lt.(iang+jang)

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  864

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  852

                                          xin(864) = xin(864) + dxij*xin(852)
                                          yin(864) = yin(864) + dyij*yin(852)
                                          zin(864) = zin(864) + dzij*zin(852)

                                          ! i3 = i4 =  852
                                          ! nn = nn-1 =    4

                                          ! i4 = in(nn)+km =  840

                                          xin(852) = xin(852) + dxij*xin(840)
                                          yin(852) = yin(852) + dyij*yin(840)
                                          zin(852) = zin(852) + dzij*zin(840)

                                          ! i3 = i4 =  840
                                          ! nn = nn-1 =    3

                                          ! end do

                                          ! min = min + 1

                                          ! nn = (iang+jang) =    5

                                          ! i3 = i5 + km =  864

                                          ! do while nn.gt.min

                                          ! i4 = in(nn)+km =  852

                                          xin(864) = xin(864) + dxij*xin(852)
                                          yin(864) = yin(864) + dyij*yin(852)
                                          zin(864) = zin(864) + dzij*zin(852)

                                          ! i3 = i4 =  852
                                          ! nn = nn-1 =    4

                                          ! end do

                                          ! min = min + 1

                                          ! end do

                                          ! i3 = km + i1 + (kang+1)*(lang+1) =  744

                                          ! do nj = 1,    2

                                          ! i4 = i3 =  744

                                          ! do ni = 1,    3

                                          xin(744) = xin(768) + dxij*xin(732)
                                          yin(744) = yin(768) + dyij*yin(732)
                                          zin(744) = zin(768) + dzij*zin(732)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  780

                                          ! ni =    2

                                          xin(780) = xin(804) + dxij*xin(768)
                                          yin(780) = yin(804) + dyij*yin(768)
                                          zin(780) = zin(804) + dzij*zin(768)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  816

                                          ! ni =    3

                                          xin(816) = xin(840) + dxij*xin(804)
                                          yin(816) = yin(840) + dyij*yin(804)
                                          zin(816) = zin(840) + dzij*zin(804)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  852

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  756

                                          ! nj =    2

                                          ! i4 = i3 =  756

                                          ! do ni = 1,    3

                                          xin(756) = xin(780) + dxij*xin(744)
                                          yin(756) = yin(780) + dyij*yin(744)
                                          zin(756) = zin(780) + dzij*zin(744)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  792

                                          ! ni =    2

                                          xin(792) = xin(816) + dxij*xin(780)
                                          yin(792) = yin(816) + dyij*yin(780)
                                          zin(792) = zin(816) + dzij*zin(780)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  828

                                          ! ni =    3

                                          xin(828) = xin(852) + dxij*xin(816)
                                          yin(828) = yin(852) + dyij*yin(816)
                                          zin(828) = zin(852) + dzij*zin(816)

                                          ! i4 = i4 + (jang+1)*(kang+1)*(lang+1) =  864

                                          ! ni =    4

                                          ! end do

                                          ! i3 = i3 + (kang+1)*(lang+1) =  768

                                          ! nj =    3

                                          ! end do

                                          ! nm = nm + 1 =    6

                                          ! end do

                                          ! ----- I(NI,NJ,NK,NL) -----

                                          ! i5 = kn(kang+lang+1) =   11

                                          ! iaa = i1 =  721

                                          ! ni = 0

                                          ! do while ni.le.iang

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  732

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  731

                                          xin(732) = xin(732) + dxkl*xin(731)
                                          yin(732) = yin(732) + dykl*yin(731)
                                          zin(732) = zin(732) + dzkl*zin(731)

                                          ! i3 = i4 =  731
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  730

                                          xin(731) = xin(731) + dxkl*xin(730)
                                          yin(731) = yin(731) + dykl*yin(730)
                                          zin(731) = zin(731) + dzkl*zin(730)

                                          ! i3 = i4 =  730
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  732

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  731

                                          xin(732) = xin(732) + dxkl*xin(731)
                                          yin(732) = yin(732) + dykl*yin(731)
                                          zin(732) = zin(732) + dzkl*zin(731)

                                          ! i3 = i4 =  731
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  722

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  722

                                          ! do nk = 1,    3

                                          xin(722) = xin(724) + dxkl*xin(721)
                                          yin(722) = yin(724) + dykl*yin(721)
                                          zin(722) = zin(724) + dzkl*zin(721)
                                          ! i4 = i4 + lang+1 =  725

                                          ! nk =    2

                                          xin(725) = xin(727) + dxkl*xin(724)
                                          yin(725) = yin(727) + dykl*yin(724)
                                          zin(725) = zin(727) + dzkl*zin(724)
                                          ! i4 = i4 + lang+1 =  728

                                          ! nk =    3

                                          xin(728) = xin(730) + dxkl*xin(727)
                                          yin(728) = yin(730) + dykl*yin(727)
                                          zin(728) = zin(730) + dzkl*zin(727)
                                          ! i4 = i4 + lang+1 =  731

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  723

                                          ! nl =    2

                                          ! i4 = i3 =  723

                                          ! do nk = 1,    3

                                          xin(723) = xin(725) + dxkl*xin(722)
                                          yin(723) = yin(725) + dykl*yin(722)
                                          zin(723) = zin(725) + dzkl*zin(722)
                                          ! i4 = i4 + lang+1 =  726

                                          ! nk =    2

                                          xin(726) = xin(728) + dxkl*xin(725)
                                          yin(726) = yin(728) + dykl*yin(725)
                                          zin(726) = zin(728) + dzkl*zin(725)
                                          ! i4 = i4 + lang+1 =  729

                                          ! nk =    3

                                          xin(729) = xin(731) + dxkl*xin(728)
                                          yin(729) = yin(731) + dykl*yin(728)
                                          zin(729) = zin(731) + dzkl*zin(728)
                                          ! i4 = i4 + lang+1 =  732

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  724

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  733

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  744

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  743

                                          xin(744) = xin(744) + dxkl*xin(743)
                                          yin(744) = yin(744) + dykl*yin(743)
                                          zin(744) = zin(744) + dzkl*zin(743)

                                          ! i3 = i4 =  743
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  742

                                          xin(743) = xin(743) + dxkl*xin(742)
                                          yin(743) = yin(743) + dykl*yin(742)
                                          zin(743) = zin(743) + dzkl*zin(742)

                                          ! i3 = i4 =  742
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  744

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  743

                                          xin(744) = xin(744) + dxkl*xin(743)
                                          yin(744) = yin(744) + dykl*yin(743)
                                          zin(744) = zin(744) + dzkl*zin(743)

                                          ! i3 = i4 =  743
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  734

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  734

                                          ! do nk = 1,    3

                                          xin(734) = xin(736) + dxkl*xin(733)
                                          yin(734) = yin(736) + dykl*yin(733)
                                          zin(734) = zin(736) + dzkl*zin(733)
                                          ! i4 = i4 + lang+1 =  737

                                          ! nk =    2

                                          xin(737) = xin(739) + dxkl*xin(736)
                                          yin(737) = yin(739) + dykl*yin(736)
                                          zin(737) = zin(739) + dzkl*zin(736)
                                          ! i4 = i4 + lang+1 =  740

                                          ! nk =    3

                                          xin(740) = xin(742) + dxkl*xin(739)
                                          yin(740) = yin(742) + dykl*yin(739)
                                          zin(740) = zin(742) + dzkl*zin(739)
                                          ! i4 = i4 + lang+1 =  743

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  735

                                          ! nl =    2

                                          ! i4 = i3 =  735

                                          ! do nk = 1,    3

                                          xin(735) = xin(737) + dxkl*xin(734)
                                          yin(735) = yin(737) + dykl*yin(734)
                                          zin(735) = zin(737) + dzkl*zin(734)
                                          ! i4 = i4 + lang+1 =  738

                                          ! nk =    2

                                          xin(738) = xin(740) + dxkl*xin(737)
                                          yin(738) = yin(740) + dykl*yin(737)
                                          zin(738) = zin(740) + dzkl*zin(737)
                                          ! i4 = i4 + lang+1 =  741

                                          ! nk =    3

                                          xin(741) = xin(743) + dxkl*xin(740)
                                          yin(741) = yin(743) + dykl*yin(740)
                                          zin(741) = zin(743) + dzkl*zin(740)
                                          ! i4 = i4 + lang+1 =  744

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  736

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  745

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  756

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  755

                                          xin(756) = xin(756) + dxkl*xin(755)
                                          yin(756) = yin(756) + dykl*yin(755)
                                          zin(756) = zin(756) + dzkl*zin(755)

                                          ! i3 = i4 =  755
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  754

                                          xin(755) = xin(755) + dxkl*xin(754)
                                          yin(755) = yin(755) + dykl*yin(754)
                                          zin(755) = zin(755) + dzkl*zin(754)

                                          ! i3 = i4 =  754
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  756

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  755

                                          xin(756) = xin(756) + dxkl*xin(755)
                                          yin(756) = yin(756) + dykl*yin(755)
                                          zin(756) = zin(756) + dzkl*zin(755)

                                          ! i3 = i4 =  755
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  746

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  746

                                          ! do nk = 1,    3

                                          xin(746) = xin(748) + dxkl*xin(745)
                                          yin(746) = yin(748) + dykl*yin(745)
                                          zin(746) = zin(748) + dzkl*zin(745)
                                          ! i4 = i4 + lang+1 =  749

                                          ! nk =    2

                                          xin(749) = xin(751) + dxkl*xin(748)
                                          yin(749) = yin(751) + dykl*yin(748)
                                          zin(749) = zin(751) + dzkl*zin(748)
                                          ! i4 = i4 + lang+1 =  752

                                          ! nk =    3

                                          xin(752) = xin(754) + dxkl*xin(751)
                                          yin(752) = yin(754) + dykl*yin(751)
                                          zin(752) = zin(754) + dzkl*zin(751)
                                          ! i4 = i4 + lang+1 =  755

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  747

                                          ! nl =    2

                                          ! i4 = i3 =  747

                                          ! do nk = 1,    3

                                          xin(747) = xin(749) + dxkl*xin(746)
                                          yin(747) = yin(749) + dykl*yin(746)
                                          zin(747) = zin(749) + dzkl*zin(746)
                                          ! i4 = i4 + lang+1 =  750

                                          ! nk =    2

                                          xin(750) = xin(752) + dxkl*xin(749)
                                          yin(750) = yin(752) + dykl*yin(749)
                                          zin(750) = zin(752) + dzkl*zin(749)
                                          ! i4 = i4 + lang+1 =  753

                                          ! nk =    3

                                          xin(753) = xin(755) + dxkl*xin(752)
                                          yin(753) = yin(755) + dykl*yin(752)
                                          zin(753) = zin(755) + dzkl*zin(752)
                                          ! i4 = i4 + lang+1 =  756

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  748

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  757

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    1

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  757

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  768

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  767

                                          xin(768) = xin(768) + dxkl*xin(767)
                                          yin(768) = yin(768) + dykl*yin(767)
                                          zin(768) = zin(768) + dzkl*zin(767)

                                          ! i3 = i4 =  767
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  766

                                          xin(767) = xin(767) + dxkl*xin(766)
                                          yin(767) = yin(767) + dykl*yin(766)
                                          zin(767) = zin(767) + dzkl*zin(766)

                                          ! i3 = i4 =  766
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  768

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  767

                                          xin(768) = xin(768) + dxkl*xin(767)
                                          yin(768) = yin(768) + dykl*yin(767)
                                          zin(768) = zin(768) + dzkl*zin(767)

                                          ! i3 = i4 =  767
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  758

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  758

                                          ! do nk = 1,    3

                                          xin(758) = xin(760) + dxkl*xin(757)
                                          yin(758) = yin(760) + dykl*yin(757)
                                          zin(758) = zin(760) + dzkl*zin(757)
                                          ! i4 = i4 + lang+1 =  761

                                          ! nk =    2

                                          xin(761) = xin(763) + dxkl*xin(760)
                                          yin(761) = yin(763) + dykl*yin(760)
                                          zin(761) = zin(763) + dzkl*zin(760)
                                          ! i4 = i4 + lang+1 =  764

                                          ! nk =    3

                                          xin(764) = xin(766) + dxkl*xin(763)
                                          yin(764) = yin(766) + dykl*yin(763)
                                          zin(764) = zin(766) + dzkl*zin(763)
                                          ! i4 = i4 + lang+1 =  767

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  759

                                          ! nl =    2

                                          ! i4 = i3 =  759

                                          ! do nk = 1,    3

                                          xin(759) = xin(761) + dxkl*xin(758)
                                          yin(759) = yin(761) + dykl*yin(758)
                                          zin(759) = zin(761) + dzkl*zin(758)
                                          ! i4 = i4 + lang+1 =  762

                                          ! nk =    2

                                          xin(762) = xin(764) + dxkl*xin(761)
                                          yin(762) = yin(764) + dykl*yin(761)
                                          zin(762) = zin(764) + dzkl*zin(761)
                                          ! i4 = i4 + lang+1 =  765

                                          ! nk =    3

                                          xin(765) = xin(767) + dxkl*xin(764)
                                          yin(765) = yin(767) + dykl*yin(764)
                                          zin(765) = zin(767) + dzkl*zin(764)
                                          ! i4 = i4 + lang+1 =  768

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  760

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  769

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  780

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  779

                                          xin(780) = xin(780) + dxkl*xin(779)
                                          yin(780) = yin(780) + dykl*yin(779)
                                          zin(780) = zin(780) + dzkl*zin(779)

                                          ! i3 = i4 =  779
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  778

                                          xin(779) = xin(779) + dxkl*xin(778)
                                          yin(779) = yin(779) + dykl*yin(778)
                                          zin(779) = zin(779) + dzkl*zin(778)

                                          ! i3 = i4 =  778
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  780

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  779

                                          xin(780) = xin(780) + dxkl*xin(779)
                                          yin(780) = yin(780) + dykl*yin(779)
                                          zin(780) = zin(780) + dzkl*zin(779)

                                          ! i3 = i4 =  779
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  770

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  770

                                          ! do nk = 1,    3

                                          xin(770) = xin(772) + dxkl*xin(769)
                                          yin(770) = yin(772) + dykl*yin(769)
                                          zin(770) = zin(772) + dzkl*zin(769)
                                          ! i4 = i4 + lang+1 =  773

                                          ! nk =    2

                                          xin(773) = xin(775) + dxkl*xin(772)
                                          yin(773) = yin(775) + dykl*yin(772)
                                          zin(773) = zin(775) + dzkl*zin(772)
                                          ! i4 = i4 + lang+1 =  776

                                          ! nk =    3

                                          xin(776) = xin(778) + dxkl*xin(775)
                                          yin(776) = yin(778) + dykl*yin(775)
                                          zin(776) = zin(778) + dzkl*zin(775)
                                          ! i4 = i4 + lang+1 =  779

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  771

                                          ! nl =    2

                                          ! i4 = i3 =  771

                                          ! do nk = 1,    3

                                          xin(771) = xin(773) + dxkl*xin(770)
                                          yin(771) = yin(773) + dykl*yin(770)
                                          zin(771) = zin(773) + dzkl*zin(770)
                                          ! i4 = i4 + lang+1 =  774

                                          ! nk =    2

                                          xin(774) = xin(776) + dxkl*xin(773)
                                          yin(774) = yin(776) + dykl*yin(773)
                                          zin(774) = zin(776) + dzkl*zin(773)
                                          ! i4 = i4 + lang+1 =  777

                                          ! nk =    3

                                          xin(777) = xin(779) + dxkl*xin(776)
                                          yin(777) = yin(779) + dykl*yin(776)
                                          zin(777) = zin(779) + dzkl*zin(776)
                                          ! i4 = i4 + lang+1 =  780

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  772

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  781

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  792

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  791

                                          xin(792) = xin(792) + dxkl*xin(791)
                                          yin(792) = yin(792) + dykl*yin(791)
                                          zin(792) = zin(792) + dzkl*zin(791)

                                          ! i3 = i4 =  791
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  790

                                          xin(791) = xin(791) + dxkl*xin(790)
                                          yin(791) = yin(791) + dykl*yin(790)
                                          zin(791) = zin(791) + dzkl*zin(790)

                                          ! i3 = i4 =  790
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  792

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  791

                                          xin(792) = xin(792) + dxkl*xin(791)
                                          yin(792) = yin(792) + dykl*yin(791)
                                          zin(792) = zin(792) + dzkl*zin(791)

                                          ! i3 = i4 =  791
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  782

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  782

                                          ! do nk = 1,    3

                                          xin(782) = xin(784) + dxkl*xin(781)
                                          yin(782) = yin(784) + dykl*yin(781)
                                          zin(782) = zin(784) + dzkl*zin(781)
                                          ! i4 = i4 + lang+1 =  785

                                          ! nk =    2

                                          xin(785) = xin(787) + dxkl*xin(784)
                                          yin(785) = yin(787) + dykl*yin(784)
                                          zin(785) = zin(787) + dzkl*zin(784)
                                          ! i4 = i4 + lang+1 =  788

                                          ! nk =    3

                                          xin(788) = xin(790) + dxkl*xin(787)
                                          yin(788) = yin(790) + dykl*yin(787)
                                          zin(788) = zin(790) + dzkl*zin(787)
                                          ! i4 = i4 + lang+1 =  791

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  783

                                          ! nl =    2

                                          ! i4 = i3 =  783

                                          ! do nk = 1,    3

                                          xin(783) = xin(785) + dxkl*xin(782)
                                          yin(783) = yin(785) + dykl*yin(782)
                                          zin(783) = zin(785) + dzkl*zin(782)
                                          ! i4 = i4 + lang+1 =  786

                                          ! nk =    2

                                          xin(786) = xin(788) + dxkl*xin(785)
                                          yin(786) = yin(788) + dykl*yin(785)
                                          zin(786) = zin(788) + dzkl*zin(785)
                                          ! i4 = i4 + lang+1 =  789

                                          ! nk =    3

                                          xin(789) = xin(791) + dxkl*xin(788)
                                          yin(789) = yin(791) + dykl*yin(788)
                                          zin(789) = zin(791) + dzkl*zin(788)
                                          ! i4 = i4 + lang+1 =  792

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  784

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  793

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    2

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  793

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  804

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  803

                                          xin(804) = xin(804) + dxkl*xin(803)
                                          yin(804) = yin(804) + dykl*yin(803)
                                          zin(804) = zin(804) + dzkl*zin(803)

                                          ! i3 = i4 =  803
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  802

                                          xin(803) = xin(803) + dxkl*xin(802)
                                          yin(803) = yin(803) + dykl*yin(802)
                                          zin(803) = zin(803) + dzkl*zin(802)

                                          ! i3 = i4 =  802
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  804

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  803

                                          xin(804) = xin(804) + dxkl*xin(803)
                                          yin(804) = yin(804) + dykl*yin(803)
                                          zin(804) = zin(804) + dzkl*zin(803)

                                          ! i3 = i4 =  803
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  794

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  794

                                          ! do nk = 1,    3

                                          xin(794) = xin(796) + dxkl*xin(793)
                                          yin(794) = yin(796) + dykl*yin(793)
                                          zin(794) = zin(796) + dzkl*zin(793)
                                          ! i4 = i4 + lang+1 =  797

                                          ! nk =    2

                                          xin(797) = xin(799) + dxkl*xin(796)
                                          yin(797) = yin(799) + dykl*yin(796)
                                          zin(797) = zin(799) + dzkl*zin(796)
                                          ! i4 = i4 + lang+1 =  800

                                          ! nk =    3

                                          xin(800) = xin(802) + dxkl*xin(799)
                                          yin(800) = yin(802) + dykl*yin(799)
                                          zin(800) = zin(802) + dzkl*zin(799)
                                          ! i4 = i4 + lang+1 =  803

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  795

                                          ! nl =    2

                                          ! i4 = i3 =  795

                                          ! do nk = 1,    3

                                          xin(795) = xin(797) + dxkl*xin(794)
                                          yin(795) = yin(797) + dykl*yin(794)
                                          zin(795) = zin(797) + dzkl*zin(794)
                                          ! i4 = i4 + lang+1 =  798

                                          ! nk =    2

                                          xin(798) = xin(800) + dxkl*xin(797)
                                          yin(798) = yin(800) + dykl*yin(797)
                                          zin(798) = zin(800) + dzkl*zin(797)
                                          ! i4 = i4 + lang+1 =  801

                                          ! nk =    3

                                          xin(801) = xin(803) + dxkl*xin(800)
                                          yin(801) = yin(803) + dykl*yin(800)
                                          zin(801) = zin(803) + dzkl*zin(800)
                                          ! i4 = i4 + lang+1 =  804

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  796

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  805

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  816

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  815

                                          xin(816) = xin(816) + dxkl*xin(815)
                                          yin(816) = yin(816) + dykl*yin(815)
                                          zin(816) = zin(816) + dzkl*zin(815)

                                          ! i3 = i4 =  815
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  814

                                          xin(815) = xin(815) + dxkl*xin(814)
                                          yin(815) = yin(815) + dykl*yin(814)
                                          zin(815) = zin(815) + dzkl*zin(814)

                                          ! i3 = i4 =  814
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  816

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  815

                                          xin(816) = xin(816) + dxkl*xin(815)
                                          yin(816) = yin(816) + dykl*yin(815)
                                          zin(816) = zin(816) + dzkl*zin(815)

                                          ! i3 = i4 =  815
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  806

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  806

                                          ! do nk = 1,    3

                                          xin(806) = xin(808) + dxkl*xin(805)
                                          yin(806) = yin(808) + dykl*yin(805)
                                          zin(806) = zin(808) + dzkl*zin(805)
                                          ! i4 = i4 + lang+1 =  809

                                          ! nk =    2

                                          xin(809) = xin(811) + dxkl*xin(808)
                                          yin(809) = yin(811) + dykl*yin(808)
                                          zin(809) = zin(811) + dzkl*zin(808)
                                          ! i4 = i4 + lang+1 =  812

                                          ! nk =    3

                                          xin(812) = xin(814) + dxkl*xin(811)
                                          yin(812) = yin(814) + dykl*yin(811)
                                          zin(812) = zin(814) + dzkl*zin(811)
                                          ! i4 = i4 + lang+1 =  815

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  807

                                          ! nl =    2

                                          ! i4 = i3 =  807

                                          ! do nk = 1,    3

                                          xin(807) = xin(809) + dxkl*xin(806)
                                          yin(807) = yin(809) + dykl*yin(806)
                                          zin(807) = zin(809) + dzkl*zin(806)
                                          ! i4 = i4 + lang+1 =  810

                                          ! nk =    2

                                          xin(810) = xin(812) + dxkl*xin(809)
                                          yin(810) = yin(812) + dykl*yin(809)
                                          zin(810) = zin(812) + dzkl*zin(809)
                                          ! i4 = i4 + lang+1 =  813

                                          ! nk =    3

                                          xin(813) = xin(815) + dxkl*xin(812)
                                          yin(813) = yin(815) + dykl*yin(812)
                                          zin(813) = zin(815) + dzkl*zin(812)
                                          ! i4 = i4 + lang+1 =  816

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  808

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  817

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  828

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  827

                                          xin(828) = xin(828) + dxkl*xin(827)
                                          yin(828) = yin(828) + dykl*yin(827)
                                          zin(828) = zin(828) + dzkl*zin(827)

                                          ! i3 = i4 =  827
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  826

                                          xin(827) = xin(827) + dxkl*xin(826)
                                          yin(827) = yin(827) + dykl*yin(826)
                                          zin(827) = zin(827) + dzkl*zin(826)

                                          ! i3 = i4 =  826
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  828

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  827

                                          xin(828) = xin(828) + dxkl*xin(827)
                                          yin(828) = yin(828) + dykl*yin(827)
                                          zin(828) = zin(828) + dzkl*zin(827)

                                          ! i3 = i4 =  827
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  818

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  818

                                          ! do nk = 1,    3

                                          xin(818) = xin(820) + dxkl*xin(817)
                                          yin(818) = yin(820) + dykl*yin(817)
                                          zin(818) = zin(820) + dzkl*zin(817)
                                          ! i4 = i4 + lang+1 =  821

                                          ! nk =    2

                                          xin(821) = xin(823) + dxkl*xin(820)
                                          yin(821) = yin(823) + dykl*yin(820)
                                          zin(821) = zin(823) + dzkl*zin(820)
                                          ! i4 = i4 + lang+1 =  824

                                          ! nk =    3

                                          xin(824) = xin(826) + dxkl*xin(823)
                                          yin(824) = yin(826) + dykl*yin(823)
                                          zin(824) = zin(826) + dzkl*zin(823)
                                          ! i4 = i4 + lang+1 =  827

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  819

                                          ! nl =    2

                                          ! i4 = i3 =  819

                                          ! do nk = 1,    3

                                          xin(819) = xin(821) + dxkl*xin(818)
                                          yin(819) = yin(821) + dykl*yin(818)
                                          zin(819) = zin(821) + dzkl*zin(818)
                                          ! i4 = i4 + lang+1 =  822

                                          ! nk =    2

                                          xin(822) = xin(824) + dxkl*xin(821)
                                          yin(822) = yin(824) + dykl*yin(821)
                                          zin(822) = zin(824) + dzkl*zin(821)
                                          ! i4 = i4 + lang+1 =  825

                                          ! nk =    3

                                          xin(825) = xin(827) + dxkl*xin(824)
                                          yin(825) = yin(827) + dykl*yin(824)
                                          zin(825) = zin(827) + dzkl*zin(824)
                                          ! i4 = i4 + lang+1 =  828

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  820

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  829

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    3

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  829

                                          ! nj = 0

                                          ! ib = iaa

                                          ! do while nj.le.jang

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  840

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  839

                                          xin(840) = xin(840) + dxkl*xin(839)
                                          yin(840) = yin(840) + dykl*yin(839)
                                          zin(840) = zin(840) + dzkl*zin(839)

                                          ! i3 = i4 =  839
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  838

                                          xin(839) = xin(839) + dxkl*xin(838)
                                          yin(839) = yin(839) + dykl*yin(838)
                                          zin(839) = zin(839) + dzkl*zin(838)

                                          ! i3 = i4 =  838
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  840

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  839

                                          xin(840) = xin(840) + dxkl*xin(839)
                                          yin(840) = yin(840) + dykl*yin(839)
                                          zin(840) = zin(840) + dzkl*zin(839)

                                          ! i3 = i4 =  839
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  830

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  830

                                          ! do nk = 1,    3

                                          xin(830) = xin(832) + dxkl*xin(829)
                                          yin(830) = yin(832) + dykl*yin(829)
                                          zin(830) = zin(832) + dzkl*zin(829)
                                          ! i4 = i4 + lang+1 =  833

                                          ! nk =    2

                                          xin(833) = xin(835) + dxkl*xin(832)
                                          yin(833) = yin(835) + dykl*yin(832)
                                          zin(833) = zin(835) + dzkl*zin(832)
                                          ! i4 = i4 + lang+1 =  836

                                          ! nk =    3

                                          xin(836) = xin(838) + dxkl*xin(835)
                                          yin(836) = yin(838) + dykl*yin(835)
                                          zin(836) = zin(838) + dzkl*zin(835)
                                          ! i4 = i4 + lang+1 =  839

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  831

                                          ! nl =    2

                                          ! i4 = i3 =  831

                                          ! do nk = 1,    3

                                          xin(831) = xin(833) + dxkl*xin(830)
                                          yin(831) = yin(833) + dykl*yin(830)
                                          zin(831) = zin(833) + dzkl*zin(830)
                                          ! i4 = i4 + lang+1 =  834

                                          ! nk =    2

                                          xin(834) = xin(836) + dxkl*xin(833)
                                          yin(834) = yin(836) + dykl*yin(833)
                                          zin(834) = zin(836) + dzkl*zin(833)
                                          ! i4 = i4 + lang+1 =  837

                                          ! nk =    3

                                          xin(837) = xin(839) + dxkl*xin(836)
                                          yin(837) = yin(839) + dykl*yin(836)
                                          zin(837) = zin(839) + dzkl*zin(836)
                                          ! i4 = i4 + lang+1 =  840

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  832

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  841

                                          ! nj = nj + 1 =    1

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  852

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  851

                                          xin(852) = xin(852) + dxkl*xin(851)
                                          yin(852) = yin(852) + dykl*yin(851)
                                          zin(852) = zin(852) + dzkl*zin(851)

                                          ! i3 = i4 =  851
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  850

                                          xin(851) = xin(851) + dxkl*xin(850)
                                          yin(851) = yin(851) + dykl*yin(850)
                                          zin(851) = zin(851) + dzkl*zin(850)

                                          ! i3 = i4 =  850
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  852

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  851

                                          xin(852) = xin(852) + dxkl*xin(851)
                                          yin(852) = yin(852) + dykl*yin(851)
                                          zin(852) = zin(852) + dzkl*zin(851)

                                          ! i3 = i4 =  851
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  842

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  842

                                          ! do nk = 1,    3

                                          xin(842) = xin(844) + dxkl*xin(841)
                                          yin(842) = yin(844) + dykl*yin(841)
                                          zin(842) = zin(844) + dzkl*zin(841)
                                          ! i4 = i4 + lang+1 =  845

                                          ! nk =    2

                                          xin(845) = xin(847) + dxkl*xin(844)
                                          yin(845) = yin(847) + dykl*yin(844)
                                          zin(845) = zin(847) + dzkl*zin(844)
                                          ! i4 = i4 + lang+1 =  848

                                          ! nk =    3

                                          xin(848) = xin(850) + dxkl*xin(847)
                                          yin(848) = yin(850) + dykl*yin(847)
                                          zin(848) = zin(850) + dzkl*zin(847)
                                          ! i4 = i4 + lang+1 =  851

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  843

                                          ! nl =    2

                                          ! i4 = i3 =  843

                                          ! do nk = 1,    3

                                          xin(843) = xin(845) + dxkl*xin(842)
                                          yin(843) = yin(845) + dykl*yin(842)
                                          zin(843) = zin(845) + dzkl*zin(842)
                                          ! i4 = i4 + lang+1 =  846

                                          ! nk =    2

                                          xin(846) = xin(848) + dxkl*xin(845)
                                          yin(846) = yin(848) + dykl*yin(845)
                                          zin(846) = zin(848) + dzkl*zin(845)
                                          ! i4 = i4 + lang+1 =  849

                                          ! nk =    3

                                          xin(849) = xin(851) + dxkl*xin(848)
                                          yin(849) = yin(851) + dykl*yin(848)
                                          zin(849) = zin(851) + dzkl*zin(848)
                                          ! i4 = i4 + lang+1 =  852

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  844

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  853

                                          ! nj = nj + 1 =    2

                                          ! min = kang

                                          ! do while min.lt.(kang+lang)

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  864

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  863

                                          xin(864) = xin(864) + dxkl*xin(863)
                                          yin(864) = yin(864) + dykl*yin(863)
                                          zin(864) = zin(864) + dzkl*zin(863)

                                          ! i3 = i4 =  863
                                          ! nm = nm -1 =    4

                                          ! i4 = ib+kn(nm) =  862

                                          xin(863) = xin(863) + dxkl*xin(862)
                                          yin(863) = yin(863) + dykl*yin(862)
                                          zin(863) = zin(863) + dzkl*zin(862)

                                          ! i3 = i4 =  862
                                          ! nm = nm -1 =    3

                                          ! end do

                                          ! min = min + 1 =    4

                                          ! nm = (kang+lang) =    5

                                          ! i3 = ib+i5 =  864

                                          ! do while nm.gt.min

                                          ! i4 = ib+kn(nm) =  863

                                          xin(864) = xin(864) + dxkl*xin(863)
                                          yin(864) = yin(864) + dykl*yin(863)
                                          zin(864) = zin(864) + dzkl*zin(863)

                                          ! i3 = i4 =  863
                                          ! nm = nm -1 =    4

                                          ! end do

                                          ! min = min + 1 =    5

                                          ! end do

                                          ! i3 = ib + 1 =  854

                                          ! do nl = 1,    2

                                          ! i4 = i3 =  854

                                          ! do nk = 1,    3

                                          xin(854) = xin(856) + dxkl*xin(853)
                                          yin(854) = yin(856) + dykl*yin(853)
                                          zin(854) = zin(856) + dzkl*zin(853)
                                          ! i4 = i4 + lang+1 =  857

                                          ! nk =    2

                                          xin(857) = xin(859) + dxkl*xin(856)
                                          yin(857) = yin(859) + dykl*yin(856)
                                          zin(857) = zin(859) + dzkl*zin(856)
                                          ! i4 = i4 + lang+1 =  860

                                          ! nk =    3

                                          xin(860) = xin(862) + dxkl*xin(859)
                                          yin(860) = yin(862) + dykl*yin(859)
                                          zin(860) = zin(862) + dzkl*zin(859)
                                          ! i4 = i4 + lang+1 =  863

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  855

                                          ! nl =    2

                                          ! i4 = i3 =  855

                                          ! do nk = 1,    3

                                          xin(855) = xin(857) + dxkl*xin(854)
                                          yin(855) = yin(857) + dykl*yin(854)
                                          zin(855) = zin(857) + dzkl*zin(854)
                                          ! i4 = i4 + lang+1 =  858

                                          ! nk =    2

                                          xin(858) = xin(860) + dxkl*xin(857)
                                          yin(858) = yin(860) + dykl*yin(857)
                                          zin(858) = zin(860) + dzkl*zin(857)
                                          ! i4 = i4 + lang+1 =  861

                                          ! nk =    3

                                          xin(861) = xin(863) + dxkl*xin(860)
                                          yin(861) = yin(863) + dykl*yin(860)
                                          zin(861) = zin(863) + dzkl*zin(860)
                                          ! i4 = i4 + lang+1 =  864

                                          ! nk =    4

                                          ! end do

                                          ! i3 = i3 + 1 =  856

                                          ! nl =    3

                                          ! end do

                                          ! ib = ib + (kang+1)*(lang+1) =  865

                                          ! nj = nj + 1 =    3

                                          ! end do

                                          ! ni = ni + 1 =    4

                                          ! iaa = iaa + (jang+1)*(kang+1)*(lang+1) =  865

                                          ! end do

                                          ! *** Now root =    7

                                          ! mm = mm + (iang+1)*(jang+2)*(kang+1)*(lang+1) =  864

                                          !                     --- END XYZINT ---

                                          !                       --- FORMS ---
                                          ! Form final integrals adding the 2D auxiliaries over all roots

                                          j = 1

                                          do n = 1, 3600! loop over all integrals

                                            l = n - 60*(j - 1) ! index for the ket cartesian pair

                                            mx = ijx(j) + klx(l)
                                            my = ijy(j) + kly(l)
                                            mz = ijz(j) + klz(l)

                                            eri_value(n) = eri_value(n) + d23bra(j)*d23ket(l)* &
                                                           (xin(mx)*yin(my)*zin(mz) & ! root  1
                                                            + xin(mx + 144)*yin(my + 144)*zin(mz + 144) & ! root  2
                                                            + xin(mx + 288)*yin(my + 288)*zin(mz + 288) & ! root  3
                                                            + xin(mx + 432)*yin(my + 432)*zin(mz + 432) & ! root  4
                                                            + xin(mx + 576)*yin(my + 576)*zin(mz + 576) & ! root  5
                                                            + xin(mx + 720)*yin(my + 720)*zin(mz + 720)) ! root  6

                                            j = int(n/60) + 1 ! index for the next bra cartesian pair

                                          end do

                                          !                     --- END FORMS ---

                                        end do ! ij primitve loop

                                      end do ! kl primitve loop

                                      !                     --- DIRFCK_RHF ---
                                      !          Compute Fock matrix elements from 2EIs

                                      maxl = 6
                                      same = (ish .eq. ksh) .and. (jsh .eq. lsh)

                                      loci = res%atom_loc(ish) - 1
                                      locj = res%atom_loc(jsh) - 1
                                      lock = res%atom_loc(ksh) - 1
                                      locl = res%atom_loc(lsh) - 1

                                      nij = 0

                                      do i = 1, 10 ! # of cartesians in i

                                        ii1 = i + loci
                                        ip = (i - 1)*360 ! Stride between functions in i

                                        do j = 1, 6 ! # of cartesians in j

                                          nij = nij + 1

                                          maxl2 = maxl

                                          jj1 = j + locj
                                          i2 = ii1
                                          j2 = jj1
                                          if (ii1 .lt. jj1) then ! Sort <ij|
                                            i2 = jj1
                                            j2 = ii1
                                          end if

                                          ijp = (j - 1)*60 + ip ! Add stride between functions in j

                                          nkl = nij

                                          do k = 1, 10 ! # of cartesians in k

                                            kk1 = k + lock

                                            ijkp = (k - 1)*6 + ijp ! Add stride between functions in k

                                            if (same) then ! Account for non-unique permutations
                                              itmp = min(maxl2 - 1 + 1, nkl)
                                              if (itmp .eq. 0) exit ! Move to the next j iteration
                                              maxl2 = 1 + itmp - 1
                                              nkl = nkl - itmp
                                            end if

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

                                    end do ! kkll - or single loop for Do Concurrent

                                end do ! iijj
                                !$omp end target teams distribute parallel do


                                deallocate (n23bra)
                                deallocate (xint23bra)

                                end subroutine int3232
                                end submodule
