!! Shell pair precomputation for efficient integral evaluation
module liberi_shell_pair
  !! Precomputes and stores shell-pair data for all angular momentum combinations.
  !!
  !! This module is responsible for the critical preprocessing step that
  !! significantly accelerates integral evaluation. For each unique shell pair
  !! (ss, sp, pp, sd, pd, dd, sf, pf, df, ff), it precomputes:
  !!
  !! - Combined exponents and their inverses
  !! - Contraction coefficient products
  !! - Schwarz integral estimates for screening
  !! - Location/offset arrays for accessing primitive pair data
  !!
  !! ### Usage
  !!
  !! Call `shell_pair(basis, pairs)` with a populated basis_t. This must be
  !! done before calling the integral driver.
  !!
  !! ```fortran
  !! type(basis_t) :: basis
  !! type(shell_pair_container_t) :: pairs
  !! ! ... populate basis ...
  !! call shell_pair(basis, pairs)  ! Precompute all shell-pair data
  !! call compute_integrals(pairs, density, fock)  ! Now ready to compute
  !! ```
  !!
  !! The shell-pair data is returned in the `pairs` container, which holds
  !! `shell_pair_t` instances for each angular momentum combination
  !! (e.g., `pairs%ss`, `pairs%sp`, `pairs%pp`, etc.).
  use parameters, only: cux, twopi_5_2, sqrt3
  use liberi_types, only: int64, dp, basis_t, shell_pair_container_t, eri_resources_t
  use boys, only: boys_grid_zero, exponent_grid
  implicit none
  private
  public :: shell_pair, shell_pair_deallocate

contains

  subroutine shell_pair(basis, pairs)
    !! Main entry point for shell-pair precomputation.
    !!
    !! Classifies shells by angular momentum and computes shell-pair data
    !! for all relevant combinations based on which shell types are present
    !! in the basis set.
    implicit none
    type(basis_t), intent(in) :: basis
      !! Basis set data
    type(shell_pair_container_t), intent(inout) :: pairs
      !! Container for all shell-pair types (resources%n_rank/n_size must be set)
    real(dp), allocatable :: cmax(:)

    call setup_shell_coordinates(basis, pairs%resources)
    call classify_shells(basis, cmax, pairs%resources)

    if (pairs%resources%n_s_shl > 0) call compute_ss_pairs(basis, cmax, pairs)

    if (pairs%resources%n_p_shl > 0) then
      call compute_sp_pairs(basis, cmax, pairs)
      call compute_pp_pairs(basis, cmax, pairs)
    end if

    if (pairs%resources%n_d_shl > 0) then
      call compute_sd_pairs(basis, cmax, pairs)
      call compute_pd_pairs(basis, cmax, pairs)
      call compute_dd_pairs(basis, cmax, pairs)
    end if

    if (allocated(cmax)) deallocate (cmax)

    if (pairs%resources%n_f_shl > 0) then
      call compute_sf_pairs(basis, pairs)
      call compute_pf_pairs(basis, pairs)
      call compute_df_pairs(basis, pairs)
      call compute_ff_pairs(basis, pairs)
    end if

    ! Populate remaining resources from basis
    pairs%resources%contr_num = basis%contr_num(1:basis%nsh)
    pairs%resources%atom_loc = basis%atom_loc(1:basis%nsh)
    pairs%resources%num_bas = basis%num_bas

#ifdef ENABLE_GPU
    ! Map resources to GPU
    !$omp target enter data &
    !$omp map(to: pairs%resources, pairs%resources%coord_sh, &
    !$omp   pairs%resources%ia, pairs%resources%contr_num, &
    !$omp   pairs%resources%atom_loc)

    if (pairs%resources%n_s_shl > 0) then
      !$omp target enter data map(to: pairs%resources%i_s_shl)
    end if
    if (pairs%resources%n_p_shl > 0) then
      !$omp target enter data map(to: pairs%resources%i_p_shl)
    end if
    if (pairs%resources%n_d_shl > 0) then
      !$omp target enter data map(to: pairs%resources%i_d_shl)
    end if
    if (pairs%resources%n_f_shl > 0) then
      !$omp target enter data map(to: pairs%resources%i_f_shl)
    end if

    ! Map Boys function grids to GPU (static data, persists until deallocate)
    !$omp target enter data map(to: boys_grid_zero, exponent_grid)
#endif

  end subroutine shell_pair

  subroutine setup_shell_coordinates(basis, res)
    !! Sets up shell coordinates and triangular index array.
    implicit none
    type(basis_t), intent(in) :: basis
    type(eri_resources_t), intent(inout) :: res
    real(dp) :: coord(3, basis%natoms)
    integer(kind=int64) :: ii, i, iatom

    if (.not. allocated(res%coord_sh)) allocate (res%coord_sh(basis%nsh, 3))
    if (.not. allocated(res%ia)) allocate (res%ia(basis%num_bas))

    do iatom = 1, basis%natoms
      coord(1, iatom) = basis%coords(3*(iatom - 1) + 1)
      coord(2, iatom) = basis%coords(3*(iatom - 1) + 2)
      coord(3, iatom) = basis%coords(3*(iatom - 1) + 3)
    end do
    do ii = 1, basis%nsh
      i = basis%atom_num(ii)
      res%coord_sh(ii, 1) = coord(1, i)
      res%coord_sh(ii, 2) = coord(2, i)
      res%coord_sh(ii, 3) = coord(3, i)
    end do

    do ii = 1, basis%num_bas
      res%ia(ii) = (ii*ii - ii)/2
    end do

  end subroutine setup_shell_coordinates

  subroutine classify_shells(basis, cmax, res)
    !! Classifies shells by angular momentum and sorts by contraction level.
    implicit none
    type(basis_t), intent(in) :: basis
    real(dp), allocatable, intent(out) :: cmax(:)
    type(eri_resources_t), intent(inout) :: res

    integer(kind=int64), allocatable :: buffer(:)
    integer(kind=int64) :: isc(10), ipc(10), idc(10), ifc(10)
    integer(kind=int64) :: isc_buffer(10), ipc_buff(10)
    integer(kind=int64) :: idc_buff(10), ifc_buff(10)
    integer(kind=int64) :: nsc, npc, ndc, nfc
    integer(kind=int64) :: ii, jj, ish, jsh, iang, ii_buff
    integer(kind=int64) :: l, n

    res%n_s_shl = 0
    res%n_p_shl = 0
    res%n_d_shl = 0
    res%n_f_shl = 0
    do ii = 1, basis%nsh
      iang = basis%ang_mom(ii)
      if (iang == 1) then
        res%n_s_shl = res%n_s_shl + 1
      else if (iang == 2) then
        res%n_p_shl = res%n_p_shl + 1
      else if (iang == 3) then
        res%n_d_shl = res%n_d_shl + 1
      else if (iang == 4) then
        res%n_f_shl = res%n_f_shl + 1
      end if
    end do

    if (.not. allocated(res%i_s_shl)) allocate (res%i_s_shl(1:res%n_s_shl))
    if (.not. allocated(res%i_p_shl)) allocate (res%i_p_shl(1:res%n_p_shl))
    if (.not. allocated(res%i_d_shl)) allocate (res%i_d_shl(1:res%n_d_shl))
    if (.not. allocated(res%i_f_shl)) allocate (res%i_f_shl(1:res%n_f_shl))

    res%n_s_shl = 0
    res%n_p_shl = 0
    res%n_d_shl = 0
    res%n_f_shl = 0
    do ii = 1, basis%nsh
      iang = basis%ang_mom(ii)
      if (iang == 1) then
        res%n_s_shl = res%n_s_shl + 1
        res%i_s_shl(res%n_s_shl) = ii
      else if (iang == 2) then
        res%n_p_shl = res%n_p_shl + 1
        res%i_p_shl(res%n_p_shl) = ii
      else if (iang == 3) then
        res%n_d_shl = res%n_d_shl + 1
        res%i_d_shl(res%n_d_shl) = ii
      else if (iang == 4) then
        res%n_f_shl = res%n_f_shl + 1
        res%i_f_shl(res%n_f_shl) = ii
      end if
    end do

    ! Get contraction level of shells
    nsc = 0
    npc = 0
    ndc = 0
    nfc = 0
    do ii = maxval(basis%contr_num(1:basis%nsh)), minval(basis%contr_num(1:basis%nsh)), -1
      do ish = 1, res%n_s_shl
        if (basis%contr_num(res%i_s_shl(ish)) == ii) then
          nsc = nsc + 1
          isc(nsc) = ii
          exit
        end if
      end do
      do ish = 1, res%n_p_shl
        if (basis%contr_num(res%i_p_shl(ish)) == ii) then
          npc = npc + 1
          ipc(npc) = ii
          exit
        end if
      end do
      do ish = 1, res%n_d_shl
        if (basis%contr_num(res%i_d_shl(ish)) == ii) then
          ndc = ndc + 1
          idc(ndc) = ii
          exit
        end if
      end do
      do ish = 1, res%n_f_shl
        if (basis%contr_num(res%i_f_shl(ish)) == ii) then
          nfc = nfc + 1
          ifc(nfc) = ii
          exit
        end if
      end do
    end do

    ! Sort shell arrays by contraction level
    if (.not. allocated(buffer)) allocate (buffer(res%n_s_shl))
    ii_buff = 1
    do ii = 1, nsc
      isc_buffer(ii) = ii_buff
      do ish = 1, res%n_s_shl
        if (basis%contr_num(res%i_s_shl(ish)) == isc(ii)) then
          buffer(ii_buff) = res%i_s_shl(ish)
          ii_buff = ii_buff + 1
        end if
      end do
    end do
    res%i_s_shl = buffer
    if (allocated(buffer)) deallocate (buffer)

    if (.not. allocated(buffer)) allocate (buffer(res%n_p_shl))
    ii_buff = 1
    do ii = 1, npc
      ipc_buff(ii) = ii_buff
      do ish = 1, res%n_p_shl
        if (basis%contr_num(res%i_p_shl(ish)) == ipc(ii)) then
          buffer(ii_buff) = res%i_p_shl(ish)
          ii_buff = ii_buff + 1
        end if
      end do
    end do
    res%i_p_shl = buffer
    if (allocated(buffer)) deallocate (buffer)

    if (.not. allocated(buffer)) allocate (buffer(res%n_d_shl))
    ii_buff = 1
    do ii = 1, ndc
      idc_buff(ii) = ii_buff
      do ish = 1, res%n_d_shl
        if (basis%contr_num(res%i_d_shl(ish)) == idc(ii)) then
          buffer(ii_buff) = res%i_d_shl(ish)
          ii_buff = ii_buff + 1
        end if
      end do
    end do
    res%i_d_shl = buffer
    if (allocated(buffer)) deallocate (buffer)

    if (.not. allocated(buffer)) allocate (buffer(res%n_f_shl))
    ii_buff = 1
    do ii = 1, nfc
      ifc_buff(ii) = ii_buff
      do ish = 1, res%n_f_shl
        if (basis%contr_num(res%i_f_shl(ish)) == ifc(ii)) then
          buffer(ii_buff) = res%i_f_shl(ish)
          ii_buff = ii_buff + 1
        end if
      end do
    end do
    res%i_f_shl = buffer
    if (allocated(buffer)) deallocate (buffer)

    ! Find maximum contraction coefficient for each primitive
    if (.not. allocated(cmax)) allocate (cmax(basis%nsh*basis%mxgtot))
    do ii = 1, basis%nsh
      l = basis%sh_loc(ii)
      n = l + basis%contr_num(ii) - 1
      do jj = l, n
        cmax(jj) = max(abs(basis%contr_coef_s(jj)), &
                       abs(basis%contr_coef_p(jj)), &
                       abs(basis%contr_coef_d(jj)), &
                       abs(basis%contr_coef_f(jj)))
      end do
    end do

    ! Sort shell arrays in descending order
    do ii = 1, res%n_s_shl
      do jj = 1, ii
        ish = res%i_s_shl(ii)
        jsh = res%i_s_shl(jj)
        if (ish < jsh) then
          res%i_s_shl(ii) = jsh
          res%i_s_shl(jj) = ish
        end if
      end do
    end do

    do ii = 1, res%n_p_shl
      do jj = 1, ii
        ish = res%i_p_shl(ii)
        jsh = res%i_p_shl(jj)
        if (ish < jsh) then
          res%i_p_shl(ii) = jsh
          res%i_p_shl(jj) = ish
        end if
      end do
    end do

    do ii = 1, res%n_d_shl
      do jj = 1, ii
        ish = res%i_d_shl(ii)
        jsh = res%i_d_shl(jj)
        if (ish < jsh) then
          res%i_d_shl(ii) = jsh
          res%i_d_shl(jj) = ish
        end if
      end do
    end do

    do ii = 1, res%n_f_shl
      do jj = 1, ii
        ish = res%i_f_shl(ii)
        jsh = res%i_f_shl(jj)
        if (ish < jsh) then
          res%i_f_shl(ii) = jsh
          res%i_f_shl(jj) = ish
        end if
      end do
    end do

  end subroutine classify_shells

  subroutine compute_ss_pairs(basis, cmax, pairs)
    !! Computes ss shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_ss
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)

    allocate (pairs%ss%xints(pairs%resources%n_s_shl*(pairs%resources%n_s_shl + 1)/2))
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_s_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%ss%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%ss%pair_loc(pairs%resources%n_s_shl*(pairs%resources%n_s_shl + 1)/2 + 1))
    pairs%ss%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_s_shl(jj)
        pairs%ss%pair_loc(iijj + 1) = pairs%ss%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_ss = pairs%ss%pair_loc(pairs%resources%n_s_shl*(pairs%resources%n_s_shl + 1)/2 + 1)
    pairs%ss%n_pairs = n_ss
    allocate (pairs%ss%sq(n_ss))
    allocate (ismlp(n_ss))
    allocate (pairs%ss%t_expon_ab(n_ss))
    allocate (pairs%ss%expon_a(n_ss))
    allocate (pairs%ss%expon_b(n_ss))
    allocate (pairs%ss%ismlp(n_ss))
    allocate (pairs%ss%d_coeff(n_ss))
    allocate (pairs%ss%d_coeff_alt(n_ss))
    allocate (pairs%ss%t_inverse_expon_ab(n_ss))
    allocate (pairs%ss%t_alpha(n_ss))
    allocate (pairs%ss%t_beta(n_ss))

    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, ii

        iijj = jj + ii*(ii - 1)/2

        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_s_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            pairs%ss%expon_a(pairs%ss%pair_loc(iijj) + ji) = buff(1)
            pairs%ss%expon_b(pairs%ss%pair_loc(iijj) + ji) = buff(2)
            pairs%ss%t_expon_ab(pairs%ss%pair_loc(iijj) + ji) = buff(3)
            pairs%ss%t_inverse_expon_ab(pairs%ss%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%ss%t_alpha(pairs%ss%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%ss%t_beta(pairs%ss%pair_loc(iijj) + ji) = pairs%ss%t_alpha(pairs%ss%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%ss%sq(pairs%ss%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_s(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)
            pairs%ss%d_coeff(pairs%ss%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_s(basis%sh_loc(jsh) + j - 1)

            pairs%ss%d_coeff_alt(pairs%ss%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_s(basis%sh_loc(jsh) + j - 1)

          end do
        end do
        pairs%ss%ismlp(pairs%ss%pair_loc(iijj) + 1:pairs%ss%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do

    if (allocated(ismlp)) deallocate (ismlp)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%ss, pairs%ss%d_coeff, pairs%ss%d_coeff_alt, pairs%ss%expon_a, pairs%ss%expon_b, pairs%ss%ismlp, pairs%ss%pair_loc, pairs%ss%sq, pairs%ss%t_alpha, pairs%ss%t_beta, pairs%ss%t_expon_ab, pairs%ss%t_inverse_expon_ab, pairs%ss%xints)
#endif

  end subroutine compute_ss_pairs

  subroutine compute_sp_pairs(basis, cmax, pairs)
    !! Computes sp shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_sp
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)

    allocate (pairs%sp%xints(pairs%resources%n_s_shl*pairs%resources%n_p_shl))
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_p_shl
        iijj = jj + (ii - 1)*pairs%resources%n_p_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%sp%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%sp%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_p_shl + 1))
    pairs%sp%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_p_shl
        iijj = jj + (ii - 1)*pairs%resources%n_p_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)
        pairs%sp%pair_loc(iijj + 1) = pairs%sp%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_sp = pairs%sp%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_p_shl + 1)
    pairs%sp%n_pairs = n_sp
    allocate (pairs%sp%sq(n_sp))
    allocate (ismlp(n_sp))
    allocate (pairs%sp%expon_a(n_sp))
    allocate (pairs%sp%expon_b(n_sp))
    allocate (pairs%sp%ismlp(n_sp))
    allocate (pairs%sp%d_coeff(n_sp))
    allocate (pairs%sp%d_coeff_alt(n_sp))
    allocate (pairs%sp%t_expon_ab(n_sp))
    allocate (pairs%sp%t_inverse_expon_ab(n_sp))
    allocate (pairs%sp%t_alpha(n_sp))
    allocate (pairs%sp%t_beta(n_sp))

    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_p_shl

        iijj = jj + (ii - 1)*pairs%resources%n_p_shl

        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            pairs%sp%expon_a(pairs%sp%pair_loc(iijj) + ji) = buff(1)
            pairs%sp%expon_b(pairs%sp%pair_loc(iijj) + ji) = buff(2)
            pairs%sp%t_expon_ab(pairs%sp%pair_loc(iijj) + ji) = buff(3)
            pairs%sp%t_inverse_expon_ab(pairs%sp%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%sp%t_alpha(pairs%sp%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%sp%t_beta(pairs%sp%pair_loc(iijj) + ji) = pairs%sp%t_alpha(pairs%sp%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%sp%sq(pairs%sp%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)
            pairs%sp%d_coeff(pairs%sp%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

            pairs%sp%d_coeff_alt(pairs%sp%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

          end do
        end do
        pairs%sp%ismlp(pairs%sp%pair_loc(iijj) + 1:pairs%sp%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do
    if (allocated(ismlp)) deallocate (ismlp)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%sp,  pairs%sp%d_coeff, pairs%sp%d_coeff_alt, pairs%sp%expon_a, pairs%sp%expon_b, pairs%sp%ismlp, pairs%sp%pair_loc, pairs%sp%sq, pairs%sp%t_alpha, pairs%sp%t_beta, pairs%sp%t_expon_ab, pairs%sp%t_inverse_expon_ab, pairs%sp%xints)
#endif

  end subroutine compute_sp_pairs

  subroutine compute_pp_pairs(basis, cmax, pairs)
    !! Computes pp shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, lgmax, n_pp
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)
    logical :: iandj, double
    ! Rys-specific locals (not stored in shell_pair_t)
    integer(kind=int64) :: n_pp_rys
    integer(kind=int64), allocatable :: locpp_rys(:)
    real(dp), allocatable :: expon_a_pp_rys(:), expon_b_pp_rys(:)
    real(dp), allocatable :: t_expon_ab_pp_rys(:)
    real(dp), allocatable :: d11p_pp_rys(:, :)

    allocate (pairs%pp%xints(pairs%resources%n_p_shl*(pairs%resources%n_p_shl + 1)/2))
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%pp%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%pp%pair_loc(pairs%resources%n_p_shl*(pairs%resources%n_p_shl + 1)/2 + 1))
    pairs%pp%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)

        pairs%pp%pair_loc(iijj + 1) = pairs%pp%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    allocate (locpp_rys(pairs%resources%n_p_shl*(pairs%resources%n_p_shl + 1)/2 + 1))
    locpp_rys(1) = 0
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)

        iandj = ish == jsh
        ji = 0
        lgmax = basis%contr_num(jsh)
        do i = 1, basis%contr_num(ish)
          if (iandj) lgmax = i
          do j = 1, lgmax
            ji = ji + 1
          end do
        end do
        locpp_rys(iijj + 1) = locpp_rys(iijj) + ji
      end do
    end do

    n_pp = pairs%pp%pair_loc(pairs%resources%n_p_shl*(pairs%resources%n_p_shl + 1)/2 + 1)
    pairs%pp%n_pairs = n_pp
    allocate (pairs%pp%sq(n_pp))
    allocate (pairs%pp%expon_a(n_pp))
    allocate (pairs%pp%expon_b(n_pp))
    allocate (ismlp(n_pp))
    allocate (pairs%pp%t_expon_ab(n_pp))
    allocate (pairs%pp%ismlp(n_pp))
    allocate (pairs%pp%d_coeff(n_pp))
    allocate (pairs%pp%d_coeff_alt(n_pp))
    allocate (pairs%pp%t_inverse_expon_ab(n_pp))
    allocate (pairs%pp%t_alpha(n_pp))
    allocate (pairs%pp%t_beta(n_pp))

    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, ii

        iijj = jj + ii*(ii - 1)/2

        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        lgmax = basis%contr_num(jsh)
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            pairs%pp%expon_a(pairs%pp%pair_loc(iijj) + ji) = buff(1)
            pairs%pp%expon_b(pairs%pp%pair_loc(iijj) + ji) = buff(2)
            pairs%pp%t_expon_ab(pairs%pp%pair_loc(iijj) + ji) = buff(3)
            pairs%pp%t_inverse_expon_ab(pairs%pp%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%pp%t_alpha(pairs%pp%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%pp%t_beta(pairs%pp%pair_loc(iijj) + ji) = pairs%pp%t_alpha(pairs%pp%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%pp%sq(pairs%pp%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)
            pairs%pp%d_coeff(pairs%pp%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

            pairs%pp%d_coeff_alt(pairs%pp%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

          end do
        end do
        pairs%pp%ismlp(pairs%pp%pair_loc(iijj) + 1:pairs%pp%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do

    if (allocated(ismlp)) deallocate (ismlp)

    n_pp_rys = locpp_rys(pairs%resources%n_p_shl*(pairs%resources%n_p_shl + 1)/2 + 1)
    allocate (expon_a_pp_rys(n_pp_rys))
    allocate (expon_b_pp_rys(n_pp_rys))
    allocate (t_expon_ab_pp_rys(n_pp_rys))
    allocate (d11p_pp_rys(3, n_pp_rys))

    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, ii

        iijj = jj + ii*(ii - 1)/2

        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_p_shl(jj)

        iandj = ish == jsh
        ji = 0
        lgmax = basis%contr_num(jsh)
        do i = 1, basis%contr_num(ish)

          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          if (iandj) lgmax = i

          do j = 1, lgmax

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            expon_a_pp_rys(locpp_rys(iijj) + ji) = buff(1)  ! exponent a
            expon_b_pp_rys(locpp_rys(iijj) + ji) = buff(2)  ! exponent b
            t_expon_ab_pp_rys(locpp_rys(iijj) + ji) = buff(3)  ! exponent_a_b

            double = iandj .and. (i /= j)

            d11p_pp_rys(1, locpp_rys(iijj) + ji) = buff(4)* &
                                                   basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                   basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)
            if (double) then

              d11p_pp_rys(1, locpp_rys(iijj) + ji) = buff(4)* &
                                                     basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                     basis%contr_coef_p(basis%sh_loc(jsh) + j - 1) + &
                                                     buff(4)* &
                                                     basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                     basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

            end if

            d11p_pp_rys(2, locpp_rys(iijj) + ji) = sqrt3*buff(4)* &
                                                   basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                   basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)
            d11p_pp_rys(3, locpp_rys(iijj) + ji) = sqrt3*sqrt3*buff(4)* &
                                                   basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                   basis%contr_coef_p(basis%sh_loc(jsh) + j - 1)

          end do
        end do
      end do
    end do

    if (allocated(locpp_rys)) deallocate (locpp_rys)
    if (allocated(expon_a_pp_rys)) deallocate (expon_a_pp_rys)
    if (allocated(expon_b_pp_rys)) deallocate (expon_b_pp_rys)
    if (allocated(t_expon_ab_pp_rys)) deallocate (t_expon_ab_pp_rys)
    if (allocated(d11p_pp_rys)) deallocate (d11p_pp_rys)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%pp, pairs%pp%d_coeff, pairs%pp%d_coeff_alt, pairs%pp%expon_a, pairs%pp%expon_b, pairs%pp%ismlp, pairs%pp%pair_loc, pairs%pp%sq, pairs%pp%t_alpha, pairs%pp%t_beta, pairs%pp%t_expon_ab, pairs%pp%t_inverse_expon_ab, pairs%pp%xints)
#endif

  end subroutine compute_pp_pairs

  subroutine compute_sd_pairs(basis, cmax, pairs)
    !! Computes sd shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_sd
    integer(kind=int64) :: jmax, ij_loop, iii, jjj
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)
    logical :: iandj
    ! Rys-specific local
    real(dp), allocatable :: d02p_sd_rys(:, :)

    allocate (pairs%sd%xints(pairs%resources%n_s_shl*pairs%resources%n_d_shl))
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_d_shl
        iijj = jj + (ii - 1)*pairs%resources%n_d_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%sd%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%sd%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_d_shl + 1))
    pairs%sd%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_d_shl
        iijj = jj + (ii - 1)*pairs%resources%n_d_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        pairs%sd%pair_loc(iijj + 1) = pairs%sd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_sd = pairs%sd%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_d_shl + 1)
    pairs%sd%n_pairs = n_sd
    allocate (pairs%sd%sq(n_sd))
    allocate (ismlp(n_sd))
    allocate (pairs%sd%expon_a(n_sd))
    allocate (pairs%sd%expon_b(n_sd))
    allocate (pairs%sd%ismlp(n_sd))
    allocate (pairs%sd%d_coeff(n_sd))
    allocate (pairs%sd%d_coeff_alt(n_sd))
    allocate (d02p_sd_rys(2, n_sd))
    allocate (pairs%sd%t_expon_ab(n_sd))
    allocate (pairs%sd%t_inverse_expon_ab(n_sd))
    allocate (pairs%sd%t_alpha(n_sd))
    allocate (pairs%sd%t_beta(n_sd))

    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_d_shl

        iijj = jj + (ii - 1)*pairs%resources%n_d_shl

        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            pairs%sd%expon_a(pairs%sd%pair_loc(iijj) + ji) = buff(1)
            pairs%sd%expon_b(pairs%sd%pair_loc(iijj) + ji) = buff(2)
            pairs%sd%t_expon_ab(pairs%sd%pair_loc(iijj) + ji) = buff(3)
            pairs%sd%t_inverse_expon_ab(pairs%sd%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%sd%t_alpha(pairs%sd%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%sd%t_beta(pairs%sd%pair_loc(iijj) + ji) = pairs%sd%t_alpha(pairs%sd%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%sd%sq(pairs%sd%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)

            pairs%sd%d_coeff(pairs%sd%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            pairs%sd%d_coeff_alt(pairs%sd%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            iandj = ish == jsh

            jmax = 1
            ij_loop = 0
            do iii = 1, 6
              if (iandj) jmax = iii
              do jjj = 1, jmax

                ij_loop = ij_loop + 1
                d02p_sd_rys(1, pairs%sd%pair_loc(iijj) + ji) = buff(4)* &
                                                               basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                               basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

                d02p_sd_rys(2, pairs%sd%pair_loc(iijj) + ji) = sqrt3*buff(4)* &
                                                               basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                               basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

              end do
            end do

          end do
        end do
        pairs%sd%ismlp(pairs%sd%pair_loc(iijj) + 1:pairs%sd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do

    if (allocated(ismlp)) deallocate (ismlp)
    if (allocated(d02p_sd_rys)) deallocate (d02p_sd_rys)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%sd, pairs%sd%d_coeff, pairs%sd%d_coeff_alt, pairs%sd%expon_a, pairs%sd%expon_b, pairs%sd%ismlp, pairs%sd%pair_loc, pairs%sd%sq, pairs%sd%t_alpha, pairs%sd%t_beta, pairs%sd%t_expon_ab, pairs%sd%t_inverse_expon_ab, pairs%sd%xints)
#endif

  end subroutine compute_sd_pairs

  subroutine compute_pd_pairs(basis, cmax, pairs)
    !! Computes pd shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_pd
    integer(kind=int64) :: jmax, ij_loop, iii, jjj
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)
    logical :: iandj
    ! Rys-specific local
    real(dp), allocatable :: d12p_pd_rys(:, :)

    allocate (pairs%pd%xints(pairs%resources%n_p_shl*pairs%resources%n_d_shl))
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_d_shl
        iijj = jj + (ii - 1)*pairs%resources%n_d_shl
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%pd%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do
    allocate (pairs%pd%pair_loc(pairs%resources%n_p_shl*pairs%resources%n_d_shl + 1))
    pairs%pd%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_d_shl
        iijj = jj + (ii - 1)*pairs%resources%n_d_shl
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        pairs%pd%pair_loc(iijj + 1) = pairs%pd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do
    n_pd = pairs%pd%pair_loc(pairs%resources%n_p_shl*pairs%resources%n_d_shl + 1)
    pairs%pd%n_pairs = n_pd
    allocate (pairs%pd%sq(n_pd))
    allocate (ismlp(n_pd))
    allocate (pairs%pd%expon_a(n_pd))
    allocate (pairs%pd%expon_b(n_pd))
    allocate (pairs%pd%ismlp(n_pd))
    allocate (pairs%pd%d_coeff(n_pd))
    allocate (pairs%pd%d_coeff_alt(n_pd))
    allocate (d12p_pd_rys(2, n_pd))
    allocate (pairs%pd%t_expon_ab(n_pd))
    allocate (pairs%pd%t_inverse_expon_ab(n_pd))
    allocate (pairs%pd%t_alpha(n_pd))
    allocate (pairs%pd%t_beta(n_pd))

    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_d_shl

        iijj = jj + (ii - 1)*pairs%resources%n_d_shl

        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse exp_a_b

            pairs%pd%expon_a(pairs%pd%pair_loc(iijj) + ji) = buff(1)
            pairs%pd%expon_b(pairs%pd%pair_loc(iijj) + ji) = buff(2)
            pairs%pd%t_expon_ab(pairs%pd%pair_loc(iijj) + ji) = buff(3)
            pairs%pd%t_inverse_expon_ab(pairs%pd%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%pd%t_alpha(pairs%pd%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%pd%t_beta(pairs%pd%pair_loc(iijj) + ji) = pairs%pd%t_alpha(pairs%pd%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%pd%sq(pairs%pd%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)
            pairs%pd%d_coeff(pairs%pd%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            pairs%pd%d_coeff_alt(pairs%pd%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            iandj = ish == jsh

            jmax = 3
            ij_loop = 0
            do iii = 1, 6
              if (iandj) jmax = iii
              do jjj = 1, jmax

                ij_loop = ij_loop + 1
                d12p_pd_rys(1, pairs%pd%pair_loc(iijj) + ji) = buff(4)* &
                                                               basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                               basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

                d12p_pd_rys(2, pairs%pd%pair_loc(iijj) + ji) = sqrt3*buff(4)* &
                                                               basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                               basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

              end do
            end do

          end do
        end do
        pairs%pd%ismlp(pairs%pd%pair_loc(iijj) + 1:pairs%pd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do

    if (allocated(ismlp)) deallocate (ismlp)
    if (allocated(d12p_pd_rys)) deallocate (d12p_pd_rys)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%pd, pairs%pd%d_coeff, pairs%pd%d_coeff_alt, pairs%pd%expon_a, pairs%pd%expon_b, pairs%pd%ismlp, pairs%pd%pair_loc, pairs%pd%sq, pairs%pd%t_alpha, pairs%pd%t_beta, pairs%pd%t_expon_ab, pairs%pd%t_inverse_expon_ab, pairs%pd%xints)
#endif

  end subroutine compute_pd_pairs

  subroutine compute_dd_pairs(basis, cmax, pairs)
    !! Computes dd shell pair data.
    implicit none
    real(dp), intent(in) :: cmax(:)
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_dd
    integer(kind=int64) :: bra_loop
    real(dp) :: buff(20), r12, rab
    integer(kind=int64), allocatable :: ismlp(:)
    ! Rys-specific locals
    integer(kind=int64), allocatable :: locdd_rys(:)
    real(dp), allocatable :: d22p_dd_rys(:, :)

    allocate (pairs%dd%xints(pairs%resources%n_d_shl*(pairs%resources%n_d_shl + 1)/2))
    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%dd%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%dd%pair_loc(pairs%resources%n_d_shl*(pairs%resources%n_d_shl + 1)/2 + 1))
    allocate (locdd_rys(pairs%resources%n_d_shl*(pairs%resources%n_d_shl + 1)/2 + 1))
    pairs%dd%pair_loc(1) = 0
    locdd_rys(1) = 0
    bra_loop = 0
    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, ii
        bra_loop = bra_loop + 1
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)
        pairs%dd%pair_loc(iijj + 1) = pairs%dd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
        locdd_rys(iijj + 1) = locdd_rys(iijj) + basis%contr_num(ish)*basis%contr_num(jsh) + bra_loop
      end do
    end do

    n_dd = pairs%dd%pair_loc(pairs%resources%n_d_shl*(pairs%resources%n_d_shl + 1)/2 + 1)
    pairs%dd%n_pairs = n_dd
    allocate (pairs%dd%sq(n_dd))
    allocate (pairs%dd%expon_a(n_dd))
    allocate (pairs%dd%expon_b(n_dd))
    allocate (ismlp(n_dd))
    allocate (pairs%dd%t_expon_ab(n_dd))
    allocate (pairs%dd%ismlp(n_dd))
    allocate (pairs%dd%d_coeff(n_dd))
    allocate (pairs%dd%d_coeff_alt(n_dd))
    allocate (d22p_dd_rys(3, n_dd))
    allocate (pairs%dd%t_inverse_expon_ab(n_dd))
    allocate (pairs%dd%t_alpha(n_dd))
    allocate (pairs%dd%t_beta(n_dd))

    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, ii

        iijj = jj + ii*(ii - 1)/2

        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_d_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                         ! inverse of exponent_ab

            pairs%dd%expon_a(pairs%dd%pair_loc(iijj) + ji) = buff(1)
            pairs%dd%expon_b(pairs%dd%pair_loc(iijj) + ji) = buff(2)
            pairs%dd%t_expon_ab(pairs%dd%pair_loc(iijj) + ji) = buff(3)
            pairs%dd%t_inverse_expon_ab(pairs%dd%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(6) = 1.0_dp - buff(5)                                         ! 1 - exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2

            pairs%dd%t_alpha(pairs%dd%pair_loc(iijj) + ji) = buff(6)*rab
            pairs%dd%t_beta(pairs%dd%pair_loc(iijj) + ji) = pairs%dd%t_alpha(pairs%dd%pair_loc(iijj) + ji) - rab

            if (buff(8) > cux) then
              ismlp(ji) = 2
              cycle
            end if

            buff(9) = buff(4)*exp(-buff(8))
            buff(10) = buff(9)*cmax(basis%sh_loc(ish) + i - 1)*cmax(basis%sh_loc(jsh) + j - 1)

            pairs%dd%sq(pairs%dd%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            ismlp(ji) = 0

            buff(9) = twopi_5_2*buff(9)
            pairs%dd%d_coeff(pairs%dd%pair_loc(iijj) + ji) = buff(9)* &
                                                             basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                             basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            pairs%dd%d_coeff_alt(pairs%dd%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

            d22p_dd_rys(1, pairs%dd%pair_loc(iijj) + ji) = buff(4)* &
                                                           basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                           basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)
            d22p_dd_rys(2, pairs%dd%pair_loc(iijj) + ji) = sqrt3*buff(4)* &
                                                           basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                           basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)
            d22p_dd_rys(3, pairs%dd%pair_loc(iijj) + ji) = sqrt3*sqrt3*buff(4)* &
                                                           basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                           basis%contr_coef_d(basis%sh_loc(jsh) + j - 1)

          end do
        end do
        pairs%dd%ismlp(pairs%dd%pair_loc(iijj) + 1:pairs%dd%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)) = ismlp(1:basis%contr_num(ish)*basis%contr_num(jsh))
      end do
    end do
    if (allocated(ismlp)) deallocate (ismlp)
    if (allocated(locdd_rys)) deallocate (locdd_rys)
    if (allocated(d22p_dd_rys)) deallocate (d22p_dd_rys)

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%dd, pairs%dd%d_coeff, pairs%dd%d_coeff_alt, pairs%dd%expon_a, pairs%dd%expon_b, pairs%dd%ismlp, pairs%dd%pair_loc, pairs%dd%sq, pairs%dd%t_alpha, pairs%dd%t_beta, pairs%dd%t_expon_ab, pairs%dd%t_inverse_expon_ab, pairs%dd%xints)
#endif

  end subroutine compute_dd_pairs

  subroutine compute_sf_pairs(basis, pairs)
    !! Computes sf shell pair data.
    implicit none
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_sf
    real(dp) :: buff(20), r12, rab

    allocate (pairs%sf%xints(pairs%resources%n_s_shl*pairs%resources%n_f_shl))
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%sf%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%sf%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_f_shl + 1))
    pairs%sf%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        pairs%sf%pair_loc(iijj + 1) = pairs%sf%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_sf = pairs%sf%pair_loc(pairs%resources%n_s_shl*pairs%resources%n_f_shl + 1)
    pairs%sf%n_pairs = n_sf
    allocate (pairs%sf%sq(n_sf))
    allocate (pairs%sf%t_expon_ab(n_sf))
    allocate (pairs%sf%expon_a(n_sf))
    allocate (pairs%sf%expon_b(n_sf))
    allocate (pairs%sf%d_coeff_alt(n_sf))
    allocate (pairs%sf%t_inverse_expon_ab(n_sf))

    do ii = 1, pairs%resources%n_s_shl
      do jj = 1, pairs%resources%n_f_shl

        iijj = jj + (ii - 1)*pairs%resources%n_f_shl

        ish = pairs%resources%i_s_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                        ! inverse exp_a_b

            pairs%sf%expon_a(pairs%sf%pair_loc(iijj) + ji) = buff(1)
            pairs%sf%expon_b(pairs%sf%pair_loc(iijj) + ji) = buff(2)
            pairs%sf%t_expon_ab(pairs%sf%pair_loc(iijj) + ji) = buff(3)
            pairs%sf%t_inverse_expon_ab(pairs%sf%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                  ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2
            buff(9) = buff(4)*exp(-buff(8))                         ! Intermediate for ERIC

            pairs%sf%sq(pairs%sf%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

            pairs%sf%d_coeff_alt(pairs%sf%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_s(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

          end do
        end do
      end do
    end do

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%sf, pairs%sf%d_coeff, pairs%sf%d_coeff_alt, pairs%sf%expon_a, pairs%sf%expon_b, pairs%sf%ismlp, pairs%sf%pair_loc, pairs%sf%sq, pairs%sf%t_alpha, pairs%sf%t_beta, pairs%sf%t_expon_ab, pairs%sf%t_inverse_expon_ab, pairs%sf%xints)
#endif

  end subroutine compute_sf_pairs

  subroutine compute_pf_pairs(basis, pairs)
    !! Computes pf shell pair data.
    implicit none
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_pf
    real(dp) :: buff(20), r12, rab

    allocate (pairs%pf%xints(pairs%resources%n_p_shl*pairs%resources%n_f_shl))
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%pf%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%pf%pair_loc(pairs%resources%n_p_shl*pairs%resources%n_f_shl + 1))
    pairs%pf%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        pairs%pf%pair_loc(iijj + 1) = pairs%pf%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_pf = pairs%pf%pair_loc(pairs%resources%n_p_shl*pairs%resources%n_f_shl + 1)
    pairs%pf%n_pairs = n_pf
    allocate (pairs%pf%sq(n_pf))
    allocate (pairs%pf%t_expon_ab(n_pf))
    allocate (pairs%pf%expon_a(n_pf))
    allocate (pairs%pf%expon_b(n_pf))
    allocate (pairs%pf%d_coeff_alt(n_pf))
    allocate (pairs%pf%t_inverse_expon_ab(n_pf))

    do ii = 1, pairs%resources%n_p_shl
      do jj = 1, pairs%resources%n_f_shl

        iijj = jj + (ii - 1)*pairs%resources%n_f_shl

        ish = pairs%resources%i_p_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                        ! inverse exp_a_b

            pairs%pf%expon_a(pairs%pf%pair_loc(iijj) + ji) = buff(1)
            pairs%pf%expon_b(pairs%pf%pair_loc(iijj) + ji) = buff(2)
            pairs%pf%t_expon_ab(pairs%pf%pair_loc(iijj) + ji) = buff(3)
            pairs%pf%t_inverse_expon_ab(pairs%pf%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                 ! exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                  ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                         ! reduced exponent·|A-B|**2
            buff(9) = buff(4)*exp(-buff(8))                         ! Intermediate for ERIC

            pairs%pf%sq(pairs%pf%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

            pairs%pf%d_coeff_alt(pairs%pf%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_p(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

          end do
        end do
      end do
    end do

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%pf, pairs%pf%d_coeff, pairs%pf%d_coeff_alt, pairs%pf%expon_a, pairs%pf%expon_b, pairs%pf%ismlp, pairs%pf%pair_loc, pairs%pf%sq, pairs%pf%t_alpha, pairs%pf%t_beta, pairs%pf%t_expon_ab, pairs%pf%t_inverse_expon_ab, pairs%pf%xints)
#endif

  end subroutine compute_pf_pairs

  subroutine compute_df_pairs(basis, pairs)
    !! Computes df shell pair data.
    implicit none
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, n_df
    real(dp) :: buff(20), r12, rab

    allocate (pairs%df%xints(pairs%resources%n_d_shl*pairs%resources%n_f_shl))
    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%df%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%df%pair_loc(pairs%resources%n_d_shl*pairs%resources%n_f_shl + 1))
    pairs%df%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, pairs%resources%n_f_shl
        iijj = jj + (ii - 1)*pairs%resources%n_f_shl
        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        pairs%df%pair_loc(iijj + 1) = pairs%df%pair_loc(iijj) + basis%contr_num(ish)*basis%contr_num(jsh)
      end do
    end do

    n_df = pairs%df%pair_loc(pairs%resources%n_d_shl*pairs%resources%n_f_shl + 1)
    pairs%df%n_pairs = n_df
    allocate (pairs%df%sq(n_df))
    allocate (pairs%df%t_expon_ab(n_df))
    allocate (pairs%df%expon_a(n_df))
    allocate (pairs%df%expon_b(n_df))
    allocate (pairs%df%d_coeff_alt(n_df))
    allocate (pairs%df%t_inverse_expon_ab(n_df))

    do ii = 1, pairs%resources%n_d_shl
      do jj = 1, pairs%resources%n_f_shl

        iijj = jj + (ii - 1)*pairs%resources%n_f_shl

        ish = pairs%resources%i_d_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, basis%contr_num(jsh)

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                        ! inverse exp_a_b

            pairs%df%expon_a(pairs%df%pair_loc(iijj) + ji) = buff(1)
            pairs%df%expon_b(pairs%df%pair_loc(iijj) + ji) = buff(2)
            pairs%df%t_expon_ab(pairs%df%pair_loc(iijj) + ji) = buff(3)
            pairs%df%t_inverse_expon_ab(pairs%df%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                ! exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                        ! reduced exponent·|A-B|**2
            buff(9) = buff(4)*exp(-buff(8))                        ! Intermediate for ERIC

            pairs%df%sq(pairs%df%pair_loc(iijj) + ji) = buff(9)* &
                                                        basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                        basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

            pairs%df%d_coeff_alt(pairs%df%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_d(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

          end do
        end do
      end do
    end do

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%df, pairs%df%d_coeff, pairs%df%d_coeff_alt, pairs%df%expon_a, pairs%df%expon_b, pairs%df%ismlp, pairs%df%pair_loc, pairs%df%sq, pairs%df%t_alpha, pairs%df%t_beta, pairs%df%t_expon_ab, pairs%df%t_inverse_expon_ab, pairs%df%xints)
#endif

  end subroutine compute_df_pairs

  subroutine compute_ff_pairs(basis, pairs)
    !! Computes ff shell pair data.
    implicit none
    type(shell_pair_container_t), intent(inout) :: pairs
    type(basis_t), intent(in) :: basis

    integer(kind=int64) :: ii, jj, iijj, ish, jsh, ijsh, ji, i, j, jmax, n_ff
    real(dp) :: buff(20), r12, rab
    logical :: iandj

    allocate (pairs%ff%xints(pairs%resources%n_f_shl*(pairs%resources%n_f_shl + 1)/2))
    do ii = 1, pairs%resources%n_f_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_f_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        ijsh = pairs%resources%ia(ish) + jsh
        if (jsh > ish) ijsh = pairs%resources%ia(jsh) + ish
        pairs%ff%xints(iijj) = basis%schwrz_int(ijsh)
      end do
    end do

    allocate (pairs%ff%pair_loc(pairs%resources%n_f_shl*(pairs%resources%n_f_shl + 1)/2 + 1))
    pairs%ff%pair_loc(1) = 0
    do ii = 1, pairs%resources%n_f_shl
      do jj = 1, ii
        iijj = jj + ii*(ii - 1)/2
        ish = pairs%resources%i_f_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)
        ji = basis%contr_num(ish)*basis%contr_num(jsh)
        pairs%ff%pair_loc(iijj + 1) = pairs%ff%pair_loc(iijj) + ji
      end do
    end do

    n_ff = pairs%ff%pair_loc(pairs%resources%n_f_shl*(pairs%resources%n_f_shl + 1)/2 + 1)
    pairs%ff%n_pairs = n_ff
    allocate (pairs%ff%t_expon_ab(n_ff))
    allocate (pairs%ff%expon_a(n_ff))
    allocate (pairs%ff%expon_b(n_ff))
    allocate (pairs%ff%d_coeff_alt(n_ff))
    allocate (pairs%ff%t_inverse_expon_ab(n_ff))

    do ii = 1, pairs%resources%n_f_shl
      do jj = 1, ii

        iijj = jj + ii*(ii - 1)/2

        ish = pairs%resources%i_f_shl(ii)
        jsh = pairs%resources%i_f_shl(jj)

        r12 = 0.0_dp
        r12 = r12 + (pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))*(pairs%resources%coord_sh(jsh, 1) - pairs%resources%coord_sh(ish, 1))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))*(pairs%resources%coord_sh(jsh, 2) - pairs%resources%coord_sh(ish, 2))
        r12 = r12 + (pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))*(pairs%resources%coord_sh(jsh, 3) - pairs%resources%coord_sh(ish, 3))

        rab = 0.0_dp
        if (r12 /= 0.0_dp) then
          rab = sqrt(r12)
        end if

        ji = 0
        iandj = ish == jsh
        jmax = basis%contr_num(jsh)
        do i = 1, basis%contr_num(ish)
          buff(1) = basis%exponents(basis%sh_loc(ish) + i - 1)                 ! exponent_a
          do j = 1, jmax

            ji = ji + 1

            buff(2) = basis%exponents(basis%sh_loc(jsh) + j - 1)         ! exponent_b
            buff(3) = buff(1) + buff(2)                                 ! exponent_a_b
            buff(4) = 1.0_dp/buff(3)                                        ! inverse exp_a_b

            pairs%ff%expon_a(pairs%ff%pair_loc(iijj) + ji) = buff(1)
            pairs%ff%expon_b(pairs%ff%pair_loc(iijj) + ji) = buff(2)
            pairs%ff%t_expon_ab(pairs%ff%pair_loc(iijj) + ji) = buff(3)
            pairs%ff%t_inverse_expon_ab(pairs%ff%pair_loc(iijj) + ji) = buff(4)*0.5_dp

            buff(5) = buff(1)*buff(4)                                ! exp_a/exp_a+b
            buff(7) = buff(5)*buff(2)                                 ! exp_a*exp_b/exp_a+b
            buff(8) = r12*buff(7)                                        ! reduced exponent·|A-B|**2
            buff(9) = buff(4)*exp(-buff(8))                        ! Intermediate for K in ERIC

            pairs%ff%d_coeff_alt(pairs%ff%pair_loc(iijj) + ji) = buff(4)* &
                                                                 basis%contr_coef_f(basis%sh_loc(ish) + i - 1)* &
                                                                 basis%contr_coef_f(basis%sh_loc(jsh) + j - 1)

          end do
        end do
      end do
    end do

#ifdef ENABLE_GPU
    !$omp target enter data &
    !$omp map(to: pairs%ff, pairs%ff%d_coeff, pairs%ff%d_coeff_alt, pairs%ff%expon_a, pairs%ff%expon_b, pairs%ff%ismlp, pairs%ff%pair_loc, pairs%ff%sq, pairs%ff%t_alpha, pairs%ff%t_beta, pairs%ff%t_expon_ab, pairs%ff%t_inverse_expon_ab, pairs%ff%xints)
#endif

  end subroutine compute_ff_pairs

  subroutine shell_pair_deallocate(pairs)
    !! Deallocates shell pair data and unmaps from device.
    !! Call this when shell pair data is no longer needed.
    implicit none
    type(shell_pair_container_t), intent(inout) :: pairs

#ifdef ENABLE_GPU
    if (pairs%resources%n_s_shl > 0) then
      !$omp target exit data &
      !$omp map(delete: pairs%ss%d_coeff, pairs%ss%d_coeff_alt, pairs%ss%expon_a, pairs%ss%expon_b, pairs%ss%ismlp, pairs%ss%pair_loc, pairs%ss%sq, pairs%ss%t_alpha, pairs%ss%t_beta, pairs%ss%t_expon_ab, pairs%ss%t_inverse_expon_ab, pairs%ss%xints)
      !$omp target exit data map(delete: pairs%resources%i_s_shl)
    end if

    if (pairs%resources%n_p_shl > 0) then
      !$omp target exit data &
      !$omp map(delete: pairs%sp%d_coeff, pairs%sp%d_coeff_alt, pairs%sp%expon_a, pairs%sp%expon_b, pairs%sp%ismlp, pairs%sp%pair_loc, pairs%sp%sq, pairs%sp%t_alpha, pairs%sp%t_beta, pairs%sp%t_expon_ab, pairs%sp%t_inverse_expon_ab, pairs%sp%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%pp%d_coeff, pairs%pp%d_coeff_alt, pairs%pp%expon_a, pairs%pp%expon_b, pairs%pp%ismlp, pairs%pp%pair_loc, pairs%pp%sq, pairs%pp%t_alpha, pairs%pp%t_beta, pairs%pp%t_expon_ab, pairs%pp%t_inverse_expon_ab, pairs%pp%xints)
      !$omp target exit data map(delete: pairs%resources%i_p_shl)
    end if

    if (pairs%resources%n_d_shl > 0) then
      !$omp target exit data &
      !$omp map(delete: pairs%sd%d_coeff, pairs%sd%d_coeff_alt, pairs%sd%expon_a, pairs%sd%expon_b, pairs%sd%ismlp, pairs%sd%pair_loc, pairs%sd%sq, pairs%sd%t_alpha, pairs%sd%t_beta, pairs%sd%t_expon_ab, pairs%sd%t_inverse_expon_ab, pairs%sd%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%pd%d_coeff, pairs%pd%d_coeff_alt, pairs%pd%expon_a, pairs%pd%expon_b, pairs%pd%ismlp, pairs%pd%pair_loc, pairs%pd%sq, pairs%pd%t_alpha, pairs%pd%t_beta, pairs%pd%t_expon_ab, pairs%pd%t_inverse_expon_ab, pairs%pd%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%dd%d_coeff, pairs%dd%d_coeff_alt, pairs%dd%expon_a, pairs%dd%expon_b, pairs%dd%ismlp, pairs%dd%pair_loc, pairs%dd%sq, pairs%dd%t_alpha, pairs%dd%t_beta, pairs%dd%t_expon_ab, pairs%dd%t_inverse_expon_ab, pairs%dd%xints)
      !$omp target exit data map(delete: pairs%resources%i_d_shl)
    end if

#ifdef ENABLE_F
    if (pairs%resources%n_f_shl > 0) then
      !$omp target exit data &
      !$omp map(delete: pairs%sf%d_coeff, pairs%sf%d_coeff_alt, pairs%sf%expon_a, pairs%sf%expon_b, pairs%sf%ismlp, pairs%sf%pair_loc, pairs%sf%sq, pairs%sf%t_alpha, pairs%sf%t_beta, pairs%sf%t_expon_ab, pairs%sf%t_inverse_expon_ab, pairs%sf%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%pf%d_coeff, pairs%pf%d_coeff_alt, pairs%pf%expon_a, pairs%pf%expon_b, pairs%pf%ismlp, pairs%pf%pair_loc, pairs%pf%sq, pairs%pf%t_alpha, pairs%pf%t_beta, pairs%pf%t_expon_ab, pairs%pf%t_inverse_expon_ab, pairs%pf%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%df%d_coeff, pairs%df%d_coeff_alt, pairs%df%expon_a, pairs%df%expon_b, pairs%df%ismlp, pairs%df%pair_loc, pairs%df%sq, pairs%df%t_alpha, pairs%df%t_beta, pairs%df%t_expon_ab, pairs%df%t_inverse_expon_ab, pairs%df%xints)
      !$omp target exit data &
      !$omp map(delete: pairs%ff%d_coeff, pairs%ff%d_coeff_alt, pairs%ff%expon_a, pairs%ff%expon_b, pairs%ff%ismlp, pairs%ff%pair_loc, pairs%ff%sq, pairs%ff%t_alpha, pairs%ff%t_beta, pairs%ff%t_expon_ab, pairs%ff%t_inverse_expon_ab, pairs%ff%xints)
      !$omp target exit data map(delete: pairs%resources%i_f_shl)
    end if
#endif

    ! Unmap Boys function grids from device
    !$omp target exit data map(delete: boys_grid_zero, exponent_grid)

    ! Unmap resources from device
    !$omp target exit data &
    !$omp map(delete: pairs%resources%coord_sh, pairs%resources%ia, &
    !$omp   pairs%resources%contr_num, pairs%resources%atom_loc, &
    !$omp   pairs%resources)
#endif

  end subroutine shell_pair_deallocate

end module liberi_shell_pair
