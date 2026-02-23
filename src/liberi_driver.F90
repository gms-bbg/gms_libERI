!! Top-level integral dispatch routine
module liberi_integral_driver
  !! Orchestrates the evaluation of all two-electron integrals.
  !!
  !! This module provides the main entry point for computing two-electron
  !! repulsion integrals (ERIs) and building the Fock matrix. It dispatches
  !! to the appropriate integral kernels based on the angular momentum
  !! combinations present in the basis set.
  !!
  !! ### Kernel Types
  !!
  !! The driver calls three types of integral kernels:
  !!
  !! 1. **Rotated-axis kernels** (`rot_axis_kernels`): For s, p, d shell
  !!    combinations using the Head-Gordon & Pople method
  !! 2. **ERIC kernels** (`eric_kernels`): McMurchie-Davidson style for
  !!    f-function quartets (requires `ENABLE_F`)
  !! 3. **Rys kernels** (`rys_kernels`): Rys quadrature for higher angular
  !!    momentum combinations (requires `ENABLE_F`)
  !!
  !! ### Kernel Naming Convention
  !!
  !! Kernels are named `intABCD` where A, B, C, D are angular momentum indices:
  !!
  !! - 0 = s-type
  !! - 1 = p-type
  !! - 2 = d-type
  !! - 3 = f-type
  !!
  !! For example, `int0111` computes (sp|pp) integrals.
  use liberi_types, only: dp, shell_pair_container_t, eri_resources_t
  use rot_axis_kernels, only: int0000, int0001, int0011, int0101, int0111, int1111, &
                              int0002, int0012, int0102, int0112, int0022, &
                              int0122, int1102, int0202, &
                              int2211, int2121, int2120, int2111, int2220, int2221, int2222gen
#ifdef ENABLE_F
  use eric_kernels, only: int0030, int0031, int0032, int1030
  use rys_kernels, only: int3333, int3330, int3011, int3020, int2130, int2230, int3030, &
                         int3110, int3111, int3120, int3121, int2231, &
                         int3130, int3131, int3210, int3331, int3211, &
                         int3220, int3221, int3222, int3230, int3231, &
                         int3232, int3300, int3310, int3311, int3320, &
                         int3321, int3322, int3332
#endif
  implicit none
  private

  public :: compute_integrals

contains

  subroutine compute_integrals(pairs, density, fock)
    !! Compute all two-electron integrals and accumulate into Fock matrix.
    !!
    !! This is the clean API with no global dependencies.
    !! density and fock are passed directly; resources come from pairs%resources.
    type(shell_pair_container_t), intent(in) :: pairs
    real(dp), intent(in) :: density(:)
    real(dp), intent(inout) :: fock(:)

    associate (res => pairs%resources)

!-----------------RHF kernels---------------
      fock = 0.0_dp
#ifdef ENABLE_GPU
      !$omp target enter data map(to: density, fock)
#endif

      if (res%n_s_shl .gt. 0) then
        call int0000(pairs%ss, density, fock, res)
      end if
      if (res%n_p_shl .gt. 0) then

        call int0001(pairs%ss, pairs%sp, density, fock, res)
        call int0011(pairs%ss, pairs%pp, density, fock, res)
        call int0111(pairs%sp, pairs%pp, density, fock, res)
        call int0101(pairs%sp, density, fock, res)
        call int1111(pairs%pp, density, fock, res)
      end if
      if (res%n_d_shl .gt. 0) then
        call int0002(pairs%ss, pairs%sd, density, fock, res)
        call int0012(pairs%ss, pairs%pd, density, fock, res)
        call int0102(pairs%sp, pairs%sd, density, fock, res)
        call int0112(pairs%sp, pairs%pd, density, fock, res)
        call int0022(pairs%ss, pairs%dd, density, fock, res)
        call int0122(pairs%sp, pairs%dd, density, fock, res)
        call int1102(pairs%pp, pairs%sd, density, fock, res)
        call int0202(pairs%sd, density, fock, res)
        !*********************************
        !*                               *
        !* GENERATED RYS KERNELS FOLLOW  *
        !*                               *
        !*********************************

        call int2211(pairs%dd, pairs%pp, density, fock, res)
        call int2121(pairs%pd, density, fock, res)
        call int2120(pairs%pd, pairs%sd, density, fock, res)
        call int2111(pairs%pd, pairs%pp, density, fock, res)
        call int2220(pairs%dd, pairs%sd, density, fock, res)
        call int2221(pairs%dd, pairs%pd, density, fock, res)
        call int2222gen(pairs%dd, density, fock, res)
      end if
#ifdef ENABLE_F

      if (res%n_f_shl .gt. 0) then

        !*********************************
        !*                               *
        !* GENERATED RYS KERNELS PART 1  *
        !*                               *
        !*********************************

        call int3333(pairs%ff, density, fock, res)
        !*********************************
        !*                               *
        !*      ERIC KERNELS FOLLOW      *
        !*                               *
        !*********************************

        call int0030(pairs%ss, pairs%sf, density, fock, res)
        call int0031(pairs%ss, pairs%pf, density, fock, res)
        call int0032(pairs%ss, pairs%df, density, fock, res)
        call int1030(pairs%sp, pairs%sf, density, fock, res)
        !*********************************
        !*                               *
        !* GENERATED RYS KERNELS FOLLOW  *
        !*                               *
        !*********************************

        call int3330(pairs%ff, pairs%sf, density, fock, res)
        call int3011(pairs%sf, pairs%pp, density, fock, res)
        call int3020(pairs%sf, pairs%sd, density, fock, res)
        call int2130(pairs%pd, pairs%sf, density, fock, res)
        call int2230(pairs%dd, pairs%sf, density, fock, res)
        call int3030(pairs%sf, density, fock, res)
        call int3110(pairs%pf, pairs%sp, density, fock, res)
        call int3111(pairs%pf, pairs%pp, density, fock, res)
        call int3120(pairs%pf, pairs%sd, density, fock, res)
        call int3121(pairs%pf, pairs%pd, density, fock, res)
        call int2231(pairs%dd, pairs%pf, density, fock, res)
        call int3130(pairs%pf, pairs%sf, density, fock, res)
        call int3131(pairs%pf, density, fock, res)
        call int3210(pairs%df, pairs%sp, density, fock, res)
        call int3331(pairs%ff, pairs%pf, density, fock, res)
        call int3211(pairs%df, pairs%pp, density, fock, res)
        call int3220(pairs%df, pairs%sd, density, fock, res)
        call int3221(pairs%df, pairs%pd, density, fock, res)
        call int3222(pairs%df, pairs%dd, density, fock, res)
        call int3230(pairs%df, pairs%sf, density, fock, res)
        call int3231(pairs%df, pairs%pf, density, fock, res)
        call int3232(pairs%df, density, fock, res)
        call int3300(pairs%ff, pairs%ss, density, fock, res)
        call int3310(pairs%ff, pairs%sp, density, fock, res)
        call int3311(pairs%ff, pairs%pp, density, fock, res)
        call int3320(pairs%ff, pairs%sd, density, fock, res)
        call int3321(pairs%ff, pairs%pd, density, fock, res)
        call int3322(pairs%ff, pairs%dd, density, fock, res)
        call int3332(pairs%ff, pairs%df, density, fock, res)
      end if
#endif

#ifdef ENABLE_GPU
      !$omp target exit data map(from: fock) map(delete: density)
#endif

    end associate
  end subroutine compute_integrals

end module liberi_integral_driver
