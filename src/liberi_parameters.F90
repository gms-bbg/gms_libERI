!! Numerical parameters and constants for integral evaluation
module parameters
  !! Contains derived mathematical constants and screening thresholds
  !! used throughout the integral kernels.
  !!
  !! This module provides:
  !!
  !! - Pre-computed pi-based constants for integral normalization
  !! - Square roots used in angular momentum normalization factors
  !! - Schwarz screening cutoff for integral prescreening
  use liberi_types, only: dp, int64
  use liberi_constants, only: M_PI

  real(dp), parameter :: tenm12 = 1.0e-12_dp
    !! Small value threshold (10^-12)

  real(dp), parameter :: pi4 = M_PI/4.0_dp
    !! Pi divided by 4

  real(dp), parameter :: twopi_5_2 = 2.0_dp*(M_PI**2)*sqrt(M_PI)
    !! 2 * pi^(5/2), used in integral prefactors
  real(dp), parameter :: sqrt2_pi_5_4 = sqrt(2.0_dp)*(M_PI**1.25_dp)
    !! sqrt(2) * pi^(5/4)
  real(dp), parameter :: pi_1_4_div_sqrt2 = (M_PI**0.25_dp)/sqrt(2.0_dp)
    !! pi^(1/4) / sqrt(2)

  real(dp), parameter :: sqrt3 = sqrt(3.0_dp)
    !! Square root of 3, used in d-function normalization
  real(dp), parameter :: sqrt5 = sqrt(5.0_dp)
    !! Square root of 5, used in f-function normalization
  real(dp), parameter :: sqrt7 = sqrt(7.0_dp)
    !! Square root of 7, used in f-function normalization
  real(dp), parameter :: sqrt15 = sqrt(15.0_dp)
    !! Square root of 15, used in angular momentum normalization
  real(dp), parameter :: sqrt_pi_div_2 = sqrt(M_PI)/2.0_dp
    !! sqrt(pi) / 2

  real(dp), parameter :: cutoff_schwarz = 5.00000000000000018e-11_dp
    !! Schwarz screening cutoff for integral prescreening.
    !! Integrals with Schwarz bound below this value are skipped.

  real(dp), parameter :: acycut = 1.0e-10_dp
    !! Asymptotic cutoff threshold
  real(dp), parameter :: cux = 50.0_dp
    !! Upper bound for Boys function argument before asymptotic expansion

end module parameters
