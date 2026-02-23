!! Mathematical constants for libERI
module liberi_constants
  !! Provides fundamental mathematical constants used in integral evaluation.
  use liberi_types, only: dp
  implicit none
  private

  real(dp), parameter, public :: M_PI = 3.14159265358979323846_dp
    !! Pi to double precision (3.14159265358979...)

end module liberi_constants
