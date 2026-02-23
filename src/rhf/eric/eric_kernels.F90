module eric_kernels
  use liberi_types, only: dp, shell_pair_t, eri_resources_t
  implicit none
  interface
        module subroutine int0030(ss_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0031(ss_pair, pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0032(ss_pair, df_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, df_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int1030(sp_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int1031(sp_pair, pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int1130(pp_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pp_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2030(sd_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sd_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
  end interface
end module eric_kernels
