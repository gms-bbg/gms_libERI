module rot_axis_kernels
  use liberi_types, only: dp, shell_pair_t, eri_resources_t
  implicit none
  interface
    module subroutine int0000(ss_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0001(ss_pair, sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0002(ss_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0011(ss_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0012(ss_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0022(ss_pair, dd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ss_pair, dd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0101(sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0102(sp_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0111(sp_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0112(sp_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0122(sp_pair, dd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sp_pair, dd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int0202(sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int1102(pp_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pp_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int1111(pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2111(pd_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pd_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2120(pd_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pd_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2121(pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2211(dd_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2220(dd_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2221(dd_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2222gen(dd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
  end interface
end module rot_axis_kernels
