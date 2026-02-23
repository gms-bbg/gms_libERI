module rys_kernels
  use liberi_types, only: dp, shell_pair_t, eri_resources_t
  implicit none
  interface
        module subroutine int2130(pd_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pd_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2230(dd_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int2231(dd_pair, pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: dd_pair, pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3000(sf_pair, ss_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sf_pair, ss_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3010(sf_pair, sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sf_pair, sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3011(sf_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sf_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3020(sf_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sf_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3030(sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3100(pf_pair, ss_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, ss_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3110(pf_pair, sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3111(pf_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3120(pf_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3121(pf_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3130(pf_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3131(pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3200(df_pair, ss_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, ss_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3210(df_pair, sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3211(df_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3220(df_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3221(df_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3222(df_pair, dd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, dd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3230(df_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3231(df_pair, pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair, pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3232(df_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: df_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3300(ff_pair, ss_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, ss_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3310(ff_pair, sp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, sp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3311(ff_pair, pp_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, pp_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3320(ff_pair, sd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, sd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3321(ff_pair, pd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, pd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3322(ff_pair, dd_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, dd_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3330(ff_pair, sf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, sf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3331(ff_pair, pf_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, pf_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3332(ff_pair, df_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair, df_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
        module subroutine int3333(ff_pair, density, fock, res)
      type(shell_pair_t), intent(in) :: ff_pair
      real(dp), intent(in) :: density(:)
      real(dp), intent(inout) :: fock(:)
      type(eri_resources_t), intent(in) :: res
    end subroutine
  end interface
end module rys_kernels
