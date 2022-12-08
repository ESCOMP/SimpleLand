module lnd2atmMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle lnd2atm mapping
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use clm_varpar           , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use clm_varcon           , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl           , only : iulog
  use decompMod            , only : bounds_type
  use subgridAveMod        , only : p2g, c2g 
  use lnd2atmType          , only : lnd2atm_type
  use atm2lndType          , only : atm2lnd_type
  use EnergyFluxType       , only : energyflux_type
  use SurfaceAlbedoType    , only : surfalb_type
  use TemperatureType      , only : temperature_type
  use WaterstateType       , only : waterstate_type
  use glcBehaviorMod       , only : glc_behavior_type
  use ColumnType           , only : col
  use LandunitType         , only : lun
  use GridcellType         , only : grc                
  use landunit_varcon      , only : istice_mec
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm_minimal
  !
  ! !PRIVATE MEMBER FUNCTIONS:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine lnd2atm_minimal(bounds, &
      waterstate_inst, surfalb_inst, energyflux_inst, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_inst component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !USES:
    use clm_varcon, only : sb
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstate_type) , intent(in)    :: waterstate_inst
    type(surfalb_type)    , intent(in)    :: surfalb_inst
    type(energyflux_type) , intent(in)    :: energyflux_inst
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    real(r8), parameter :: amC   = 12.0_r8          ! Atomic mass number for Carbon
    real(r8), parameter :: amO   = 16.0_r8          ! Atomic mass number for Oxygen
    real(r8), parameter :: amCO2 = amC + 2.0_r8*amO ! Atomic mass number for CO2
    ! The following converts g of C to kg of CO2
    real(r8), parameter :: convertgC2kgCO2 = 1.0e-3_r8 * (amCO2/amC)
    !------------------------------------------------------------------------

    call c2g(bounds, &
         waterstate_inst%h2osno_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%h2osno_grc(g) = lnd2atm_inst%h2osno_grc(g)/1000._r8
    end do

    call p2g(bounds, numrad, &
         surfalb_inst%albd_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albd_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, numrad, &
         surfalb_inst%albi_patch (bounds%begp:bounds%endp, :), &
         lnd2atm_inst%albi_grc   (bounds%begg:bounds%endg, :), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    call p2g(bounds, &
         energyflux_inst%eflx_lwrad_out_patch (bounds%begp:bounds%endp), &
         lnd2atm_inst%eflx_lwrad_out_grc      (bounds%begg:bounds%endg), &
         p2c_scale_type='unity', c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    ! TODO slevis SLIM: Eliminating next assignment returns error:
    ! Longwave down is <= 0, while eliminating and then changing the
    ! initialization of t_rad_grc to tfrz in lnd2atmType works but changes
    ! answers throughout, including in MML variables. See
    ! /glade/scratch/slevis/ERS_D_Ld60.f19_g16.H_MML_2000_CAM5.cheyenne_gnu.clm-global_uniform_g16_SOM.C.20221207_111523_9gqwvp/ERS_D_Ld60.f19_g16.H_MML_2000_CAM5.cheyenne_gnu.clm-global_uniform_g16_SOM.C.20221207_111523_9gqwvp.clm2.h0.0001-03-02-00000.nc.cprnc.out
    do g = bounds%begg,bounds%endg
       lnd2atm_inst%t_rad_grc(g) = sqrt(sqrt(lnd2atm_inst%eflx_lwrad_out_grc(g)/sb))
    end do

  end subroutine lnd2atm_minimal

end module lnd2atmMod
