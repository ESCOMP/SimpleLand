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
      waterstate_inst, lnd2atm_inst)
    !
    ! !DESCRIPTION:
    ! Compute clm_l2a_inst component of gridcell derived type. This routine computes
    ! the bare minimum of components necessary to get the first step of a
    ! run started.
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in)    :: bounds  
    type(waterstate_type) , intent(in)    :: waterstate_inst
    type(lnd2atm_type)    , intent(inout) :: lnd2atm_inst 
    !
    ! !LOCAL VARIABLES:
    integer :: g                                    ! index
    !------------------------------------------------------------------------

    call c2g(bounds, &
         waterstate_inst%h2osno_col (bounds%begc:bounds%endc), &
         lnd2atm_inst%h2osno_grc    (bounds%begg:bounds%endg), &
         c2l_scale_type= 'urbanf', l2g_scale_type='unity')

    do g = bounds%begg,bounds%endg
       lnd2atm_inst%h2osno_grc(g) = lnd2atm_inst%h2osno_grc(g)/1000._r8
    end do

  end subroutine lnd2atm_minimal

end module lnd2atmMod
