module AerosolMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use shr_infnan_mod   , only : nan => shr_infnan_nan, assignment(=)
  use abortutils       , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC DATA MEMBERS:
  real(r8), public, parameter :: snw_rds_min = 54.526_r8          ! minimum allowed snow effective radius (also cold "fresh snow" value) [microns]
  real(r8), public            :: fresh_snw_rds_max = 204.526_r8   ! maximum warm fresh snow effective radius [microns]
  !
  type, public :: aerosol_type

  end type aerosol_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

end module AerosolMod
