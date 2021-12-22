module CNDriverMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use decompMod                       , only : bounds_type
  use perf_mod                        , only : t_startf, t_stopf
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNDriverInit         ! Ecosystem dynamics: initialization
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNDriverInit(bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds      
    character(len=*)                       , intent(in)    :: NLFilename     ! Namelist filename
    !-----------------------------------------------------------------------
    
  end subroutine CNDriverInit

end  module CNDriverMod
