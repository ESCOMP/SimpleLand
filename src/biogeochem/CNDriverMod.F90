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
  subroutine CNDriverInit(bounds, NLFilename, cnfire_method)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    use CNPhenologyMod              , only : CNPhenologyInit
    use CNFireMethodMod             , only : cnfire_method_type
    use SoilBiogeochemCompetitionMod, only : SoilBiogeochemCompetitionInit
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)    :: bounds      
    character(len=*)                       , intent(in)    :: NLFilename     ! Namelist filename
    class(cnfire_method_type)              , intent(inout) :: cnfire_method 
    !-----------------------------------------------------------------------
    call SoilBiogeochemCompetitionInit(bounds)
    call CNPhenologyInit(bounds)
    call cnfire_method%CNFireInit(bounds, NLFilename)
    
  end subroutine CNDriverInit

end  module CNDriverMod
