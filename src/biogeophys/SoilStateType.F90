module SoilStateType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use clm_varpar      , only : nlevgrnd
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: soilstate_type

     ! hydraulic properties
     real(r8), pointer :: watsat_col(:,:)  ! col volumetric soil water at saturation (porosity) 

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 

  end type soilstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(soilstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    allocate(this%watsat_col(begc:endc,nlevgrnd)); this%watsat_col(:,:) = nan

  end subroutine InitAllocate

end module SoilStateType
