module CropType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing variables needed for the crop model
  !
  ! TODO(wjs, 2014-08-05) Move more crop-specific variables into here - many are
  ! currently in CNVegStateType
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use spmdMod             , only : masterproc
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use clm_varcon          , only : spval
  use clm_varctl          , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC DATA TYPES:
  !
  ! Crop state variables structure
  type, public :: crop_type

     ! Note that cropplant and harvdate could be 2D to facilitate rotation
     integer , pointer :: nyrs_crop_active_patch  (:)   ! number of years this crop patch has been active (0 for non-crop patches)
     logical , pointer :: croplive_patch          (:)   ! patch Flag, true if planted, not harvested
     logical , pointer :: cropplant_patch         (:)   ! patch Flag, true if planted
     integer , pointer :: harvdate_patch          (:)   ! patch harvest date
     real(r8), pointer :: fertnitro_patch         (:)   ! patch fertilizer nitrogen
     real(r8), pointer :: gddplant_patch          (:)   ! patch accum gdd past planting date for crop       (ddays)
     real(r8), pointer :: gddtsoi_patch           (:)   ! patch growing degree-days from planting (top two soil layers) (ddays)
     real(r8), pointer :: vf_patch                (:)   ! patch vernalization factor for cereal
     real(r8), pointer :: cphase_patch            (:)   ! phenology phase
     real(r8), pointer :: latbaset_patch          (:)   ! Latitude vary baset for gddplant (degree C)
     character(len=20) :: baset_mapping
     real(r8) :: baset_latvary_intercept
     real(r8) :: baset_latvary_slope

   contains
     ! Public routines
     procedure, public  :: Init               ! Initialize the crop type
     procedure, public  :: InitAccVars

     ! Private routines
     procedure, private :: InitAllocate 
     procedure, private, nopass :: checkDates

  end type crop_type

  character(len=*), parameter, private :: baset_map_constant = 'constant'
  character(len=*), parameter, private :: baset_map_latvary  = 'varytropicsbylat'
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------
    
    call this%InitAllocate(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    ! !USES:
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%nyrs_crop_active_patch(begp:endp)) ; this%nyrs_crop_active_patch(:) = 0
    allocate(this%croplive_patch (begp:endp)) ; this%croplive_patch (:) = .false.
    allocate(this%cropplant_patch(begp:endp)) ; this%cropplant_patch(:) = .false.
    allocate(this%harvdate_patch (begp:endp)) ; this%harvdate_patch (:) = huge(1) 
    allocate(this%fertnitro_patch (begp:endp)) ; this%fertnitro_patch (:) = spval
    allocate(this%gddplant_patch (begp:endp)) ; this%gddplant_patch (:) = spval
    allocate(this%gddtsoi_patch  (begp:endp)) ; this%gddtsoi_patch  (:) = spval
    allocate(this%vf_patch       (begp:endp)) ; this%vf_patch       (:) = 0.0_r8
    allocate(this%cphase_patch   (begp:endp)) ; this%cphase_patch   (:) = 0.0_r8
    allocate(this%latbaset_patch (begp:endp)) ; this%latbaset_patch (:) = spval

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES:
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    
    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level patch field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg=" allocation error for rbufslp"//&
            errMsg(sourcefile, __LINE__))
    endif

    nstep = get_nstep()

    call extract_accum_field ('GDDPLANT', rbufslp, nstep) 
    this%gddplant_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('GDDTSOI', rbufslp, nstep) 
    this%gddtsoi_patch(begp:endp)  = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine checkDates( )
    !
    ! !DESCRIPTION: 
    ! Make sure the dates are compatible. The date given to startup the model
    ! and the date on the restart file must be the same although years can be
    ! different. The dates need to be checked when the restart file is being
    ! read in for a startup or branch case (they are NOT allowed to be different
    ! for a restart case).
    !
    ! For the prognostic crop model the date of planting is tracked and growing
    ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
    ! messes up these bits of saved information.
    !
    ! !ARGUMENTS:
    use clm_time_manager, only : get_driver_start_ymd, get_start_date
    use clm_varctl      , only : iulog
    use clm_varctl      , only : nsrest, nsrBranch, nsrStartup
    !
    ! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !-----------------------------------------------------------------------
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( masterproc ) &
            write(iulog,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(iulog,formDate) 'Start date is: ', styr, stmon_day/100, &
               (stmon_day - stmon_day/100)
          call endrun(msg=' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
               ' and the date on the restart file needs to match (years can be different)'//&
               errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine checkDates

end module CropType

