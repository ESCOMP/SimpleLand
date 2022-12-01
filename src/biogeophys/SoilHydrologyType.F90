Module SoilHydrologyType

  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use decompMod             , only : bounds_type
  use clm_varpar            , only : nlevgrnd, nlayer, nlayert, nlevsoi 
  use clm_varcon            , only : spval
  use clm_varctl            , only : iulog
  use LandunitType          , only : lun                
  use ColumnType            , only : col                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  type, public :: soilhydrology_type

     ! NON-VIC
     real(r8), pointer :: zwt_col           (:)     ! col water table depth
     real(r8), pointer :: wa_col            (:)     ! col water in the unconfined aquifer (mm)

     ! VIC 
     real(r8), pointer :: porosity_col      (:,:)   ! col VIC porosity (1-bulk_density/soil_density)
     real(r8), pointer :: vic_clm_fract_col (:,:,:) ! col VIC fraction of VIC layers in CLM layers 
     real(r8), pointer :: depth_col         (:,:)   ! col VIC layer depth of upper layer  
     real(r8), pointer :: expt_col          (:,:)   ! col VIC pore-size distribution related paramter(Q12) 

   contains

     ! Public routines
     procedure, public  :: Init
     procedure, public  :: Restart

     ! Private routines
     procedure, private :: InitAllocate

  end type soilhydrology_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)

    class(soilhydrology_type) :: this
    type(bounds_type), intent(in)    :: bounds  
    character(len=*), intent(in) :: NLFilename

    call this%InitAllocate(bounds) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%zwt_col           (begc:endc))                 ; this%zwt_col           (:)     = nan

    allocate(this%wa_col            (begc:endc))                 ; this%wa_col            (:)     = nan

    allocate(this%depth_col         (begc:endc,nlayert))         ; this%depth_col         (:,:)   = nan
    allocate(this%porosity_col      (begc:endc,nlayer))          ; this%porosity_col      (:,:)   = nan
    allocate(this%vic_clm_fract_col (begc:endc,nlayer, nlevsoi)) ; this%vic_clm_fract_col (:,:,:) = nan
    allocate(this%expt_col          (begc:endc,nlayer))          ; this%expt_col          (:,:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(soilhydrology_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='WA', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water in the unconfined aquifer', units='mm', &
         interpinic_flag='interp', readvar=readvar, data=this%wa_col)

    call restartvar(ncid=ncid, flag=flag, varname='ZWT', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='water table depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%zwt_col)

  end subroutine Restart

end Module SoilHydrologyType
