module WaterstateType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module variables for hydrology
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use LandunitType   , only : lun                
  use ColumnType     , only : col                
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type
     real(r8), pointer :: snow_depth_col         (:)   ! col snow height of snow covered area (m)
     real(r8), pointer :: h2osno_col             (:)   ! col snow water (mm H2O)
     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)

   contains

     procedure          :: Init         
     procedure          :: Restart      
     procedure, private :: InitAllocate 
     procedure, private :: InitCold

  end type waterstate_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, h2osno_input_col, snow_depth_input_col)

    class(waterstate_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(inout) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: snow_depth_input_col(bounds%begc:)

    call this%InitAllocate(bounds) 
    call this%InitCold(bounds, h2osno_input_col, snow_depth_input_col)

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
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%h2osno_col             (begc:endc))                     ; this%h2osno_col             (:)   = nan   
    allocate(this%q_ref2m_patch          (begp:endp))                     ; this%q_ref2m_patch          (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, h2osno_input_col, snow_depth_input_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use ncdio_pio       , only : ncd_io
    !
    ! !ARGUMENTS:
    class(waterstate_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer :: c
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_input_col) == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))

    ! The first three arrays are initialized from the input argument
    do c = bounds%begc,bounds%endc
       this%h2osno_col(c)             = h2osno_input_col(c)
       this%snow_depth_col(c)         = snow_depth_input_col(c)
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_io, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(waterstate_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer  :: c,l,j
    logical  :: readvar
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='H2OSNO', xtype=ncd_double,  &
         dim1name='column', &
         long_name='snow water', units='kg/m2', &
         interpinic_flag='interp', readvar=readvar, data=this%h2osno_col)

    call restartvar(ncid=ncid, flag=flag, varname='SNOW_DEPTH', xtype=ncd_double,  & 
         dim1name='column', &
         long_name='snow depth', units='m', &
         interpinic_flag='interp', readvar=readvar, data=this%snow_depth_col) 

  end subroutine Restart

end module WaterstateType
