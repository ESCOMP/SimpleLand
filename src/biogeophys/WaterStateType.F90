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
  use clm_varpar     , only : nlevgrnd, nlevurb, nlevsno   
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
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)

     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)

     ! Balance Checks
     real(r8), pointer :: endwb_col              (:)   ! water mass end of the time step

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
  subroutine Init(this, bounds, &
       h2osno_input_col, snow_depth_input_col, t_soisno_col)

    class(waterstate_type)            :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(inout) :: h2osno_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: snow_depth_input_col(bounds%begc:)
    real(r8)          , intent(inout) :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)

#ifdef __PGI
# if __PGIC__ == 14 && __PGIC_MINOR__ == 7
    ! COMPILER_BUG(bja, 2015-04, pgi 14.7-?) occurs at: call this%InitCold(...)
    ! PGF90-F-0000-Internal compiler error. normalize_forall_array: non-conformable
    ! not sure why this fixes things....
    real(r8), allocatable :: workaround_for_pgi_internal_compiler_error(:)
# endif
#endif

    call this%InitAllocate(bounds) 

    call this%InitCold(bounds, h2osno_input_col, snow_depth_input_col, t_soisno_col)

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
    integer :: begl, endl
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl
    begg = bounds%begg; endg= bounds%endg

    allocate(this%snow_depth_col         (begc:endc))                     ; this%snow_depth_col         (:)   = nan
    allocate(this%h2osno_col             (begc:endc))                     ; this%h2osno_col             (:)   = nan   
    allocate(this%tws_grc                (begg:endg))                     ; this%tws_grc                (:)   = nan

    allocate(this%q_ref2m_patch          (begp:endp))                     ; this%q_ref2m_patch          (:)   = nan

    allocate(this%endwb_col              (begc:endc))                     ; this%endwb_col              (:)   = nan
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       h2osno_input_col, snow_depth_input_col, t_soisno_col)
    !
    ! !DESCRIPTION:
    ! Initialize time constant variables and cold start conditions
    !
    ! !USES:
    use shr_log_mod     , only : errMsg => shr_log_errMsg
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_TKFRZ
    use clm_varpar      , only : nlevsoi, nlevgrnd, nlevsno, nlevurb
    use landunit_varcon , only : istwet, istsoil, istcrop, istice_mec
    use column_varcon   , only : icol_road_perv
    use column_varcon   , only : icol_road_imperv
    use clm_varcon      , only : denice, denh2o, spval, bdsno
    use clm_varcon      , only : tfrz, spval
    use spmdMod         , only : masterproc
    use abortutils      , only : endrun
    use fileutils       , only : getfil
    use ncdio_pio       , only : file_desc_t, ncd_io
    !
    ! !ARGUMENTS:
    class(waterstate_type)                :: this
    type(bounds_type)     , intent(in)    :: bounds
    real(r8)              , intent(in)    :: h2osno_input_col(bounds%begc:)
    real(r8)              , intent(in)    :: snow_depth_input_col(bounds%begc:)
    ! volumetric soil water at saturation (porosity)
    real(r8)              , intent(in)    :: t_soisno_col(bounds%begc:, -nlevsno+1:) ! col soil temperature (Kelvin)
    !
    ! !LOCAL VARIABLES:
    integer            :: p,c,j,l,g,lev,nlevs
    real(r8)           :: d
    type(file_desc_t)  :: ncid
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(h2osno_input_col)     == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(snow_depth_input_col) == (/bounds%endc/))          , errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno_col)         == (/bounds%endc,nlevgrnd/)) , errMsg(sourcefile, __LINE__))

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
    use spmdMod          , only : masterproc
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
    integer  :: c,l,j,nlevs
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
