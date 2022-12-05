module SurfaceAlbedoType

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : numrad
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC DATA MEMBERS:
  type, public :: surfalb_type

     real(r8), pointer :: albd_patch(:,:)  ! patch surface albedo (direct)   (numrad)
     real(r8), pointer :: albi_patch(:,:)  ! patch surface albedo (diffuse)  (numrad)

   contains

     procedure, public  :: Init         
     procedure, private :: InitAllocate 
     procedure, private :: InitCold     
     procedure, public  :: Restart      

  end type surfalb_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate module variables and data structures
    !
    ! !USES:
    use shr_infnan_mod, only: nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%albd_patch(begp:endp,numrad)); this%albd_patch(:,:) = nan
    allocate(this%albi_patch(begp:endp,numrad)); this%albi_patch(:,:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! Initialize module surface albedos to reasonable values
    !
    ! !ARGUMENTS:
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%albd_patch     (begp:endp, :) = 0.2_r8
    this%albi_patch     (begp:endp, :) = 0.2_r8

  end subroutine InitCold
   
  !---------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use decompMod  , only : bounds_type
    use ncdio_pio  , only : file_desc_t, ncd_double
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(surfalb_type)               :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !---------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='albd', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (direct) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albd_patch)

    call restartvar(ncid=ncid, flag=flag, varname='albi', xtype=ncd_double,  &
         dim1name='pft', dim2name='numrad', switchdim=.true., &
         long_name='surface albedo (diffuse) (0 to 1)', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%albi_patch)

  end subroutine Restart

end module SurfaceAlbedoType
