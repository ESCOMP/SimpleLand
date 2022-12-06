module EnergyFluxType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! Energy flux data structure
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use ColumnType     , only : col                
  use PatchType      , only : patch                
  !
  implicit none
  save
  private
  !
  type, public :: energyflux_type

     ! Fluxes
     real(r8), pointer :: eflx_lwrad_out_patch    (:)   ! patch emitted infrared (longwave) radiation (W/m**2)

   contains

     procedure, public  :: Init            ! Public initialization method
     procedure, private :: InitAllocate    ! initialize/allocate
     procedure, private :: InitCold        ! initialize for cold start
     procedure, public  :: Restart         ! setup restart fields

  end type energyflux_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    !    Allocate and initialize the data type and setup history, and initialize for cold-start.
    ! !USES:
    implicit none
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    call this%InitAllocate ( bounds )
    call this%InitCold ( bounds, t_grnd_col) 

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize and allocate data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    allocate( this%eflx_lwrad_out_patch    (begp:endp))             ; this%eflx_lwrad_out_patch    (:)   = nan

  end subroutine InitAllocate
    
  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, t_grnd_col)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod    , only : r8 => shr_kind_r8
    use clm_varcon      , only : sb
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type)         :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: t_grnd_col( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(t_grnd_col) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    ! Patches
    do p = bounds%begp, bounds%endp 
       c = patch%column(p)
       this%eflx_lwrad_out_patch(p) = sb * (t_grnd_col(c))**4
    end do

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t, ncd_double
    use restUtilMod
    implicit none
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !-----------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='EFLX_LWRAD_OUT', xtype=ncd_double,  & 
         dim1name='pft', &
         long_name='emitted infrared (longwave) radiation', units='watt/m^2', &
         interpinic_flag='interp', readvar=readvar, data=this%eflx_lwrad_out_patch)

  end subroutine Restart

end module EnergyFluxType
