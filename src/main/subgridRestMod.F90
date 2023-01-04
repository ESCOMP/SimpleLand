module subgridRestMod

#include "shr_assert.h"

  !------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use abortutils         , only : endrun
  use decompMod          , only : bounds_type, ldecomp
  use domainMod          , only : ldomain
  use clm_varcon         , only : nameg
  use pio                , only : file_desc_t
  use ncdio_pio          , only : ncd_int, ncd_double
  use GridcellType       , only : grc                
  use restUtilMod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: subgridRestWrite              ! handle restart writes of subgrid variables

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: subgridRest_write_only     ! handle restart of subgrid variables that only need to be written, not read

  ! !PRIVATE TYPES:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine subgridRestWrite(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart writes (and defines) of subgrid variables
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define or write data
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'subgridRestWrite'
    !-----------------------------------------------------------------------

    call subgridRest_write_only(bounds, ncid, flag)

  end subroutine subgridRestWrite

  !-----------------------------------------------------------------------
  subroutine subgridRest_write_only(bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart for variables that only need to be written, not read. This applies
    ! to variables that are time-constant and are only put on the restart file for the
    ! sake of having some additional metadata there.
    !
    ! Note that 'active' flags appear in this routine: they don't need to be read because
    ! they can be computed using other info on the restart file (particularly subgrid
    ! weights).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds ! bounds
    type(file_desc_t), intent(inout) :: ncid   ! netCDF dataset id
    character(len=*) , intent(in)    :: flag   ! flag to determine if define, write or read data
    !
    ! !LOCAL VARIABLES:
    integer :: g,l,c,p,i             ! indices
    logical :: readvar               ! temporary
    real(r8), pointer :: rgarr(:)    ! temporary
    integer , pointer :: igarr(:)    ! temporary

    real(r8), pointer :: temp2d_r(:,:) ! temporary for multi-level variables
    integer , pointer :: temp2d_i(:,:) ! temporary for multi-level variables

    character(len=*), parameter :: subname = 'subgridRest_write_only'
    !-----------------------------------------------------------------------
    
    !------------------------------------------------------------------
    ! Write gridcell info
    !------------------------------------------------------------------

    allocate(rgarr(bounds%begg:bounds%endg), igarr(bounds%begg:bounds%endg))

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lon', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell longitude', units='degrees_east',         &
         interpinic_flag='skip', readvar=readvar, data=grc%londeg)

    call restartvar(ncid=ncid, flag=flag, varname='grid1d_lat', xtype=ncd_double, &
         dim1name='gridcell',                                          &
         long_name='gridcell latitude', units='degrees_north',         &
         interpinic_flag='skip', readvar=readvar, data=grc%latdeg)

    do g=bounds%begg,bounds%endg
       igarr(g)= mod(ldecomp%gdc2glo(g)-1,ldomain%ni) + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_ixy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d longitude index of corresponding gridcell',     &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    do g=bounds%begg,bounds%endg
       igarr(g)= (ldecomp%gdc2glo(g) - 1)/ldomain%ni + 1
    enddo
    call restartvar(ncid=ncid, flag=flag, varname='grid1d_jxy', xtype=ncd_int,    &
         dim1name='gridcell',                                          &
         long_name='2d latitude index of corresponding gridcell',      &
         interpinic_flag='skip', readvar=readvar, data=igarr)

    deallocate(rgarr,igarr)

  end subroutine subgridRest_write_only

end module subgridRestMod
