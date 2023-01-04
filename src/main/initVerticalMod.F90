module initVerticalMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize vertical components of column datatype
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_infnan_mod    , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_abort
  use decompMod         , only : bounds_type
  use spmdMod           , only : masterproc
  use clm_varpar        , only : nlevgrnd
  use clm_varpar        , only : nlevsoi
  use clm_varctl        , only : iulog
  use clm_varcon        , only : zsoi, dzsoi, zisoi, spval
  use abortUtils        , only : endrun    
  use ncdio_pio
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initVertical
  ! !PRIVATE MEMBER FUNCTIONS:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine initVertical(bounds)
    !
    ! !ARGUMENTS:
    type(bounds_type)   , intent(in)    :: bounds
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,g,i,j,lev     ! indices 
    real(r8)              :: scalez = 0.025_r8 ! Soil layer thickness discretization (m)

    !------------------------------------------------------------------------

    ! --------------------------------------------------------------------
    ! Define layer structure for soil
    ! Vertical profile of snow is not initialized here - but below
    ! --------------------------------------------------------------------
    
    ! Soil layers and interfaces (assumed same for all non-lake patches)
    ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil
       do j = 1, nlevgrnd
          zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
       enddo

       dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
       do j = 2,nlevgrnd-1
          dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
       enddo
       dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)

       zisoi(0) = 0._r8
       do j = 1, nlevgrnd-1
          zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
       enddo
       zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8*dzsoi(nlevgrnd)

    if (masterproc) then
       write(iulog, *) 'zsoi', zsoi(:) 
       write(iulog, *) 'zisoi: ', zisoi(:)
       write(iulog, *) 'dzsoi: ', dzsoi(:)
    end if

  end subroutine initVertical

end module initVerticalMod
