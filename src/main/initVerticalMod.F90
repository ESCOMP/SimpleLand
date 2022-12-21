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
  use clm_varpar        , only : nlevsno, nlevgrnd, nlevlak
  use clm_varpar        , only : nlevsoi, nlevsoifl, nlevurb 
  use clm_varctl        , only : fsurdat, iulog
  use clm_varcon        , only : zsoi, dzsoi, zisoi, spval, ispval, grlnd 
  use column_varcon     , only : icol_roof, icol_sunwall, icol_shadewall
  use landunit_varcon   , only : istdlak
  use fileutils         , only : getfil
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use glcBehaviorMod    , only : glc_behavior_type
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
  subroutine initVertical(bounds, glc_behavior, thick_wall, thick_roof)
    !
    ! !ARGUMENTS:
    type(bounds_type)   , intent(in)    :: bounds
    type(glc_behavior_type), intent(in) :: glc_behavior
    real(r8)            , intent(in)    :: thick_wall(bounds%begl:)
    real(r8)            , intent(in)    :: thick_roof(bounds%begl:)
    !
    ! LOCAL VARAIBLES:
    integer               :: c,l,g,i,j,lev     ! indices 
    type(file_desc_t)     :: ncid              ! netcdf id
    logical               :: readvar 
    integer               :: dimid             ! dimension id
    character(len=256)    :: locfn             ! local filename
    real(r8)              :: slope0            ! temporary
    real(r8)              :: slopebeta         ! temporary
    real(r8)              :: slopemax          ! temporary
    integer               :: ier               ! error status
    real(r8)              :: scalez = 0.025_r8 ! Soil layer thickness discretization (m)
    real(r8), allocatable :: zurb_wall(:,:)    ! wall (layer node depth)
    real(r8), allocatable :: zurb_roof(:,:)    ! roof (layer node depth)
    real(r8), allocatable :: dzurb_wall(:,:)   ! wall (layer thickness)
    real(r8), allocatable :: dzurb_roof(:,:)   ! roof (layer thickness)
    real(r8), allocatable :: ziurb_wall(:,:)   ! wall (layer interface)
    real(r8), allocatable :: ziurb_roof(:,:)   ! roof (layer interface)
    real(r8)              :: depthratio        ! ratio of lake depth to standard deep lake depth 
    integer               :: begc, endc
    integer               :: begl, endl

    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begl = bounds%begl; endl= bounds%endl

    SHR_ASSERT_ALL((ubound(thick_wall)  == (/endl/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(thick_roof)  == (/endl/)), errMsg(sourcefile, __LINE__))

    ! Open surface dataset to read in data below 

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile (ncid, locfn, 0)

    ! --------------------------------------------------------------------
    ! Define layer structure for soil, lakes, urban walls and roof 
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

    if (nlevurb > 0) then
       allocate(zurb_wall(bounds%begl:bounds%endl,nlevurb),    &
            zurb_roof(bounds%begl:bounds%endl,nlevurb),    &
            dzurb_wall(bounds%begl:bounds%endl,nlevurb),   &
            dzurb_roof(bounds%begl:bounds%endl,nlevurb),   &
            ziurb_wall(bounds%begl:bounds%endl,0:nlevurb), &
            ziurb_roof(bounds%begl:bounds%endl,0:nlevurb), &
            stat=ier)
       if (ier /= 0) then
          call shr_sys_abort(' ERROR allocation error for '//&
               'zurb_wall,zurb_roof,dzurb_wall,dzurb_roof,ziurb_wall,ziurb_roof'//&
               errMsg(sourcefile, __LINE__))
       end if
    end if

    ! Column level initialization for urban wall and roof layers and interfaces
    do l = bounds%begl,bounds%endl

       ! "0" refers to urban wall/roof surface and "nlevsoi" refers to urban wall/roof bottom
       if (lun%urbpoi(l)) then
             do j = 1, nlevurb
                zurb_wall(l,j) = (j-0.5)*(thick_wall(l)/float(nlevurb))  !node depths
             end do
             do j = 1, nlevurb
                zurb_roof(l,j) = (j-0.5)*(thick_roof(l)/float(nlevurb))  !node depths
             end do

             dzurb_roof(l,1) = 0.5*(zurb_roof(l,1)+zurb_roof(l,2))    !thickness b/n two interfaces
             do j = 2,nlevurb-1
                dzurb_roof(l,j)= 0.5*(zurb_roof(l,j+1)-zurb_roof(l,j-1)) 
             enddo
             dzurb_roof(l,nlevurb) = zurb_roof(l,nlevurb)-zurb_roof(l,nlevurb-1)

             dzurb_wall(l,1) = 0.5*(zurb_wall(l,1)+zurb_wall(l,2))    !thickness b/n two interfaces
             do j = 2,nlevurb-1
                dzurb_wall(l,j)= 0.5*(zurb_wall(l,j+1)-zurb_wall(l,j-1)) 
             enddo
             dzurb_wall(l,nlevurb) = zurb_wall(l,nlevurb)-zurb_wall(l,nlevurb-1)

             ziurb_wall(l,0) = 0.
             do j = 1, nlevurb-1
                ziurb_wall(l,j) = 0.5*(zurb_wall(l,j)+zurb_wall(l,j+1))          !interface depths
             enddo
             ziurb_wall(l,nlevurb) = zurb_wall(l,nlevurb) + 0.5*dzurb_wall(l,nlevurb)

             ziurb_roof(l,0) = 0.
             do j = 1, nlevurb-1
                ziurb_roof(l,j) = 0.5*(zurb_roof(l,j)+zurb_roof(l,j+1))          !interface depths
             enddo
             ziurb_roof(l,nlevurb) = zurb_roof(l,nlevurb) + 0.5*dzurb_roof(l,nlevurb)
       end if
    end do

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)

       if (lun%urbpoi(l)) then
          if (col%itype(c)==icol_sunwall .or. col%itype(c)==icol_shadewall) then
             col%z(c,1:nlevurb)  = zurb_wall(l,1:nlevurb)
             col%zi(c,0:nlevurb) = ziurb_wall(l,0:nlevurb)
             col%dz(c,1:nlevurb) = dzurb_wall(l,1:nlevurb)
             if (nlevurb < nlevgrnd) then
                col%z(c,nlevurb+1:nlevgrnd)  = spval
                col%zi(c,nlevurb+1:nlevgrnd) = spval
                col%dz(c,nlevurb+1:nlevgrnd) = spval
             end if
          else if (col%itype(c)==icol_roof) then
             col%z(c,1:nlevurb)  = zurb_roof(l,1:nlevurb)
             col%zi(c,0:nlevurb) = ziurb_roof(l,0:nlevurb)
             col%dz(c,1:nlevurb) = dzurb_roof(l,1:nlevurb)
             if (nlevurb < nlevgrnd) then
                col%z(c,nlevurb+1:nlevgrnd)  = spval
                col%zi(c,nlevurb+1:nlevgrnd) = spval
                col%dz(c,nlevurb+1:nlevgrnd) = spval
             end if
          else
             col%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
             col%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
             col%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
          end if
       else if (lun%itype(l) /= istdlak) then
          col%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          col%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          col%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
       end if
    end do

    if (nlevurb > 0) then
       deallocate(zurb_wall, zurb_roof, dzurb_wall, dzurb_roof, ziurb_wall, ziurb_roof)
    end if

    do c = bounds%begc,bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istdlak) then
          col%z(c,1:nlevgrnd)  = zsoi(1:nlevgrnd)
          col%zi(c,0:nlevgrnd) = zisoi(0:nlevgrnd)
          col%dz(c,1:nlevgrnd) = dzsoi(1:nlevgrnd)
       end if
    end do

    call ncd_pio_closefile(ncid)

  end subroutine initVertical

end module initVerticalMod
