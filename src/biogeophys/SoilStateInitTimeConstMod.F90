module SoilStateInitTimeConstMod

  !------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Set hydraulic and thermal properties 
  !
  ! !USES
  use SoilStateType , only : soilstate_type
  use LandunitType  , only : lun
  use ColumnType    , only : col
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: SoilStateInitTimeConst
  !
  ! !PRIVATE DATA:
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------
  !
contains

  !-----------------------------------------------------------------------
  subroutine SoilStateInitTimeConst(bounds, soilstate_inst) 
    !
    ! !USES:
    use shr_kind_mod        , only : r8 => shr_kind_r8
    use decompMod           , only : bounds_type
    use landunit_varcon     , only : istdlak, istwet, istice_mec
    use column_varcon       , only : icol_road_perv, icol_road_imperv
    use clm_varpar          , only : nlevgrnd
    use clm_varcon          , only : spval
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds  
    type(soilstate_type) , intent(inout) :: soilstate_inst
    !
    ! !LOCAL VARIABLES:
    integer            :: lev, c, l
    integer            :: begc, endc
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    do c = begc, endc
       l = col%landunit(c)
       if (lun%itype(l)==istwet .or. lun%itype(l)==istice_mec) then

          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = spval
          end do

       else if (lun%urbpoi(l) .and. (col%itype(c) /= icol_road_perv) .and. (col%itype(c) /= icol_road_imperv) )then

          ! Urban Roof, sunwall, shadewall properties set to special value
          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = spval
          end do

       else

          do lev = 1,nlevgrnd
             if (lun%itype(l) /= istdlak) then  ! soil columns of both urban and non-urban types
                soilstate_inst%watsat_col(c,lev) = 1._r8
             end if
          end do
       end if
    end do
    do c = begc, endc
       l = col%landunit(c)
       if (lun%itype(l)==istdlak) then
          do lev = 1,nlevgrnd
             soilstate_inst%watsat_col(c,lev) = 0.489_r8
          end do
       endif
    end do

  end subroutine SoilStateInitTimeConst

end module SoilStateInitTimeConstMod
