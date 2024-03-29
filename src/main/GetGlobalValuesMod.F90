module GetGlobalValuesMod

  !-----------------------------------------------------------------------
  ! Obtain and Write Global Index information
  !-----------------------------------------------------------------------
  implicit none
  private

  ! PUBLIC MEMBER FUNCTIONS:

  public :: GetGlobalIndex
  public :: GetGlobalWrite

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  integer function GetGlobalIndex(decomp_index, clmlevel)

    !----------------------------------------------------------------
    ! Description
    ! Determine global index space value for target point at given clmlevel
    !
    ! Uses:
    use shr_log_mod, only: errMsg => shr_log_errMsg
    use decompMod  , only: bounds_type, get_clmlevel_gsmap, get_proc_bounds
    use spmdMod    , only: iam
    use clm_varcon , only: nameg
    use clm_varctl , only: iulog
    use mct_mod    , only: mct_gsMap, mct_gsMap_orderedPoints
    use shr_sys_mod, only: shr_sys_abort
    !
    ! Arguments 
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    type(bounds_type)             :: bounds_proc   ! processor bounds
    type(mct_gsMap),pointer       :: gsmap         ! global seg map
    integer, pointer,dimension(:) :: gsmap_ordered ! gsmap ordered points
    integer                       :: beg_index     ! beginning proc index for clmlevel
    !----------------------------------------------------------------

    call get_proc_bounds(bounds_proc)

    if (trim(clmlevel) == nameg) then
       beg_index = bounds_proc%begg
    else
       call shr_sys_abort('clmlevel of '//trim(clmlevel)//' not supported' // &
            errmsg(sourcefile, __LINE__))
    end if

    call get_clmlevel_gsmap(clmlevel=trim(clmlevel), gsmap=gsmap)
    call mct_gsMap_orderedPoints(gsmap, iam, gsmap_ordered)
    GetGlobalIndex = gsmap_ordered(decomp_index - beg_index + 1)
    deallocate(gsmap_ordered)

  end function GetGlobalIndex

  !-----------------------------------------------------------------------
  subroutine GetGlobalWrite(decomp_index, clmlevel)

    !-----------------------------------------------------------------------
    ! Description:
    ! Write global index information for input local indices
    !
    use shr_sys_mod  , only : shr_sys_flush
    use shr_sys_mod  , only : shr_sys_abort
    use shr_log_mod  , only : errMsg => shr_log_errMsg
    use clm_varctl   , only : iulog
    use clm_varcon   , only : nameg
    use GridcellType , only : grc                
    !
    ! Arguments:
    integer          , intent(in) :: decomp_index
    character(len=*) , intent(in) :: clmlevel
    !
    ! Local Variables:
    integer :: igrc
    !-----------------------------------------------------------------------

    if (trim(clmlevel) == nameg) then

       igrc = decomp_index
       write(iulog,*)'local  gridcell index = ',igrc
       write(iulog,*)'global gridcell index = ',GetGlobalIndex(decomp_index=igrc, clmlevel=nameg)
       write(iulog,*)'gridcell longitude    = ',grc%londeg(igrc)
       write(iulog,*)'gridcell latitude     = ',grc%latdeg(igrc)

    else		       
       call shr_sys_abort('clmlevel '//trim(clmlevel)//'not supported '//errmsg(sourcefile, __LINE__))

    end if

    call shr_sys_flush(iulog)

  end subroutine GetGlobalWrite

end module GetGlobalValuesMod
