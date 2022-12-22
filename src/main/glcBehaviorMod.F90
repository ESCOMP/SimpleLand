module glcBehaviorMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Determines a number of aspects of the behavior of glacier_mec classes in each grid
  ! cell.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use clm_varctl     , only : iulog
  use landunit_varcon, only : istice_mec
  use clm_instur     , only : wt_lunit, wt_glc_mec
  use decompMod      , only : bounds_type
  use filterColMod   , only : filter_col_type
  use ColumnType     , only : col

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: glc_behavior_type
     private

   contains

     ! ------------------------------------------------------------------------
     ! Public routines
     ! ------------------------------------------------------------------------

     ! get number of subgrid units in glc_mec landunit on one grid cell
     procedure, public  :: get_num_glc_mec_subgrid

     ! returns true if memory should be allocated for the given glc_mec column, and its
     ! weight on the landunit
     procedure, public  :: glc_mec_col_exists

  end type glc_behavior_type

  ! !PRIVATE MEMBER DATA:

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !-----------------------------------------------------------------------
  subroutine get_num_glc_mec_subgrid(this, gi, npatches, ncols, nlunits)
    !
    ! !DESCRIPTION:
    ! Get number of subgrid units in glc_mec landunit on one grid cell
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    integer , intent(in)  :: gi       ! grid cell index
    integer , intent(out) :: npatches ! number of glacier_mec patches in this grid cell
    integer , intent(out) :: ncols    ! number of glacier_mec columns in this grid cell
    integer , intent(out) :: nlunits  ! number of glacier_mec landunits in this grid cell
    !
    ! !LOCAL VARIABLES:
    integer  :: m  ! loop index
    logical  :: col_exists
    real(r8) :: col_wt_lunit

    character(len=*), parameter :: subname = 'get_num_glc_mec_subgrid'
    !-----------------------------------------------------------------------

    ncols = 0

    do m = 1, 10
       call this%glc_mec_col_exists(gi = gi, elev_class = m, &
            exists = col_exists, col_wt_lunit = col_wt_lunit)
       if (col_exists) then
          ncols = ncols + 1
       end if
    end do

    if (ncols > 0) then
       npatches = ncols
       nlunits = 1
    else
       npatches = 0
       nlunits = 0
    end if
  
  end subroutine get_num_glc_mec_subgrid

  !-----------------------------------------------------------------------
  subroutine glc_mec_col_exists(this, gi, elev_class, exists, col_wt_lunit)
    !
    ! !DESCRIPTION:
    ! For the given glc_mec column, with elevation class index elev_class, in grid cell
    ! gi: sets exists to true if memory should be allocated for this column, and sets
    ! col_wt_lunit to the column's weight on the icemec landunit.
    !
    ! If exists is false, then col_wt_lunit is arbitrary and should be ignored.
    !
    ! !ARGUMENTS:
    class(glc_behavior_type), intent(in) :: this
    integer,  intent(in)  :: gi           ! grid cell index
    integer,  intent(in)  :: elev_class   ! elevation class index
    logical,  intent(out) :: exists       ! whether memory should be allocated for this column
    real(r8), intent(out) :: col_wt_lunit ! column's weight on the icemec landunit
    !
    ! !LOCAL VARIABLES:
    integer :: atm_elev_class ! elevation class corresponding to atmosphere topographic height
    integer :: err_code

    character(len=*), parameter :: subname = 'glc_mec_col_exists'
    !-----------------------------------------------------------------------

    ! Set default outputs
    exists = .false.
    col_wt_lunit = wt_glc_mec(gi, elev_class)

       if (wt_lunit(gi, istice_mec) > 0.0_r8 .and. &
            wt_glc_mec(gi, elev_class) > 0.0_r8) then
          ! If the landunit has non-zero weight on the grid cell, and this column has
          ! non-zero weight on the landunit...
          exists = .true.
       end if

  end subroutine glc_mec_col_exists

end module glcBehaviorMod
