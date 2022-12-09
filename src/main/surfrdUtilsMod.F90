module surfrdUtilsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains utility methods that can be used when reading surface datasets or similar
  ! datasets (such as the landuse_timeseries dataset)
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varctl   , only : iulog
  use abortutils   , only : endrun
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use spmdMod      , only : masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: check_sums_equal_1  ! Confirm that sum(arr(n,:)) == 1 for all n
  public :: renormalize         ! Renormalize an array

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  !-----------------------------------------------------------------------

contains
  
  !-----------------------------------------------------------------------
  subroutine check_sums_equal_1(arr, lb, name, caller, ier)
    !
    ! !DESCRIPTION:
    ! Confirm that sum(arr(n,:)) == 1 for all n. If this isn't true for any n, abort with a message.
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lb           ! lower bound of the first dimension of arr
    real(r8)        , intent(in) :: arr(lb:,:)   ! array to check
    character(len=*), intent(in) :: name         ! name of array
    character(len=*), intent(in) :: caller       ! identifier of caller, for more meaningful error messages
    integer, optional, intent(out):: ier         ! Return an error code rather than abort
    !
    ! !LOCAL VARIABLES:
    logical :: found
    integer :: nl
    integer :: nindx
    real(r8), parameter :: eps = 1.e-13_r8
    !-----------------------------------------------------------------------

    if( present(ier) ) ier = 0
    found = .false.

    do nl = lbound(arr, 1), ubound(arr, 1)
       if (abs(sum(arr(nl,:)) - 1._r8) > eps) then
          found = .true.
          nindx = nl
          exit
       end if
    end do

    if (found) then
       write(iulog,*) trim(caller), ' ERROR: sum of ', trim(name), ' not 1.0 at nl=', nindx
       write(iulog,*) 'sum is: ', sum(arr(nindx,:))
       if( present(ier) ) then
          ier = -10
       else
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine check_sums_equal_1

  !-----------------------------------------------------------------------
  subroutine renormalize(arr, lb, normal)
    !
    ! !DESCRIPTION:
    ! Re normalize an array so that it sums to the input value
    !
    ! !ARGUMENTS:
    integer         , intent(in) :: lb            ! lower bound of the first dimension of arr
    real(r8)        , intent(inout) :: arr(lb:,:) ! array to check
    real(r8)        , intent(in) :: normal        ! normal to sum to
    !
    ! !LOCAL VARIABLES:
    integer :: nl        ! Array index
    real(r8) :: arr_sum  ! sum of array
    real(r8) :: ratio    ! ratio to multiply by
    !-----------------------------------------------------------------------

    do nl = lbound(arr, 1), ubound(arr, 1)
       arr_sum = sum(arr(nl,:))
       if ( arr_sum /= 0.0_r8 )then
          ratio     = normal / arr_sum
          arr(nl,:) = arr(nl,:) * ratio
       end if
    end do

  end subroutine renormalize

end module surfrdUtilsMod
