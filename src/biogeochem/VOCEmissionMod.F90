module VOCEmissionMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Volatile organic compound emission
  !
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  !
  implicit none
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC TYPES:
  type, public :: vocemis_type
  end type vocemis_type
  !
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

end module VOCEmissionMod


