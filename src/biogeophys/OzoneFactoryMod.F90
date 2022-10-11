module OzoneFactoryMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Factory to create an instance of ozone_base_type. This module figures out the
  ! particular type to return.
  !
  ! !USES:
  use decompMod   , only : bounds_type

  implicit none
  save
  private

  !
  ! !PUBLIC ROUTINES:
  public :: create_and_init_ozone_type  ! create an object of class ozone_base_type

contains

  !-----------------------------------------------------------------------
  function create_and_init_ozone_type(bounds) result(ozone)
    !
    ! !DESCRIPTION:
    ! Create and initialize an object of ozone_base_type, and return this object
    !
    ! !USES:
    use OzoneBaseMod , only : ozone_base_type
    use OzoneOffMod  , only : ozone_off_type
    use OzoneMod     , only : ozone_type
    !
    ! !ARGUMENTS:
    class(ozone_base_type), allocatable :: ozone  ! function result
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'create_and_init_ozone_type'
    !-----------------------------------------------------------------------

    allocate(ozone, source = ozone_off_type())

    call ozone%Init(bounds)
    
  end function create_and_init_ozone_type

end module OzoneFactoryMod
